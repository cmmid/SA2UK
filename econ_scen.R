suppressPackageStartupMessages({
    require(data.table)
    require(countrycode)
})

.debug <- c(".", "NGA", "baseline")
#.debug <- c("~/Dropbox/Covid-WHO-vax", "00714")
.args <- if (interactive()) sprintf(c(
    "covid_other_costs.csv",
    "covid_vac_costs_per_dose.csv",
    "daly_scenarios.csv",
    if(.debug[3] == "baseline") "%s/outputs/acc_scen/%s_0.rds" else c("%s/outputs/acc_scen", "%s/outputs/acc_scen/%s_baseline.rds"),
    .debug[2],
    if(.debug[3] == "baseline") "%s/outputs/econ_scen/%s_baseline.rds" else "%s/outputs/econ_scen/%s.rds"
),.debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

isbaseline <- grepl("_baseline\\.rds$", tail(.args,1))
tariso3 <- "NGA"

scen.dt <- rbind(
  data.table(vax_mech = "none", vax_eff = 0, coverage = 0),
  data.table(expand.grid(
    vax_mech = c("infection", "disease"),
    vax_eff = c(.5, .9),
    coverage = c(.25, .50)
  ))
)[, id := (1:.N)-1 ]

totaldoses <- function(
  coverage, ages, pop = readRDS("./inputs/pops/NGA.rds")$pop[[1]]$size[ages]
) round(pop*coverage)

dt <- (if (isbaseline) {
    readRDS(.args[4])
} else {
  rbindlist(lapply(list.files(.args[4], "_[1-8]\\.rds$", full.names = TRUE), readRDS))
})[
  iyear > 0 & compartment != "R", .(id = epi_id, sampleId = sample, anni_year = iyear, outcome = compartment, age = group, value )
]

# scn[, .(id, vax_delay, strategy_str, doses_per_day)
wide.dt <- dcast(dt, id + age + sampleId + anni_year ~ outcome, value.var = "value")

wide.dt[scen.dt, on=.(id), doses := fifelse(between(age, 4, 16) & anni_year == 1, totaldoses(coverage, age), 0)]

rm(dt)

#' SETUP ECON DATA

#' vector of doses

# by perspective
othercosts <- dcast(
  fread(.args[1])[,
    iso3 := countrycode(country, "country.name", "iso3c")
  ][tariso3 == iso3 & perspective == "health_system"],
  perspective ~ name, value.var = "cost"
)

vac_cost.dt <- fread(.args[2])[,
  iso3 := countrycode(country, "country.name", "iso3c")
][tariso3 == iso3][scenario == "blended"][, cost_vac_dose := as.numeric(gsub("\\$", "", cost_vac_dose)) ]

dalys.dt <- fread(.args[3])[tariso3 == iso3]
dalys.dt[, age := age_cat ]
dalys.dt$age_cat <- NULL

#' econ scenarios: 3x perspective, 3x vac costs

econscns.dt <- data.table(expand.grid(
    perspective = othercosts[, unique(perspective)],
    vac_price = vac_cost.dt[, unique(vac_price)],
    daly_scenario = dalys.dt[, unique(daly_scenario)],
    disc.costs = dalys.dt[, max(disc_rate)],
    disc.dalys = dalys.dt[, unique(disc_rate)]
))[, econ_id := 1:.N ]

econ_digestor <- function(epi.dt, dalys.dt, econ_pars){
    
    # divide one-off annual/daily costs across age categories
    age_cats <- epi.dt[, max(age)]
    econ_pars[["cost_hs_day_erm"]] <- econ_pars[["cost_hs_day_erm"]] / age_cats
    econ_pars[["cost_hs_one_erm"]] <- econ_pars[["cost_hs_one_erm"]] / age_cats
    econ_pars[["cost_hs_day_comms"]] <- econ_pars[["cost_hs_day_comms"]] / age_cats
    
    epi.dt[, costs := with(econ_pars,
        # one-off / daily health system response costs
        365 * (cost_hs_day_erm + cost_hs_day_comms) + 
        cost_hs_one_erm +
        # 10 % of symptomatic cases tested, 7 contacts per tested case
        # TODO: how do these assumptions scale with prevalence
        cases * 0.1 * 7 * (cost_hs_per_traced + cost_hs_per_quarantined) +
        # testing costs, 11.31 tests per hospitalised case
        (nonicu_i + icu_i) * 11.31 * cost_hs_per_test +
        # testing of 10% of non-hospitalised cases
        (cases - (nonicu_i + icu_i)) * 
        0.1 * 11.31 * cost_hs_per_test +
        # daily cost of treatment on general ward
        (nonicu_p) * cost_hs_day_treat_general +
        # daily cost of critical care
        icu_p * cost_hs_day_treat_critical +
        # one-off cost of treating 10% of non-hospitalised cases at home
        (cases - (nonicu_i+icu_i)) * 
        0.1 * cost_hs_treat_home +
        # cost of death to the health system
        death_o * cost_hs_per_death +
        # household: cost of death
        death_o * (cost_hh_death_funeral + cost_hh_death_income) +
        # household: medical + non-med costs while in general ward
        (nonicu_p) *
        (cost_hh_treat_general_med_per_day + cost_hh_treat_general_non_med_per_day) +
        # household: medical + non-med costs while in icu
        icu_p * (cost_hh_treat_critical_med_per_day + cost_hh_treat_critical_non_med_per_day) +
        # household: medical + non-med costs for non-hospitalised cases
        # TODO: update with prevalent cases or alternative assumption about
        # treatment duration at home
        (cases - (nonicu_i+icu_i)) * 
        7 * (cost_hh_treat_home_med_per_day + cost_hh_treat_home_non_med_per_day) +
        # household: individual and caregiver lost income
        cases * (cost_hh_individual_income_per_case + cost_hh_caregiver_income_per_case)
    )][,
       costs := with(econ_pars, costs + doses*cost_vac_dose)
    ][,
       costs := with(econ_pars, (1/(1 + disc.costs)^(anni_year - 1)) * costs)
    ]
    
    epi.dt[
        dalys.dt,
        dalys := with(econ_pars,
            death_o * dalys_death + # dalys per death
            cases * dalys_case + # dalys per case
            # dalys per hospitalised case in general ward
            nonicu_i * dalys_hospital +
            # dalys per icu admissions that survive
            (icu_i - death_o) * dalys_icu
        ),
        on=.(age)
    ][, #' TODO revisit approach?
      dalys := with(econ_pars, (1/(1 + disc.dalys)^(anni_year-1)) * dalys)
    ]
    
    agg.dt <- epi.dt[,
        .(costs = sum(costs), dalys = sum(dalys)),
        by=.(id, sampleId, anni_year)
    ]
    
    agg.dt[
        order(anni_year),
        c("ccosts", "cdalys") := .(cumsum(costs), cumsum(dalys)),
        by=.(id, sampleId)
    ]
    
    return(agg.dt)
}

ret.dt <- econscns.dt[,{
    es <- as.list(.SD)
    econ_pars <- c(
        as.list(othercosts[perspective == es$perspective]),
        as.list(vac_cost.dt[vac_price == es$vac_price, .(cost_vac_dose)]),
        es[c("disc.costs","disc.dalys")]
    )
    econ_digestor(copy(wide.dt), dalys.dt[daly_scenario == es$daly_scenario & disc_rate == es$disc.dalys], econ_pars)
},by=econ_id]

if (!isbaseline) {
    #' join the relevant baseline
    basescn <- 0
    qtile <- function(
        v, ps = c(lo95=0.025, lo50=0.25, md=0.5, hi50=0.75, hi95=0.975),
        withMean = c("mn", NA),
        fmt = "%s",
        na.rm = TRUE
    ) {
        qs <- quantile(v, probs = ps, na.rm = na.rm)
        names(qs) <- sprintf(fmt, names(ps))
        if (!is.na(withMean[1])) {
            mn <- mean(v)
            names(mn) <- sprintf(fmt, withMean[1])
            qs <- c(qs, mn)
        }
        as.list(qs)
    }
    
    # if we aren't looking at one of the base scenarios
    both.dt <- { 
        ref <- readRDS(gsub("(\\w+)\\.rds", "\\1_baseline.rds", tail(.args, 1)))[id == basescn][, .(econ_id, sampleId, anni_year, costs, dalys, ccosts, cdalys)]
        inc.dt <- copy(ret.dt)[ref, on=.(econ_id, sampleId, anni_year), .(
            econ_id, id, sampleId, anni_year,
            costs = costs - i.costs, # -cost == savings
            dalys = i.dalys - dalys, # +dalys == dalys gained
            ccosts = ccosts - i.ccosts,
            cdalys = i.cdalys - cdalys 
        )][, view := "incremental" ][, icer := ccosts / cdalys ]
        
        both.dt <- melt(rbind(
            ret.dt[, view := "raw" ], inc.dt, fill = TRUE
        ), id.vars = c("econ_id","id","sampleId", "anni_year", "view"))
    }
    
    res.dt <- dcast(melt(
        both.dt[,qtile(value),by=.(econ_id, id, anni_year, view, variable)],
        id.vars = c("econ_id","id", "anni_year", "view", "variable"), variable.name = "qtile"
    ), qtile + id + econ_id + anni_year + view ~ variable)
    saveRDS(res.dt, tail(.args, 1))
} else {
    saveRDS(ret.dt, tail(.args, 1))
}
