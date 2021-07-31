# update to econ.R with total costs disaggregated into coi and vax costs
suppressPackageStartupMessages({
    require(data.table)
})

.debug <- c("~/Dropbox/Covid-WHO-vax", "baseline")
#.debug <- c("~/Dropbox/Covid-WHO-vax", "03778")
.args <- if (interactive()) sprintf(c(
    "covid_other_costs.csv",
    "covid_vac_costs_per_dose.csv",
    "daly_scenarios.csv",
    "%s/outputs/config.rds",
    ifelse(.debug[2]=="baseline","%s/outputs/sim","%s/outputs/sim/%s.rds"), # for the baseline, will combine several
    "%s/outputs/econ/%s.rds"
),.debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)


#' INGEST EPI DATA

isbaseline <- grepl("baseline.rds$", tail(.args, 1))

scn <- if (isbaseline) {
    readRDS(.args[4])[strategy == "none"]
} else {
    readRDS(.args[4])[id == as.integer(gsub("(\\d+)\\.rds$","\\1",basename(tail(.args, 1))))]
}

epi.fls <- if (isbaseline) {
    scns <- scn[, sprintf("(%s)", paste(sprintf("%05i",id), collapse = "|")) ]
    list.files(tail(.args, 2)[1], scns, full.names = TRUE)
} else tail(.args, 2)[1]

dt <- rbindlist(lapply(epi.fls, function(fn) readRDS(fn)[order(anni_year),.(anni_year = anni_year[-1], value = diff(value)), by=.(sampleId, age, outcome)][,fn := fn ]))
dt[, id := as.integer(gsub("(\\d+)\\.rds$","\\1", basename(fn))) ]
dt$fn <- NULL

wide.dt <- dcast(dt, id + age + sampleId + anni_year ~ outcome, value.var = "value")[scn[, .(id, vax_delay, strategy_str, doses_per_day)], on=.(id)]

rm(dt)

#' SETUP ECON DATA

#' vector of doses

increasing <- !is.na(scn$increasing[1]) & scn$increasing[1]
doses_yr_1 <- (if(increasing) { sum(c(1,4,6,8))/4 } else { 1 })*365  # year 1
doses_routine <- (if(increasing) { 8 } else { 1 })*365 
doses_per_anniversary <- c(doses_yr_1, rep(doses_routine, scn$horizon[1]))

if (is.na(scn$strategy_str[1])) {
  doses_per_anniversary <- doses_per_anniversary*0
} else if (scn$strategy_str[1] == 0) {
  # no change
} else {
  yrs <- floor(scn$strategy_str[1]/365)
  partial <- scn$strategy_str[1]/365 - floor(scn$strategy_str[1]/365)
  mul <- rep(0, length(doses_per_anniversary))
  mul[1:yrs] <- 1
  if (yrs < length(mul)) mul[yrs+1] <- partial
  doses_per_anniversary <- doses_per_anniversary * mul
}

# by perspective
othercosts <- dcast(fread(.args[1]), perspective ~ name, value.var = "cost")

vac_cost.dt <- fread(.args[2])[scenario == "campaign"]
dalys.dt <- fread(.args[3])
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
    
    epi.dt[, coi_costs := with(econ_pars,
        # one-off / daily health system response costs
        365 * (cost_hs_day_erm + cost_hs_day_comms) + 
        cost_hs_one_erm +
        # 10 % of symptomatic cases tested, 7 contacts per tested case
        # TODO: how do these assumptions scale with prevalence
        cases * 0.1 * 7 * (cost_hs_per_traced + cost_hs_per_quarantined) +
        # testing costs, 11.31 tests per hospitalised case
        (non_icu_severe_i + non_icu_critical_i) * 11.31 * cost_hs_per_test +
        # testing of 10% of non-hospitalised cases
        (cases - (non_icu_severe_i + non_icu_critical_i)) * 
        0.1 * 11.31 * cost_hs_per_test +
        # daily cost of treatment on general ward
        (non_icu_severe_p + non_icu_critical_p) * cost_hs_day_treat_general +
        # daily cost of critical care
        icu_critical_p * cost_hs_day_treat_critical +
        # one-off cost of treating 10% of non-hospitalised cases at home
        (cases - (non_icu_severe_i + non_icu_critical_i)) * 
        0.1 * cost_hs_treat_home +
        # cost of death to the health system
        death_o * cost_hs_per_death +
        # household: cost of death
        death_o * (cost_hh_death_funeral + cost_hh_death_income) +
        # household: medical + non-med costs while in general ward
        (non_icu_severe_p + non_icu_critical_p) *
        (cost_hh_treat_general_med_per_day + cost_hh_treat_general_non_med_per_day) +
        # household: medical + non-med costs while in icu
        icu_critical_p * (cost_hh_treat_critical_med_per_day + cost_hh_treat_critical_non_med_per_day) +
        # household: medical + non-med costs for non-hospitalised cases
        # TODO: update with prevalent cases or alternative assumption about
        # treatment duration at home
        (cases - (non_icu_severe_i + non_icu_critical_i)) * 
        7 * (cost_hh_treat_home_med_per_day + cost_hh_treat_home_non_med_per_day) +
        # household: individual and caregiver lost income
        cases * (cost_hh_individual_income_per_case + cost_hh_caregiver_income_per_case)
    )][,
       vax_costs := with(econ_pars, fifelse(
           is.na(vax_delay),
           0,
           doses_per_day*doses_per_anniversary[anni_year]/age_cats *
               cost_vac_dose * fifelse(vax_delay == 0, 1, 2)
       ))
    ][,
      costs := coi_costs + vax_costs
    ][,.(
      vax_costs = with(econ_pars, (1/(1 + disc.costs)^(anni_year - 1)) * vax_costs),
      coi_costs = with(econ_pars, (1/(1 + disc.costs)^(anni_year - 1)) * coi_costs),
      costs = with(econ_pars, (1/(1 + disc.costs)^(anni_year - 1)) * costs)
    )
    ]
    
    epi.dt[
        dalys.dt,
        dalys := with(econ_pars,
            death_o * dalys_death + # dalys per death
            cases * dalys_case + # dalys per case
            # dalys per hospitalised case in general ward
            non_icu_severe_i * dalys_hospital +
            # dalys per icu admissions that survive
            (icu_critical_i - death_o) * dalys_icu
        ),
        on=.(age)
    ][, #' TODO revisit approach?
      dalys := with(econ_pars, (1/(1 + disc.dalys)^(anni_year-1)) * dalys)
    ]
    
    agg.dt <- epi.dt[,
        .(costs = sum(costs), coi_costs = sum(coi_costs), vax_costs = sum(vax_costs), dalys = sum(dalys)),
        by=.(id, sampleId, anni_year)
    ]
    
    agg.dt[
        order(anni_year),
        c("ccosts","ccoi_costs","cvax_costs","cdalys") := .(cumsum(costs), cumsum(coi_costs), cumsum(vax_costs), cumsum(dalys)),
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
    this <- as.list(scn[,.(nat_imm_dur_days, start_timing)])
    basescn <- readRDS(.args[4])[strategy == "none" & nat_imm_dur_days == this$nat_imm_dur_days & start_timing == this$start_timing, id]
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
    both.dt <- { if (basescn != scn$id) {
        ref <- readRDS(gsub("\\d+\\.rds", "baseline.rds", tail(.args, 1)))[id == basescn][, .(econ_id, sampleId, anni_year, costs, coi_costs, vax_costs, dalys, ccosts, ccoi_costs, cvax_costs, cdalys)]
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
    } else melt(ret.dt[, view := "raw" ][, icer := NA_real_ ], id.vars = c("econ_id","id","sampleId", "anni_year", "view")) } 
    
    res.dt <- dcast(melt(
        both.dt[,qtile(value),by=.(econ_id, id, anni_year, view, variable)],
        id.vars = c("econ_id","id", "anni_year", "view", "variable"), variable.name = "qtile"
    ), qtile + id + econ_id + anni_year + view ~ variable)
    saveRDS(res.dt, tail(.args, 1))
} else {
    saveRDS(ret.dt, tail(.args, 1))
}
