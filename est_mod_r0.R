suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenarios/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/covidm_fit_yu.qs",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  .debug[2],
  "../covidm",
  "%s/mod_scenarios/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

# identify country / scenario
scenario <- readRDS(.args[1])
tariso <- tail(.args, 3)[1]

load("NGM.rda")

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

suppressPackageStartupMessages({
  source(file.path(cm_path,"R","covidm.R"))
})

yu_fits <- qread(.args[3])[order(ll)]
yu_fits[, eqs := (1:.N)/.N ]
#' using the median yu fits
medyu <- yu_fits[which.max(eqs > 0.5)]
yref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("y_",colnames(medyu))]))
uref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("u_",colnames(medyu))]))
ys <- rep(yref[1, ], each = 2)
us <- rep(uref[1, ], each = 2)

params <- readRDS(.args[2])
  
intros <- readRDS(.args[5])[iso3 == tariso][, intro.day := as.integer(date - date[1]) ][, Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))]
urbfrac <- readRDS(.args[6])[iso3 == tariso, value / 100]
  
params$pop[[1]]$seed_times <- intros
params$pop[[1]]$size <- round(params$pop[[1]]$size*urbfrac)
# params$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))

params$pop <- lapply(
  params$pop,
  function(x){
    x$y <- ys
    x$u <- us
    return(x)
  }
)

Rts <- readRDS(.args[4])[era != "transition"]
  
run_options <- melt(
  Rts[era == "pre"],
  measure.vars = c("lo.lo","lo","med","hi","hi.hi"),
  value.name = "r0"
)[, model_seed := 1234L ]
  
refR0 <- cm_ngm(params)$R0

params_back <- params
  
# scenario[, waning_dur := NA_integer_ ]
# scenario[!(scen_type == "unmitigated"), waning_dur := 90 ]
  
allbind <- data.table(
  run=integer(0),
  group = factor(),
  value = integer(),
  AR = numeric(),
  r_id = integer()
)
  
iv_data <- scenario[scen_id != 1 & !is.na(self_iso)][order(trigger_type)]
for(i in 1:nrow(run_options)){
  params <- params_back
  # only run until recalc point
  params$time1 <- iv_data[, max(end_day, na.rm = TRUE)]
  
  #' adjust r0 to that in current sample
  #' TODO use sample y / u?
  target_R0 <- run_options[i, r0]
  uf <- target_R0 / refR0
  params$pop <- lapply(
    params$pop,
    function(x){
      x$u <- x$u * uf
      return(x)
    }
  )
      
  # generic interventions
  for (j in 1:nrow(iv_data[population == -1])) {
    #pars <- as.list(iv_data[population == -1][j])
    cons <- with(as.list(iv_data[population == -1][j]), {
      list(1-c(home, work, school, other), c(1,1,1,1))
    })
    tms <- with(as.list(iv_data[population == -1][j]), {
      as.Date(params$date0)+c(start_day, ifelse(is.finite(end_day),end_day,params$time1))
    })
    si <- with(as.list(iv_data[population == -1][j]), {
      list(rep(1-self_iso, 16), rep(1, 16))
    })
    params$schedule[[length(params$schedule)+1]] <- list(
      parameter = "contact",
      pops = numeric(),
      mode = "multiply",
      values = cons,
      times = tms
    )
    params$schedule[[length(params$schedule)+1]] <- list(
      parameter = "fIs",
      pops = numeric(),
      mode = "multiply",
      values = si,
      times = tms
    )
    
    
  }

  #' run the model
  sim <- cm_simulate(
    params, 1,
    model_seed = run_options[i, model_seed]
  )$dynamics[compartment == "R"][
    order(t), .(value = value[.N]), keyby=.(run, group)
  ]
  
  sim[order(run, group), AR := value / params$pop[[1]]$size ]
  
  sim[, r_id := i ]
  
  allbind <- rbind(allbind, sim)
  
}

#' now that we have susceptible fraction after initial interventions,
#' need to fit reductions that best adjust R

modR <- melt(Rts[era == "modification"], id.vars = "era", measure.vars = 2:6)
tarRs <- modR[, value]
baseRs <- run_options[, r0]

#' these are the age specific reductions in susceptibility,
#' due to baseline adjustment to R0 & associated AR
ufs <- mapply(
  function(br, sfrac) br/refR0 * sfrac,
  br = baseRs,
  sfrac = lapply(1:5, function(ri) 1 - allbind[r_id == ri, AR]),
  SIMPLIFY = FALSE
)

#' these are the baseline reductions  
vs <- scenario[!is.na(self_iso) & scen_id != 1][,c(home, work, school, other, self_iso)]

#' fitting space to model space
#' del is a single factor for changes to all parameters
#' del can range from -inf to +inf
#' del < 0 => decreasing the reduction toward 0%
#' del > 0 => increasing the reduction toward 100%
#' del == 0 => no change
mvs <- function(del) {
  fac <- ifelse(del >= 0, 1, 0) - sign(del)*exp(-abs(del))
  return(if (del > 0) {
    vs*(1-fac) + fac
  } else {
    vs*fac
  })
}

fn <- function(del) {
  rs <- mvs(del)
  sum((sapply(ufs, function(uf)
    cm_ngm(
      params_back,
      uf,
      contact_reductions = rs[-5],
      fIs_reductions = rs[5]
    )$R0
  ) - tarRs)^2)
}

#' determine fitting space optimal del    
odel <- optimize(fn, c(-20, 20))$minimum

#' update scenarios
scenario[
  is.na(self_iso) & scen_id != 1 & !is.infinite(end_day),
  c("home","work","school","other", "self_iso") := as.list(mvs(odel))
]

allbind[, era := "pre" ]
  
iv_data <- scenario[scen_id != 1 & !is.na(self_iso)][order(trigger_type)]
for(i in 1:nrow(run_options)){
  
  params <- params_back
  # only run until recalc point
  params$time1 <- iv_data[, max(end_day, na.rm = TRUE)]
      
  #adjust r0 to that in current sample
  target_R0 <- run_options[i, r0]
  uf <- target_R0 / refR0
  params$pop <- lapply(
    params$pop,
    function(x){
      x$u <- x$u * uf
      return(x)
    }
  )
      
  cons <- with(iv_data[order(start_day)], mapply(
    function(h,w,s,o) list(1-c(h,w,s,o)),
    h=home, w=work, s=school, o=other
  ))
  si <- lapply(iv_data[order(start_day)]$self_iso, function(si) {
    rep(1-si, 16)
  })
  tms <- as.Date(params$date0) + iv_data[order(start_day)]$start_day
  params$schedule[[length(params$schedule)+1]] <- list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    values = cons,
    times = tms
  )
  params$schedule[[length(params$schedule)+1]] <- list(
    parameter = "fIs",
    pops = numeric(),
    mode = "multiply",
    values = si,
    times = tms
  )
  #run the model
  sim <- cm_simulate(
    params, 1,
    model_seed = run_options[i, model_seed]
  )$dynamics[compartment == "R"][
    order(t), .(value = value[.N]), keyby=.(run, group)
  ]
      
  sim[order(run, group), AR := value / params$pop[[1]]$size ]
  
  sim[, r_id := i ][, era := "variant" ]
      
  allbind <- rbind(allbind, sim)

}

#' now we have the fully depleted population
#' + observed Reff
#' 
#' option 1: full cross protective - what multiple of uf would be required to observe
#' this Reff, given depletion & modification era reductions?
#' that multiple => increase of R0

ufs <- mapply(
  function(br, sfrac) br/refR0 * sfrac,
  br = baseRs,
  sfrac = lapply(1:5, function(ri) 1 - allbind[r_id == ri, AR]),
  SIMPLIFY = FALSE
)

iv_data <- scenario[scen_id != 1 & !is.na(self_iso)][order(start_day)][.N]

varR <- melt(Rts[era == "variant"], id.vars = "era", measure.vars = 2:6)$value

Rmultipliers_with_depletion <- mapply(function(uf, R) with(
  iv_data,
  R/cm_ngm(params_back, R0_multiplier = uf, contact_reductions = c(home,work,school,other), fIs_reductions = self_iso)$R0
), uf = ufs, R = varR)

ufs_non_depl <- mapply(
  function(br, sfrac) br/refR0 * sfrac,
  br = baseRs,
  sfrac = 1,
  SIMPLIFY = FALSE
)

Rmultipliers_non_depletion <- mapply(function(uf, R) with(
  iv_data,
  R/cm_ngm(params_back, R0_multiplier = uf, contact_reductions = c(home,work,school,other), fIs_reductions = self_iso)$R0
), uf = ufs_non_depl, R = varR)

res <- data.table(
  model = c(rep("cross-protected", 5), rep("susceptible", 5)),
  Rfactor = c(Rmultipliers_with_depletion, Rmultipliers_non_depletion)
)

saveRDS(res, tail(.args, 1))