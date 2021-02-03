

suppressPackageStartupMessages({
  require(data.table)
  require(qs)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenarios/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/covidm_fit_yu.qs",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/relaxation/%s.rds",
  .debug[2],
  "../covidm",
  "%s/outputs/projections/%s.qs"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

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

Rts <- readRDS(.args[4])[era %in% c("pre", "post"), .(date, med, era)]

run_options <- melt(
  Rts[era == "pre"],
  measure.vars = c("med"),
  value.name = "r0"
)[, model_seed := 1234L ]

refR0 <- cm_ngm(params)$R0

params_back <- params

allbind <- data.table(
  run=integer(0),
  group = factor(),
  value = integer(),
  AR = numeric(),
  r_id = integer()
)

Rtrelax <- readRDS(.args[7])

startrelax <- Rtrelax[order(date)][1, date]

iv_data <- scenario[scen_id != 1 & !is.na(self_iso)][order(trigger_type)]

genschedule <- function(curschedule, iv) {
  cons <- with(iv[order(start_day)], mapply(
    function(h,w,s,o) list(1-c(h,w,s,o)),
    h=home, w=work, s=school, o=other
  ))
  si <- lapply(iv[order(start_day)]$self_iso, function(si) {
    rep(1-si, 16)
  })
  tms <- as.Date(params$date0) + iv[order(start_day)]$start_day
  curschedule[[length(curschedule)+1]] <- list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    values = cons,
    times = tms
  )
  curschedule[[length(curschedule)+1]] <- list(
    parameter = "fIs",
    pops = numeric(),
    mode = "multiply",
    values = si,
    times = tms
  )
  curschedule
}

for(i in 1:nrow(run_options)){
  params <- params_back
  # only run until recalc point
  params$time1 <- startrelax
  
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
  
  params$schedule <- genschedule(params$schedule, iv_data)
  
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

#' TODO ris <- 1:5?
ris <- 1
baseRs <- run_options[, r0]
ufs <- mapply(
  function(br, sfrac) br/refR0 * sfrac,
  br = baseRs,
  sfrac = lapply(ris, function(ri) 1 - allbind[r_id == ri, AR]),
  SIMPLIFY = FALSE
)

#' these are the baseline reductions  
vs <- scenario[!is.na(self_iso) & scen_id != 1][,c(home, work, school, other, self_iso)]

cm_ngm(params_back, ufs[[1]], vs[-5], vs[5])$R0

modR <- melt(Rts[era == "modification"], id.vars = "era", measure.vars = 2:6)
tarRs <- modR[, value]


#' these are the age specific reductions in susceptibility,
#' due to baseline adjustment to R0 & associated AR
ufs <- mapply(
  function(br, sfrac) br/refR0 * sfrac,
  br = baseRs,
  sfrac = lapply(1:5, function(ri) 1 - allbind[r_id == ri, AR]),
  SIMPLIFY = FALSE
)


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
  ) - rev(tarRs))^2)
}

#' determine fitting space optimal del    
odel <- optimize(fn, c(-20, 20))$minimum

#' update scenarios
scenario[
  is.na(self_iso) & scen_id != 1 & !is.infinite(end_day),
  c("home","work","school","other", "self_iso") := as.list(mvs(odel))
]

allbind <- data.table(
  run = integer(), t = integer(),
  group = factor(), compartment = factor(),
  value = numeric(), r_id = integer()
)

iv_data <- scenario[scen_id != 1 & !is.na(self_iso)][order(trigger_type)]
for(i in 1:nrow(run_options)){
  
  params <- params_back
  # only run until recalc point
  params$time1 <- iv_data[, max(end_day, na.rm = TRUE)+60]
  
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
  )$dynamics
  
  allbind <- rbind(sim[
    compartment %in% c("cases", "subclinical", "R"),
    .(value),
    keyby=.(run, t, group, compartment)
  ][, r_id := i ], allbind)
  
}

#' @examples 
#' require(ggplot2)
#' ggplot(allbind[compartment == "cases"]) + aes(t, value) + facet_grid(group ~ ., scales = "free_y") + geom_line()

qsave(allbind, tail(.args, 1))
