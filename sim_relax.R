
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenarios/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/yuqs/%s.rds",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "../covidm",
  "%s/outputs/projections/%s.rds"
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
  source(file.path(cm_path, "R", "covidm.R"))
})

params <- readRDS(.args[2])
yuref <- readRDS(.args[3])[order(eqs)]
qs.inds <- with(yuref, c(which.max(eqs >= 0.25),which.max(eqs >= 0.5),which.max(eqs >= 0.75)))
yuuse <- yuref[qs.inds][, variable := c("lo","med","hi") ][,-c("trial","chain","lp","ll","mult","size")]

intros.dt <- readRDS(.args[5])[iso3 == tariso]
day0 <- as.Date(intros.dt[, min(date)])
intros <- intros.dt[, intro.day := as.integer(date - date[1]) ][, Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))]

urbfrac <- readRDS(.args[6])[iso3 == tariso, value / 100]

params$date0 <- day0
params$pop[[1]]$seed_times <- intros
params$pop[[1]]$size <- round(params$pop[[1]]$size*urbfrac)
params$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))

Rts <- readRDS(.args[4])[era != "transition"]

run_options <- melt(
  Rts[era == "pre"],
  id.vars = "era",
  measure.vars = c("lo", "med", "hi"),
  value.name = "r0"
)[, model_seed := 1234L ][yuuse, on=.(variable) ][,
  umul := r0 / baseR
]

timings <- readRDS(.args[7])
#' TODO update timings?
startrelax <- as.integer(as.Date("2020-05-01") - day0)
endsim <- as.integer(timings[era == "relaxation", end] - day0) + 30

startpost <- as.integer(timings[era == "transition", start[1]] - day0)

params$time1 <- endsim

# scenario[, waning_dur := NA_integer_ ]
# scenario[!(scen_type == "unmitigated"), waning_dur := 90 ]

iv_data <- scenario[scen_id != 1 & !is.na(self_iso)]
tms <- as.Date(params$date0) + startpost:iv_data[order(start_day)]$start_day[1]
relaxtms <- as.Date(params$date0) + startrelax:endsim
redfrac <- 1:length(tms)/length(tms)

params_back <- params

sim <- data.table(); for(i in 1:nrow(run_options)){

  params <- params_back

  #' only one intervention
  cons <- with(iv_data[i], mapply(
    function(h,w,s,o) { lapply(redfrac, function(rf) c(1-c(h,w,s,o)))  },
    h=home, w=work, s=school, o=other, SIMPLIFY = FALSE
  ))[[1]]
  
  si <- lapply(iv_data[i]$self_iso, function(si) {
    lapply(redfrac, function(rf) rep((1-si), 16))
  })[[1]]
  
  # relaxrate <- c(0.001,0.005,0.006)[i]
  # relaxfact <- pmax(1-seq(0,by=relaxrate,length.out = length(relaxcons)),0)
  
  tier2 <- as.Date("2020-08-15")+7
  relaxrate <- c(0.03, 0.035, 0.075)[i]
  relaxfact <- 1-(1+exp(-relaxrate*as.numeric(relaxtms-tier2)))^-1
  relaxfact <- 0.5*(relaxfact-0.5)/(relaxfact[1]-0.5)+0.5
  #' @example 
  #' plot(1:length(relaxfact), relaxfact)
  
  relaxcons <- lapply(
    relaxfact, function(rf) 1-(1-cons[[length(cons)]])*rf
  )
  
  relaxsi <- lapply(
    relaxfact, function(rf) 1-(1-si[[length(si)]])*rf
  )
  
  params$schedule[[length(params$schedule)+1]] <- list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    values = c(cons, relaxcons),
    times = c(tms, relaxtms)
  )
  params$schedule[[length(params$schedule)+1]] <- list(
    parameter = "fIs",
    pops = numeric(),
    mode = "multiply",
    values = c(si, relaxsi),
    times = c(tms, relaxtms)
  )
  
  #' adjust r0 to that in current sample
  uf <- run_options[i, umul*rep(as.numeric(.SD), each = 2), .SDcols = grep("^u_", names(run_options))]
  ys <- run_options[i, rep(as.numeric(.SD), each = 2), .SDcols = grep("^y_", names(run_options))]
  
  params$pop <- lapply(
    params$pop,
    function(x){
      x$u <- x$u * uf
      x$y <- ys
      return(x)
    }
  )
  
  #run the model
  sim <- rbind(
    cm_simulate(
      params, 1,
      model_seed = run_options[i, model_seed]
    )$dynamics[compartment %in% c("R","cases","death_o")][, r_id := i],
    sim, fill = TRUE
  )
  
}

res <- sim[, .(date = t + day0, group, compartment, value, q = c("lo","md","hi")[r_id])]

#' @examples 
#' ggplot(sim[compartment == "cases", .(value = sum(value)), by=.(r_id, date = t+params_back$date0)]) +
#'   aes(date, value, color = factor(r_id), group = r_id) +
#'   geom_line() +
#'   geom_line(
#'     aes(date, cases),
#'     data = readRDS("~/Dropbox/SA2UK/inputs/epi_data.rds")[iso3 == "ZAF"],
#'     color = "black", inherit.aes = FALSE
#'   ) +
#'   scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#'   scale_y_log10() +
#'   theme_minimal()

saveRDS(res, tail(.args, 1))

