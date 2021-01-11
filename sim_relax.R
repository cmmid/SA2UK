
suppressPackageStartupMessages({
  require(data.table)
  require(optimization)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/params/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "../covidm",
  "%s/outputs/projections/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 3)[1]

fits <- readRDS(.args[1])
params <- readRDS(.args[2])
intros.dt <- readRDS(.args[3])[iso3 == tariso]
urbfrac <- readRDS(.args[4])[iso3 == tariso, value / 100]
timings <- readRDS(.args[5])

day0 <- as.Date(intros.dt[, min(date)])
intros <- intros.dt[,
  intro.day := as.integer(date - date[1])
][, Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))]

params$date0 <- day0
params$pop[[1]]$seed_times <- intros
params$pop[[1]]$size <- round(params$pop[[1]]$size*urbfrac)
params$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))

startrelax <- as.integer(timings[era == "relaxation", start] - day0)
endrelax <- as.integer(timings[era == "relaxation", end] - day0)
endsim <- as.integer(timings[era == "variant", start] - day0)

startpost <- as.integer(timings[era == "transition", start[1]] - day0)

params$time1 <- endsim

tms <- day0 + startpost
relaxtms <- day0 + startrelax:endrelax

tier2 <- as.Date("2020-08-15")

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

scheduler <- function(large, small, symp, k, shft) {
  cons <- list(1-c(0, small, large, small))
  si <- list(rep(1-symp, 16))
  
  relaxfact <- 1-(1+exp(-k*as.numeric(relaxtms-tier2-shft)))^-1
  relaxfact <- (1-relaxfact[1]) + relaxfact
  relaxcons <- lapply(relaxfact, function(rf) 1-c(0, small, large, small)*rf)
  relaxsi <- lapply(relaxfact, function(rf) 1-rep(symp, 16)*rf)
  
  list(
    list(
      parameter = "contact",
      pops = numeric(),
      mode = "multiply",
      values = c(cons, relaxcons),
      times = c(tms, relaxtms)
    ),
    list(
      parameter = "fIs",
      pops = numeric(),
      mode = "multiply",
      values = c(si, relaxsi),
      times = c(tms, relaxtms)
    )
  )
}

sims <- rbindlist(lapply(1:nrow(fits), function(i) with(as.list(fits[i,.(large, small, sympt, k, shft, asc)]), {
  us <- rep(fits[i, as.numeric(.SD)*umod, .SDcols = grep("^u_",names(fits))], each = 2)
  ys <- rep(fits[i, as.numeric(.SD), .SDcols = grep("^y_",names(fits))], each = 2)
  testpop <- params; testpop$pop[[1]]$y <- ys
  testpop$pop[[1]]$u <- testpop$pop[[1]]$u*us
  testpop$schedule <- scheduler(large, small, sympt, k, shft)
  res <- cm_simulate(
    testpop, 1,
    model_seed = 42L
  )$dynamics[compartment %in% c("cases","death_o","R")][, sample := i ]
})))

res <- sims[,
  .(
    sample, date = t + day0,
    group, compartment,
    value
  )
]

#' @examples 
#' ggplot(res[compartment == "cases"][fits[, .(sample, asc)], on=.(sample)][, .(asc.value = sum(value)*asc, value = sum(value)), by=.(sample, date)]) +
#'   aes(date, asc.value, group = sample) +
#'   geom_line(alpha = 0.1) +
#'   geom_line(aes(y=value), alpha = 0.1, color = "red") +
#'   geom_line(
#'     aes(date, cases),
#'     data = readRDS("~/Dropbox/SA2UK/inputs/epi_data.rds")[iso3 == "ZAF"],
#'     color = "black", inherit.aes = FALSE
#'   ) +
#'   scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#'   scale_y_log10() +
#'   theme_minimal()

saveRDS(res, tail(.args, 1))

