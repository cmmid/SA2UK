

suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/params/%s_0001.rds",
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
intros.dt <- readRDS(.args[3])[iso3 == tariso]
urbfrac <- readRDS(.args[4])[iso3 == tariso, value / 100]
timings <- readRDS(.args[5])

day0 <- as.Date(intros.dt[, min(date)])
intros <- intros.dt[,
                    intro.day := as.integer(date - date[1])
][, .(t=Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))), by=sample ]

popsetup <- function(basep, urbanfraction, day0) {
  basep$date0 <- day0
  basep$pop[[1]]$size <- round(basep$pop[[1]]$size*urbanfraction)
  basep$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))
  basep
}

base <- popsetup(readRDS(.args[2]), urbfrac, day0)

startrelax <- as.integer(timings[era == "relaxation", start] - day0)
endrelax <- as.integer(timings[era == "relaxation", end] - day0)
endsim <- as.integer(timings[era == "pre" & period == 3, start] - day0)

startpost <- as.integer(timings[era == "transition", start[1]] - day0)

base$time1 <- endsim

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

sims <- fits[,{
  us <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^u_",names(.SD))], each = 2)*umod
  ys <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^y_",names(.SD))], each = 2)
  testpop <- base;
  testpop$pop[[1]]$y <- ys
  testpop$pop[[1]]$u <- testpop$pop[[1]]$u*us
  sid <- sample
  testpop$pop[[1]]$seed_times <- intros[sample == sid, t]
  testpop$schedule <- scheduler(large, small, sympt, k, shft)
  res <- cm_simulate(
    testpop, 1, model_seed = 42L
  )$dynamics[compartment %in% c("cases","death_o","R")]
}, by=sample]

res <- sims[,
            .(
              sample, date = t + day0,
              group, compartment,
              value
            )
]

#' @examples 
#' comparison <- res[compartment == "cases" & between(date, tarwindow[1], tarwindow[2]), .(value = sum(value)), by=.(sample, date)][sample == 1]
#' ggplot(res[compartment == "cases"][fits[, .(sample, asc)], on=.(sample)][, .(asc.value = sum(value)*asc, value = sum(value)), by=.(sample, date)]) +
#'   aes(date, asc.value, group = sample) +
#'   geom_line(alpha = 0.1) +
#'   geom_line(aes(y=value), alpha = 0.1, color = "red") +
#'   geom_line(
#'     aes(date, cases),
#'     data = readRDS("~/Dropbox/Covid_LMIC/All_Africa_paper/inputs/epi_data.rds")[iso3 == .debug[2]],
#'     color = "black", inherit.aes = FALSE
#'   ) +
#'   annotate("rect", xmin=as.Date("2020-09-01"), xmax =as.Date("2020-10-01"), ymin = 0.01, ymax = Inf, alpha = 0.2, fill = "dodgerblue") +
#'   scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#'   scale_y_log10() +
#'   theme_minimal()

saveRDS(res, tail(.args, 1))