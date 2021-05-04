suppressPackageStartupMessages({
  require(data.table)
  require(doParallel)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper","PAK")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/pops/%s.rds",
  "%s/outputs/params/%s_consolidated.rds",
  "%s/inputs/mobility.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/outputs/adj_data.rds",
  .debug[2], # PAK
  "../covidm",
  "%s/outputs/params/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

fits.dt <- readRDS(.args[2])

tariso <- tail(.args, 3)[1]

mob <- readRDS(.args[3])[iso3 == tariso]
timings <- readRDS(.args[4])
tarwindow <- timings[era == "relaxation", start]
tarwindow[2] <- timings[era == "pre" & period == 3, start]

case.dt <- readRDS(.args[6])[
  (iso3 == tariso),
  .(date, croll = frollmean(cases, 7, align = "center"))
]

intros.dt <- readRDS(.args[5])[iso3 == tariso]

day0 <- as.Date(intros.dt[, min(date)])

contact_schedule <- with(mob[date >= day0], mapply(
  function(work, other, school, home = 1) c(home, work, other, school),
  work = workr, other = otherr, school = school_multiplier,
  SIMPLIFY = FALSE
))

intros <- intros.dt[,
  intro.day := as.integer(date - date[1])
][, .(t=Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))), by=.(sid=sample) ]

popsetup <- function(basep, day0) {
  basep$date0 <- day0
  basep$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))
  basep
}

params <- popsetup(readRDS(.args[1]), day0)
params$schedule <- list(
  list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    values = contact_schedule,
    times = day0 + 0:(length(contact_schedule)-1)
  )
)

#' TODO: fix warning here; providing correct value, however
startpost <- as.integer(timings[era == "transition" & period == 1, start[1]] - day0)
startrelax <- as.integer(timings[period == 2, start[1]] - day0)

tart <- as.numeric(tarwindow - day0)
case.slc <- case.dt[between(date, tarwindow[1], tarwindow[2]), round(croll)]
params$time1 <- tarwindow[2]

#' TODO extract fIs model into it's own thing
fIs_amp <- function(
  case0, k, # fit elements
  reff, # sampling element: value of function at css[1]
  css = case.slc, # data element
  baseline = (reff*(1+exp(-k*(css[1]-case0))) - 1)/exp(-k*(css[1]-case0)) # entailed remaining coefficient
) (1-baseline)/(1+exp(-k*(css-case0))) + baseline

#   case0 = 10^mean(range(log10(case.slc)))

add_fIs <- function(
  case0, k, # fit elements
  fIs_reduction_at_post, # value at post-intervention
  model_t = startpost, d0 = startrelax
) {
  reds <- fIs_amp(case0, k, fIs_reduction_at_post)
  vals <- c(list(rep(1-fIs_reduction_at_post, 16)), lapply(reds, function(d) rep(1-d, 16)))
  tms <- d0 + 0:(length(vals)-2)
  # vals <- c(list(rep(1-fIs_reduction_at_post, 16)), list(rep(0, 16)))
  return(list(
    parameter = "fIs",
    pops = numeric(),
    mode = "multiply",
    values = vals,
    times = day0 + c(model_t, tms)
    #times = day0 + c(model_t, model_t+30)
  ))
}

sim_step <- function(p, comps = "cases") cm_simulate(p, 1, model_seed = 42L)$dynamics[
  compartment %in% comps
]

sim_consolidate <- function(dt, grp = c("t")) dt[, .(value = sum(value)), by = grp]

us <- function(sdt.row, umod = sdt.row[1, umod]) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt.row))], each = 2)*umod
}

ys <- function(sdt.row) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt.row))], each = 2)
}

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

.cl <- makeCluster(getDTthreads())
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, { 
  require(data.table)
  require(optimization)
})

span <- nrow(fits.dt)
#' span <- 2

est <- rbindlist(parLapply(.cl, X = 1:span, function(i) {
  suppressPackageStartupMessages({
    source(file.path(cm_path, "R", "covidm.R"))
  })
  sdt <- fits.dt[i]
  testpop <- params;
  testpop$pop[[1]]$y <- ys(sdt)
  testpop$pop[[1]]$u <- us(sdt)
  testpop$pop[[1]]$seed_times <- intros[i == sid, t]
  testpop$schedule[[2]] <- add_fIs(sdt$case0, sdt$k, sdt$sympt)
  sim_consolidate(sim_step(testpop))[order(t), .(sample = i, date = t + day0, rv = value, value = value*sdt$asc, asc = sdt$asc)]
}))

parplot <- ggplot(est) +
 aes(date, rv, group = sample) +
 geom_line(aes(color="sim. total"), alpha = 0.2) +
 geom_line(aes(y=value, color="sim. ascertained"), alpha = 0.2) +
 geom_line(
    aes(date, croll, color = "observed"),
    data = case.dt,
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = as.Date(c("2020-03-01",NA))) +
 scale_color_manual(
   name=NULL,
   breaks = c("sim. total", "sim. ascertained", "observed"),
   values = c("firebrick", "goldenrod", "black")
 ) +
 scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
 scale_y_log10(
   "Cases",
   breaks = 10^c(0,3,6),
   minor_breaks = 10^(1:6),
   labels = scales::label_number_si()
  ) +
 theme_minimal()

ggsave(tail(.args, 1), parplot, height = 6.5, width = 15, units = "in", dpi = 900)
