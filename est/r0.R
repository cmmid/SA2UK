suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
  require(EpiNow2)
  require(future)
})

.debug <- c("analysis", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/ins/adj_data.rds",
  "%s/ins/sampling.json",
  "%s/gen/intervention_timing/%s.rds",
  "%s/gen/pops/%s.rds",
  "%s/gen/yuqs/%s.rds", # not used currently
  .debug[2],
  "%s/outputs/r0/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

sampler <- read_json(.args[2])

smps <- sampler$samplen * sampler$rtsamplemul
crs <- fcoalesce(
  as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")),
  getDTthreads()
)

tariso <- tail(.args, 2)[1]

case.dt <- readRDS(.args[1])[iso3 == tariso][, .(date, confirm = cases )]
fill.case <- case.dt[
  case.dt[, .(date = seq(min(date),max(date), by="day"))],
  on=.(date),
  .(date, confirm = fifelse(is.na(confirm), 0, confirm))
]

lims.dt <- readRDS(.args[3])
pop <- readRDS(.args[4])

#' TODO this approximation relies on
#'  - asymptomatic == presymptomatic + symptomatic infector waiting time distributions
#'  - (which comes from the average of those distros being equal + infectiousness for pre- and symptomatic being equal)
#'  - not true (in the model world) once fIs modified (by intervention)
#'  - not true (ibid) for *distribution*, given dIp + dIs variance, etc will differ from dIa variance
#'  - probably not true in the real world (fIp != fIs, because age specific fIX effects, etc)
#'  
#' Alternatives:
#'  - use NGM information to approximate the correct distribution of age & symptom combinations of *infectors*,
#'  weighted by contribution to creating *infectees*; downside: some gnarly bayes-theorem work to do here to get it right?
#'  - run the simulation, and extract times? not really setup for that - don't have the info about how long an
#'  infector has been infectious  

#' drawing infection generation times:
#'  - assumes non-depleting / infinite infectee population
#'  - therefore, for time of infection: each interval of infectivity has same
#'  hazard for producing infections
#'  - however, the relative distribution of those intervals is influenced by the infectious duration

p_of_having_an_interval <- 1-cumsum(c(0, pop$pop[[1]]$dIa))
p_infection_in_interval <- head(p_of_having_an_interval/sum(p_of_having_an_interval),-1)

gi_sample <- sample(
  seq(0,by=pop$time_step,length.out = length(p_infection_in_interval)),
  size = 1e4, replace = T,
  prob = p_infection_in_interval
) + sample(
  seq(0,by=pop$time_step,length.out = length(pop$pop[[1]]$dE)),
  size = 1e4, replace = T,
  prob = pop$pop[[1]]$dE
)

warning("...sampling generation interval.")

plan(multicore, workers = crs)

generation_time <- bootstrapped_dist_fit(
  gi_sample, dist = "gamma", samples = 4000,
  max_value = (length(pop$pop[[1]]$dIa)+length(pop$pop[[1]]$dE)-2)*pop$time_step,
  verbose = TRUE
)

warning("...sampling incubation period.")

incubation_period <- estimate_delay(
  sample(length(pop$pop[[1]]$dE), 1e4, replace = TRUE, prob = pop$pop[[1]]$dE)*pop$time_step
)

plan(sequential)

est.window <- 30

est.qs <- unique(c(pnorm(seq(-1,0,by=0.25)), pnorm(seq(0,1,by=0.25))))

Rtcalc <- function(case.dt, gp = NULL, rt = rt_opts()) estimate_infections(
  reported_cases = case.dt,
  generation_time = generation_time,
  delays = delay_opts(incubation_period),
  rt = rt,
  stan = stan_opts(
    samples = smps*2,
    warmup = 200, 
    cores = crs,
    control = list(adapt_delta = 0.99, max_treedepth = 20)
  ),
  gp = gp,
  verbose = TRUE,
  CrIs = est.qs,
  horizon = 0
)

processRt <- function(
  rt, keep.start, keep.end,
  era.labels, tarvar = "R"
) rt$samples[variable == tarvar, .(value, variable), by=.(sample, date)][
  between(date, keep.start, keep.end)
][, era := eval(era.labels), by= sample ]

results <- list()

for (grpi in lims.dt[era != "censor", sort(unique(period))]) {
  sublims <- lims.dt[period == grpi]
  daterange <- sublims[, range(c(start, end))]
  incslice <- fill.case[between(date, daterange[1]-est.window, daterange[2]+est.window)]
  breakbased <- sublims[era == "transition", .N]
  if (breakbased) {
    incslice[, era := "tail" ]
    for (e in c("post", "transition", "pre")) {
      incslice[date <= sublims[era == e]$end, era := e  ]
    }
    incslice[date < sublims[era == e]$start, era := "censor" ]
    incslice[, breakpoint := era %in% c("censor", "transition", "tail") ]
    incslice[era == "post", breakpoint := c(TRUE, rep(FALSE, .N-1))]
    results[[grpi]] <- processRt(
      Rtcalc(incslice),
      sublims[era == "pre", end], sublims[era == "post", start],
      expression(c("pre",rep("transition",.N-2),"post"))
    )[, period := grpi ]
  } else {
    results[[grpi]] <- processRt(
      Rtcalc(incslice, rt = NULL, gp = gp_opts()),
      sublims[, min(start)], sublims[, max(end)],
      era.labels = "relaxation", tarvar = "infections"
    )[, period := grpi ]
  }
}

ret <- rbindlist(results)

saveRDS(ret, tail(.args, 1))
