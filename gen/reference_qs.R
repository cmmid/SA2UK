
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
})

.debug <- c("analysis", "ETH")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/covidm_fit_yu.qs",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/inputs/mobility.rds",
  .debug[2],
  "covidm",
  "%s/inputs/yuqs/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

SAMPLESIZE <- 2000

yu <- qread(.args[1])[sample(.N, SAMPLESIZE, replace = T)]
yu[, sample_id := 1:.N ]
pop <- readRDS(.args[2])

timings <- readRDS(.args[3])

tariso <- tail(.args, 3)[1]

mob.dt <- readRDS(.args[4])[iso3 == tariso]

#' MAGIC NUMBER: average delay from infection -> reporting
rep_delay <- 7 

c_reductions <- timings[(period == 1) & (era == "pre")][
  mob.dt, on=.(iso3), allow.cartesian = TRUE
][
  between(date, start-rep_delay, end-rep_delay)
][,
  1-c(
    home = 1,
    work = prod(work_multiplier)^(1/.N),
    other = prod(other_multiplier)^(1/.N),
    school = mean(school_multiplier)
  )
]

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

source(file.path(cm_path, "R", "covidm.R"))

uids <- rep(grep("^u_", colnames(yu)), each = 2)
yids <- rep(grep("^y_", colnames(yu)), each = 2)

qs.dt <- yu[, {
  umod <- as.numeric(.SD[1,])[1:16]
  ymod <- as.numeric(.SD[1,])[17:32]
  ngm <- cm_eigen_ngm(
    pop,
    uval = umod,
    yval = ymod,
    contact_reductions = c_reductions
  )
  si <- cm_generation_interval(
    pop$pop[[1]], yval = ymod,
    eigen_vector = ngm$ss, time_step = pop$time_step
  )
  .(baseR = ngm$R0, si = si)
}, by=sample_id, .SDcols = c(uids, yids)]

ret <- qs.dt[order(baseR), eqs := (1:.N)/.N ][yu, on=.(sample_id), nomatch = 0]

saveRDS(ret, tail(.args, 1))
