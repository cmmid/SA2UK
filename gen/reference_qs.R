
suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
  require(qs)
})

.debug <- c("analysis", "NGA")
.args <- if (interactive()) sprintf(c(
  "%s/ins/covidm_fit_yu.qs",
  "%s/ins/sampling.json",
  "%s/gen/pops/%s.rds",
  "%s/gen/intervention_timing/%s.rds",
  "%s/gen/mobility.rds",
  .debug[2],
  "covidm",
  "%s/gen/yuqs/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

sampler <- read_json(.args[2])

SAMPLESIZE <- sampler$samplen

yu <- qread(.args[1])[
  sample(.N, SAMPLESIZE, replace = T)
][, sample_id := 1:.N ]

pop <- readRDS(.args[3])
timings <- readRDS(.args[4])

tariso <- tail(.args, 3)[1]

mob.dt <- readRDS(.args[5])[iso3 == tariso]

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

#' TODO parallelize this?
#' annoying bit is exporting the relevant covidm parts vs loading it in parallel vs how much time savings accrue
qs.dt <- yu[, {
  mod <- as.numeric(.SD[1,])
  ngm <- cm_eigen_ngm(
    pop,
    uval = mod[1:16], yval = mod[17:32],
    contact_reductions = c_reductions
  )
  ngmfrac <- ngm$frac[seq(1,15,by=2)] + ngm$frac[seq(2,16,by=2)]
  names(ngmfrac) <- sprintf("o_%i0", (1:length(ngmfrac))-1)
  c(.(baseR = ngm$R0), as.list(ngmfrac))
}, by=sample_id, .SDcols = c(uids, yids)]

ret <- qs.dt[order(baseR), eqs := (1:.N)/.N ][yu, on=.(sample_id), nomatch = 0]

saveRDS(ret, tail(.args, 1))
