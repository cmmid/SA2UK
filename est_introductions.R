suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
 # "%s/outputs/Rg.rds", #' estimated
  "%s/outputs/r0/%s.rds", #' estimated
  "%s/inputs/populations.rds", #' assembled from other inputs, no estimation
  "%s/inputs/pops/%s.rds", #' assembled from other inputs, no estimation
  "%s/inputs/epi_data.rds", #' cleaned input
  "ene-ifr.csv", #' input sourced from mbevands estimate
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "%s/outputs/introductions/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

Rl <- readRDS(.args[1])[era == "pre", med]

window_start <- readRDS(.args[6])[era == "pre", end]

tariso <- tail(.args, 2)[1]
pop <- readRDS(.args[2])[iso3 == tariso]
pars <- readRDS(.args[3])

load("NGM.rda")

#' the stable age distribution of infection in the exponential growth period == the principle eigenvector
agedist <- cm_ss_age_distro(cm_ngm(pars)$ss) 
#' generation interval depends on age distro + model infectious period durations
gen_int <- cm_generation_time(pars)

ref <- readRDS(.args[4])[iso3 == tariso, .SD][date >= window_start][1:(which.max(deaths > 0)+round(gen_int))]

ifr <- fread(.args[5])$ifr/100
# MAGIC NUMBER WARNING
pop_vs_ifr <- 2 #' pop age categories per IFR age category

ifr.ext <- rep(ifr, each = pop_vs_ifr)
#' ASSERT: ifr.ext has more categories than agedist / pars$pop[[1]]$group_names
#' so need to expand last category by population there, but reference population
pd <- pop[length(agedist):.N, { 
  tmp <- pop
  tmp[.N-1] <- tmp[.N-1]+tmp[.N]
  tmp[-.N]
}]

pd <- pd/sum(pd)
agedist.ext <- c(agedist[-length(agedist)], agedist[length(agedist)] * pd)
exp.ifr <- sum(ifr.ext*agedist.ext)
#' per actual death
exp.infs <- 1/exp.ifr

idur <- sum(pars$pop[[1]]$dIa*seq(0,by=0.25,length.out = length(pars$pop[[1]]$dIa)))

#' MAGIC NUMBER WARNING
inf2death_dur <- 22

slc <- ref[which.max(deaths > 0):.N][1:(gen_int/2)][deaths > 0]

if (slc[.N, date] - inf2death_dur < window_start) {
  gens <- cumsum(Rl^(0:10))
  ind <- which.max(gens > exp.infs)-1
  missing <- exp.infs/sum(Rl^(0:(ind-1)))
  infs <- missing*Rl^(ind-1)
} else stop("not implemented")

#' MAGIC NUMBER WARNING
death_underreporting <- 1 #' assume initially no under-reporting 2.7

#' infections per OBSERVED death
tot.infs <- infs*death_underreporting

intros <- slc[,.(
  infections = ceiling(tot.infs*deaths)
), by=.(continent, iso3, date = date - inf2death_dur)]
# }

saveRDS(intros, tail(.args, 1))