suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/yuqs/%s.rds", #' estimated
  "%s/outputs/r0/%s.rds", #' estimated
  "%s/inputs/populations.rds", #' assembled from other inputs, no estimation
  "%s/inputs/pops/%s.rds", #' assembled from other inputs, no estimation
  "%s/inputs/epi_data.rds", #' cleaned input
  "ene-ifr.csv", #' input sourced from mbevands estimate
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "%s/outputs/introductions/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

yuref <- readRDS(.args[1])[order(eqs)][which.max(eqs >= .5)]
us <- rep(yuref[, as.numeric(.SD), .SDcols = grep("^u_", colnames(yuref))], each = 2)
ys <- rep(yuref[, as.numeric(.SD), .SDcols = grep("^y_", colnames(yuref))], each = 2)
Rl <- readRDS(.args[2])[, median(pre)]
window_start <- readRDS(.args[7])[era == "pre", end]

tariso <- tail(.args, 2)[1]
pop <- readRDS(.args[3])[iso3 == tariso]
pars <- readRDS(.args[4])

load("NGM.rda")

#' the stable age distribution of infection in the exponential growth period
#' == the principle eigenvector
#' not affected by rescaling R0
refngm <- cm_ngm(pars, us, ymod = ys)

agedist <- cm_ss_age_distro(refngm$ss) 
#' generation interval depends on age distro + model infectious period durations
gen_int <- cm_generation_time(pars, ymod = ys, ngm = refngm)

ref <- readRDS(.args[5])[iso3 == tariso, .SD][date >= window_start][1:(which.max(deaths > 0)+round(gen_int))]

ifr <- fread(.args[6])$ifr/100
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
exp.ifr <- sum((ifr.ext*agedist.ext))

# if exp.ifr purely from initial random exposure:
subus <- us[seq(2,by=2,length.out=8)]
subus[9:10] <- subus[8]
alt.exp.ifr <- sum(ifr*subus/sum(subus))
mix.exp.ifr <- sqrt(exp.ifr*alt.exp.ifr)

#' assume early transmission is a mixture
exp.infs <- 1/mix.exp.ifr

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