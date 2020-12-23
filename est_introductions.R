suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/covidLMIC", "UGA")
.args <- if (interactive()) sprintf(c(
 # "%s/outputs/Rg.rds", #' estimated
  "%s/outputs/r0/%s.rds", #' estimated
  "%s/inputs/populations.rds", #' assembled from other inputs, no estimation
  "%s/inputs/pops/%s.rds", #' assembled from other inputs, no estimation
  "%s/inputs/ecdc_data.rds", #' cleaned input
  "ene-ifr.csv", #' input sourced from mbevands estimate
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "%s/outputs/introductions/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

#' TODO replace w/ .args 1 and 2
Rg <- 2
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

#' going to back project introductions to earlier date range
#' if the are report cases *prior* to that date,
#' use reported cases in the introduction window and proposed infections to
#' create inflation factor
#' inflate the first serial interval worth of cases by that factor

first.death.date.offset <- ref[
  which.max(deaths > 0), {
    tmp <- date - inf2death_dur
    gen <- 0
    while (tmp > window_start) {
      tmp <- tmp - gen_int
      gen <- gen + 1
    }
    .(date = tmp, gen = gen, offset = round(inf2death_dur + gen*gen_int))
  }
]

#' per actual death
est.infs <- function(
  Rg, Rl,
  tar = exp.infs,
  gen = ceiling(log(tar)/log(Rg)),
  #' MAGIC NUMBER WARNING
  inf.dur = idur/gen_int, #' in generation intervals
  genback = first.death.date.offset$gen[1]
) {
  It <- Rg^(0:gen)
  Ct <- Reduce(function(l, r) l*Rl+r, x = It, accumulate = TRUE)
  tots <- cumsum(Ct)
  infs.ind <- which.max(tots > tar)
  #' this is first gen that exceeds target - is target closer to previous gen or this one?
  # if ( (tar - tots[infs.ind-1]) < (tots[infs.ind]-tar)) infs.ind <- infs.ind - 1
  Ct[infs.ind - genback]
}

#' MAGIC NUMBER WARNING
death_underreporting <- 2.7

#' infections per OBSERVED death
tot.infs <- est.infs(Rg, Rl)*death_underreporting

# if (ref[date < first.death.date.offset, sum(cases) > 0]) {
#   infl.fact <- (tot.infs*ref[which.max(deaths > 0):.N][1:round(gen_int), sum(deaths)])/ref[between(date, first.death.date.offset, first.death.date.offset+gen_int), sum(cases)]
#   intros <- ref[which.max(cases > 0):.N][1:gen_int][, infections := round(cases*infl.fact) ][, .(infections), keyby=.(continent, iso3, date) ]
# } else {
  slc <- ref[which.max(deaths > 0):.N][1:gen_int][deaths > 0]
  intros <- slc[,.(
    i.date = date - first.death.date.offset$offset + 0:(gen_int-1)-ceiling(gen_int/2),
    infections = rep(ceiling(deaths*tot.infs/gen_int), gen_int)
  ), by=.(continent, iso3, date)][, .(infections = sum(infections)), keyby=.(continent, iso3, date=trunc(i.date, "days"))]
# }

saveRDS(intros, tail(.args, 1))