suppressPackageStartupMessages({
  require(data.table)
  require(wpp2019)
  require(countrycode)
})

.debug <- c("analysis", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/gen/intervention_timing/%s.rds",
  "%s/gen/pops/%s.rds", #' assembled from other inputs, no estimation
  "%s/gen/yuqs/%s.rds", #' estimated elsewhere
  "%s/est/sample/%s.rds", #' estimated
  "%s/ins/adj_data.rds", #' cleaned input
  "ene-ifr.csv", #' input sourced from mbevands estimate
  .debug[2],
  "covidm",
  "%s/est/introductions/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2
source(file.path(cm_path, "R", "covidm.R"))

tariso <- tail(.args, 3)[1]

data(popM)[countrycode(country_code)]

pop <- readRDS(.args[2])
samples.dt <- readRDS(.args[4])[period == 1]

ifr <- fread(.args[6])$ifr/100
pop_vs_ifr <- 2 #' pop age categories per IFR age category
ifr.ext <- rep(ifr, each = pop_vs_ifr)
#' however, will need to consolidate upper categories
consol.inds <- length(pop$pop[[1]]$size):length(ifr.ext)

samples.dt[, {
  #' umod unnecessary here - ss unaffected by R0
  us <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^u_", colnames(yuref))], each = 2)
  ys <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^y_", colnames(yuref))], each = 2)
  #' ASSUMING: transmission largely at steady state by the time the infections which lead
  #' to detected deaths has occurred
  ss <- cm_eigen_ngm(pop, uval = us, yval = ys)$ss
  #' determine the expected aggregate IFR

}, by=.(sample)]

Rs <- readRDS(.args[4])[era == "pre" & period == 1]
window_start <- readRDS(.args[7])[era == "pre" & period == 1, end]

pars <- readRDS(.args[4])

#' the stable age distribution of infection in the exponential growth period
#' == the principle eigenvector
#' not affected by rescaling R0
refngm <- cm_ngm(pars, us, ymod = ys)

agedist <- cm_ss_age_distro(refngm$ss) 
#' generation interval depends on age distro + model infectious period durations
gen_int <- cm_generation_time(pars, ymod = ys, ngm = refngm)

ref <- readRDS(.args[5])[
  iso3 == tariso, .SD
][1:(which.max(deaths > 0)+round(gen_int))]

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

#' MAGIC NUMBER WARNING
death_underreporting <- 2.7 #' assume initially no under-reporting vs later estimate of ~2.7 for SA

if (slc[.N, date] - inf2death_dur < window_start) {
  expansion <- Rs[,{
    gens <- cumsum(value^(0:10))
    ind <- which.max(gens > exp.infs)-1
    missing <- exp.infs/sum(value^(0:(ind-1)))
    infs <- missing*value^(ind-1)
    #' infections per OBSERVED death
    .(tot.infs = infs*death_underreporting)
  }, by=sample]
} else stop("not implemented")

intros <- expansion[, slc[,.(
  infections = ceiling(tot.infs*deaths)
), by=.(continent, iso3, date = date - inf2death_dur)], by=sample]

saveRDS(intros, tail(.args, 1))