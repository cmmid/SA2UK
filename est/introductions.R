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

#' step 1: get IFR
ifr <- fread(.args[6])$ifr/100
#' because the input doesn't match the breakdown in model,
#' we have to:
#'  - expand it (from 10 year windows to 5 year windows)
pop_vs_ifr <- 2 #' pop age categories per IFR age category
ifr.ext <- rep(ifr, each = pop_vs_ifr)
#'  - consolidate the upper end
#'  - the upper-upper end is outside the last population window,
#'   so need to get the relevant pop from elsewhere (wpp2019)
popMF <- 
  as.data.table(get(data(popM)))[countrycode(country_code, "iso3n", "iso3c") == tariso, as.numeric(`2020`)] +
  as.data.table(get(data(popF)))[countrycode(country_code, "iso3n", "iso3c") == tariso, as.numeric(`2020`)]
if (length(ifr.ext) < length(popMF)) {
#' of course this population data is *also* mismatched
  popMF[length(ifr.ext)] <- sum(popMF[length(ifr.ext):length(popMF)])
  popMF <- popMF[1:length(ifr.ext)]
}

pop <- readRDS(.args[2])
#' now, consolidate ifr to match population categories
consol.inds <- length(pop$pop[[1]]$size):length(ifr.ext)
ifr.ext[consol.inds[1]] <- sum(popMF[consol.inds]*ifr.ext[consol.inds])/sum(popMF[consol.inds])
ifr.ext <- ifr.ext[1:length(pop$pop[[1]]$size)]
#' that gives us IFR by age, which we can then weight according
#' to steady-state distribution of infections in each generation
#' 
#' now need mean generation interval

p_of_having_an_interval <- 1-cumsum(c(0, pop$pop[[1]]$dIa))
p_infection_in_interval <- head(p_of_having_an_interval/sum(p_of_having_an_interval),-1)

gen_int.mu <- (
  sum(pop$pop[[1]]$dE*((1:length(pop$pop[[1]]$dE))-1)) + sum(p_infection_in_interval*((1:length(p_infection_in_interval))-1))
)*pop$time_step

#' get relevant cases and deaths
ref <- readRDS(.args[5])[
  iso3 == tariso, .SD
][1:(which.max(deaths > 0)+round(gen_int.mu))]

#' get samples
samples.dt <- readRDS(.args[4])[period == 1]

window_start <- readRDS(.args[1])[era == "pre" & period == 1, end]
slc <- ref[which.max(deaths > 0):.N][1:(gen_int.mu/2)][deaths > 0]

#' MAGIC NUMBERS WARNING
inf2death_dur <- 22
death_underreporting <- 2.7 #' assume initially no under-reporting vs later estimate of ~2.7 for SA

if (slc[.N, date] - inf2death_dur < window_start) {
  expansion <- samples.dt[, {
    yuref <- colnames(.SD)
    #' umod unnecessary here - ss unaffected by R0
    us <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^u_", yuref)], each = 2)
    ys <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^y_", yuref)], each = 2)
    #' ASSUMING: transmission largely at steady state by the time the infections which lead
    #' to detected deaths has occurred
    ss <- cm_eigen_ngm(pop, uval = us, yval = ys)$ss
    #' determine the expected aggregate IFR
    ifr.mu <- sum(ifr.ext*ss)/sum(ss)
    #' if instead, weighted purely by proportion of population / susceptibility
    # sus.weight <- (pop$pop[[1]]$size*us)/sum(pop$pop[[1]]$size*us)
    # ifr.mu2 <- sum(ifr.ext*sus.weight)
    R0 <- pre
    mu.infections <- (1/ifr.mu)
    gens <- cumsum(R0^(0:10))
    # how many generations until enough infections to (eventually) have a death?
    ind <- which.max(gens > mu.infections)-1
    missing <- mu.infections/sum(R0^(0:(ind-1)))
    infs <- missing*R0^(ind-1)
    .(tot.infs = infs*death_underreporting)
  }, by=.(sample)]
} else stop("not implemented")

intros <- expansion[, slc[,.(
  infections = ceiling(tot.infs*deaths)
), by=.(continent, iso3, date = date - inf2death_dur)], by=sample]

saveRDS(intros, tail(.args, 1))