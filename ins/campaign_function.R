
.debug <- file.path("analysis", "ins")
.args <- if (interactive()) file.path(
  .debug[1], "campaign_function.rds"
) else commandArgs(trailingOnly = TRUE)

campf <- function(
  target_id = NULL, recurring_id = NULL, camp_dur_id = NULL,
  covidmpars = NULL, tvax, v2delay = NA, horizon = 10
) {
  total_doses <- list(
    over60 = function(cmpars, coverage = 0.9) {
      doses <- cmpars$pop[[1]]$size
      doses[-(13:16)] <- 0
      return(round(doses*coverage))
    },
    plus15to60at20 = function(cmpars, coverage = 0.2) {
      basedoses <- total_doses[[1]](cmpars)
      adddoses <- cmpars$pop[[1]]$size
      adddoses[-(4:12)] <- 0
      return(basedoses + round(adddoses*coverage))
    },
    plus15to60at60 = function(cmpars, coverage = 0.6) total_doses[[2]](cmpars, coverage),
    from15 = function(cmpars, ...) {
      reftot = sum(total_doses[[2]](cmpars))
      totpop = sum(cmpars$pop[[1]]$size[-(1:3)])
      coverage = reftot/totpop
      doses <- round(cmpars$pop[[1]]$size*coverage)
      doses[1:3] <- 0
      return(doses)
    }
  )
  #' the resulting schedule will need to be modified:
  #'  - first, shift all times by t_vax, when vaccination starts relative
  #'  to t0
  #'  - optionally: duplicate schedule item, update parameter to v2, and offset
  #'  by between dose delay
  campaign_duration <- round(365*c(0.5, 1))
  periods <- c(2, 5)
  yrs <- 10
  recurrence <- list(
    everyother = function(duration, doses, tv, v2 = v2delay, period = periods[1]*365, horizon = yrs*365) {
      basevalues <- list(round(doses/duration), rep(0, length(doses)))
      basetimes <- c(0, duration)
      reps <- floor(horizon/period)
      vs <- rep(basevalues, reps)
      ts <- basetimes + rep(seq(0,reps-1)*period, each = length(basetimes))
      vaxschedule <- list(list(
        parameter = 'v', pops = 0, mode = 'assign',
        times = ts + tv,
        values = vs
      ))
      if (!is.na(v2)) {
        vaxschedule[[2]] <- vaxschedule[[1]]
        vaxschedule[[2]]$parameter <- 'v2'
        vaxschedule[[2]]$times <- vaxschedule[[2]]$times + v2
      }
      return(vaxschedule)
    },
    every5 = function(duration, doses, tv, v2 = v2delay, period = periods[2]*365, horizon = yrs*365) recurrence[[1]](
      duration, doses, tv, v2, period, horizon
    )
  )
#  browser()
  if (is.null(target_id)) {
    return(list(
      doses_per_campaign = if (
        !is.null(covidmpars)
      ) sapply(1:length(total_doses), function(i) sum(total_doses[[i]](covidmpars))) else names(total_doses),
      number_of_campaigns = horizon/periods, # how many recurrences
      campaign_durations = campaign_duration,
      horizon = horizon
    ))
  } else {
    return(
      recurrence[[recurring_id]](
        campaign_duration[camp_dur_id],
        total_doses[[target_id]](covidmpars),
        tvax, v2delay
      )
    )
  }
}

saveRDS(campf, tail(.args, 1))