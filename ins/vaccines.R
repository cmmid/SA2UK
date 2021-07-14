
require(data.table)

.debug <- file.path("analysis", "ins")
.args <- if (interactive()) sprintf(c(
  file.path("%s", "vaccines.rds")
), .debug[1]) else commandArgs(trailingOnly = TRUE)

#' placeholder for efficacy type
#' if efficacy is all-or-0, then the model has people => V1 and become immune to infection OR NOT
#' that is, there is currently no direct support for a version where people have immunizing prevention
#' of infection OR some enhanced reduction in disease. Would need distinct V routes (all-or-0 route, leaky route)
vaccinef <- function(
  product = NULL,
  vax_exp_dur_days,
  covidmpop,
  eff_type = "leaky"
) {
  if (is.null(product)) return(c("pfizer", "az", "j&j"))

  vaccine <- { data.table(
    label = c(
      "pfizer", "pfizer",
      "az", "az",
      "j&j"
    ),
    dose = c(
      1:2, 
      1:2,
      1
    ),
    #' baseline efficacy against symptoms
    eff.inf = c(
      0.7, 0.85,
      0.67, 0.68,
      0.5
    ),
    #' additional efficacy, conditional on infection, against symptoms
    eff.sym = c(
      0, .267,
      0, .267,
      0.5
    ),
    #' additional efficacy, conditional on infection (not symptoms?), against severe disease
    eff.sev = c(
      0.483, 0.333,
      0.483, 0.333,
      0.75
    )
  )[label == product] }
  
  if (eff_type == "leaky") {
    newpop <- covidmpop
    newpop$pop[[1]]$ev = newpop$pop[[1]]$ev2 = rep(1, newpop$pop[[1]]$n_groups)
    if (missing(vax_exp_dur_days) | is.infinite(vax_exp_dur_days[1])) {
      newpop$pop[[1]]$wv = newpop$pop[[1]]$wv2 = rep(0, newpop$pop[[1]]$n_groups)
    } else {
      newpop$pop[[1]]$wv = newpop$pop[[1]]$wv2 = rep(vax_exp_dur_days[1], newpop$pop[[1]]$n_groups)
      if (length(vax_exp_dur_days) > 1) newpop$pop[[1]]$wv2 = rep(vax_exp_dur_days[2], newpop$pop[[1]]$n_groups)
    }
    
    oldprocesses <- newpop$processes
    #' copy the E processes
    #' we'll make duplicates of these for Ev (and Ev2, if necessary)
    Eprocesses <- Filter(function(p) p$source == "E", oldprocesses)
    
    #' manage vaccine-related changes to uv, yv by multiplying it from
    #' t = 0.
    #' CppChanges (i.e., using the posterior sample inputs) occurs
    #' prior to constructing the simulation, with schedule-based changes
    #' happening *during* construction
    #' this also allows setting different vaccine efficacy for after switch
    #' to new variants
    
    effschedule <- function(
      par = "uv",
      ts = 0,
      efficacy = 0,
      grps = newpop$pop[[1]]$n_groups
    ) list(
      parameter = par, pops = 0, mode = 'multiply',
      times = ts,
      values = lapply(efficacy, function(e) rep(1-e, grps))
    )
  
    effprocess <- function(p, dose.dt, src = "Ev") {
      p$source <- src
      if (any(p$names == "death")) { # death process
        p$prob["death", ] <- p$prob["death", ] * (1 - dose.dt$eff.sev)
        p$prob["null", ] <- 1 - p$prob["death", ]
      } else if (any(p$names == "to_icu")) {
        p$prob[c("to_icu", "to_nonicu"), ] <- p$prob[c("to_icu", "to_nonicu"), ] * (1 - dose.dt$eff.sev)
        p$prob["null", ] <- 1 - colSums(p$prob[c("to_icu", "to_nonicu"), ])
      } else stop("don't know how to deal with process %s!", paste(p$names, collapse = "|"))
      p
    }
    
    #' based on number of doses,
    #' TODO: and number of variants...?
    #' note, apparently no way to change processes on a schedule
    d1 <- vaccine[dose == 1,]
    vaxeffmods <- list(
      effschedule(eff = d1$eff.inf), # efficacy against infection reduces uv
      effschedule(par = "yv", eff = d1$eff.sym) # efficacy against symptoms reduces yv
    )
    # and efficacy against severe disease reduces hospitalization and death probabilities
    vaxxedprocesses <- lapply(Eprocesses, effprocess, dose.dt = d1)
    
    if (dim(vaccine)[1] > 1) {
      d2 <- vaccine[dose == 2,]
      vaxeffmods <- c(vaxeffmods, list(
        effschedule(par = "uv2", eff = d2$eff.inf), # efficacy against infection reduces uv
        effschedule(par = "yv2", eff = d2$eff.sym) # efficacy against symptoms reduces yv
      ))
      vaxxedprocesses <- c(vaxxedprocesses, lapply(Eprocesses, effprocess, dose.dt = d2, src = "Ev2"))
    }
    
    newpop$processes <- c(vaxxedprocesses, newpop$processes)
    newpop$schedule <- c(newpop$schedule, vaxeffmods)
    
    return(newpop)
  } else stop(sprintf("efficacy type '%s' is not supported."))

}

saveRDS(vaccinef, tail(.args, 1))