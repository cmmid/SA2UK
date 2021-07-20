
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "analysis"
.args <- if (interactive()) sprintf(c(
  "%s/ins/epi_data.rds", # this is the JHU
  "%s/est/adj_data.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

raw.dt <- readRDS(.args[1])

#' want to impute daily incidence for iso3s where there is an
#' apparent shift from daily reporting to some other sporadic
#' interval
#' 
#' for those intervals, asserting 0s are actually NAs
#' 
#' going to impute replacements for those 0s and reductions
#' to the values at the end points

#' want: iso3s with "apparent sporadic reporting intervals" (ASRI)
#' define ASRI as
#'  - an interval of days
#'  - starting and ending with a value above a incidence threshold
#'  - with all intermediate values == 0
#'  - and the length of the intervening period longer than some duration threshold
#'  - OR for even a duration of a single non-reporting day (or other minimum time unit),
#'  a sufficiently extreme day
#'
#' propose: incidence threshold >= 10
#' propose: duration threshold >= 2 days
#' propose: extreme incidince: 10x basic incidence

#' @param series, integer vector which begins with a non-zero value
#' @param inc.threshold, integer > 0, the incidence threshold for ends of the ASRI
#' @param dur.threshold, integer > 0, the duration threshold for ASRI
#' @param ext.threshold, integer > inc.threshold, the incidence threshold at which to ignore the duration threshold
#'   i.e., when the incidence is this high, even single day intermediate zeros are considered asri
detect.asri <- function(
  series,
  inc.threshold = 10,
  dur.threshold = 2,
  ext.threshold = inc.threshold*10
) {
  lv <- rle(series)
  len <- length(lv$values)
  #' can't address tail 0 series (no end value to distribute back)
  #' so going to replace those with NAs
  if (lv$values[len] == 0) {
    lv$values <- lv$values[-len]
    lastl <- lv$lengths[len]
    lv$lengths <- lv$lengths[-len]
  }
  
  #' all the 0 intervals
  intervening <- which(lv$values == 0)
  validstarts <- which(
    (lv$lengths[intervening - 1] == 1) & (#' non-runs of values
      # above the duration threshold & lower incidence threshold
      (lv$lengths[intervening] > dur.threshold & lv$values[intervening - 1] >= inc.threshold) |
      # OR above the extreme incidence threshold
      lv$values[intervening - 1] >= ext.threshold
    )  
  )
  #' subset to valid starts, then check for valid ends
  intervening <- intervening[validstarts]
  validends <- which(
    (lv$lengths[intervening + 1] == 1) & (
      (lv$lengths[intervening] > dur.threshold & lv$values[intervening + 1] >= inc.threshold) |
      lv$values[intervening + 1] >= ext.threshold
    )
  )
  intervening <- intervening[validends]
  
  #' change the RLE object to be about the intervals
  #' initially, default everything to not an interval (FALSE)
  lv$values[1:length(lv$values)] <- FALSE
  #' set the remaining zero intervals w/ their start + end as ASRI
  lv$values[c(intervening-1, intervening, intervening+1)] <- TRUE
  #' if necessary, add back final interval
  if (len > length(lv$values)) {
    lv$values[len] <- NA
    lv$lengths[len] <- lastl
  }
  #' return the series (rather than RLE of the series)
  inverse.rle(lv)
}

#' for each iso, look at the time series subset from first case
#' and evaluate ASRI
asri.dt <- raw.dt[,
  if (sum(cases) != 0) .SD[which.max(cases > 0):.N,
    .(date, asri = detect.asri(cases))
  ] else .SD[,.(date, asri = NA_integer_)], by=iso3
]

#' if there is a short gap in asri, link up series; basically assumes
#' that short intervals between large asri intervals should be considered
#' potential asri as well. will only end up affecting 0 run -> non zero ends
#' in that period 
#' 
#' @param asri.series series of 0 or 1, no NAs.
#' @param merge.threshold numeric between 0 and 1; if the gap between two asri's
#'   is less than this percentage of new total length, merge asri's
#' @return the modified series 
fill.asri <- function(asri.series, merge.threshold = 0.25) {
  lv <- rle(asri.series)
  from1 <- which.max(lv$values)
  while(from1 < (length(lv$values)-1)) {
    tot <- sum(lv$lengths[from1:(from1+2)])
    frac <- lv$lengths[from1+1]/tot
    if (frac < merge.threshold) {
      lv$lengths[from1] <- tot
      lv$lengths <- lv$lengths[-(from1 + (1:2))]
      lv$values <- lv$values[-(from1 + (1:2))]
    } else { #' if not merging, go to the next window
      if (from1 != length(lv$lengths)) {
        from1 <- from1 + which.max(lv$values[-(1:from1)])
        if (from1 == 1) from1 == length(lv$lengths)
      }
    }
  }
  inverse.rle(lv)
}

asri.dt[, merge.asri := asri ]
asri.dt[!is.na(merge.asri), merge.asri := fill.asri(merge.asri), by=iso3 ]
#' run backwards as well
asri.dt[!is.na(merge.asri), merge.asri := rev(fill.asri(rev(merge.asri))), by=iso3 ]

raw.dt[asri.dt, on=.(iso3, date), asri := merge.asri ]

#' TODO recapitulate ASRI in terms of geometric mean / standard deviation,
#' and sigma-based departures?

#' @examples 
#' has_asri <- raw.dt[continent == "Africa"][,any(asri == 1, na.rm = TRUE),by=iso3][V1 == TRUE][order(iso3), iso3]
#' for (is in has_asri) {
#'   print(ggplot(raw.dt[iso3 == is][which.max(cases > 0):.N]) +
#'   aes(date, cases, color = c("normal","asri")[asri+1]) +
#'   geom_point() + ggtitle(is) +
#'   theme_minimal() + coord_cartesian(xlim = as.Date(c("2020-03-01", NA))) +
#'   scale_y_log10() +
#'   scale_x_date(name = NULL, date_breaks = "month", minor_breaks = NULL, date_labels = "%b %y"))
#' }

#' now have identified asri 
#'  - need to replace 0 run -> non-zero
#'  - with ... run -> lower non-zero
#'  
#' question is what to use for ...?
#' 
#' 1. cannot reduce final value to 0
#' 2. require that the imputed run is total preserving
#' 3. enforce a process assumption: constant daily multiplier r over that period
#'   [implies geometric growth or decline (or stasis w/ r ~= 1)]
#' 4. starting from a value at the outset of the interval
#' 5. raw start values may be artificially high or low (
#'   from prior sporadic reporting, from partial reporting starting on that day, etc
#' )
#' 
#' so:
#'  - estimate a starting value from rolling average process
#'  - from that starting value, compute the daily r that results
#'  in the correct total for the period
#'  - for each day, floor to integer value, and roll leftovers to
#'  the next day
#'  - if the last day is zero, shift one case from most recent non-zero day to end
#'  
#' rolling average process:
#'  - assume in general a poisson process, convoluted with probability p of observation
#'  - give that process a daily multiplicative increase or decrease
#'  - constant lambda, p over rolling window

frollstarting <- function(x) {
  n <- length(x)
  ml <- optim(
    c(r=1, lm=1),
    function(par, obs, n) -sum(dpois(obs, par[2]*(par[1]^(1:n - 1)), log = TRUE)),
    obs = x, n = length(x),
    lower = c(1e-6, 1e-6),
    upper = c(2, max(x)*5),
    method = "L-BFGS-B"
  )
  ml$par[2]*ml$par[1]^(n-1)
}

matchfun <- function(r, dt, inc0, incf) inc0*sum(r^(1:dt)) - incf

impute.single <- function(series) {
  n <- length(series)
  inc0 <- series[1]
  incf <- series[n]
  r <- uniroot(matchfun, c(0, incf/inc0), dt = n-1, inc0 = inc0, incf = incf, tol = .Machine$double.eps^0.5)$root
  incn <- inc0*r^(1:(n-1))
  incd <- floor(incn)
  delinc <- cumsum(incn-incd)
  idel <- floor(delinc)
  adds <- c(0, diff(idel))
  fininc <- incd + adds
  del <- incf - sum(fininc)
  if (del != 0) fininc[n-1] <- fininc[n-1]+del
  fininc
}

#' @param series numeric vector, non-zero start, interspersed with non-zeros
impute.asri <- function(series) {
  lv <- rle(series)
  intervals <- which(lv$values == 0)
  iends <- cumsum(lv$lengths)[intervals]+1
  istarts <- cumsum(lv$lengths)[intervals-1]
  for (i in 1:length(istarts)) {
    series[(istarts[i]+1):iends[i]] <- impute.single(series[istarts[i]:iends[i]])
  }
  series[-1]
}

raw.dt[, adj := {
  if (any(asri == 1, na.rm = TRUE)) {
    lv <- rle(asri)
    runstarts <- cumsum(lv$lengths)[which(lv$values == 1)-1]+1
    inc0 <- numeric(length(runstarts))
    for (i in 1:length(runstarts)) inc0[i] <- frollstarting(cases[runstarts[i]-6:0])
    runends <- cumsum(lv$lengths)[which(lv$values == 1)]
    res <- cases
    for (i in 1:length(runstarts)) {
      series <- cases[runstarts[i]:runends[i]]
      series[1] <- inc0[i]
      res[(runstarts[i]+1):runends[i]] <- impute.asri(series)
    }
    # print(iso3)
    res
  } else cases
}, by=iso3 ]

#' @examples 
#' for (is in has_asri) {
#'   print(ggplot(raw.dt[iso3 == is][which.max(cases > 0):.N]) +
#'   aes(date, cases, color = c("normal","asri")[asri+1], alpha = c("normal","asri")[asri+1]) +
#'   geom_point(aes(y=adj, alpha = "normal"), data = function(dt) dt[asri == 1], pch = 17) +
#'   geom_point() +
#'   geom_line(data = function(dt) dt[, .(date, cases = frollmean(cases, 7, align = "center"), asri = 0)]) +
#'   ggtitle(is) +
#'   theme_minimal() + coord_cartesian(xlim = as.Date(c("2020-03-01", NA))) +
#'   scale_color_discrete(NULL) +
#'   scale_alpha_manual(NULL, values = c(asri = 0.2, normal = 1)) +
#'   scale_y_log10() +
#'   scale_x_date(name = NULL, date_breaks = "month", minor_breaks = NULL, date_labels = "%b %y"))
#' }

saveRDS(raw.dt, tail(.args, 1))
