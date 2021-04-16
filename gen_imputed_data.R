
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/outputs/epi_mod_data.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

raw.dt <- readRDS(.args[1])

#' want: iso3s with "apparent sporadic reporting intervals"
#' define ASRI as
#'  - an interval of days
#'  - starting and ending with a value above a incidence threshold
#'  - with all intermediate values == 0
#'  - and the length of the intervening period longer than some duration threshold
#'
#' propose: incidence threshold > 10
#' propose: duration threshold >= 2 days

#' @param series, integer vector which begins with a non-zero value
#' @param inc.threshold, integer > 0, the incidence threshold for ends of the ASRI
#' @param dur.threshold, integer > 0, the duration threshold for ASRI
detect.asri <- function(series, inc.threshold = 10, dur.threshold = 2) {
  lv <- rle(series)
  len <- length(lv$values)
  if (lv$values[len] == 0) {
    lv$values <- lv$values[-len]
    lastl <- lv$lengths[len]
    lv$lengths <- lv$lengths[-len]
  }
  
  intervening <- which(lv$values == 0 & lv$lengths >= dur.threshold)
  validstarts <- which((lv$values[intervening - 1] >= inc.threshold) & (lv$lengths[intervening - 1] == 1))
  intervening <- intervening[validstarts]
  validends <- which((lv$values[intervening + 1] >= inc.threshold) & (lv$lengths[intervening + 1] == 1))
  intervening <- intervening[validends]
  
  lv$values[1:length(lv$values)] <- FALSE
  lv$values[c(intervening-1, intervening, intervening+1)] <- TRUE
  if (len > length(lv$values)) {
    lv$values[len] <- NA
    lv$lengths[len] <- lastl
  }
  inverse.rle(lv)
}

asri.dt <- raw.dt[,.SD[which.max(cases > 0):.N, .(date, asri = detect.asri(cases))], by=iso3]

#' if there is a short gap in asri, link up series
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

raw.dt[asri.dt, on=.(iso3, date), asri := merge.asri ]

matchfun <- function(r, dt, inc0, incf) inc0*sum(r^(1:dt)) - incf

#' @param series integer vector; 
#' @param asri integer vector; the asri series
impute.asri <- function(series, asri) {
  # find the first asri block
  alv <- rle(asri)
  block <- which.max(alv$values == 0)
  
  
  lv <- rle(series)
  dts <- lv$lengths[seq(2, length(lv$lengths), by = 2)]+1
  lv$values[1] <- initial
  res <- series
  wins <- cumsum(lv$lengths)+1
  for (step in 1:length(dts)) {
    inc0 <- lv$values[(step-1)*2+1]
    incf <- lv$values[(step-1)*2+3]
    dt <- dts[step]
    r <- uniroot(matchfun, c(0, 5), dt = dt, inc0 = inc0, incf = incf, tol = .Machine$double.eps^0.5)$root
    incn <- round(inc0*r^(1:dt))
    del <- incf - sum(incn)
    if (del < 0) incn[1:-del] <- incn[1:-del]-1
    if (del > 0) incn[dt:(dt-del+1)] <- incn[dt:(dt-del+1)]+1
    if (incn[dt] == 0) {
      warning(sprintf("mismatch on %i", step))
    } else {
      res[wins[(step-1)*2+1]:wins[(step-1)*2+2]] <- incn
      lv$values[(step-1)*2+3] <- incn[dt]
    }
  }
  res
}

raw.dt[!is.na(asri), i.cases := {
  if (sum(asri)) {
    
  } else cases
}, by=iso3]
