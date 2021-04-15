
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
#' propose: duration threshold >= 3 days (exceeds weekend)

#' @param series, integer vector which begins with a non-zero value
#' @param inc.threshold, integer > 0, the incidence threshold for ends of the ASRI
#' @param dur.threshold, integer > 0, the duration threshold for ASRI
detect.asri <- function(series, inc.threshold = 10, dur.threshold = 3) {
  lv <- rle(series)
  len <- length(lv$values)
  if (lv$values[len] == 0) {
    lv$values <- lv$values[-len]
    lastl <- lv$lengths[len]
    lv$lengths <- lv$lengths[-len]
  }
  
  intervening <- which(lv$values == 0 & lv$lengths > dur.threshold)
  validstarts <- which((lv$values[intervening - 1] > inc.threshold) & (lv$lengths[intervening - 1] == 1))
  intervening <- intervening[validstarts]
  validends <- which((lv$values[intervening + 1] > inc.threshold) & (lv$lengths[intervening + 1] == 1))
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

raw.dt[asri.dt, on=.(iso3, date), asri := asri ]
raw.dt[, i.cases := cases ]
raw.dt[asri == 1 & i.cases == 0, i.cases := NA_integer_ ]

raw.dt