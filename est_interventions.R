suppressPackageStartupMessages({
  require(data.table)
  require(lubridate)
})

.debug <- "~/Dropbox/covidLMIC"
.args <- if (interactive()) sprintf(c(
  "%s/outputs/ox_si_timing.csv",
  "%s/inputs/isos/all.iso",
  "timing_adjust.csv",
  "%s/outputs/interventions.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

tarisos <- readLines(.args[2])

ref <- fread(.args[1], select = c("iso3","date"))[iso3 %in% tarisos]

updates <- fread(tail(.args, 2)[1], header = FALSE, col.names = c("iso3","intervention","censor","deaths","modification"))
rmv <- updates[is.na(intervention), iso3]

res <- merge(ref[!(iso3 %in% rmv)], updates[!(iso3 %in% rmv)], by="iso3", all = TRUE)
res[is.na(intervention), c("intervention","censor","deaths") := .(0, 0, 0)]

# TODO: move over intervention start timing analysis
# dayeff: is the estimate of the day the intervention becomes effective

# ref <- melt(
#   fread(.args[1], drop = "CountryName"),
#   id.vars = "CountryCode"
# )[CountryCode != "", .(
#   si=value
# ), keyby=.(
#   iso3=CountryCode,
#   date=as.Date(variable,"%d%b%Y")
# )]
# 
# mxs <- ref[date < as.Date("2020-06-01"),
#   .SD[which.max(si)],
#   by=iso3
# ]
# 
# xs <- ref[
#   mxs[,.(iso3,si)], on=.(iso3)
# ][,.SD[which.max(si > 2*i.si/3)], by=iso3]

#' @examples 
#' require(ggplot2)
#' require(ggrepel)
#' filt <- expression(iso3 %in% c("NGA","RWA","ZAF","ZWE"))
#' filt <- expression(1:.N)
#' plotter <- function(ref, mxs, xs, filt) ggplot(
#'   ref[eval(filt)]
#' ) + aes(date, si, color=iso3) +
#'   geom_line() + theme_minimal() +
#'   geom_point(data=mxs[eval(filt)], color = "red") +
#'   geom_point(data=xs[eval(filt)], color = "blue") +
#'   scale_color_discrete(guide = "none") +
#'   scale_x_date(
#'     date_breaks = "month",
#'     date_minor_breaks = "week",
#'     date_labels = "%b"
#'   )

saveRDS(res, tail(.args, 1))