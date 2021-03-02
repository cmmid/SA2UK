suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "~/Dropbox/SA2UK"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

target <- tail(.args, 1)
jhurl <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
casesurl <- sprintf("%s/time_series_covid19_confirmed_global.csv", jhurl)
deathsurl <- sprintf("%s/time_series_covid19_deaths_global.csv", jhurl)

fetch <- function(url, vn) melt(fread(url)[
  `Country/Region` %in% c("South Africa", "United Kingdom") &
    `Province/State` == ""
][, -c(1,3,4) ], id.vars = "Country/Region", variable.name = "date", value.name = vn)

#' fetch ECDC data; requires network connection
cases.dt <- fetch(casesurl, "cases")
deaths.dt <- fetch(deathsurl, "deaths")

res <- cases.dt[deaths.dt, on=.(`Country/Region`, date)]
res[, date := as.Date(date, format = "%m/%d/%y") ]

#' select the columns of interest; order by key columns
final <- res[order(date),
  .(date, cases = c(cases[1], diff(cases)), deaths = c(deaths[1], diff(deaths))),
  keyby=.(
    continent = fifelse(`Country/Region`=="South Africa","Africa","Europe"),
    iso3 = fifelse(`Country/Region`=="South Africa","ZAF","GBR")
  )
]

#' @examples 
#' ggplot(final) + aes(date) + facet_grid(iso3 ~ .) +
#'   geom_line(aes(y=cases)) + theme_minimal() + scale_y_log10() +
#'   scale_x_date(name=NULL, date_breaks = "month", date_labels = "%b")

#' TODO capture stderr in makefile?
#' alert functions for next steps
neg.warn <- function(
  isos,
  msg = "the following countries have negative cases:\n%s"
) warning(sprintf(
  msg,
  paste(isos, collapse = "\n")
))

smooth.warn <- function(width) warning(
  sprintf(
    "attempting to eliminate negative case counts by averaging surrounding +/- %i days.",
    width
  )
)

#' if there are any countries with negative case counts,
#' we assume those are reporting corrections to be smoothed out into
#' earlier cases
if (final[, any(cases < 0)]) {
  fix.isos <- final[, unique(iso3[cases < 0])]
  unrepairable <- final[iso3 %in% fix.isos, any(cumsum(cases) < 0), by=iso3][V1 == TRUE, iso3]
  final <- final[!(iso3 %in% unrepairable)]
  neg.warn(fix.isos)
  if (length(unrepairable) == length(fix.isos)) {
    warning("and all have total cases < 0 at some point; all will be excluded.")
  } else {
    if (length(unrepairable)) {
      neg.warn(unrepairable, "and the following have total cases < 0 at some point and will be excluded:\n%s")
    }
    fix.isos <- setdiff(fix.isos, unrepairable)
    repair <- final[iso3 %in% fix.isos][, corrected := cases ]
    while(repair[, any(corrected < 0)]) {
      repair[iso3 %in% fix.isos, corrected := if(any(corrected < 0)) {
        ind <- which.max(corrected < 0)
        del <- corrected[ind]
        prop <- corrected[1:(ind-1)]/sum(corrected[1:(ind-1)])
        #' slightly emphasize large case counts
        relprop <- prop^1.1/sum(prop^1.1)
        #' TODO: rounding slightly changes total case count;
        #' possible to (in|de)flate del to get exact match?
        c(round(relprop*del)+corrected[1:(ind-1)], 0, corrected[(ind+1):.N])
      } else corrected, by=iso3 ]
    }
    final <- setkeyv(rbind(
      final[!(iso3 %in% fix.isos)],
      repair[, .(continent, iso3, date, cases = corrected, deaths)]
    ), key(final))
  }
}

#' @examples 
#' require(ggplot2); require(ggrepel)
#' repairp <- ggplot(repair) + aes(date) + 
#'   facet_grid(iso3 ~ ., scales = "free_y") +
#'   geom_line(aes(y=cases, alpha="original")) + 
#'   geom_line(aes(y=corrected, alpha="repaired")) +
#'   geom_label_repel(
#'     aes(y=10, label=cases),
#'     data = function(dt) dt[cases < 0],
#'     color = "red",
#'     fill = alpha(c("white"), 0.1),
#'     size = 2, fontface = "bold",
#'     show.legend = FALSE,
#'     direction = "y", max.iter = 1e4,
#'     segment.alpha = 0.2,
#'     label.size = NA, label.padding = .05
#'   ) +
#'   coord_cartesian(ylim = c(1, NA), expand = FALSE) +
#'   scale_alpha_manual("Series", values = c(original = 0.4, repaired = 1)) +
#'   scale_x_date(NULL, date_breaks = "months", date_minor_breaks = "weeks", date_labels = "%b") +
#'   scale_y_log10(minor_breaks = NULL, labels = scales::label_number_si()) +
#'   theme_minimal() +
#'   theme(legend.position = c(0, 1), legend.justification = c(0, 1))
#' ggsave("ecdc_repair.png", repairp, width = 6, height = 8, units = "in", dpi = 600)

saveRDS(final, target)
