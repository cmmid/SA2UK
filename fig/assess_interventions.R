suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/adj_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/phylo.rds",
  .debug[2],
  "%s/outputs/intervention_timing/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[1])[iso3 == tariso]

eras <- readRDS(.args[2])

phylo.frac <- readRDS(.args[3])

#' for visualization convenience, approximating this phylo frac
#' change as universal
#' so transition last date to begining of period 3

newend <- eras[period == 3 & era == "pre", start ]
phylo.frac[, date := newend + -(.N-1):0 ]

outcomes[, var.frac := 0 ]
outcomes[phylo.frac[!is.na(binop)], on=.(date), var.frac := binop ]
outcomes[!between(date, phylo.frac[, min(date)], phylo.frac[, max(date)]), var.frac := NA ]

death.delay <- 21

outcomes[, del.var.frac := 0 ]
outcomes[phylo.frac[!is.na(binop), .(date = date + 21, binop)], on=.(date), del.var.frac := binop ]
outcomes[!between(date, phylo.frac[, min(date)+21], phylo.frac[, max(date)+21]), del.var.frac := NA ]

roll.window <- 7

#' set to adjusted cases
outcomes[, cases := adj ]

outcomes[order(date), rollcases := frollmean(cases, roll.window), by=iso3 ]
outcomes[order(date), rolldeaths := frollmean(deaths, roll.window), by=iso3 ]
outcomes[, rollvarcases := rollcases * var.frac ]
outcomes[, rollvardeaths := rolldeaths * del.var.frac ]
outcomes$cases <- as.numeric(outcomes$cases)
outcomes$deaths <- as.numeric(outcomes$deaths)

mlt <- melt(
  outcomes, id.vars = "date", measure.vars = c("cases", "deaths", "rollcases", "rolldeaths", "rollvarcases", "rollvardeaths")
)[!is.na(value)]

mlt[, outcome := fifelse(grepl("case", variable), "cases", "deaths")]
mlt[, measure := fifelse(grepl("roll", variable), "rolling", "raw")]
mlt[, variant := fifelse(grepl("var", variable), "variant", "all")]

p <- force(ggplot(mlt[measure == "raw" | (value > 0.1)]) +
             aes(date, y= value, color = variant, alpha = measure, linetype = outcome) +
             geom_line() +
             geom_rect(
               aes(
                 ymin = 0, ymax = Inf,
                 xmin=start-0.5, xmax=end+0.5,
                 fill = era
               ), data = eras[!(era %in% c("censor","transition"))], inherit.aes = FALSE,
               alpha = 0.2
             ) +
             theme_minimal(base_size = 14) +
             scale_y_log10(sprintf("Incidence", roll.window), breaks = 10^(0:5), labels = scales::label_number_si()) +
             scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
             scale_color_manual(name = NULL, labels = c(all = "all", variant = "est. variant"), values = c(all="black", variant = "red")) +
             scale_linetype_manual(name = NULL, values = c(cases="solid", deaths = "longdash")) +
             scale_alpha_manual(name = NULL, values = c(raw=0.5, rolling = 1)) +
             scale_fill_manual(
               name = NULL,
               breaks=c("pre","post","relaxation","variant"),
               labels=c(pre="pre-intervention",post="post-intervention",relaxation="relaxation",variant="emergent variant"),
               values = c(pre="firebrick", post="dodgerblue",relaxation="goldenrod",variant="red")
             ) +
             coord_cartesian(ylim = c(1, NA), xlim = as.Date(c("2020-03-01", NA)), expand = FALSE))

ggsave(tail(.args, 1), p, height = 6.5, width = 15, units = "in", dpi = 900)
