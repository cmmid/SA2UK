
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
  require(ggplot2)
  require(ggrepel)
  require(patchwork)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/projections/%s.qs",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/epi_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/scenarios/%s.rds",
  .debug[2],
  "%s/outputs/projections/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <-tail(.args, 2)[1]

day0 <- readRDS(.args[3])[iso3 == tariso, date[1]]

urbfrac <- readRDS(.args[6])[iso3 == tariso, value/100 ]

cases <- melt(
  readRDS(.args[4])[iso3 == tariso],
  id.vars = c("date"),
  measure.vars = c("cases", "deaths"),
  variable.name = "compartment"
)
#cases[order(date), cvalue := cumsum(value), by=compartment ]
cases[, c("scen_id", "run") := .(0, 3) ]

lims.dt <- readRDS(.args[5])

dt <- qread(.args[1])[
  compartment %in% c("cases", "subclinical", "R")
]

dt[compartment == "R", compartment := "exposed" ]
dt[, compartment := factor(compartment, levels = c("cases","subclinical","exposed"))]

pars <- readRDS(.args[2])

pop <- pars$pop[[1]]
capita <- data.table(
  pop=urbfrac*c(
    # pop$size,
    # sum(pop$size),
    c(sum(pop$size[1:4]), sum(pop$size[5:8]), sum(pop$size[9:16]))
  ),
  age=c(
    #pop$group_names, "all",
    "youth", "20-40", ">40"
  )
)

dt[, date := day0 + t ]

allage.dt <- rbind(dt[,
  .(value = sum(value), scen_id = 2),
  keyby=.(run = r_id, compartment, date)
], cases)

byage.group.dt <- dt[compartment == "exposed",
  .(value = sum(value), scen_id = 2),
  keyby=.(
    run = r_id, compartment, date,
    age = fifelse(
      as.numeric(group) < 5, "youth",
      fifelse(as.numeric(group) > 8, ">40", "20-40")
    )
  )
]

allage.dt[order(date), cvalue := cumsum(value), by=.(scen_id, run, compartment) ]
allage.dt[compartment == "exposed", cvalue := value ]

r0refs <- readRDS(.args[7])
#' TODO mod scenario info?
scenrefs <- readRDS(.args[8])[scen_id > 1]

qraw <- c("lo.lo","lo","med","hi","hi.hi")

imms.mlt <- setkeyv(
  melt(r0refs[era %in% c("pre", "post")], id.vars = "era", measure.vars = qraw, variable.name = "run")[, run := as.integer(factor(run, levels = qraw, ordered = TRUE))],
  c("era", "run")
)

imms.mlt[, threshold := 1 - 1/value ]
imms.mlt[era == "pre", scen_id := 1 ]
imms.mlt[era == "post", scen_id := 2 ]

qs <- c('ll','lo','md','hi','hh')

# colvals <- c("black","firebrick",rep("dodgerblue", 7))

crd <- function(xlim = as.Date(c("2020-02-01","2021-01-31")), ...) do.call(
  coord_cartesian, c(as.list(environment()), list(...))
)

p.inc <- ggplot(allage.dt[
  compartment %in% c("cases") & (scen_id == 0 | (date < max(date)-60))
][between(run, 2, 4)]) +
#  facet_grid(compartment ~ ., scales = "free_y") +
  aes(
    date, value, color = factor(scen_id),
    linetype = qs[run], group = interaction(scen_id, run),
    alpha = factor(scen_id)
  ) +
  geom_rect(
    aes(
      ymin = 0.1, ymax = Inf,
      xmin=start-0.5, xmax=end+0.5,
      fill = era
    ), data = eras[!(era %in% c("censor","transition"))], inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line() +
  geom_line(data = function(dt) dt[
    scen_id == 0,.(date, value = frollmean(value, 7)), keyby=.(scen_id, run, compartment)
  ], alpha = 1) +
  theme_minimal() +
  theme(
    panel.spacing.y = unit(1, "line")
  ) +
  scale_x_date(name = NULL, date_breaks = "months", date_labels = "%b", date_minor_breaks = "weeks") +
  scale_y_log10(
    expression("daily incidence ("*log[10]*" scale)"), breaks = function(lims) unique(c(1, scales::log10_trans()$breaks(lims))),
    labels = scales::label_number_si()
  ) + scale_color_manual(
    name=NULL,
    breaks = c(0, 2),
    labels = c("reported","calibrated\ninterventions"),
    values = c("black","dodgerblue"),
    guide = "legend"
  ) +
  scale_linetype_manual(
    "quantile",
    breaks = c("md","lo","ll"),
    labels = c(ll="2.5-97.5%",lo="25-75%",md="50%"),
    values = c(md="solid", lo="dashed", hi="dashed", ll="dotted", hh="dotted")
  ) +
  scale_fill_manual(
    name = NULL,
    breaks=c("pre","post","modification","variant"),
    labels=c(pre="pre-intervention",post="post-intervention",modification="relaxation",variant="emergent variant"),
    values = c(pre="firebrick", post="dodgerblue",modification="goldenrod",variant="red")
  ) +
  crd(ylim=c(10,NA), expand = FALSE) +
  scale_alpha_discrete(guide = "none")
# +
#   theme(
#     legend.position = c(1,0), legend.justification = c(1,0)
#   )

c.attack <- byage.group.dt[
  compartment == "exposed"
][order(date)][capita, on=.(age)][, .(date, value = value/pop), keyby=.(scen_id, run, age)]

en <- c.attack[, max(date)]

att.p <- ggplot(c.attack[between(run, 2, 4)]) + aes(
  date, value, color = age,
  linetype = qs[run], group = interaction(age, run)
) + geom_line() + 
  theme_minimal() +
  crd(ylim = c(0, 0.5), xlim = as.Date(c("2020-06-01","2020-10-01")), expand = FALSE) +
  scale_x_date(name = NULL, date_breaks = "months", date_labels = "%b", date_minor_breaks = "weeks") +
  scale_y_continuous(
    "cumulative attack proportion\nof urban population"
  ) +
  scale_linetype_manual(
    "quantile",
    breaks = c("md","lo","ll"),
    labels = c(ll="2.5-97.5%",lo="25-75%",md="50%"),
    values = c(md="solid", lo="dashed", hi="dashed", ll="dotted", hh="dotted"),
    guide = "none"
  )

res <- ((((p.inc + eventop + theme(strip.text = element_blank())) + plin)) / 
  (rep.p + att.p)) + plot_layout(guides = "collect")

ggsave(tail(.args, 1), res, width = 10, height = 6, units = "in")

# exp.dt <- dcast(allage.dt[
#   compartment == "exposed", .(
#     scen_id, date, run, value = value / capita[group=="all", pop*urbfrac]
#   )
# ], scen_id + date ~ qs[run])
# 
# p.pre <- ggplot(exp.dt) +
#   aes(date, color = factor(scen_id), fill = factor(scen_id)) +
#   geom_ribbon(aes(ymin=ll, ymax=hh, color = NULL), alpha = 0.1) +
#   geom_ribbon(aes(ymin=lo, ymax=hi, color = NULL), alpha = 0.2) +
#   geom_line(aes(y=md)) +
#   theme_minimal() +
#   crd(ylim = c(0, 1), expand = FALSE) +
#   scale_x_date(NULL, date_breaks = "months", date_labels = "%b") +
#   scale_y_continuous("Prevalence") +
#   scale_color_discrete(guide = "legend", aesthetics = c("color", "fill"))
#   
# ggsave(tail(.args, 1), p, width = 9, height = 6, units = "in")