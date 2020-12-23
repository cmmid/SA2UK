
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
  require(ggplot2)
  require(ggrepel)
  require(patchwork)
})

.debug <- c("~/Dropbox/covidLMIC", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/projections/%s.qs",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/ecdc_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/mod_scenarios/%s.rds",
  .debug[2],
  sprintf("%s.png",.debug[2])
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

smoother <- function(dt, n = 7, align = "center") {
  setkeyv(rbind(
    dt[, smoothed := FALSE ],
    dt[order(date),
       .(value = frollmean(value, n, align), smoothed = TRUE),
      by = setdiff(key(dt),"date")
    ]
  ), key(dt))
}

lims.dt <- readRDS(.args[5])

scenkey <- scens <- rbind(
  data.table(work = NA, school = NA, other = NA),
  data.table(work = NA, school = NA, other = NA),
  data.table(
    expand.grid(
    #    home = c("none", "some"),
    work = c("small", "large"),
    school = c("small", "large"),
    other = c("small", "large"),
    #    symptrans = c("none","some"),
    # qtile = 1:length(preR0),
    stringsAsFactors = FALSE
  )
)[!(work == "small" & school == "small" & other == "small")]
)

scenkey[, scen_id := (1:.N)-1 ]
scenkey[scen_id == 0, label := "observed" ]
scenkey[scen_id == 1, label := "unmitigated" ]
scenkey[scen_id > 1, label := sprintf(
  "%s%s%s",
  fifelse(work=="large","W","w"),
  fifelse(school=="large","S","s"),
  fifelse(other=="large","O","o")
)]

dt <- qread(.args[1])[
  compartment %in% c("cases", "death_o", "R")
]

dt[compartment == "death_o", compartment := "deaths" ]
dt[compartment == "R", compartment := "exposed" ]
dt[, compartment := factor(compartment, levels = c("cases","deaths","exposed"))]

pars <- readRDS(.args[2])

pop <- pars$pop[[1]]
capita <- data.table(
  pop=c(pop$size, sum(pop$size))*urbfrac, group=c(pop$group_names, "all")
)

dt[, date := day0 + t ]

allage.dt <- rbind(dt[,
  .(value = sum(value)),
  keyby=.(scen_id, run, compartment, date)
], cases)

allage.dt[order(date), cvalue := cumsum(value), by=.(scen_id, run, compartment) ]
allage.dt[compartment == "exposed", cvalue := value ]

r0refs <- readRDS(.args[7])
scenrefs <- readRDS(.args[8])[scen_id > 1]

qraw <- c("lo.lo","lo","med","hi","hi.hi")

imms.mlt <- setkeyv(
  melt(r0refs[era %in% c("pre", "post")], id.vars = "era", measure.vars = qraw, variable.name = "run")[, run := as.integer(factor(run, levels = qraw, ordered = TRUE))],
  c("era", "run")
)

imms.mlt[, threshold := 1 - 1/value ]
imms.mlt[era == "pre", scen_id := 1 ]
imms.mlt[era == "post", scen_id := 2 ]

load("NGM.rda")

refR0 <- cm_ngm(pars)$R0

imms.shft <- imms.mlt[era == "pre"][, era := "modification" ][, {
  uf <- value / refR0
  scenrefs[is.infinite(end_day), {
    .(value = cm_ngm(pars, uf, contact_reductions = c(home, work, school, other), self_iso)$R0)
  },by=scen_id]
}, by=run][, threshold := 1 - 1/value ]

#' drop a line for intervention turn on
#' area for ignore window
#' area for end of fitting period

qs <- c('ll','lo','md','hi','hh')

collbls <- scenkey[,{ 
  ret <- label
  names(ret) <- scen_id
  ret
}]

colvals <- c("black","firebrick",rep("dodgerblue", 7))

crd <- function(xlim = as.Date(c("2020-03-01","2020-12-31")), ...) do.call(
  coord_cartesian, c(as.list(environment()), list(...))
)

eventop <- geom_blank(
  aes(date, value),
  data = function(dt) dt[,.SD[which.max(value), .(date, value = 10^ceiling(log10(value)))],by=compartment],
  inherit.aes = FALSE
)

p.inc <- ggplot(allage.dt[
  compartment %in% c("cases", "deaths")
][between(run, 2, 4)]) +
  facet_wrap(compartment ~ ., nrow = 2, scales = "free_y") +
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
    ), data = lims.dt, inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_text_repel(
    aes(x=date, y=Inf, label = era, angle = 90),
    data = lims.dt[era %in% c("pre","post","modification"), .(date = end, era)],
    inherit.aes = FALSE, hjust = 1
  ) +
  geom_line(data = function(dt) dt[scen_id != 1]) +
  geom_line(data = function(dt) dt[scen_id == 1]) +
  geom_line(data = function(dt) dt[
    scen_id == 0,.(date, value = frollmean(value, 7, align = "center")), keyby=.(scen_id, run, compartment)
  ], alpha = 1) +
#  geom_line(aes(color = "observed", linetype = "md", group = NULL), cases) +
  theme_minimal() +
  theme(
    panel.spacing.y = unit(1, "line")
  ) +
  scale_x_date(name = NULL, date_breaks = "months", date_labels = "%b", date_minor_breaks = "weeks") +
  scale_y_log10(
    expression("daily incidence ("*log[10]*" scale)"), breaks = function(lims) unique(c(1, scales::log10_trans()$breaks(lims))),
    labels = scales::label_number_si()
  ) +
  scale_color_manual(
    name=NULL,
    breaks = c(0, 1, 2),
    labels = c("reported","unmitigated\nscenario","calibrated\nintervention\nscenarios"),
    values = colvals,
    guide = "legend"
  ) +
  scale_linetype_manual(
    "quantile",
    breaks = c("md","lo","ll"),
    labels = c(ll="2.5-97.5%",lo="25-75%",md="50%"),
    values = c(md="solid", lo="dashed", hi="dashed", ll="dotted", hh="dotted")
  ) + 
  scale_fill_manual(
    guide = "none",
    values = c(censor = "grey", pre = "firebrick", transition = "grey", post = "dodgerblue", modification = "goldenrod")
  ) +
  crd(ylim=c(1,NA), expand = FALSE) + {
    nonrefscens <- allage.dt[,length(unique(scen_id))-2]
    scale_alpha_manual(
      guide = "none",
      values = c(0.5,1,rep((1/nonrefscens)^(2/3), nonrefscens))
    )
  }

#' @examples 
#' p <- p.inc + theme(legend.position = c(.95, .5), legend.justification = c(1,0.5), legend.background = element_rect(fill="white", color = "white"))
#' ggsave("incidence.png", p, width = 8, height = 5, units = "in")

plin <- p.inc + 
  scale_y_continuous(
    "daily incidence",
    labels = scales::label_number_si(),
    position = "right"
  )

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