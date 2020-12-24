suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv",
  "%s/inputs/epi_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "fig1.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

rsa.dt <- fread(.args[1])[Country == "South Africa"][order(`Collection Data`)]

plot.dt <- rsa.dt[,
                  .(.N, iso3c = "ZAF"),
                  by = .(date=round(`Collection Data`, digits = "weeks"), newvariant = `Clade`=="20C")
]

plot.dt[, total := sum(N), by=.(date, iso3c) ]

plot.dt[,
        c("lo95","hi95") := 
          as.data.table(t(mapply(
            function(x, n, p=x/n) binom.test(x, n, p, conf.level = .95)$conf.int,
            x = N, n = total
          )))
][,
  c("lo50","hi50") := 
    as.data.table(t(mapply(
      function(x, n, p=x/n) binom.test(x, n, p, conf.level = .50)$conf.int,
      x = N, n = total
    )))
]

p.phylo <- ggplot(plot.dt[newvariant == TRUE]) + aes(date) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.1, fill = "red") +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.2, fill = "red") +
  geom_line(aes(y=N/total), color = "red") +
  scale_x_date(
    name = NULL,
    date_breaks = "months", date_minor_breaks = "weeks",
    date_labels = "%b"
  ) +
  scale_y_continuous("501Y.V2 Fraction") +
  coord_cartesian(ylim = c(0, 1), xlim = c(as.Date("2020-10-01")-5, as.Date("2020-12-01")+5), expand = FALSE) +
  theme_minimal()

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[2])[iso3 == tariso]

eras <- readRDS(.args[3])

p <- ggplot(outcomes) +
  aes(date) +
  geom_line(aes(y=cases, linetype="cases")) +
  geom_line(aes(y=deaths, linetype="deaths")) +
  geom_rect(
    aes(
      ymin = 0.1, ymax = Inf,
      xmin=start-0.5, xmax=end+0.5,
      fill = era
    ), data = eras[!(era %in% c("censor","transition"))], inherit.aes = FALSE,
    alpha = 0.2
  ) +
  theme_minimal() +
  scale_y_log10("incidence", breaks = 10^(0:5), labels = scales::label_number_si()) +
  scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
  scale_fill_manual(
    name = NULL,
    breaks=c("pre","post","modification","variant"),
    labels=c(pre="pre-intervention",post="post-intervention",modification="relaxation",variant="emergent variant"),
    values = c(pre="firebrick", post="dodgerblue",modification="goldenrod",variant="red")
  ) +
  scale_linetype_discrete(name=NULL) +
  coord_cartesian(ylim = c(1, NA), xlim = as.Date(c("2020-02-01","2021-01-01")), expand = FALSE)


p.inc <- ggplot(allage.dt[
  compartment %in% c("cases") & (scen_id == 0 | (date < max(date)-60))
][between(run, 2, 4)]) +
  #  facet_grid(compartment ~ ., scales = "free_y") +
  aes(
    date, value, color = factor(scen_id),
    linetype = qs[run], group = interaction(scen_id, run),
    alpha = factor(scen_id)
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
  crd(ylim=c(10,NA), expand = FALSE) +
  scale_alpha_discrete(guide = "none") +
  theme(
    legend.position = c(1,0), legend.justification = c(1,0)
  )

c.attack <- byage.group.dt[
  compartment == "exposed"
][order(date)][capita, on=.(age)][, .(date, value = value/pop), keyby=.(scen_id, run, age)]

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
  ) + theme(
    legend.position = c(0,1), legend.justification = c(0, 1)
  )

p.inc / (att.p + p.phylo) & theme(text = element_text(size = 8))
