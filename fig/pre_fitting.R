suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.debug <- c("analysis", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/ins/adj_data.rds",
  "%s/gen/intervention_timing/%s.rds",
  "%s/gen/pops/%s.rds",
  "%s/gen/yuqs/%s.rds",
  "%s/gen/mobility.rds",
  "%s/est/r0/%s.rds",
  .debug[2],
  "%s/figs/prefitting/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[1])[iso3 == tariso]
eras <- readRDS(.args[2])[period == 1 & era != "transition"]
#' for comparing r0 estimates, want a reference case series
#' essentially want to draw a line (in log space) through that series
#' to do that, find the appropriate anchor point
#' == geometric mean, center of the estimation period?
#' then draw with slope R/gi through that point
ref.infs <- outcomes[
  eras, on=.(iso3), allow.cartesian = TRUE
][
  between(date, start-4, end+4),
  .(date, adj, era)
][-(1:(which.max(adj > 0)-1))]
wind <- ref.infs[order(date), {
  enc <- rle(adj)
  max(c(enc$lengths[enc$values == 0],0))
}, by=era][, max(V1)]+1
ref.infs[order(date), rm := frollmean(adj, wind, align = "left"), by=era]
for (i in 1:(wind-1)) {
  ref.infs[c(which(adj==0),which(adj==0)+1), adj := rep(rm[which(adj==0)],2)]
}
ref.infs[, gm := exp(frollmean(log(adj), 7, align = "center")), by=era ]
anchors <- ref.infs[
  !is.na(gm),
  .(date=date[.N/2], gm=prod(gm)^(1/.N), extent = .N), by=era
][,
  .(date, gm, mod=(1:extent)-floor(1+extent/2)), by=era
][, .(date=date+mod, gm, mod), by=era ]

pop <- readRDS(.args[3])

p_of_having_an_interval <- 1-cumsum(c(0, pop$pop[[1]]$dIa))
p_infection_in_interval <- head(p_of_having_an_interval/sum(p_of_having_an_interval),-1)

mean_gi <- sum(seq(0,by=pop$time_step,length.out = length(p_infection_in_interval))*p_infection_in_interval) +
  sum(seq(0,by=pop$time_step,length.out = length(pop$pop[[1]]$dE))*pop$pop[[1]]$dE)

mobility <- readRDS(.args[5])[tariso == iso3]

r0 <- readRDS(.args[6])[(period == 1 & era != "transition") | variable == "infections"]
infections.dt <- r0[variable == "infections"]
interventionrs <- dcast(r0[period == 1][, r := value^(1/mean_gi) ], sample + era ~ ., value.var = "r")
rplot.dt <- interventionrs[anchors, on=.(era), allow.cartesian = TRUE][,
  .(sample, era, date, value = gm*`.`^mod)
][order(sample, era, date)]

p <- ggplot(outcomes) +
  aes(date) +
  geom_line(
    aes(y=adj, color="report")
  ) +
  geom_line(
    aes(y = value, color="Rt", group = sample),
    data = infections.dt, alpha = 0.01
  ) +
  geom_line(
    aes(y = value, color=era, group = interaction(era, sample)),
    data = rplot.dt, alpha = 0.01
  ) +
  scale_color_manual(
    name="Case Measures",
    breaks = c("report", "pre", "post", "Rt"),
    labels = c(report="reported (adj.)", pre="pre-intervention R0", post="post-intervention R0", Rt="Rt-based infections"),
    values = c(report="black", pre="firebrick", post="dodgerblue", Rt="goldenrod")
  ) +
  scale_y_log10("Cases (log scale)")

mob.plot <- melt(mobility, id.vars = "date", measure.vars = c("school_multiplier","workr","otherr"))

p.mob <- ggplot(mob.plot) + aes(date, value, color=variable) +
  geom_line() +
  scale_color_discrete(
    "Contact Multipliers",
    labels = c(school_multiplier="Schools", workr="Work", otherr="Non-Home Other")
  )

merge.p <- (p / p.mob) + plot_layout(guides = "collect") & theme_minimal() &
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "left"
  ) &
  scale_x_date(NULL, date_breaks = "months", date_labels = "%b %y") &
  coord_cartesian(xlim = as.Date(c("2020-02-01","2021-07-01")), expand = FALSE)

ggsave(tail(.args, 1), merge.p, height = 6.5, width = 15, units = "in", dpi = 900)
