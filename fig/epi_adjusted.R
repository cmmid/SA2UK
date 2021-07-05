suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("analysis", "ZMB")
.args <- if (interactive()) sprintf(c(
  "%s/ins/adj_data.rds",
  .debug[2],
  "%s/fig/adjusted/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[1])[iso3 == tariso]

outcomes[is.na(asri), asri := 0]

roll.window <- 7

cats <- c("normal","asri")

p <- ggplot(outcomes) +
  aes(date, cases, color = factor(cats[asri+1], levels = cats), alpha = factor(cats[asri+1], levels = cats)) +
  geom_point(aes(y=adj, alpha = NULL), data = function(dt) dt[asri == 1], pch = 17, alpha = 1) +
  geom_point() +
  geom_line(data = function(dt) dt[, .(date, cases = frollmean(cases, roll.window, align = "center"), asri = 0)]) +
  ggtitle(tariso) +
  theme_minimal() + coord_cartesian(xlim = as.Date(c("2020-03-01", NA))) +
  scale_color_discrete(NULL, drop = F) +
  scale_alpha_manual(NULL, values = c(asri = 0.2, normal = 1), drop = F) +
  scale_y_log10() +
  scale_x_date(
    name = NULL,
    date_breaks = "month", date_minor_breaks = "week",
    date_labels = "%b %y"
  )

#ggsave(tail(.args, 1), p(mlt[outcome == "cases" & date >= "2020-03-01"], 1600), width = 16, height = 8, dpi = 300)
ggsave(tail(.args, 1),
  p,
  width = 16, height = 8, dpi = 300
)  
