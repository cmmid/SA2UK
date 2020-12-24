suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- "~/Dropbox/SA2UK"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv",
  "comparison.png"
), .debug) else commandArgs(trailingOnly = TRUE)

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

p <- ggplot(plot.dt[newvariant == TRUE]) + aes(date) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.1) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.2) +
  geom_line(aes(y=N/total)) +
  scale_x_date(
    name = NULL,
    date_breaks = "months", date_minor_breaks = "weeks",
    date_labels = "%b"
  ) +
  scale_y_continuous("501Y.V2 Fraction") +
  coord_cartesian(ylim = c(0, 1), xlim = c(as.Date("2020-10-01")-5, as.Date("2020-12-01")+5), expand = FALSE) +
  theme_minimal()

saveRDS(p, gsub("png$","rds", tail(.args, 1)))

ggsave(tail(.args, 1), p, width = 6, height = 3, units = "in", dpi = 300)