
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/covidLMIC")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/r0", "summary_r0.png"
), .debug) else commandArgs(trailingOnly = TRUE)

r0s <- list.files(.args[1])

extractR0 <- function(fn, pth = .args[1]) {
  readRDS(file.path(pth, fn))[era %in% c("pre","post")][, iso3 := gsub("^.*(.{3})\\.rds$", "\\1", fn)]
}

allR0s <- rbindlist(lapply(r0s, extractR0))
allR0s[, era := factor(era, levels = c("pre","post"), ordered = TRUE)]

p <- ggplot(allR0s) + aes(
  date, med,
  shape = era, color = era, fill = era
) +
  facet_grid(. ~ era) +
  geom_point() +
  theme_minimal() +
  scale_x_date(
    date_breaks = "months", date_minor_breaks = "weeks",
    date_labels = "%b"
  ) +
  scale_shape_discrete(
    palette = function(i) c(pre=24, post=25),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(pre="firebrick", post="dodgerblue"),
    aesthetics = c("color","fill"),
    guide = "none"
  )
