suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/cfrs.rds",
  "%s/outputs/figs/ccfr.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

dt <- readRDS(.args[1])
plot1.dt <- dt[ver == "cCFR" & date >= "2020-09-01" & province != "all"]
plot2.dt <- dt[ver == "cCFR" & date >= "2020-09-01" & province == "all"]

levels(plot1.dt$province) <- list(
                     "Eastern Cape" = "EC",  
                     "Free State" = "FS",
                     "Gauteng" = "GP",
                     "KwaZulu-Natal" = "KZN",
                     "Limpopo" = "LP",
                     "Mpumalanga" = "MP",
                     "Northern Cape" = "NC",
                     "North West" = "NW",
                     "Western Cape" = "WC",
                     "South Africa" = "all")



cfr.p1 <- ggplot(plot1.dt) + aes(date, md) +
    facet_wrap(~province, scale = "free", ncol = 3, nrow = 3) +  
    geom_line(aes(color = ver)) +
    geom_ribbon(aes(fill = ver, ymin = lo, ymax = hi), alpha = 0.2) +
    coord_cartesian(ylim = c(0, .075),
                    ##xlim = as.Date(c("2020-04-01", NA)),
                    expand = FALSE) +
    scale_y_continuous("Case Fatality Ratio (CFR)",
                       breaks = c(0,0.025,0.05,0.075),
                       labels = scales::percent) +
    scale_x_date(name = "", 
                 date_breaks = "months",
                 date_minor_breaks = "weeks",
                 date_labels = "%b") +
    scale_color_manual(NULL,
                       breaks = c("cCFR"),
                       labels = c(cCFR = ""),
                       values = c(cCFR = "darkorchid4"),
                       aesthetics = c("color", "fill")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") 
    #labs(tag = "")

cfr.p2 <- ggplot(plot2.dt) + aes(date, md) +
    geom_line(aes(color = ver)) +
    geom_ribbon(aes(fill = ver, ymin = lo, ymax = hi), alpha = 0.2) +
    coord_cartesian(ylim = c(0, .075),
                    ##xlim = as.Date(c("2020-04-01", NA)),
                    expand = FALSE) +
    scale_y_continuous("Case Fatality Ratio (CFR)",
                       breaks = c(0,0.025,0.05,0.075),
                       labels = scales::percent) +
    scale_x_date(name = "", 
                 date_breaks = "months",
                 date_minor_breaks = "weeks",
                 date_labels = "%b") +
    scale_color_manual(NULL,
                       breaks = c("cCFR"),
                       labels = c(cCFR = ""),
                       values = c(cCFR = "darkorchid4"),
                       aesthetics = c("color", "fill")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) + 
    labs(title = "South Africa")
 

library(patchwork)
layout <- "AAAAAA
           AAAAAA
           AAAAAA
           BBBBBB"

(cfr.p1 + cfr.p2 + plot_layout(design = layout))

saveRDS(cfr.p, tail(.args, 1))
