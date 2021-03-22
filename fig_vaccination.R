
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/params/%s_0001.rds",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "../covidm",
  "%s/outputs/projections/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 3)[1]

fits <- readRDS(.args[1])
intros.dt <- readRDS(.args[3])[iso3 == tariso]
urbfrac <- readRDS(.args[4])[iso3 == tariso, value / 100]
timings <- readRDS(.args[5])

day0 <- as.Date(intros.dt[, min(date)])
intros <- intros.dt[,
  intro.day := as.integer(date - date[1])
][, .(t=Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))), by=sample ]

popsetup <- function(basep, urbanfraction, day0) {
  basep$date0 <- day0
  basep$pop[[1]]$size <- round(basep$pop[[1]]$size*urbanfraction)
  basep$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))
  basep
}

base <- popsetup(readRDS(.args[2]), urbfrac, day0)

startrelax <- as.integer(timings[era == "relaxation", start] - day0)
endrelax <- as.integer(timings[era == "relaxation", end] - day0)
endsim <- as.integer(timings[era == "pre" & period == 3, start] - day0)

startpost <- as.integer(timings[era == "transition", start[1]] - day0)

base$time1 <- endsim

tms <- day0 + startpost
relaxtms <- day0 + startrelax:endrelax

tier2 <- as.Date("2020-08-15")

load("NGM.rda")

res1 <- readRDS("example_0001.rds")
res2 <- readRDS("example_0002.rds")

averted <- res1[res2, on=.(sample, t, group, compartment), .(sample, t, group, compartment, del = value - i.value)]

pl2.dt <- averted[, .(del = sum(del)), by=.(sample, date = t + day0, compartment)][compartment %in% c("cases","death_o")]
pl2.dt[order(date), cdel := cumsum(del), by=.(sample, compartment)]

qp2 <- pl2.dt[between(date, "2021-03-01","2021-09-01"),{
  qs <- quantile(cdel, probs = c(0.25, 0.5, 0.75))
  names(qs) <- c("lo","md","hi")
  as.list(qs)
}, by=.(date, compartment)
]

p2 <- ggplot(qp2) +
  aes(date) +
  facet_wrap(compartment ~ ., nrow=2, scale = "free_y", labeller = labeller(compartment = c(cases="cases",death_o="deaths"))) +
  geom_ribbon(aes(ymax=hi, ymin=lo), alpha = 0.2) +
  geom_line(aes(y=md)) +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  scale_x_date(NULL,
    date_breaks = "3 months", date_minor_breaks = "month", date_labels = "%b %Y"
  ) +
  scale_y_continuous("Cumulative Averted Outcomes")

ggsave("thing2.png", p2, width = 8, height = 8, dpi = 300)

res <- rbind(res1, res2)
comparison <- res[compartment == "cases", .(value = sum(value)), by=.(sample, date, intervention_id)]

pl1.dt <- comparison[fits[, .(sample, asc)], on=.(sample), .(sample, date, intervention_id, asc.value = value*asc, value)]

pl <- pl1.dt[,.(date, v = value), by=.(sample, intervention_id)][,{
  qs <- quantile(v, probs = c(.25, 0.5, .75))
  names(qs) <- c("lo","md","hi")
  as.list(qs)
}, by=.(date, intervention_id)]

psim <- ggplot(pl1.dt[between(date,"2020-03-01","2021-04-01") & intervention_id == 1]) +
  aes(date) +
#  geom_ribbon(aes(ymax=hi, ymin=lo, fill=c("none","vax")[intervention_id]), alpha = 0.2) +
  geom_line(aes(y=value, color=c("none","vax")[intervention_id], group = sample)) +
  geom_point(
    aes(date, cases),
    data = readRDS("~/Dropbox/Covid_LMIC/All_Africa_paper/inputs/epi_data.rds")[iso3 == .debug[2]][date >= "2020-03-01"],
    color = "black", inherit.aes = FALSE, alpha = 0.2
  ) +
  geom_line(
    aes(date, cases),
    data = readRDS("~/Dropbox/Covid_LMIC/All_Africa_paper/inputs/epi_data.rds")[iso3 == .debug[2]][date >= "2020-03-01",.(date, cases = frollmean(cases, 7))],
    color = "black", inherit.aes = FALSE
  ) +
  coord_cartesian(
    xlim = as.Date(c("2020-03-01", "2021-04-01")),
    expand = F, clip = "off"
  ) +
  scale_x_date(NULL,
               date_breaks = "3 months", date_minor_breaks = "month", date_labels = "%b %Y"
  ) +
  scale_y_log10("Case Incidence") +
  scale_color_manual(
    name=NULL,
    guide = "none", #guide_legend(override.aes = list(alpha = 1)),
    values = c("firebrick","dodgerblue"), labels = c("No Vaccination", "90% VE, 4k doses per day")) +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1))

ggsave("thing.png", psim, width = 16, height = 8, dpi = 300)

#' @examples 

saveRDS(res, tail(.args, 1))

