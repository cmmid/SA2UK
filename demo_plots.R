
require(data.table)
require(ggplot2)

figrt <- file.path("~", "Dropbox", "Covid-WHO-vax", "nigeria")

scen.dt <- rbind(
  data.table(vax_mech = "none", vax_eff = 0, coverage = 0),
  data.table(expand.grid(
    vax_mech = c("infection", "disease"),
    vax_eff = c(.5, .9),
    coverage = c(.25, .90)
  ))
)[, id := (1:.N)-1 ]

episrc <- file.path(figrt, "outputs", "acc_scen")
fig.path <- function(fn) file.path(figrt, "fig", fn)

epi.dt <- rbindlist(lapply(list.files(episrc, full.names = TRUE), readRDS))
epi.ref <- epi.dt[epi_id == 0, .SD, .SDcols = -c('epi_id') ]
epi.int <- epi.dt[epi_id != 0][epi.ref, on=.(sample, iyear, compartment, group)]
epi.int[, averted := i.value - value ]

sx <- function(name="Years since start of vaccination", ...) scale_x_continuous(name=name, ...)

thm <- theme_minimal(base_size = 14)

epi.p <- ggplot(
  epi.dt[ iyear > 0 & compartment == "cases", .(value = sum(value)), by=.(epi_id, sample, iyear, compartment) ][scen.dt, on=.(epi_id = id)]
) + aes(
  iyear, value, color=vax_mech,
  group = interaction(vax_mech, sample)
) + facet_grid(vax_eff ~ coverage, labeller = labeller(
  vax_eff = function(n) sprintf("Efficacy = %i%%", as.integer(as.numeric(n)*100)),
  coverage = function(n) sprintf("Coverage = %i%%", as.integer(as.numeric(n)*100))
)) +
  geom_line(alpha = 0.2) +
  thm +
  sx() +
  scale_y_continuous("Cases", labels = scales::label_number_si()) +
  scale_color_manual("Vaccine\nMechanism",
    values = c(none="black", disease="green", infection="blue"),
    guide = guide_legend(override.aes = list(alpha=1))
  )

ggsave(fig.path("cases_slide.png"), epi.p, width = 15, height = 8, units = "in")

epi.del <- ggplot(
  epi.int[ iyear >= 1 & compartment == "cases", .(value = sum(averted)), by=.(epi_id, sample, iyear, compartment) ][scen.dt[id!=0], on=.(epi_id = id)]
) + aes(
  iyear, value, color=vax_mech,
  group = interaction(vax_mech, sample)
) + facet_grid(vax_eff ~ coverage, labeller = labeller(
  vax_eff = function(n) sprintf("Efficacy = %i%%", as.integer(as.numeric(n)*100)),
  coverage = function(n) sprintf("Coverage = %i%%", as.integer(as.numeric(n)*100))
)) +
  geom_line(alpha = 0.2) +
  thm +
  sx() +
  scale_y_continuous("Cases Averted In Year X", labels = scales::label_number_si()) +
  scale_color_manual("Vaccine\nMechanism",
                     values = c(disease="green", infection="blue"),
                     drop = TRUE,
                     guide = guide_legend(override.aes = list(alpha=1))
  )

ggsave(fig.path("cases_averted.png"), epi.del, width = 15, height = 8, units = "in")

epi.cdel <- ggplot(
  epi.int[ iyear >= 1 & compartment == "cases", .(
    value = sum(averted)
  ), by=.(epi_id, sample, iyear, compartment) ][,
    .(iyear, value = cumsum(value)),
    by=.(epi_id, sample, compartment)
  ][scen.dt[id!=0], on=.(epi_id = id)]
) + aes(
  iyear, value, color=vax_mech,
  group = interaction(vax_mech, sample)
) + facet_grid(vax_eff ~ coverage, labeller = labeller(
  vax_eff = function(n) sprintf("Efficacy = %i%%", as.integer(as.numeric(n)*100)),
  coverage = function(n) sprintf("Coverage = %i%%", as.integer(as.numeric(n)*100))
)) +
  geom_line(alpha = 0.2) +
  thm +
  sx() +
  scale_y_continuous("Cumulative Cases Averted", labels = scales::label_number_si()) +
  scale_color_manual("Vaccine\nMechanism",
    values = c(disease="green", infection="blue"),
    drop = TRUE,
    guide = guide_legend(override.aes = list(alpha=1))
  )

ggsave(fig.path("ccases_averted.png"), epi.cdel, width = 15, height = 8, units = "in")

src <- file.path(figrt, "outputs", "econ_scen", "NGA.rds")

dt <- readRDS(src)[view == "incremental" & econ_id == 1]

r.dt <- dcast(
  melt(dt,
       id.vars = c("id", "econ_id", "anni_year", "qtile"),
       measure.vars = c("costs", "dalys", "ccosts", "cdalys", "icer")
  ), 
  id + econ_id + anni_year + variable ~ qtile 
)[scen.dt, on=.(id)]

fulle.dt <- dcast(
  melt(readRDS(src)[view == "incremental"],
       id.vars = c("id", "econ_id", "anni_year", "qtile"),
       measure.vars = c("costs", "dalys", "ccosts", "cdalys", "icer")
  ), 
  id + econ_id + anni_year + variable ~ qtile 
)[scen.dt[id != 0], on=.(id)]

# note: only health_system econ perspective

# build up slide show from slides
# fill notes to slide
ccost.p <- ggplot(r.dt[variable == "ccosts"]) + aes(anni_year) +
  geom_ribbon(aes(ymin=lo95, ymax=hi95, fill=vax_mech), alpha = 0.1) +
  geom_ribbon(aes(ymin=lo50, ymax=hi50, fill=vax_mech), alpha = 0.1) +
  geom_line(aes(y=md, color=vax_mech)) +
  facet_grid(vax_eff ~ coverage, labeller = labeller(
    vax_eff = function(n) sprintf("Efficacy = %i%%", as.integer(as.numeric(n)*100)),
    coverage = function(n) sprintf("Coverage = %i%%", as.integer(as.numeric(n)*100))
  )) +
  thm +
  scale_y_continuous("Cumulative Incremental Cost", labels = scales::label_number_si()) +
  sx() +
  scale_color_manual("Vaccine\nMechanism",
                     values = c(disease="green", infection="blue"),
                     aesthetics = c("color","fill")
  )

ggsave(fig.path("ccost_incremental.png"), ccost.p, width = 15, height = 8, units = "in")

cdaly.p <- ggplot(r.dt[variable == "cdalys"]) + aes(anni_year) +
  geom_ribbon(aes(ymin=lo95, ymax=hi95, fill=vax_mech), alpha = 0.1) +
  geom_ribbon(aes(ymin=lo50, ymax=hi50, fill=vax_mech), alpha = 0.1) +
  geom_line(aes(y=md, color=vax_mech)) +
  facet_grid(vax_eff ~ coverage, labeller = labeller(
    vax_eff = function(n) sprintf("Efficacy = %i%%", as.integer(as.numeric(n)*100)),
    coverage = function(n) sprintf("Coverage = %i%%", as.integer(as.numeric(n)*100))
  )) +
  thm +
  scale_y_continuous("Cumulative DALYs Averted", labels = scales::label_number_si()) +
  sx() +
  scale_color_manual("Vaccine\nMechanism",
                     values = c(disease="green", infection="blue"),
                     aesthetics = c("color","fill")
  )

ggsave(fig.path("cdalys_averted.png"), cdaly.p, width = 15, height = 8, units = "in")

icer.tab <- fulle.dt[variable == "icer" & anni_year == 5, .(
  sprintf("%.1f (50Q %.1f - %.1f, 95Q %.1f - %.1f; mean %.1f)", md, lo50, hi50, lo95, hi95, mn)
), by=.(vax_mech, vax_eff, coverage, econ_id)]

fwrite(icer.tab, fig.path("icers.csv"))
