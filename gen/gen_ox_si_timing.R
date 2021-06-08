suppressPackageStartupMessages({
  require(data.table)
  require(TTR)
})

.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/ox_si_timing.csv"
), .debug) else commandArgs(trailingOnly = TRUE)

oxsiurlpat <- "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/timeseries/%s"
#' going to use all the metrics that are ordinal & could contribute to local control
#' so exclude international travel (c8) and various E / H measures in USD
#' also exclude vaccine policy (h7), as not relevant for historical trends
#' since predicting at national level, only going to count policies when flag == 1
fls <- sprintf(oxsiurlpat, c(
  "c1_flag.csv",
  "c1_school_closing.csv",
  "c2_flag.csv",
  "c2_workplace_closing.csv",
  "c3_cancel_public_events.csv",
  "c3_flag.csv",
  "c4_flag.csv",
  "c4_restrictions_on_gatherings.csv",
  "c5_close_public_transport.csv",
  "c5_flag.csv",
  "c6_flag.csv",
  "c6_stay_at_home_requirements.csv",
  "c7_flag.csv",
  "c7_movementrestrictions.csv",
  "c8_internationaltravel.csv",
  "e1_flag.csv",
  "e1_income_support.csv",
  "e2_debtrelief.csv",
  "h1_flag.csv",
  "h1_public_information_campaigns.csv",
  "h2_testing_policy.csv",
  "h3_contact_tracing.csv",
  "h6_facial_coverings.csv",
  "h6_flag.csv"
))

readox <- function(url) {
  res <- melt.data.table(fread(url, drop = c(1, 3)), id.vars = "country_code")
  res[, date := as.Date(variable, "%d%b%Y")]
  res[, .(iso3=country_code, date, variable = gsub("^.+/(\\w+_\\w+)\\.csv$","\\1",url), value = value)][order(iso3, date)]
}

allox <- rbindlist(lapply(fls, readox))
allox[, base := gsub("^([^_]+)_.+$","\\1", variable) ]
allox[, isflag := grepl("flag", variable) ]
setkey(allox, iso3, date, base)

baseox <- allox[isflag == FALSE]
flagox <- allox[isflag == TRUE, .(iso3, date, base, flag = value)]
combox <- baseox[flagox]
fl
w = 1/c(
  c1=3, c2=3, c3=3, c4=4, c5=2, c6=3, c7=2, 
  h1=2,
  h6=4
)
shft <- sum(w)/(length(w)+2)

alt_si <- mlt.indicators[,
  si := rowSums(cbind(
    (1-shft)*c1_schoolclosing*w["c1"] + shft*c1_flag,
    (1-shft)*c2_workplaceclosing*w["c2"] + shft*c2_flag,
    (1-shft)*c3_cancelpublicevents*w["c3"] + shft*c3_flag,
    (1-shft)*c4_restrictionsongatherings*w["c4"] + shft*c4_flag,
    (1-shft)*c5_closepublictransport*w["c5"] + shft*c5_flag,
    (1-shft)*c6_stayathomerequirements*w["c6"] + shft*c6_flag,
    (1-shft)*c7_movementrestrictions*w["c7"] + shft*c7_flag,
    (1-shft)*h1_publicinformationcampaigns*w["h1"] + shft*h1_flag,
    h2_testingpolicy/3,
    h3_contacttracing/2,
    (1-shft)*h6_facialcoverings*w["h6"] + shft*h6_flag
  ), na.rm = TRUE)*100/(length(w) + 2)
][order(date), .(
  date, si = frollmean(si, 7, align = "right")
), keyby = .(iso3)]

alt_si[order(date), zz := ZigZag(si, 1), by=iso3]
alt_si[order(date), zz.peak := FALSE ]
alt_si[order(date), zz.del := c(NA, sign(diff(zz))), by=iso3]
alt_si[, zz.peak := {
  tmp <- zz.peak
  ind <- which.max(zz.del == -1)-1
  if (ind) {
    zzv <- zz[ind]
    ind <- which.max(si == zzv)
    tmp[ind] <- TRUE
  }
  tmp
},by=iso3]
alt_si[, abs.max := FALSE ]
alt_si[, abs.max := {
  tmp <- abs.max
  tmp[which.max(si)] <- TRUE
  tmp
},by=iso3]


threshold <- .5

res <- alt_si[, .SD[which.max(si > threshold*zz[zz.peak])], keyby=iso3]

fwrite(res, tail(.args, 1))
# timing_assess <- rbindlist(lapply(c("good","bad","short","long"), function(fn) {
#   fread(sprintf("./%s.txt", fn), header = FALSE, col.names = "iso3")[, cat := fn ]
# }))
# add_assess <- function(dt) {
#   dt[timing_assess, on=.(iso3), cat := cat ]
#   dt[is.na(cat), cat := "none" ]
#   setkeyv(dt, c(key(dt), "cat"))
# }
# add_assess(alt_si)
#' 
#' oth <- setkeyv(melt(
#'   fread("~/Dropbox/covidLMIC/inputs/raw_ox_si.csv", drop = "CountryName"),
#'   id.vars = "CountryCode"
#' )[CountryCode != ""][, .(
#'   date=as.Date(variable,"%d%b%Y"), si=value, iso3=CountryCode
#' )][, ver := "orig" ][order(date)], c("iso3","ver"))
#' add_assess(oth)
#' 
#' assemble <- function(oth, alt, thresh) {
#'   oth.mxs <- extract.mx(oth)
#'   oth.xs <- extract.xs(oth, oth.mxs)
#'   alt.mxs <- extract.mx(alt)
#'   alt.xs <- extract.xs(alt, alt.mxs, thresh)
#'   ref <- rbind(oth, alt)
#'   mxs <- rbind(oth.mxs, alt.mxs)
#'   xs <- rbind(oth.xs, alt.xs)
#'   del <- oth.xs[alt.xs, on=setdiff(key(oth.xs), "ver")][, diff := as.numeric(date - i.date)]
#'   list(ref=ref, mxs=mxs, xs=xs, del=del)
#' }
#' 
#' comp <- assemble(oth, alt_si, thresh = 2/3)
#' 
#' comp$del[cat == "short" & diff > -2]
#' 
#' ggplot(comp$del[!(cat %in% c("bad"))]) + aes(diff) + facet_grid(cat ~ .) + geom_histogram()
#' 
#' #' @examples 
#' require(ggplot2)
#' 
#' plotter <- function(oth, alt, filt, delonly = FALSE, thresh = 2/3) with(assemble(oth, alt, thresh), {
#' 
#'   otherstuff <- if(delonly) { list() } else list(
#'     geom_line(),
#'     geom_point(data=mxs[eval(filt)]),
#'     geom_text_repel(aes(label=iso3), data = function(dt) dt[date == "2020-06-01" & ver == "orig"])
#'   )
#'   
#' return(ggplot(
#'   ref[eval(filt)]
#' ) + aes(date, si, color=iso3, linetype = ver, shape = ver, group=interaction(iso3, ver)) +
#'   geom_point(data=xs[eval(filt)]) +
#'   geom_segment(
#'     aes(x=date, xend=i.date, y=si, yend=i.si, color = iso3),
#'     data = del[eval(filt)],
#'     arrow = arrow(angle = 15, length = unit(0.1,"inches")),
#'     inherit.aes = FALSE
#'   ) +
#'   otherstuff +
#'   theme_minimal() +
#'   coord_cartesian(xlim = as.Date(c("2020-02-01","2020-06-01")), ylim = c(25, 100), expand = FALSE) +
#'   scale_color_discrete(guide = "none") +
#'   scale_x_date(
#'     date_breaks = "month",
#'     date_minor_breaks = "week",
#'     date_labels = "%b"
#'   ))
#' })
#' 
#' # scheme needs to not change good isos
#' good_isos <- fread("./good.txt", header = F)[,V1]
#' # slightly change short isos
#' short_isos <- fread("./short.txt", header = F)[,V1]
#' # really change long isos
#' long_isos <- fread("./long.txt", header = F)[,V1]
#' 
#' newthresh <- 3/4
#' 
#' all <- assemble(oth, alt_si, thresh = newthresh)
#' 
#' plotter(oth, alt_si, expression(iso3 %in% good_isos), delonly = TRUE, thresh = newthresh)
#' plotter(oth, alt_si, expression(iso3 == "BGR"))
#' # diff > 0 => intervention starts earlier by this amount
#' # diff < 0 => intervention starts later by this amount
#' all$del[
#'   iso3 %in% good_isos,
#'   .(diff = as.integer(date-i.date), iso3)
#' ][diff < 0]
#' 
#' plotter(oth, alt_si, expression(iso3 %in% short_isos))
#' # doesn't move enough / wrong direction or moves too much
#' all$del[
#'   iso3 %in% long_isos, .(diff = as.integer(date-i.date), iso3)
#' ][, sum((diff >= -7))/.N]
#' 
#' plotter(oth, alt_si, expression(iso3 %in% long_isos))
#' all$del[
#'   iso3 %in% long_isos, .(diff = as.integer(date-i.date), iso3)
#' ][(diff >= -7)]
#' 
#' # maybe do delay
#' urbanization <- readRDS("~/workspaces/epinow2_hpcarray/intros/urban.rds")
#' urbanization[iso3 %in% good_isos, lbl := "good"]
#' urbanization[iso3 %in% short_isos, lbl := "short"]
#' urbanization[iso3 %in% long_isos, lbl := "long"]
#' 
#' pop <- readRDS("~/Dropbox/covidLMIC/inputs/populations.rds")[,
#'  .(pop = log10(sum(pop))), keyby=iso3
#' ]
#' pop[iso3 %in% good_isos, lbl := "good"]
#' pop[iso3 %in% short_isos, lbl := "short"]
#' pop[iso3 %in% long_isos, lbl := "long"]
#' 
#' bth <- urbanization[pop, on=.(iso3, lbl)][!is.na(lbl)]
#' diffs <- assemble(oth, alt_si, newthresh)$del[, .(diff = as.numeric(date - i.date), pk = i.si), keyby=iso3]
#' thing <- bth[diffs, on=.(iso3)][!is.na(lbl)]
#' 
#' ggplot(thing) +
#'   aes(x=pop+log10(value/100), y = diff, color = pk) + geom_point() +
#'   facet_grid(lbl ~ .)
#' 
#' 
#' checkerr <- function(oth, alt, dt, thresh = newthresh) {
#'   diffs <- assemble(oth, alt_si, thresh)$del[, .(diff = as.numeric(date - i.date)), keyby=iso3]
#'   diffs[iso3 %in% good_isos, lbl := "good"]
#'   diffs[iso3 %in% short_isos, lbl := "short"]
#'   diffs[iso3 %in% long_isos, lbl := "long"]
#'   
#'   ggplot(dt[!is.na(lbl)]) +
#'     aes(x=value, fill = lbl) + geom_histogram() +
#'     facet_grid(lbl ~ .)
#'   
#' }
#' 
#' checkerr(oth, alt, urbanization)
#' checkerr(oth, alt, pop)
#' 
#' ggplot(
#'   melt(alt_si[cat != "bad"], id.vars = c("iso3","date","ver","cat"))
#' ) + facet_grid(variable ~ .) + aes(date, value, group = iso3, color = cat) +
#'   geom_line(alpha = 0.1) 
#' 
#' min.pk <- alt_si[!(iso3 %in% bad.isos), max(si, na.rm = TRUE), by=iso3][, min(V1)]
#' 
#' by.pk <- alt_si[!(iso3 %in% bad.isos),.SD[which.max(si >= min.pk)], by=iso3][, ver := "bymin" ]
#' by.pk[iso3 %in% good_isos, lbl := "good"]
#' by.pk[iso3 %in% short_isos, lbl := "short"]
#' by.pk[iso3 %in% long_isos, lbl := "long"]
#' by.pk[is.na(lbl), lbl := "other"]
#' 
#' thing <- by.pk[oth.xs, on=.(iso3)][!is.na(ver), .(diff = date - i.date, lbl, iso3)]
#' 
#' ggplot() + aes(date, si, shape=ver) +
#'   geom_point(data = oth.xs) +
#'   geom_point(data = by.pk) +
#'   geom_segment(
#'     aes(xend=i.date, yend=i.si, color = lbl),
#'     by.pk[oth.xs, on=.(iso3)][!is.na(ver)],
#'     arrow = arrow(15, unit(0.1, "inches"))
#'   ) +
#'   theme_minimal()
#' 
#' extract.rmx <- function(dt) dt[
#'   date < as.Date("2020-06-01"),
#'   .SD[which.max(rolling_si)],
#'   keyby=setdiff(key(dt),"date")
#'   ]
#' 
#' extract.rxs <- function(dt,
#'                         mxs = extract.rmx(dt)[,.SD,.SDcols = c(setdiff(key(dt),"date"), "rolling_si")],
#'                         threshold = 3/4
#' ) setkeyv(dt[
#'   mxs, on=setdiff(key(dt), "date")
#'   ][,.SD[which.max(rolling_si > threshold*i.rolling_si)], keyby=setdiff(key(dt), "date")][,.(iso3, ver, date, rolling_si)], key(dt))
#' 
#' thing <- extract.rxs(alt_si)
#' 
#' del <- oth.xs[thing, on=.(iso3)][,.(diff = date - i.date, iso3)]
#' del[timing_assess, on=.(iso3), cat := cat ]
#' 
#' ggplot(del[!is.na(cat)][!(cat %in% c("bad","none"))]) + aes(diff) + facet_grid(cat ~ .) + geom_histogram()
