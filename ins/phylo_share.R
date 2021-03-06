suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "~/Dropbox/SA2UK"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv",
  "%s/outputs/phylo.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

rsa.dt <- fread(.args[1])[Country == "South Africa"][order(`Collection Data`)]

rolling.int <- 7
var.label <- "501Y.V2"

#' this choice seems to make minor difference, but unclear which is correct
red.dt <- rsa.dt[,
  .N,
  by = .(
   date=`Collection Data`,
#   variant = fifelse(`Clade`=="20C", var.label, "other"),
   variant = fifelse(`Pangolin lineage`=="B.1.351", var.label, "other")
  )
][date > (as.Date("2020-10-01")-rolling.int) ]

fill.dt <- red.dt[, .(date = as.Date(rep(min(date):max(date), 2), origin = "1970-01-01"))]
fill.dt[, variant := rep(c(var.label, "other"), each = .N/2)]

plot.dt <- dcast(
  red.dt[
    fill.dt, on=.(date, variant),
    .(date, variant, N=fifelse(is.na(N), 0, N))
  ],
  date ~ variant, value.var = "N"
)

plot.dt[, total := as.integer(rowSums(.SD)), .SDcols = -c("date") ]
plot.dt[, rolling.tot := as.integer(frollsum(total, rolling.int)) ]
plot.dt[, rolling.var := as.integer(frollsum(get(var.label), rolling.int)) ]

bino <- function(ci, pos, tot) as.data.table(t(mapply(
  function(x, n, p=x/n) binom.test(x, n, p, conf.level = ci)$conf.int,
  x = pos, n = tot
)))

plot.dt[, binop := rolling.var/rolling.tot ]
plot.dt[!is.na(rolling.var),
        c("lo95","hi95") := bino(.95, rolling.var, rolling.tot)
]
plot.dt[!is.na(rolling.var),
        c("lo50","hi50") := bino(.50, rolling.var, rolling.tot)
]

saveRDS(plot.dt, tail(.args, 1))
