
suppressPackageStartupMessages({
  require(data.table)
  require(optimization)
})

.debug <- c("~/Dropbox/SA2UK","ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/r0/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/yuqs/%s.rds",
  "%s/outputs/fits/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

R0ref <- readRDS(.args[1])
pars <- readRDS(.args[2])

load("NGM.rda")

#' TODO remove magic numbers for quantile targets
tarqs <- c("lo","med","hi")
extractRt <- function(dt, cols = tarqs) melt(
  dt[,.SD, .SDcols = cols],
  measure.vars = cols, variable.name = "q"
)
#' TODO handle exceptions
preR0 <- extractRt(R0ref[era == "pre"])
postR0 <- extractRt(R0ref[era == "post"])
yuref <- readRDS(.args[3])
inds <- yuref[order(eqs)][with(yuref[order(eqs)], c(which.max(eqs >= .25), which.max(eqs >= .5), which.max(eqs >= .75)))]

inds[, q := tarqs ]
inds[preR0, on = .(q), umul := value / baseR ]
inds[postR0, on = .(q), tarR := value ]

ucols <- rep(grep("^u_", names(inds)), each = 2)
ycols <- rep(grep("^y_", names(inds)), each = 2)

#' TODO constrain to be similar?
ret <- inds[,{
  us <- as.numeric(.SD[, 1:16, with = F])*umul
  ys <- as.numeric(.SD[, -(1:16), with = F])
  res <- optim_sa(
    function(ps) { 
      lrg <- ps[1]; sml <- ps[2]; symp <- ps[3]
      if ((lrg < sml) | (lrg < symp)) NA_real_ else (cm_ngm(
        pars, us, contact_reductions = c(home=0, work=sml, school=lrg, other=sml),
        fIs_reductions = symp, ymod = ys
    )$R0-tarR)^2 },
    c(lrg=0.7, sml=0.5, symp=0.25),
    lower = c(0.2, 0.2, 0.2), upper = c(0.8, 0.8, 0.8),
    control = list(nlimit = 100)
  )
  names(res$par) <- c("largered", "smallred", "sympred")
  as.list(res$par)
}, by=q, .SDcols = c(ucols, ycols)]

#' @examples 
checker <- function(lrg, sml, symp, us, ys) cm_ngm(
  pars, us, contact_reductions = c(home=0, work=sml, school=lrg, other=sml),
  fIs_reductions = symp, ymod = ys
)$R0

ret[, era := "post" ][, date := R0ref[era == "post"][1, date]][,
  c("work", "school", "other") := .("small", "large", "small")
][, Rp := {
  checker(
    largered, smallred, sympred,
    rep(as.numeric(inds[q == qtar, grep("^u_", names(inds)), with = F]), each = 2)*inds[q == qtar, umul],
    rep(as.numeric(inds[q == qtar, grep("^y_", names(inds)), with = F]), each = 2)
  )
}, by=.(qtar = q)]

saveRDS(ret, tail(.args, 1))
