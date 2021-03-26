
.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper"
.args <- if (interactive()) sprintf(c(
  "matrices.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

targets <- c("home", "others", "school", "work")

maturlfmt <- "https://github.com/kieshaprem/synthetic-contact-matrices/raw/master/output/syntheticcontactmatrices2020/urban/contact_%s_urban.rdata"

res <- lapply(targets, function(tar) get(load(url(sprintf(maturlfmt, tar)))))
names(res) <- targets

isos <- names(res[[1]])

ret <- lapply(isos, function(iso) {
  mats <- lapply(names(res), function(nm) res[[nm]][[iso]])
  names(mats) <- names(res)
  mats
})
names(ret) <- isos

saveRDS(ret, tail(.args, 1))