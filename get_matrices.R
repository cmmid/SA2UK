
suppressPackageStartupMessages({
  require(socialmixr)
})

.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper"
.args <- if (interactive()) sprintf(c(
  "%s/matrices.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

targets <- c("home", "work", "others", "school")

maturlfmt <- "https://github.com/kieshaprem/synthetic-contact-matrices/raw/master/output/syntheticcontactmatrices2020/urban/contact_%s_urban.rdata"

res <- lapply(targets, function(tar) get(load(url(sprintf(maturlfmt, tar)))))
#' for consistency w/ covidm convention
targets[3] <- "other"
names(res) <- targets

isos <- names(res[[1]])

ages <- list(levels(socialmixr::limits_to_agegroups(integer(), seq(0,75,by=5))))
ages[[2]] <- ages[[1]]

ret <- lapply(isos, function(iso) {
  mats <- lapply(names(res), function(nm) { 
    m <- res[[nm]][[iso]]
    dimnames(m) <- ages
    m
  })
  names(mats) <- names(res)
  mats
})
names(ret) <- isos

#' one substitution for all-Africa analysis
ret[["SOM"]] <- ret[["DJI"]]

saveRDS(ret, tail(.args, 1))