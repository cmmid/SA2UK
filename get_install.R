
.args <- if (interactive()) c(
  "rpack.txt", "../covidm", ".install"
) else commandArgs(trailingOnly = TRUE)

pcks <- readLines(.args[1])
need <- !sapply(pcks, require, character.only = TRUE)
if (length(pcks[need])) install.packages(pcks[need])
if (dim(installed.packages()[pcks,])[1] == length(pcks)) {
  writeLines("COMPLETE\n", tail(.args, 1))
} else stop("Unable to install all packages.")

cm_path <- .args[2]

#' force build of covidm
#' all subsequent uses should have cm_force_shared = T
cm_force_rebuild = T;
cm_build_verbose = T;
cm_force_shared = F;
cm_version = 2;
source(file.path(cm_path, "R", "covidm.R"))