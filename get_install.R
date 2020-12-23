
.args <- if (interactive()) c(
  "rpack.txt", ".install"
) else commandArgs(trailingOnly = TRUE)

pcks <- readLines(.args[1])
need <- !sapply(pcks, require, character.only = TRUE)
if (length(pcks[need])) install.packages(pcks[need])
if (dim(installed.packages()[pcks,])[1] == length(pcks)) {
  writeLines("COMPLETE\n", tail(.args, 1))
} else stop("Unable to install all packages.")
