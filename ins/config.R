
require(data.table) #' only for getting available threads
require(jsonlite)

.args <- if (interactive()) c(
  "4000",
  "sampling.json"
) else commandArgs(trailingOnly = TRUE)

write_json(
  x = list(
    samplen = as.integer(.args[1]),
    cores = getDTthreads()
  ),
  path = tail(.args, 1),
  pretty = TRUE, auto_unbox = TRUE
)
