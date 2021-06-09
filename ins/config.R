
require(data.table) #' only for getting available threads
require(jsonlite)

.debug <- "analysis"
.args <- if (interactive()) c(
  "200", # low for debugging
  file.path(.debug, "ins", "sampling.json")
) else commandArgs(trailingOnly = TRUE)

write_json(
  x = list(
    samplen = as.integer(.args[1]),
    rtsamplemul = 10,
    cores = getDTthreads()
  ),
  path = tail(.args, 1),
  pretty = TRUE, auto_unbox = TRUE
)
