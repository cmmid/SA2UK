
require(data.table)
require(jsonlite)

.args <- if (interactive()) c(
  "4000",
  "sampling.json"
) else commandArgs(trailingOnly = TRUE)

write_json(
  toJSON(
    list(
      samplen = as.integer(.args[1]),
      cores = getDTthreads()
    ),
    pretty = TRUE, auto_unbox = TRUE
  ),
  tail(.args, 1)
)
