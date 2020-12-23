
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
})

.debug <- c("~/Dropbox/covidLMIC", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/mod_scenarios/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/covidm_fit_yu.qs",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  .debug[2],
  "../covidm",
  "%s/outputs/projections/%s.qs"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

load("NGM.rda")

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

suppressPackageStartupMessages({
  source(file.path(cm_path, "R", "covidm.R"))
})

# identify country / scenario
scenario <- readRDS(.args[1])
tarfile <- tail(.args, 1)

yu_fits <- qread(.args[3])[order(ll)]
yu_fits[, eqs := (1:.N)/.N ]
#' using the median yu fits
medyu <- yu_fits[which.max(eqs > 0.5)]
yref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("y_",colnames(medyu))]))
uref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("u_",colnames(medyu))]))
ys <- rep(yref[1, ], each = 2)
us <- rep(uref[1, ], each = 2)
#if (.debug == "NGA") {
#  ushft <- max(us)/us
#  us <- 1/(1:length(us)+1)
#  ys <- ys / ushft
#}

params <- readRDS(.args[2])

tariso <- tail(.args, 3)[1]

intros <- readRDS(.args[5])[iso3 == tariso][, intro.day := as.integer(date - date[1]) ][, Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))]
urbfrac <- readRDS(.args[6])[iso3 == tariso, value / 100]

params$pop[[1]]$seed_times <- intros
params$pop[[1]]$size <- round(params$pop[[1]]$size*urbfrac)
# params$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))

params$pop <- lapply(
  params$pop,
  function(x){
    x$y <- ys
    x$u <- us
    return(x)
  }
)

run_options <- melt(
  readRDS(.args[4])[era == "pre"],
  measure.vars = c("lo.lo","lo","med","hi","hi.hi"),
  value.name = "r0"
)[, model_seed := 1234L ]

refR0 <- cm_ngm(params)$R0

params_back <- params

allbind <- data.table()

scenario[, waning_dur := NA_integer_ ]
#scenario[!(scen_type == "unmitigated"), waning_dur := 90 ]

prg <- txtProgressBar(max = scenario[,.N]*run_options[,.N], style = 3)
prgind <- 0

for (scenario_index in 1:max(scenario$scen_id)) {
  #' sub
  iv_data <- scenario[scen_id == scenario_index][order(trigger_type)]
  for(i in 1:nrow(run_options)){
    
    params <- params_back
    
    #adjust r0 to that in current sample
    target_R0 <- run_options[i, r0]
    uf <- target_R0 / refR0
    params$pop <- lapply(
      params$pop,
      function(x){
        x$u <- x$u * uf
        return(x)
      }
    )
    
    if (iv_data[!(scen_type == "unmitigated"), .N]) {
      
      # generic interventions
      for (j in 1:nrow(iv_data[population == -1])) {
        #pars <- as.list(iv_data[population == -1][j])
        cons <- with(as.list(iv_data[population == -1][j]), {
          list(1-c(home, work, school, other), c(1,1,1,1))
        })
        tms <- with(as.list(iv_data[population == -1][j]), {
          as.Date(params$date0)+c(start_day, ifelse(is.finite(end_day),end_day,params$time1))
        })
        si <- with(as.list(iv_data[population == -1][j]), {
          list(rep(1-self_iso, 16), rep(1, 16))
        })
        params$schedule[[length(params$schedule)+1]] <- list(
          parameter = "contact",
          pops = numeric(),
          mode = "multiply",
          values = cons,
          times = tms
        )
        params$schedule[[length(params$schedule)+1]] <- list(
          parameter = "fIs",
          pops = numeric(),
          mode = "multiply",
          values = si,
          times = tms
        )
      }
    }
    #run the model
    sim <- cm_simulate(
      params, 1,
      model_seed = run_options[i, model_seed]
    )$dynamics
    
    result <- sim[,
      .(value = sum(value)),
      keyby = .(run, t, group, compartment)
    ][, run := i ][, scen_id := scenario_index ]
    
    allbind <- rbind(allbind, result)
    prgind <- prgind + 1
    setTxtProgressBar(prg, prgind)
  }
  
}

#' @examples 
#' require(ggplot2)
#' ggplot(allbind[compartment == "cases"]) + aes(t, value) + facet_grid(group ~ ., scales = "free_y") + geom_line()

qsave(allbind, tarfile)
