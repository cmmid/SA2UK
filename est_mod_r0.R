suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/covidLMIC","PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenarios/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/covidm_fit_yu.qs",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  .debug[2],
  "../covidm",
  "%s/mod_scenarios/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

# identify country / scenario
scenario <- readRDS(.args[1])
tariso <- tail(.args, 3)[1]

if (scenario[scen_type != "unmitigated" & is.na(self_iso), .N]) {
  load("NGM.rda")
  
  # load covidm
  cm_path = tail(.args, 2)[1]
  cm_force_rebuild = F;
  cm_build_verbose = F;
  cm_force_shared = T
  cm_version = 2
  
  suppressPackageStartupMessages({
    source(file.path(cm_path,"R","covidm.R"))
  })

  yu_fits <- qread(.args[3])[order(ll)]
  yu_fits[, eqs := (1:.N)/.N ]
  #' using the median yu fits
  medyu <- yu_fits[which.max(eqs > 0.5)]
  yref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("y_",colnames(medyu))]))
  uref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("u_",colnames(medyu))]))
  ys <- rep(yref[1, ], each = 2)
  us <- rep(uref[1, ], each = 2)

  params <- readRDS(.args[2])
  
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
  
  scenario[, waning_dur := NA_integer_ ]
  # scenario[!(scen_type == "unmitigated"), waning_dur := 90 ]
  
  prg <- txtProgressBar(max = scenario[,(.N-1)/2]*run_options[,.N], style = 3)
  prgind <- 0
  
  allbind <- data.table(run=integer(0), group = factor(), value = integer(), AR = numeric(), scen_id = integer(), r_id = integer())
  
  for (scenario_index in 2:max(scenario$scen_id)) {
    #' sub
    iv_data <- scenario[scen_id == scenario_index & !is.na(self_iso)][order(trigger_type)]
    for(i in 1:nrow(run_options)){
      
      params <- params_back
      # only run until recalc point
      params$time1 <- scenario[scen_id == scenario_index & !is.na(self_iso), max(end_day, na.rm = TRUE)]
      
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
      )$dynamics[compartment == "R"][
        order(t), .(value = value[.N]), keyby=.(run, group)
      ]
      
      sim[order(run, group), AR := value / params$pop[[1]]$size ]
      
      sim[, scen_id := scenario_index ][, r_id := i ]
      
      allbind <- rbind(allbind, sim)
      prgind <- prgind + 1
      setTxtProgressBar(prg, prgind)
    }
    
  }
  
  #' now that we have susceptible fraction after initial interventions, need to fit reductions
  #' that best adjust R
  
  modR <- melt(readRDS(.args[4])[era == "modification"], id.vars = "era", measure.vars = 2:6)
  tarRs <- modR[, value]
  baseRs <- run_options[, r0]
  
  for (scenario_index in 2:max(scenario$scen_id)) {
    basescn <- scenario[!is.na(self_iso) & scen_id == scenario_index]

    ufs <- mapply(
      function(br, sfrac) br/refR0 * sfrac,
      br=baseRs,
      sfrac = lapply(1:5, function(ri) 1-allbind[scen_id == scenario_index & r_id == ri, AR]),
      SIMPLIFY = FALSE
    )
    
    vs <- basescn[,c(home, work, school, other, self_iso)]
    
    mvs <- function(del) {
      fac <- ifelse(del >= 0, 1, 0) - sign(del)*exp(-abs(del))
      return(if (del > 0) {
        vs*(1-fac) + fac
      } else {
        vs*fac
      })
    }
    
    fn <- function(del) {
      #' del can range from -inf to +inf
      #' del < 0 => decreasing the contact reduction toward 0
      #' del > 0 => increasing the contact reduction toward 1
      #' del == 0 => no change
      rs <- mvs(del)
      sum((sapply(ufs, function(uf)
        cm_ngm(
          params_back,
          uf,
          contact_reductions = rs[-5],
          fIs_reductions = rs[5]
        )$R0
      ) - tarRs)^2)
    }
    
    odel <- optimize(fn, c(-20, 20))$minimum
    
    scenario[
      is.na(self_iso) & scen_id == scenario_index,
      c("home","work","school","other", "self_iso") := as.list(mvs(odel))
    ]
    
  }
  
  saveRDS(scenario, tail(.args, 1))
  
} else {
  warning(sprintf("no scenario modification to generate for %s; making a link instead.", tariso))
  file.link(.args[1], tail(.args, 1))
}
