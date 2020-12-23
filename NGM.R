suppressPackageStartupMessages({
  require(data.table)
})

.args <- if (interactive()) c(
  "NGM.rda"
) else commandArgs(trailingOnly = TRUE)

#' TODO remove boxcar stuff - recall that it has no
#' impact on mean stay times in compartments
#' (by design - nor would any other mean preserving transform),
#' therefore introducing them has no impact on R0 calcs
cm_ngm <- function(
  cm_params,
  R0_multiplier = 1,
  contact_reductions = c(0, 0, 0, 0),
  fIs_reductions = rep(0, 16),
  infected_states = c("E","Ia","Ip","Is"),
  ts = cm_params$time_step
) {
  
  ages <- 1:max(sapply(cm_params$pop, function(p) length(p$size)))
  pops <- 1:length(cm_params$pop)

  dur <- function(dX) weighted.mean((seq_along(dX)-1)*ts, dX)
  getNms <- function(fmt, nms, p) grep(sprintf(fmt, p), nms, value = TRUE)
  v.dur <- function(dX) {
    xm <- dur(dX)
    sum(dX*((seq_along(dX)-1)*ts - xm)^2)
  }
  
  #' for boxcar compartments -- i.e. Erlang instead of exponential distro
  #' overall mean dur = k * compartment mean dur = distro mean dur
  #' overall var = k * (compartment mean dur)^2 = 1/k (distro mean dur)^2
  #' want to do 'well' on sd, while exactly matching mean dur
  #' so k = floor((distro mean dur)^2 / (distro var))
  bestk <- function(dX) floor(dur(dX)^2 / v.dur(dX))
  
  #' assumes all pops have same dE
  ks <- with(cm_params$pop[[1]], c(bestk(dE), bestk(dIa), bestk(dIp), bestk(dIs)))
  k_states <- sprintf(
    "%sk%i",
    rep(infected_states, ks),
    Reduce(c, sapply(ks, function(m) 1:m))
  )
  Ekexit <- k_states[ks[1]]
  Iaexit <- k_states[cumsum(ks)[2]]
  Ipexit <- k_states[cumsum(ks)[3]]
  Isexit <- k_states[cumsum(ks)[4]]
  
  #' TODO synthesize boxcars
  nm.dt <- data.table(expand.grid(X=k_states, A=ages, P=pops))[
    , nm := sprintf("%s_a%02i_p%02i", X, A, P)
  ]
  nms <- nm.dt$nm
  Enms <- nm.dt[X %like% "E", nm]
  Ianms <- nm.dt[X %like% "Ia", nm]
  Ipnms <- nm.dt[X %like% "Ip", nm]
  Isnms <- nm.dt[X %like% "Is", nm]
  N <- length(nms)
  
  Sigma_ij <- matrix(0, N, N, dimnames = list(to=nms, from=nms))
  T_ij <- Sigma_ij #' copy constructor
  
  #' compute transitions
  #' consistent with model, assumes all ages alike
  for (p in pops) {
    pop <- cm_params$pop[[p]]
    tarEnms <- getNms("p%02i$", Enms, p)
    tarIanms <- getNms("p%02i$", Ianms, p)
    tarIpnms <- getNms("p%02i$", Ipnms, p)
    tarIsnms <- getNms("p%02i$", Isnms, p)
    
    durE <- dur(pop$dE)
    durIa <- dur(pop$dIa)
    durIp <- dur(pop$dIp)
    durIs <- dur(pop$dIs)
    
    #' leaving X
    diag(Sigma_ij[tarEnms, tarEnms])   <- -ks[1]/durE
    diag(Sigma_ij[tarEnms, tarEnms][-1,]) <- ks[1]/durE
    diag(Sigma_ij[tarEnms, tarEnms][-1,][
      head(ages*ks[1],-1), ages*ks[1]
    ]) <- 0
    
    diag(Sigma_ij[tarIanms, tarIanms]) <- -ks[2]/durIa
    diag(Sigma_ij[tarIanms, tarIanms][-1,]) <- ks[2]/durIa
    diag(Sigma_ij[tarIanms, tarIanms][-1,][
      head(ages*ks[2],-1), ages*ks[2]
    ]) <- 0
    
    diag(Sigma_ij[tarIpnms, tarIpnms]) <- -ks[3]/durIp
    diag(Sigma_ij[tarIpnms, tarIpnms][-1,]) <- ks[3]/durIp
    diag(Sigma_ij[tarIpnms, tarIpnms][-1,][
      head(ages*ks[3],-1), ages*ks[3]
    ]) <- 0
    
    diag(Sigma_ij[tarIsnms, tarIsnms]) <- -ks[4]/durIs
    diag(Sigma_ij[tarIsnms, tarIsnms][-1,]) <- ks[4]/durIs
    diag(Sigma_ij[tarIsnms, tarIsnms][-1,][
      head(ages*ks[4],-1), ages*ks[4]
    ]) <- 0
    
    #' incubation: E => Ia | Ip
    #' n.b. tarXnms all have same age order
    diag(Sigma_ij[
      tarIanms[(ages-1)*ks[2]+1], tarEnms[ages*ks[1]]
    ]) <- (1-pop$y)*ks[1]/durE
    diag(Sigma_ij[
      tarIpnms[(ages-1)*ks[3]+1], tarEnms[ages*ks[1]]
    ]) <- pop$y*ks[1]/durE
    #' onset: Ip => Is
    diag(Sigma_ij[
      tarIsnms[(ages-1)*ks[4]+1], tarIpnms[ages*ks[3]]
    ]) <- ks[3]/durIp
  }
  
  invSigma <- solve(Sigma_ij)
  
  ordbyage <- function(Xnms, k) as.vector(t(matrix(Xnms, nrow = k)))
  
  #' transmission
  for (srcp in pops) {
    cm = Reduce('+', mapply(
      function(c, m) c * m,
      cm_params$pop[[srcp]]$contact*(1-contact_reductions),
      cm_params$pop[[srcp]]$matrices,
      SIMPLIFY = F
    ))
    srcIanms <- getNms("p%02i$", Ianms, srcp)
    srcIpnms <- getNms("p%02i$", Ipnms, srcp)
    srcIsnms <- getNms("p%02i$", Isnms, srcp)
    for (tarp in pops) {
      tarEnms <- getNms("p%02i$", Enms, tarp)
      travelrel <- cm_params$travel[tarp, srcp]
      foibase <- travelrel * (cm * cm_params$pop[[tarp]]$u) * R0_multiplier
      #' covidm uses dX = -binom(X, 1-exp(-foi dt)) => E[dX] as dt->0 = -X(1-(1-foi dt)) = -X foi dt => dot X = -foi*X
      #' so in the expected-value, linearization framework of NGM, standard foi
      #' relationship
      reordIa <- as.vector(t(matrix(srcIanms, nrow = ks[2])))
      
      T_ij[tarEnms[(ages-1)*ks[1]+1], ordbyage(srcIanms, ks[2])] <- foibase %*% diag(cm_params$pop[[srcp]]$fIa)
      T_ij[tarEnms[(ages-1)*ks[1]+1], ordbyage(srcIpnms, ks[3])] <- foibase %*% diag(cm_params$pop[[srcp]]$fIp)
      T_ij[tarEnms[(ages-1)*ks[1]+1], ordbyage(srcIsnms, ks[4])] <- foibase %*% diag(cm_params$pop[[srcp]]$fIs * (1-fIs_reductions))
    }
  }
  
  res <- eigen(-T_ij %*% invSigma)
  
  if (Im(res$values[1]) != 0) warning(
    sprintf("Non-real principle eigenvalue: %02g", Im(res$values[1]))
  )
  if (any(Im(res$vector[,1]) != 0)) warning(
    sprintf("Non-real component(s) of principle eigenvector.")
  )
  if (!(all(Re(res$vector[,1]) <= 0) | all(Re(res$vector[,1]) >= 0))) warning(
    sprintf("Negative value(s) in principle eigenvector.")
  )
  ss <- abs(Re(res$vector[,1]))
  names(ss) <- nms
  list(R0 = Re(res$values[1]), ss = ss)
}

cm_ss_age_distro <- function(eigenv_ss) {
  sts <- eigenv_ss[which(eigenv_ss != 0)]
  sts/sum(sts)
}

#' assumes single pop model currently
cm_generation_time <- function(cm_params, ...) {
  res <- cm_ngm(cm_params, ...)
  ws <- cm_ss_age_distro(res$ss)
  clinfrac <- cm_params$pop[[1]]$y
  dur <- function(dX) weighted.mean((seq_along(dX)-1)*cm_params$time_step, dX)
  durE <- dur(cm_params$pop[[1]]$dE)
  durIa <- dur(cm_params$pop[[1]]$dIa)
  durIp <- dur(cm_params$pop[[1]]$dIp)
  durIs <- dur(cm_params$pop[[1]]$dIs)
  #' each E takes E[dur E] time to become infectious
  #' the asymp fraction leads integral (foi * time asymp) infections
  ave_gen <- 0
  for (a in 1:length(ws)) {
    outws <- c(
      cm_params$pop[[1]]$fIa[a]*durIa,
      cm_params$pop[[1]]$fIp[a]*durIp,
      cm_params$pop[[1]]$fIs[a]*durIs
    )
    ave_gen <- ave_gen + ws[a]*sum((c(durIa,durIp,durIp+durIs)+durE)*outws)/sum(outws)
  }
  ave_gen
}

save(list=ls(), file = tail(.args, 1))
