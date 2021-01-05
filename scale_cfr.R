#--- template script for calculating the delay-adjusted CFR
#--- using case and death time-series data
library(patchwork)
library(ggplot)

# Set parameters for the three delay distributions
dd_mean_mid <- 13.0
dd_median_mid <- 9.1

dd_mu_mid <- log(dd_median_mid)
dd_sigma_mid <- sqrt(2 * (log(dd_mean_mid) - dd_mu_mid))

dd_mean_low <- 8.7
dd_median_low <- 6.7

dd_mu_low <- log(dd_median_low)
dd_sigma_low <- sqrt(2 * (log(dd_mean_low) - dd_mu_low))

dd_mean_high <- 20.9
dd_median_high <- 13.7

dd_mu_high <- log(dd_median_high)
dd_sigma_high<- sqrt(2 * (log(dd_mean_high) - dd_mu_low))

#--- for bootstrapping over distribution of distributions
#--- not required, as ordering changes exactly once between
#--- distributions. Therefore choosing minimum and
#--- maximum values is equivalent to bootstrapping (in 
#--- the limit of bootstrap samples)
#length_out_arg <- 100
#mean_interval <- seq(8.7, 20.9, length.out = length_out_arg)
#median_interval <- seq(6.7, 13.7, length.out = length_out_arg)

#mu_interval <- seq(log(min(median_interval)), 
                   #log(max(median_interval)), 
                   #length.out = length_out_arg)

#sigma_interval  <- seq(sqrt(2*(log(min(mean_interval)) - min(mu_interval))),
                       #sqrt(2*(log(max(mean_interval)) - max(mu_interval))),
                       #length.out = length_out_arg)

#--- testing normally distributed variates over the 
#--- 95% confidence interval of the reported delay distribution
#--- for bootstrapping
#mu_variates <- rnorm(n = 500,
                     #mean = mean(mu_interval),
                     #sd = 2*sd(mu_interval))
#sigma_variates <- rnorm(n = 500,
                        #mean = mean(sigma_interval),
                        #sd = 2*sd(sigma_interval))


#hosp_to_death_bootstrapping <- function(x) {
    #mu_ins <- rnorm(n = 1,
                    #mean = mean(mu_interval),
                    #sd = sd(mu_interval))
    #sigma_ins <- rnorm(n = 1,
                       #mean = mean(sigma_interval),
                       #sd = sd(sigma_interval))

    #plnorm(x + 1, mu_ins, sigma_ins) - plnorm(x, mu_ins, sigma_ins)

#}

##--- testing the bootstrapped delay distributions
#for(i in 1:1000) {
    #x_seq %>%
        #hosp_to_death_bootstrapping(.) %>% lines(., col = rgb(0,0,0,0.3))
#}

# Hospitalisation to death distributions
# One using the median estimate and two using the lower and
# Upper 95% confidence interval estimates
hosp_to_death_mid <- function(x, mu, sigma) {
  plnorm(x + 1, dd_mu_mid, dd_sigma_mid) - plnorm(x, dd_mu_mid, dd_sigma_mid)
}

hosp_to_death_low <- function(x, mu, sigma) {
  plnorm(x + 1, dd_mu_low, dd_sigma_low) - plnorm(x, dd_mu_low, dd_sigma_low)
}

hosp_to_death_high <- function(x, mu, sigma) {
  plnorm(x + 1, dd_mu_high, dd_sigma) - plnorm(x, dd_mu_high, dd_sigma)
}

#--- plotting three distributions to inspect them
x_seq <- seq(0, 40, 0.05)
plot(hosp_to_death_mid(x_seq), type = "l", col = "green")
lines(hosp_to_death_low(x_seq), type = "l", col = "red")
lines(hosp_to_death_high(x_seq), type = "l", col = "blue")


# function to calculate the adjusted number of 'known outcomes'
# methods from Nishiura et al. (2009)
scale_cfr_rolling <- function(data_arg,
                              delay_fun) {
  #if (smoothed == TRUE)
  #{
    #case_incidence <-  data_arg$new_cases_smoothed
    #death_incidence <- data_arg$new_deaths_smoothed
  #}
  #else 
  #{
    #case_incidence <-  data_arg$new_cases
    #death_incidence <- data_arg$new_deaths
  #}
  point_known_t <- NULL # point estimate of cases with known outcome at time tt
  
  # Sum over cases up to time tt
  for (ii in 1:length(data_arg$new_cases_smoothed)) {
    
    known_i <- 0 # number of cases with known outcome at time ii
    
    for (jj in 0:(ii - 1)) {
      known_jj <- (data_arg$new_cases_smoothed[ii - jj] * delay_fun(jj))
      known_i <- known_i + known_jj
    }

    # point estimate of known outcomes
    point_known_t <- c(point_known_t, known_i)
  }
  
  # naive CFR value
  #b_tt_series <- death_incidence / case_incidence

  # corrected CFR estimator (rolling)
  p_tt_series <- data_arg$new_deaths_smoothed / point_known_t
  
  results_df <- dplyr::tibble(
    date = data_arg$date,
    #known_t = point_known_t,
    #nCFR = b_tt_series,
    cCFR = p_tt_series)
  
  return(results_df)
  
}

# function which simplifies downloading and importing all data
case_data_function <- function() {

    owid_path <- "https://covid.ourworldindata.org/data/owid-covid-data.csv"
    case_data <- read.csv(owid_path) %>%
    dplyr::tibble() %>%
    dplyr::mutate(date = lubridate::ymd(date)) %>%
    dplyr::select(iso_code,
                  country = location,
                  date, new_cases,
                  new_deaths,
                  population)

    return(case_data)
}

# downloading all case data from JHU
case_data <- case_data_function()

# filtering for just the SA data
sa_case_data <- case_data %>%
  dplyr::filter(iso_code == "ZAF") %>%
  tidyr::drop_na()

# smoothing the new cases and new deaths data
sa_case_data_smoothed  <- sa_case_data %>%
    dplyr::mutate(new_cases_smoothed = zoo::rollmean(new_cases,
                                                     k = 21, 
                                                     na.pad = TRUE)) %>%
    dplyr::mutate(new_deaths_smoothed = zoo::rollmean(new_deaths, 
                                                      k = 21,
                                                      na.pad = TRUE)) %>%
    tidyr::drop_na()

#library(foreach)
#sa_case_data_dt <- data.table::data.table(sa_case_data_smoothed)
#cl <- parallel::makeCluster(4)
#doParallel::registerDoParallel(cl)
#res <- double(10)
##res <- 
    #system.time(foreach::foreach(i = 1:10, .combine = 'c') %dopar% {  
    #sa_case_data_dt %>%
        #scale_cfr_rolling(delay_fun = hosp_to_death_bootstrapping) 
#} )

#replicate(5,
          #sa_case_data_dt %>% 
          #scale_cfr_rolling(delay_fun = hosp_to_death_bootstrapping),
          #simplify = FALSE) %>%
#system.time()

# calculating the scaled CFR using the median, lower and upper 95% 
# of the delay distribution
cfr_rolling_smoothed_mid <- sa_case_data_smoothed %>%
    tidyr::drop_na() %>%
    scale_cfr_rolling(delay_fun = hosp_to_death_mid,
                      smoothed = TRUE) %>%
    dplyr::select(date, dCFR_mid = cCFR)

cfr_rolling_smoothed_low <- sa_case_data_smoothed %>%
    tidyr::drop_na() %>%
    scale_cfr_rolling(delay_fun = hosp_to_death_low,
                      smoothed = TRUE) %>%
    dplyr::select(date, dCFR_low = cCFR)

cfr_rolling_smoothed_high <- sa_case_data_smoothed %>%
    tidyr::drop_na() %>%
    scale_cfr_rolling(delay_fun = hosp_to_death_high,
                      smoothed = TRUE) %>%
    dplyr::select(date, dCFR_high = cCFR)

scaled_cfr_smoothed_all <- cfr_rolling_smoothed_mid %>%
    dplyr::left_join(cfr_rolling_smoothed_low) %>%
    dplyr::left_join(cfr_rolling_smoothed_high) 

scaled_cfr <- scaled_cfr_smoothed_all %>%
    dplyr::rowwise() %>%
    dplyr::mutate(dCFR_low_all = min(dCFR_mid, dCFR_low, dCFR_high)) %>%
    dplyr::mutate(dCFR_high_all = max(dCFR_mid, dCFR_low, dCFR_high))

# plotting the daily scaled CFR, as example plot only showing
# slightly adjusted new method
scaled_cfr %>%
    dplyr::filter(date > "2020-05-01") %>% 
    ggplot2::ggplot(ggplot2::aes(x = date)) +  
    ggplot2::geom_line(ggplot2::aes(y = dCFR_mid), linetype = "dashed") +  
    #ggplot2::geom_line(ggplot2::aes(y = nCFR), color = "red") +  
    ggplot2::geom_ribbon(ggplot2::aes(ymin = dCFR_low_all, 
                                      ymax = dCFR_high_all), 
                         alpha = 0.3, fill = "dodgerblue") + 
    ggplot2::ylim(0, 0.1)
