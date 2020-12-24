#--- template script for calculating the delay-adjusted CFR
#--- using case and death time-series data

# Set parameters for the delay distribution 
mean <- 13
median <- 9.1

mu <- log(median)
sigma <- sqrt(2*(log(mean) - mu))

# Hospitalisation to death distribution
hospitalisation_to_death_truncated <- function(x) {
  plnorm(x + 1, mu, sigma) - plnorm(x, mu, sigma)
}

# Function to work out correction to naive CFR
scale_cfr <- function(data_arg, delay_fun = hospitalisation_to_death_truncated)
{
  case_incidence <- data_arg$new_cases
  death_incidence <- data_arg$new_deaths
  cumulative_known_t <- NULL # cumulative cases with known outcome at time tt
  # Sum over cases up to time tt
  for(ii in 1:nrow(data_arg)){
    known_i <- 0 # number of cases with known outcome at time ii
    for(jj in 0:(ii - 1)){
      known_jj <- (case_incidence[ii - jj]*delay_fun(jj))
      known_i <- known_i + known_jj
    }
    cumulative_known_t <- c(cumulative_known_t,known_i) # Tally cumulative known
  }
  
  # naive CFR value
  b_tt <- sum(death_incidence)/sum(case_incidence) 
  # corrected CFR estimator
  p_tt <- (death_incidence/cumulative_known_t) %>% pmin(.,1)
  
  results.df <- dplyr::tibble(
    nCFR = b_tt,
    dCFR = p_tt, 
    total_deaths = sum(death_incidence),
    deaths = death_incidence,
    cum_known_t = round(cumulative_known_t),
    total_cases = sum(case_incidence))
  
  return(results.df)
}

#--- calculating delay-adjusted CFR, called cCFR in results dataframe
#--- AN EXAMPLE USING THE CONFIRMED NUMBERS FOR SOUTH AFRICA 
#--- taken from the John Hopkins database
case_data_function <- function()
{
  
  case_data <- readr::read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv") %>%
    dplyr::mutate(date = lubridate::ymd(date)) %>%
    dplyr::select(iso_code, country = location, date, new_cases, new_deaths, population)
  
  return(case_data)
}

case_data <- case_data_function()

south_africa_case_data <- case_data %>%
  dplyr::filter(iso_code == "ZAF") %>%
  tidyr::drop_na()

scaled_cfr_SA_confirmed <- scale_cfr(south_africa_case_data)
