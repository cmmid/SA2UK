#--- template script for calculating the delay-adjusted CFR
#--- using case and death time-series data
library(patchwork)
library(ggplot)

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
  
  date <- data_arg$date
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
    date = date,
    nCFR = b_tt,
    dCFR = p_tt,
    new_cases = case_incidence,
    new_deaths = death_incidence,
    total_cases = cumsum(case_incidence),
    total_deaths = cumsum(death_incidence),
    cum_known_t = round(cumulative_known_t)
  )
  
  return(results.df)
}

#--- calculating delay-adjusted CFR, called cCFR in results dataframe
#--- AN EXAMPLE USING THE CONFIRMED NUMBERS FOR SOUTH AFRICA 
#--- taken from the John Hopkins database
case_data_function <- function()
{
  
  case_data <- read.csv("https://covid.ourworldindata.org/data/owid-covid-data.csv") %>%
    dplyr::tibble() %>%
    dplyr::mutate(date = lubridate::ymd(date)) %>%
    dplyr::select(iso_code, country = location, date, new_cases, new_deaths, population)
  
  return(case_data)
}

case_data <- case_data_function()

south_africa_case_data <- case_data %>%
  dplyr::filter(iso_code == "ZAF") %>%
  tidyr::drop_na()

scaled_cfr_data <- scale_cfr(south_africa_case_data) 

scaled_cfr_SA_confirmed <- scaled_cfr_data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(nCFR = total_deaths/total_cases,
                dCFR = cum_known_t/total_cases,
                dCFR_lower_ci = binom.test(cum_known_t, total_cases)$conf.int[[1]],
                dCFR_upper_ci = binom.test(cum_known_t, total_cases)$conf.int[[2]])

scaled_cfr_SA_confirmed_longer <- scaled_cfr_SA_confirmed %>%
  tidyr::pivot_longer(cols = c(nCFR, dCFR),
                      names_to = "cfr_type", values_to = "cfr")

#--- plotting the dCFR results along with the case and death time-series data
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 1
cols = gg_color_hue(n)


cfr_subplot <- scaled_cfr_SA_confirmed %>%
  dplyr::select(-nCFR) %>%
  dplyr::mutate(date = lubridate::ymd(date)) %>%
  dplyr::filter(date > "2020-11-01") %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x = date, y = dCFR), color = cols) + 
  ggplot2::geom_line(ggplot2::aes(x = date, y = dCFR), color = cols) +
  ggplot2::geom_errorbar(ggplot2::aes(x = date, ymin = dCFR_lower_ci, ymax = dCFR_upper_ci),
                         color = "black", width = 0.01) + 
  ggplot2::labs(x = "Date", y = "Case Fatality Ratio") + 
  ggplot2::scale_y_continuous(labels = scales::percent) 

cfr_plot <- scaled_cfr_SA_confirmed %>%
  dplyr::select(-nCFR) %>%
  dplyr::mutate(date = lubridate::ymd(date)) %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x = date, y = dCFR), color = cols) + 
  ggplot2::geom_line(ggplot2::aes(x = date, y = dCFR), color = cols) +
  ggplot2::geom_errorbar(ggplot2::aes(x = date, ymin = dCFR_lower_ci, ymax = dCFR_upper_ci),
                         color = "black", width = 0.01) +
  ggplot2::labs(x = "Date", y = "Case Fatality Ratio", title = "A") + 
  ggplot2::scale_y_continuous(labels = scales::percent)

inset <- ggplotGrob(cfr_subplot)

cfr_plot_with_inset <- cfr_plot + annotation_custom(grob = inset, 
                                                    xmin = as.Date("2020-07-01"), 
                                                    xmax = as.Date("2021-01-01"),
                                                    ymin = 0.04, ymax = 0.095)

confirmed_cases_plot <- south_africa_case_data %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = new_cases)) + 
  ggplot2::geom_col(fill = "dodgerblue", alpha = 0.7, width = 1) +
  ggplot2::labs(x = "Date", y = "New cases",  title = "B") 

confirmed_deaths_plot <- south_africa_case_data %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = new_deaths)) + 
  ggplot2::geom_col(fill = "firebrick1", alpha = 0.6, width = 1) +
  ggplot2::labs(x = "Date", y = "New deaths", title = "C")

cfr_plot_final <- cfr_plot_with_inset / confirmed_cases_plot / confirmed_deaths_plot

ggplot2::ggsave(here::here("SA2UK/figure_2.png"), cfr_plot_final,
                width = 8,
                height = 12,
                dpi = 300)
