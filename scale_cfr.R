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
dd_sigma_high <- sqrt(2 * (log(dd_mean_high) - dd_mu_high))

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
  plnorm(x + 1, dd_mu_high, dd_sigma_high) - plnorm(x, dd_mu_high, dd_sigma_high)
}

# function to calculate the adjusted number of 'known outcomes'
# methods from Nishiura et al. (2009)
scale_cfr <- function(data_arg,
                      delay_fun = hosp_to_death_mid,
                      date_num) {

    case_incidence <-  data_arg$new_cases[1:date_num]
    death_incidence <- data_arg$new_deaths[1:date_num]
    cumulative_known_t <- 0 # cumulative cases with known outcome at time tt

    # Sum over cases up to time tt
    for (ii in 1:date_num) {

    known_i <- 0 # number of cases with known outcome at time ii

    for (jj in 0:(ii - 1)) {
            known_jj <- (case_incidence[ii - jj] * delay_fun(jj))
            known_i <- known_i + known_jj
    }

            # tallying cumulative known outcomes
            cumulative_known_t <- cumulative_known_t + known_i
    }

    # naive CFR value
    b_tt <- sum(death_incidence) / sum(case_incidence)

    # corrected CFR estimator
    p_tt <- sum(death_incidence) / cumulative_known_t

    results_df <- dplyr::tibble(
        date = data_arg$date[date_num],
        nCFR = b_tt,
        cCFR = p_tt,
        total_cases = sum(case_incidence),
        total_deaths = sum(death_incidence),
        cum_known_t = round(cumulative_known_t))

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

dates_num <- seq(1, nrow(south_africa_case_data), 1)

# scaling the CFR using the median, lower and upper 95% 
# confidence intervals of the delay distribution
scaled_cfr_sa_mid <- dates_num %>%
    purrr::map_dfr(~scale_cfr(sa_case_data, hosp_to_death_mid, .))

scaled_cfr_sa_low <- dates_num %>%
    purrr::map_dfr(~scale_cfr(sa_case_data, hosp_to_death_low, .)) %>%
    dplyr::select(date, dCFR_low = cCFR)

scaled_cfr_sa_high <- dates_num %>%
    purrr::map_dfr(~scale_cfr(sa_case_data, hosp_to_death_high, .)) %>%
    dplyr::select(date, dCFR_high = cCFR)

# the uncertainty range doesn't make sense for the first 
# row, given that it requires more data for the method
# to make sense
scaled_cfr_all <- scaled_cfr_sa_mid %>%
    dplyr::left_join(scaled_cfr_sa_low) %>%
    dplyr::left_join(scaled_cfr_sa_high) %>%
    dplyr::rename(dCFR_mid = cCFR) %>%
    dplyr::filter(dCFR_low < dCFR_mid & dCFR_high > dCFR_mid) %>%
    dplyr::select(date, nCFR, dCFR_mid, dCFR_low, dCFR_high,
                  dplyr::everything())


# just re-calculating the new cases and new deaths each day and
# dropping and NAs
scaled_cfr_sa <- scaled_cfr_all %>%
    dplyr::mutate(cases  = total_cases - dplyr::lag(total_cases),
                  deaths = total_deaths - dplyr::lag(total_deaths),
                  known_outcomes = cum_known_t - dplyr::lag(cum_known_t)) %>%
    tidyr::drop_na()


#--- plotting the dCFR results along with the case and death time-series data
cfr_subplot_11_data <- scaled_cfr_all %>%
  dplyr::mutate(date = lubridate::ymd(date)) %>%
  dplyr::filter(date >= "2020-11-01") %>%
  tidyr::pivot_longer(c(nCFR, dCFR_mid),
                      names_to  = "cfr_type",
                      values_to = "cfr")

cfr_subplot_12_data <- scaled_cfr_all %>%
  dplyr::mutate(date = lubridate::ymd(date)) %>%
  dplyr::filter(date >= "2020-11-01")

cfr_subplot_2_data <- scaled_cfr_all %>%
  dplyr::mutate(date = lubridate::ymd(date)) %>%
  tidyr::pivot_longer(c(nCFR, dCFR_mid),
                      names_to  = "cfr_type",
                      values_to = "cfr") %>%
  dplyr::filter(date >= "2020-05-01")

p1 <- ggplot2::ggplot(cfr_subplot_1_data,
                      ggplot2::aes(x = date)) +
  ggplot2::geom_line(ggplot2::aes(y = cfr, color = cfr_type)) +
  ggplot2::geom_ribbon(data = cfr_subplot_12_data,
                       ggplot2::aes(ymin = dCFR_low, ymax = dCFR_high),
                       fill = "dodgerblue", alpha = 0.2) +
  ggplot2::labs(x = "Date", y = "Case Fatality Ratio", color = "") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  viridis::scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7,
  labels = c("dCFR median", "nCFR"))

p2 <- ggplot2::ggplot(cfr_subplot_2_data,
                      ggplot2::aes(x = date)) +
  ggplot2::geom_line(ggplot2::aes(y = cfr, color = cfr_type)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = dCFR_low, ymax = dCFR_high),
                       fill = "dodgerblue", alpha = 0.2) +
  ggplot2::labs(x = "Date", y = "Case Fatality Ratio", color = "") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  viridis::scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7,
  labels = c("dCFR median", "nCFR"))



inset <- ggplot2::ggplotGrob(p1)
cfr_plot <- p2 + ggplot2::annotation_custom(grob = inset,
                                            xmin = as.Date("2020-09-01"),
                                            xmax = as.Date("2021-01-08"),
                                            ymin = 0.0425, ymax = 0.06)

#confirmed_cases_plot <- south_africa_case_data %>%
  #ggplot2::ggplot(ggplot2::aes(x = date, y = new_cases)) + 
  #ggplot2::geom_col(fill = "dodgerblue", alpha = 0.7, width = 1) +
  #ggplot2::labs(x = "Date", y = "New cases",  title = "B") 

#confirmed_deaths_plot <- south_africa_case_data %>%
  #ggplot2::ggplot(ggplot2::aes(x = date, y = new_deaths)) + 
  #ggplot2::geom_col(fill = "firebrick1", alpha = 0.6, width = 1) +
  #ggplot2::labs(x = "Date", y = "New deaths", title = "C")

#cfr_plot_final <- cfr_plot_with_inset / confirmed_cases_plot / confirmed_deaths_plot

# saving plot as png
ggplot2::ggsave(here::here("figure_2.png"), cfr_plot,
                width = 12,
                height = 8,
                dpi = 300)

# saving plot as pdf
ggplot2::ggsave(here::here("figure_2.pdf"), cfr_plot,
                width = 12,
                height = 8,
                dpi = 300)
