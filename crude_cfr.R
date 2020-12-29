case_data_function <- function() {

    owid_path <- "https://covid.ourworldindata.org/data/owid-covid-data.csv"
    case_data <- read.csv(owid_path) %>%
    dplyr::tibble() %>%
    dplyr::mutate(date = lubridate::ymd(date)) %>%
    dplyr::select(iso_code,
                  country = location,
                  date,
                  new_cases,
                  new_deaths,
                  population)

  return(case_data)
}

case_data <- case_data_function()

sa_case_data <- case_data %>%
  dplyr::filter(iso_code == "ZAF") %>%
  tidyr::drop_na()

sa_case_data_reduced <- sa_case_data %>%
  dplyr::mutate(cases_ma_7  = zoo::rollmean(new_cases, 7,  fill = NA)) %>%
  dplyr::select(date, new_cases, cases_ma_7)

sa_death_data_reduced <- sa_case_data %>%
  dplyr::mutate(date = date) %>%
  dplyr::mutate(deaths_ma_7 = zoo::rollmean(new_deaths, 7, fill = NA)) %>%
  dplyr::select(date, new_deaths, deaths_ma_7)

sa_death_data_lagged_reduced <- sa_case_data %>%
  dplyr::mutate(date = date - 21) %>%
  dplyr::mutate(deaths_ma_7 = zoo::rollmean(new_deaths, 7, fill = NA)) %>%
  dplyr::select(date, deaths_ma_7)

cfr_crude <- sa_case_data_reduced %>%
  dplyr::left_join(sa_death_data_reduced) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(cfr_lagged_crude = deaths_ma_7 / cases_ma_7)

cfr_lagged_crude <- sa_case_data_reduced %>%
  dplyr::left_join(sa_death_data_lagged_reduced) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(cfr_lagged_crude = deaths_ma_7 / cases_ma_7)

library(ggplot2)

p1 <- cfr_crude %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = date, y = cases_ma_7,
                         colour = "dodgerblue")) +
    ggplot2::geom_line(ggplot2::aes(x = date, y = deaths_ma_7,
                         colour = "firebrick1")) +
    ggplot2::geom_line(ggplot2::aes(x = date, y = new_cases),
            colour = "black",
            linetype = "dashed",
            alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(x = date, y = new_deaths),
                       colour = "black",
                       linetype = "dashed",
                       alpha = 0.3) +
    ggplot2::labs(x = "Date",
                  y = "Incidence",
                  colour = "",
                  title = "A") +
    ggplot2::scale_colour_manual(values = c("dodgerblue",
                                            "firebrick1",
                                            "black"),
                               labels = c("7-day moving average of cases",
                                          "7-day moving average of deaths",
                                          "Raw data"))

p2 <- ggplot2::ggplot() +
      ggplot2::geom_line(data = cfr_lagged_crude,
                         ggplot2::aes(x = date,
                                      y = cfr_lagged_crude,
                                      colour = "dodgerblue"),
                         alpha = 0.7) +
      ggplot2::geom_line(data = cfr_crude,
                         ggplot2::aes(x = date,
                                      y = cfr_lagged_crude,
                                      colour = "firebrick1"),
                         alpha = 0.7) +
      ggplot2::labs(x = "Date",
                    y = "CFR",
                    colour = "",
                    title = "B") +
     ggplot2::scale_y_continuous(labels = scales::percent) +
     ggplot2::scale_colour_manual(values = c("dodgerblue",
                                             "firebrick1"),
                               labels = c("Deaths lagged by 21 days",
                                          "nCFR"))

library(patchwork)

plot_together <- p1 / p2
plot(plot_together)

ggplot2::ggsave("figure_2_crude.png",
                plot_together,
                width = 10,
                height = 10)


ggplot2::ggsave("figure_2_crude.pdf",
                plot_together,
                width = 10,
                height = 10)
