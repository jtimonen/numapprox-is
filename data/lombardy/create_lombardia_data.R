# Adapted from: https://github.com/jriou/covid_adjusted_cfr

# Set-up
library(tidyverse)
library(lubridate)
library(readxl)
library(xtable)

# Italy contact matrix, 9 age classes
contact_matrix_italy <- function() {
  italy <- socialmixr::get_survey("https://doi.org/10.5281/zenodo.1043437")
  ages <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
  m <- socialmixr::contact_matrix(
    survey = italy,
    age.limits = ages,
    symmetric = TRUE
  )
  data.frame(contacts = as.numeric(m$matrix)) %>%
    tbl_df() %>%
    mutate(age1 = rep(1:9, 9), age2 = rep(1:9, each = 9)) %>%
    ggplot() +
    geom_tile(aes(x = age2, y = age1, fill = contacts))
  contact <- structure(m$matrix, dim = c(9, 9))
  return(contact)
}

# Controls
data_lombardia <- function() {
  day_start <- as.Date("2020-02-10")
  day_data <- as.Date("2020-02-11")
  day_max <- as.Date("2020-04-25")
  day_quarantine <- as.Date("2020-02-24")

  # Age distribution in lombardy for 9 age classes
  age2 <- c(0, 0, 10, 10, 20, 20, 30, 30, 40, 40, 50, 50, 60, 60, 70, 70, 80, 80, 80, 80, 80)
  age_dist <- read_excel("covid_adjusted_cfr/data/age_structure.xlsx") %>%
    filter(country == "Italy") %>%
    gather("age", "n", 3:23) %>%
    mutate(n = as.numeric(n), age2 = age2) %>%
    group_by(age2) %>%
    summarise(n = sum(n)) %>%
    pull(n)
  age_dist <- age_dist / sum(age_dist)

  # Population in Lombardia (source: Eurostat 2018)
  pop_t <- 10.04e6

  # Case incidence by date of symptoms in lombardy
  # https://www.epicentro.iss.it/coronavirus/bollettino/Bolletino-sorveglianza-integrata-COVID-19_28-aprile-2020_appendix.pdf
  fp <- "covid_adjusted_cfr/data/lombardy/lombardy_data_4may.csv"
  lombardy_data_4may <- read.csv(fp) %>%
    tbl_df() %>%
    mutate(date = ymd(paste("2020", month, day, sep = "-"))) %>%
    filter(date >= day_data, date <= day_max)
  incidence_cases <- pull(lombardy_data_4may, cases)

  # Deaths incidence in lombardy
  # (https://github.com/pcm-dpc/COVID-19/blob/master/dati-regioni/dpc-covid19-ita-regioni.csv)
  incidence_deaths <- pull(lombardy_data_4may, deaths)


  # Age distribution of cases
  # (https://www.epicentro.iss.it/coronavirus/bollettino/Bolletino-sorveglianza-integrata-COVID-19_28-aprile-2020_appendix.pdf)
  fp <- "covid_adjusted_cfr/data/lombardy/age_distributions_cases_deaths_4may.csv"
  age_distributions_cases_deaths_4may <- read.csv(fp) %>% tbl_df()

  cases_tmax <- pull(age_distributions_cases_deaths_4may, cases)
  prop_cases_tmax <- cases_tmax / sum(cases_tmax)

  mort_tmax <- pull(age_distributions_cases_deaths_4may, deaths)
  prop_mort_tmax <- mort_tmax / sum(mort_tmax)

  agedistr_cases <- cases_tmax
  agedistr_deaths <- mort_tmax

  # Return
  list(
    pop_t = pop_t,
    agedistr_cases = agedistr_cases,
    agedistr_deaths = agedistr_deaths,
    age_dist = age_dist,
    incidence_deaths = incidence_deaths,
    incidence_cases = incidence_cases,
    lombardy_data_4may = lombardy_data_4may
  )
}

# Download data and write to file
mat <- contact_matrix_italy()
dat <- data_lombardia()
dat$contacts <- mat
saveRDS(dat, file = "data_lombardia.rds")
