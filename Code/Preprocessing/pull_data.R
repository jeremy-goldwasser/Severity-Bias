git_directory <- system("git rev-parse --show-toplevel", intern = TRUE)
data_folder <- file.path(git_directory, "Data", "National_Data")
code_folder <- file.path(git_directory, "Code", "Analysis")
source(file.path(code_folder, "helper_functions.R"))
library(epidatr)

############### Load and save finalized counts ###############
jhu_df <- data.frame(pub_covidcast(source="jhu-csse", signals="deaths_7dav_incidence_num",
                                   geo_type = "nation", geo_values="*", 
                                   time_type="day")) 
max_jhu_date <- as.Date("2023-03-01") # Not reliable after this
idx <- (jhu_df$time_value <= max_jhu_date)
jhu_deaths <- jhu_df$value[idx]
death_dates <- jhu_df$time_value[idx]
names(jhu_deaths) <- death_dates

saveRDS(jhu_deaths, file.path(data_folder, "JHU_finalized.RData"))

hhs_df <- data.frame(pub_covidcast(source="hhs",signals= "confirmed_admissions_covid_1d_7dav", 
                                   geo_type = "nation", geo_values="*", 
                                   time_type="day"))
first_hosp_date <- as.Date("2020-07-28")
hosp_idx <- (hhs_df$time_value >= first_hosp_date) & (hhs_df$time_value <= max_jhu_date)
hosp_dates <- hhs_df$time_value[hosp_idx]
hhs_hosps <- hhs_df$value[hosp_idx]; names(hhs_hosps) <- hosp_dates

saveRDS(hhs_hosps, file.path(data_folder, "HHS_finalized.RData"))


############### Load and save real-time counts ###############

nhcs_hfrs_all <- readRDS(file.path(data_folder, "hfrs_rescaled_v2.RData"))
nhcs_hfr_dates <- as.Date(names(nhcs_hfrs_all))
d <- 75
est_weeks <- seq(first_hosp_date+d+14, min(max(hosp_dates)-4*7, max(nhcs_hfr_dates)), by=7)
n_days_after <- 2
week_st_idx <- 5
n_ests <- length(est_weeks)
est_weeks_realtime <- est_weeks[week_st_idx:n_ests]
### HHS

hhs_hosps_realtime <- c()
for (i in 1:length(est_weeks_realtime)) { #
  t <- est_weeks_realtime[i]
  if (i==1) {
    dates <- "*"
  } else {
    dates <- seq(t-6, t, by="day")
  }
  hhs_df <- data.frame(pub_covidcast(source="hhs",signals= "confirmed_admissions_covid_1d", #_7dav
                                   geo_type = "nation", geo_values="*",
                                   time_type="day", as_of=t+n_days_after, time_values=dates))
  new_hosps <- hhs_df$value; names(new_hosps) <- as.Date(hhs_df$time_value)

  # Account for bad delay distr
  if (i==1) {
    new_hosps <- new_hosps[(names(new_hosps) >= first_hosp_date) & names(new_hosps) <= t]
    dates <- seq(first_hosp_date, t, by="day")
  }

  # Extend if most recent values are unavailable
  n_new_hosps <- length(new_hosps)
  n_desired_hosps <- t - as.Date(min(names(new_hosps))) + 1
  if (n_new_hosps < n_desired_hosps) {
    new_hosps <- c(new_hosps, rep(new_hosps[length(new_hosps)], n_desired_hosps-n_new_hosps))
    names(new_hosps) <- dates
  }

  hhs_hosps_realtime <- c(hhs_hosps_realtime, new_hosps)

}
hhs_hosps_realtime <- seven_day_smoothing(hhs_hosps_realtime, centered=FALSE)
hosp_dates_realtime <- as.Date(names(hhs_hosps_realtime))
plot(hosp_dates_realtime, hhs_hosps_realtime)
saveRDS(hhs_hosps_realtime, file.path(data_folder, "HHS_real_time.RData"))


### JHU
jhu_deaths_realtime <- c()
for (i in 1:length(est_weeks_realtime)) {
  t <- est_weeks_realtime[i]
  if (i==1) {
    dates <- "*"
  } else {
    dates <- seq(t-6, t, by="day")
  }
  jhu_df <- data.frame(pub_covidcast(source="jhu-csse", signals="deaths_incidence_num",#7dav_
                                    geo_type = "nation", geo_values="*",
                                    time_type="day", as_of=t+n_days_after, time_values=dates))
  new_deaths <- jhu_df$value; names(new_deaths) <- jhu_df$time_value
  if (i==1) {
    new_deaths <- new_deaths[as.Date(names(new_deaths)) <= t]
    dates <- seq(as.Date(min(jhu_df$time_value)), t, by="day")
  }

  # Extend if most recent values are unavailable
  n_new_deaths <- length(new_deaths)
  n_desired_deaths <- t - as.Date(min(names(new_deaths))) + 1
  if (n_new_deaths < n_desired_deaths) {
    new_deaths <- c(new_deaths, rep(new_deaths[length(new_deaths)], n_desired_deaths-length(new_deaths)))

    names(new_deaths) <- dates
  }
  jhu_deaths_realtime <- c(jhu_deaths_realtime, new_deaths)
}
jhu_deaths_realtime <- seven_day_smoothing(jhu_deaths_realtime, centered=FALSE)
plot(as.Date(names(jhu_deaths_realtime)), jhu_deaths_realtime)
saveRDS(jhu_deaths_realtime, file.path(data_folder, "JHU_real_time.RData"))
