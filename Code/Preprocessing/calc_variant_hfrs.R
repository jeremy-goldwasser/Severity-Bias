library(epidatr)
library(lubridate)
library(tidyverse)

# Load data
# git_directory <- system("git rev-parse --show-toplevel", intern = TRUE)
git_directory <- here::here()
nat_data_path <- here::here(git_directory, "Data", "Real_data")
variant_path <- here::here(git_directory, "Data", "Real_data", "Variants")

###### Preprocess ######
df <- read.csv(here::here(variant_path, "seq_df_us_biweekly.csv"))

# Let other = OG variant. It shouldn't come back in 2022
df$Other[df$Date >= as.Date("2022-01-01")] <- 0
df <- df %>% rename(Original = Other)

# Normalize to proportions
all_variants <- names(df)[!names(df) %in% c("State", "Date")]
df[all_variants] <- df[all_variants]/rowSums(df[all_variants])
head(df)



###### Take only important variants ######
apply(df[all_variants], 2, max)
round(apply(df[all_variants], 2, mean),3)
thresh <-  0.02 # 0.01 includes Epsilon
variants <- all_variants[colMeans(df[,all_variants]) > thresh]
apply(df[variants], 2, max)
df <- df[c("Date", variants)]
df[variants] <- df[variants]/rowSums(df[variants])
n_variants <- length(variants)


###### Plot variant proportions over time ######
df <- data.frame(df)
df$Date <- as.Date(df$Date)

dfLong <- df %>%
  pivot_longer(cols = -Date, names_to = "Variant", values_to = "Prop")
dfLong %>%
  filter(Date >= as.Date("2020-06-01")) %>%
  filter(Date <= as.Date("2022-06-01")) %>%
  ggplot(aes(x=Date, y=Prop, color=Variant, linetype=Variant, group=Variant)) + geom_line() +
    labs(title = "US COVID Proportions by Variant",
         x = "Date",
         y = "Proportion of Cases")  +
    # scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + 
    theme_minimal()

# Get dates for each wave
original_start <- min(df$Date)
alpha_start <- min(df$Date[df$Alpha > 0.5])
delta_start <- min(df$Date[df$Delta > 0.5])
omicron_start <- min(df$Date[df$Omicron > 0.5])

###### Load national hospitalizations and deaths ######

# Deaths in the period from Sunday to the following Saturday
nchs_df <- data.frame(pub_covidcast(source="nchs-mortality",signals= "deaths_covid_incidence_num", 
                                   geo_type = "nation", geo_values="*", 
                                   time_type="week"))

# Shift, to use end-of-week
nchs_df$time_value <- nchs_df$time_value + 6

# Load hospitalizations
hhs_df <- data.frame(pub_covidcast(source="hhs",signals= "confirmed_admissions_covid_1d", 
                                   geo_type = "nation", geo_values="*", 
                                   time_type="day"))
# plot(hhs_df$time_value, hhs_df$value)
# hhs_df$time_value[hhs_df$value > 0][1:50]
# hhs_df$value[hhs_df$value > 0][1:50]
first_hosp_date <- as.Date("2020-07-14") # Basically all 0 before this
hhs_df <- hhs_df[hhs_df$time_value >= first_hosp_date, c("time_value", "value")]

# Filter dates - includes last hosps
death_sats <- nchs_df$time_value
last_hosp_date <- max(hhs_df$time_value)
sats_with_hosps <- death_sats[(death_sats >= first_hosp_date+6) & (death_sats <= last_hosp_date)]

weekly_hosps <- sapply(sats_with_hosps, function(day) {
  day <- as.Date(day)
  week <- hhs_df[hhs_df$time_value %in% seq(day-6, day, "day"), ]
  sum(week$value)
}); names(weekly_hosps) <- sats_with_hosps
plot(sats_with_hosps, weekly_hosps)

weekly_deaths <- nchs_df$value[nchs_df$time_value %in% sats_with_hosps]
names(weekly_deaths) <- sats_with_hosps


mean(weekly_deaths)/mean(weekly_hosps)

###### Calculate HFR based on dominant variant periods ######
dates <- sats_with_hosps[sats_with_hosps<=(as.Date("2023-03-11")+7)]

start_date <- max(original_start, min(dates))
calc_naive_hfr <- function(date1, date2, offset) {
  deaths <- sum(weekly_deaths[dates >= (date1+offset) & dates < (date2+offset)])
  hosps <- sum(weekly_hosps[dates >= date1 & dates < date2])
  hfr <- deaths/hosps
  return(hfr)
}
calc_naive_hfrs <- function(offset) {
  og_hfr <- calc_naive_hfr(start_date, alpha_start, offset)
  alpha_hfr <- calc_naive_hfr(alpha_start, delta_start, offset)
  delta_hfr <- calc_naive_hfr(delta_start, omicron_start, offset)
  omicron_hfr <- calc_naive_hfr(omicron_start, max(dates), offset)
  
  early_omicron_hfr <- calc_naive_hfr(omicron_start, late_omicron_start, offset)
  late_omicron_hfr <- calc_naive_hfr(late_omicron_start, max(dates), offset)
  
  variant_hfrs <- c(og_hfr, alpha_hfr, delta_hfr, omicron_hfr, early_omicron_hfr, late_omicron_hfr)
  names(variant_hfrs) <- c("Original", "Alpha", "Delta", "Omicron", "Early Omicron", "Late Omicron")
  return(variant_hfrs)
}
late_omicron_start <- as.Date("2022-04-01")
variant_hfrs0 <- calc_naive_hfrs(offset=0)
variant_hfrs14 <- calc_naive_hfrs(offset=14)
round(variant_hfrs0*100, 1)
round(variant_hfrs14*100, 1)
#         Original    Alpha   Delta   Omicron   Early_Omicron     Late_Omicron
# No lag: 17.9%       10.1%   16.8%   11.6%     16.2%             8.4%
# 2-wks:  17.6%       9.2%   18.1%   10.8%     14.4%             8.3%
variant_starts <- c(start_date, alpha_start, delta_start, omicron_start, late_omicron_start, max(dates))
names(variant_starts) <- c("Start", "Alpha", "Delta", "Omicron", "Late Omicron", "End")
df <- df[,c(1,5,2,3,4)] # Reorder
variant_list <- list(variant_hfrs0, variant_hfrs14, variant_starts, df)
names(variant_list) <- c("HFRs-0 lag", "HFRs-2wk lag", "Variant Dates", "Proportions")

saveRDS(variant_list, here::here(nat_data_path, 'variant_hfrs.RData'))
