data_folder <- here::here("Data", "Real_data")

dat <- read.csv(here::here(data_folder, "COVID-19_Hospital_Data_from_the_National_Hospital_Care_Survey.csv"))
library(tidyverse)
library(stats)

colnames(dat)
unique(dat$Measure) #Confirmed COVID-19 Deaths, Confirmed COVID-19 

set.seed(1)

df <- filter(dat, Measure=="Percent", 
             Indicator=="Confirmed COVID-19 Deaths",
             Group=="Total",
) %>% select(Start_Time, Value, Setting)
df$Start_Time <- as.Date(df$Start_Time, "%m/%d/%Y")

# Visualize HFRs of inpatient versus emergency department
ggplot(df, aes(x=Start_Time, group=Setting, color=Setting)) +
  geom_line(aes(y=Value)) + xlab("Week") + ylab("HFR")

# Only look at inpatient
df2 <- filter(df, Setting=="IP") %>% 
  select(-"Setting") %>%
  rename("HFR"="Value", "Week"="Start_Time")
df2$HFR <- df2$HFR/100
weeks <- df2$Week

raw_weekly_hfrs <- df2$HFR; names(raw_weekly_hfrs) <- weeks
saveRDS(raw_weekly_hfrs, here::here(data_folder, "HFRs_NHCS_raw.RData"))

####### Smooth with spline #######
spline_model <- smooth.spline(as.numeric(weeks), df2$HFR)
smoothed_weekly_hfrs <- predict(spline_model, as.numeric(weeks))$y
names(smoothed_weekly_hfrs) <- weeks
plot(weeks, df2$HFR, type="l", xlab="Date", ylab="HFR", main="Raw and smoothed NHCS HFRs")
lines(weeks, smoothed_weekly_hfrs, col="red")

# Convert weekly to daily
days <- seq(min(weeks), max(weeks)+6, by="day")
smoothed_daily_hfrs <- spline(weeks, smoothed_weekly_hfrs, n=length(days))$y
names(smoothed_daily_hfrs) <- days

####### Save #######
saveRDS(smoothed_daily_hfrs, here::here(data_folder, "HFRs_NHCS_smoothed.RData"))

