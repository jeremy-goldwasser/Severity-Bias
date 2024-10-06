git_directory <- system("git rev-parse --show-toplevel", intern = TRUE)
data_folder <- file.path(git_directory, "Data", "National_Data")

# data_folder <- "/Users/jeremygoldwasser/Desktop/HFR/repo/HFR/Data/National_Data"
dat <- read.csv(file.path(data_folder, "COVID-19_Hospital_Data_from_the_National_Hospital_Care_Survey.csv"))
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

head(df2)

####### Smooth with trend filtering #######
# library(genlasso)

# obj <- trendfilter(df2$HFR, ord=2)
# cv.obj <- cv.trendfilter(obj)
# plot(cv.obj)
# lambda <- cv.obj$lambda.min
# idx <- which(obj$lambda==lambda)
# idx
# smoothed_weekly_hfrs <- obj$beta[,idx]

spline_model <- smooth.spline(as.numeric(weeks), df2$HFR)
smoothed_weekly_hfrs <- predict(spline_model, as.numeric(weeks))$y
names(smoothed_weekly_hfrs) <- weeks
plot(weeks, df2$HFR, type="l")
lines(weeks, smoothed_weekly_hfrs, col="red")
# lines(weeks, smoothed_weekly_hfrs2, col="blue")


####### Convert weekly to daily with spline interpolation #######
days <- seq(min(weeks), max(weeks)+6, by="day")

smoothed_daily_hfrs <- spline(weeks, smoothed_weekly_hfrs, n=length(days))$y
names(smoothed_daily_hfrs) <- days

####### Save #######
# plot(days, smoothed_daily_hfrs, type="l")
saveRDS(smoothed_daily_hfrs, file.path(data_folder, "smoothed_hfrs.RData"))

