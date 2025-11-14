library(epidatr)
library(lubridate)
library(tidyverse)

code_folder <- here::here("Code", "Analysis")
figure_folder <- here::here("Figures", "Simulated")
source(here::here(code_folder, "helper_functions.R"))

# Load data
# git_directory <- system("git rev-parse --show-toplevel", intern = TRUE)
git_directory <- here::here()
data_folder <- here::here(git_directory, "Data", "Real_data")
variant_path <- here::here(git_directory, "Data", "Real_data", "Variants")

###### Preprocess ######
df <- read.csv(here::here(variant_path, "seq_df_us_biweekly.csv"))

# Let other = OG variant. It shouldn't come back in 2022
df$Other[df$Date >= as.Date("2022-01-01")] <- 0
df <- df %>% rename(Original = Other)

# Normalize to proportions
all_variants <- names(df)[!names(df) %in% c("State", "Date")]
df[all_variants] <- df[all_variants]/rowSums(df[all_variants])


###### Take only important variants ######
apply(df[all_variants], 2, max)
round(apply(df[all_variants], 2, mean),3)
thresh <-  0.02 # 0.01 includes Epsilon
variants <- all_variants[colMeans(df[,all_variants]) > thresh]
apply(df[variants], 2, max)
df <- df[c("Date", variants)]
df[variants] <- df[variants]/rowSums(df[variants])
n_variants <- length(variants)

head(df)
df$Date <- as.Date(df$Date)

original_start <- min(df$Date)
alpha_start <- min(df$Date[df$Alpha > 0.5])
delta_start <- min(df$Date[df$Delta > 0.5])
omicron_start <- min(df$Date[df$Omicron > 0.5])

###### Load national hospitalizations and deaths ######
jhu_deaths <- readRDS(here::here(data_folder, "JHU_finalized.RData"))
hhs_hosps <- readRDS(here::here(data_folder, "HHS_finalized.RData"))

# Subset hhs_hosps and jhu_deaths into the different variant periods established above
original_hosps <- hhs_hosps[names(hhs_hosps) < alpha_start]
alpha_hosps <- hhs_hosps[names(hhs_hosps) >= alpha_start & names(hhs_hosps) < delta_start]
delta_hosps <- hhs_hosps[names(hhs_hosps) >= delta_start & names(hhs_hosps) < omicron_start]
omicron_hosps <- hhs_hosps[names(hhs_hosps) >= omicron_start]
original_deaths <- jhu_deaths[names(jhu_deaths) < alpha_start]
alpha_deaths <- jhu_deaths[names(jhu_deaths) >= alpha_start & names(jhu_deaths) < delta_start]
delta_deaths <- jhu_deaths[names(jhu_deaths) >= delta_start & names(jhu_deaths) < omicron_start]
omicron_deaths <- jhu_deaths[names(jhu_deaths) >= omicron_start]

# Run compute_optimal_lag(HOSPS, DEATHS) for each variant period
original_lag <- compute_optimal_lag(original_hosps, original_deaths)
alpha_lag <- max(compute_optimal_lag(alpha_hosps, alpha_deaths), 5)
delta_lag <- compute_optimal_lag(delta_hosps, delta_deaths)
omicron_lag <- compute_optimal_lag(omicron_hosps, omicron_deaths)

head(df)
variant_lags <- data.frame(
  Variant = c("Original", "Alpha", "Delta", "Omicron"),
  Lag_days = c(original_lag, alpha_lag, delta_lag, omicron_lag)
)
# Compute weighted average lag for each timestep, weighting by proportion in circulation (in df)
df_lags <- df %>%
  filter(Date >= min(names(hhs_hosps)) & Date <= max(names(hhs_hosps))) %>%
  rowwise() %>%
  mutate(
    Lag_days = sum(c_across(all_of(variants)) * variant_lags$Lag_days)
  ) %>%
  select(Date, Lag_days)
# Plot over time
ggplot(df_lags, aes(x = Date, y = Lag_days)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(
    title = "Estimated Time-varying Delay from Hospitalization to Death",
    x = "Date",
    y = "Delay (days)"
  ) +
  theme_minimal()

# Alter df_lags so it's daily (smooth it)
all_dates <- data.frame(Date = seq(min(df_lags$Date), max(df_lags$Date), by = "day"))
df_lags <- all_dates %>%
  left_join(df_lags, by = "Date") %>%
  arrange(Date) %>%
  mutate(
    Lag_days = zoo::na.approx(Lag_days, rule = 2)
  )

# Add column for SD. it's the lag multiplied by 0.9 except for when lag > 30, in which case it multiplies by 0.8
df_lags <- df_lags %>%
  mutate(
    Lag_sd = ifelse(Lag_days > 30, Lag_days * 0.8, Lag_days * 0.9)
  )
# Save
head(df_lags)
saveRDS(df_lags, here::here(data_folder, "time_varying_lags.rds"))

##############

# Run make_delay_distr at all timepoints, where the first argument is the mean from df_lags, and the second arg is the sd from df_lags, and the third argument is just called d
d <- 75
# Remove the variable names from below
delay_distrs <- lapply(1:nrow(df_lags), function(i) {
  make_delay_distr(
    df_lags$Lag_days[i],
    df_lags$Lag_sd[i],
    d
  )
})
names(delay_distrs) <- df_lags$Date


compute_noiseless_deaths <- function(hosps, hfrs, delay_distr, death_dates) {
  d <- length(delay_distr) - 1
  noiseless_deaths <- sapply(death_dates, function(t) {
    trailing_dates <- seq(t-d,t, by="day")
    trailing_hosps <- hosps[names(hosps) %in% trailing_dates]
    trailing_hfrs <- hfrs[names(hfrs) %in% trailing_dates]
    exp_deaths <- sum(rev(delay_distr)*trailing_hosps*trailing_hfrs)
    exp_deaths
  })
  names(noiseless_deaths) <- death_dates
  return(noiseless_deaths)
}
# Change the above function to be able to handle time-varying delay distributions

compute_noiseless_deaths2 <- function(hosps, hfrs, delay_distrs, death_dates) {
  d <- length(delay_distrs[[1]]) - 1
  noiseless_deaths <- sapply(seq_along(death_dates), function(i) {
    t <- death_dates[i]
    delay_distr <- delay_distrs[[as.character(t)]]
    trailing_dates <- seq(t-d,t, by="day")
    trailing_hosps <- hosps[names(hosps) %in% trailing_dates]
    trailing_hfrs <- hfrs[names(hfrs) %in% trailing_dates]
    exp_deaths <- sum(rev(delay_distr)*trailing_hosps*trailing_hfrs)
    exp_deaths
  })
  names(noiseless_deaths) <- death_dates
  return(noiseless_deaths)
}
# # Example usage
# hfrs <- rep(0.15, length(hhs_hosps))
# names(hfrs) <- as.Date(names(hhs_hosps))
# death_dates <- as.Date(names(jhu_deaths))[names(jhu_deaths) >= min(names(hhs_hosps))]
# noiseless_deaths <- compute_noiseless_deaths2(hhs_hosps, hfrs, delay_distrs, death_dates)
# plot(as.Date(names(noiseless_deaths)), noiseless_deaths, type="l", col="blue", lwd=2,
#      xlab="Date", ylab="Noiseless Deaths",
#      main="Noiseless Deaths with Time-varying Delay Distribution")

hfrs_nhcs <- readRDS(here::here(data_folder, "HFRs_NHCS_rescaled.RData"))
hfr_dates <- as.Date(names(hfrs_nhcs))
# Subset hfrs_nhcs to start at same date as hhs_hosps
hfrs <- hfrs_nhcs[hfr_dates >= min(names(hhs_hosps)) & hfr_dates <= max(names(hhs_hosps))]
names(hfrs) <- hfr_dates[hfr_dates >= min(names(hhs_hosps)) & hfr_dates <= max(names(hhs_hosps))]

# Subset hhs_hosps and jhu_deaths to be within same dates as hfr_dates
hosps <- hhs_hosps[names(hhs_hosps) >= min(hfr_dates) & names(hhs_hosps) <= max(hfr_dates)]
hosp_and_hfr_dates <- as.Date(names(hosps))
deaths <- jhu_deaths[names(jhu_deaths) >= min(hosp_and_hfr_dates)+d & names(jhu_deaths) <= max(hosp_and_hfr_dates)]
death_dates <- as.Date(names(deaths))

noiseless_deaths_tv_delay <- compute_noiseless_deaths2(
  hosps,
  hfrs,
  delay_distrs[names(delay_distrs) %in% hosp_and_hfr_dates],
  death_dates
)


first_est_idx <- d+7 #d+1
est_dates <- hosp_and_hfr_dates[seq(first_est_idx, length(hosp_and_hfr_dates), by=7)]
gt_hfrs <- hfrs[names(hfrs) %in% est_dates]


# Modify this to have time-varying delay_shape (delay_distrs)
compute_conv_hfrs2 <- function(hosps, deaths, dates, delay_distrs, w=1) {
  # Only a real-time estimator. Not implemented since the beginning.
  ahfrs <- sapply(seq_along(dates), function(i) {
    date <- dates[i]
    delay_shape <- delay_distrs[[as.character(date)]]
    compute_conv_hfr(date, hosps=hosps, deaths=deaths, w=w, delay_shape=delay_shape)
  })
  names(ahfrs) <- dates
  return(ahfrs)
}

compute_lagged_hfrs2 <- function(hosps, deaths, lags, w=1, 
                                real_time=TRUE, dates=NULL, centered=NA) {
  if (is.null(dates)) {
    dates <- names(deaths)
  }
  lagged_hfrs <- sapply(seq_along(dates), function(i) {
    date <- dates[i]
    l <- lags[[as.character(date)]]
    compute_lagged_hfr(t=date, l=l, w=w, real_time=real_time, hosps=hosps,deaths=deaths, centered=centered)
  })
  names(lagged_hfrs) <- dates
  return(lagged_hfrs)
}

conv_hfrs_ws <- compute_conv_hfrs2(
  hosps,
  noiseless_deaths_tv_delay,
  est_dates,
  delay_distrs[names(delay_distrs) %in% est_dates]
)

# Repeat the below code with median lag
median_lag <- median(df_lags$Lag_days)
# Recalculate delay distributions with THIS median lag (same SD calculation as above)
delay_distr_median_lag <- make_delay_distr(
  median_lag,
  ifelse(median_lag > 30, median_lag * 0.8, median_lag * 0.9),
  d
)
conv_hfrs_median_delay <- compute_conv_hfrs(  
  hosps,
  noiseless_deaths_tv_delay,
  est_dates,
  delay_distr_median_lag,
  w=7
)

# Calculate mean lag in df_lags
mean_lag <- mean(df_lags$Lag_days)
# Recalculate delay distributions with THIS mean lag (same SD calculation as above)
delay_distr_mean_lag <- make_delay_distr(
  mean_lag,
  ifelse(mean_lag > 30, mean_lag * 0.8, mean_lag * 0.9),
  d
)
conv_hfrs_fixed_delay <- compute_conv_hfrs(
  hosps,
  noiseless_deaths_tv_delay,
  est_dates,
  delay_distr_mean_lag,
  w=7
)

lags <- df_lags %>% filter(Date %in% est_dates) %>% pull(Lag_days)
names(lags) <- est_dates

lagged_hfrs_tv <- compute_lagged_hfrs2(
  hosps=hosps,
  deaths=noiseless_deaths_tv_delay,
  lags=lags,
  w=7,
  real_time=TRUE,
  dates=est_dates,
  centered=NA
)
lagged_hfrs_tv
lagged_hfrs_fixed_lag <- compute_lagged_hfrs(
  hosps=hosps,
  deaths=noiseless_deaths_tv_delay,
  l=mean_lag,
  w=7,
  real_time=TRUE,
  dates=est_dates,
  centered=NA
)
lagged_hfrs_fixed_lag


# plot(est_dates, gt_hfrs, type="l", col="black", lwd=3,
#      xlab="Date", ylab="HFR", ylim=c(0.1,0.35),
#      main="Ratio HFRs with Fixed vs Time-varying Delay Distributions")
# lines(est_dates, conv_hfrs_ws, col="blue", lwd=2, lty=2)
# lines(est_dates, conv_hfrs_fixed_delay, col="red", lwd=2, lty=2)
# # lines(est_dates, conv_hfrs_median_delay, col="green", lwd=2)
# lines(est_dates, lagged_hfrs_tv, col="brown", lwd=2, lty=3)
# lines(est_dates, lagged_hfrs_fixed_lag, col="green", lwd=2, lty=3)
# legend("topright",
#        # legend=c("Fixed Delay (Mean Lag)", "Time-varying Delay", "Ground Truth", "Fixed Delay (Median Lag)"),
#        legend=c("GT", "Conv-TV", "Conv-Fixed", "Lagged-TV", "Lagged-Fixed"),
#        col=c("black", "red", "blue", "brown", "green"),
#        lty=c(1, 2,2,3,3), lwd=c(3, 2,2,2,2))









library(ggplot2)
library(tidyverse)

# Create data frame from the computed HFR estimates
df_plot <- data.frame(
  Date = rep(est_dates, 5),
  HFR = c(gt_hfrs, 
          conv_hfrs_ws, 
          conv_hfrs_fixed_delay, 
          lagged_hfrs_tv, 
          lagged_hfrs_fixed_lag) * 100,  # Convert to percentage
  Method = c(rep("External estimate", length(est_dates)),
             rep("Convolutional ratio", length(est_dates)),
             rep("Convolutional ratio", length(est_dates)),
             rep("Lagged ratio", length(est_dates)),
             rep("Lagged ratio", length(est_dates))),
  Delay = c(rep(NA, length(est_dates)),
            rep("Time-varying", length(est_dates)),
            rep("Fixed delay", length(est_dates)),
            rep("Time-varying", length(est_dates)),
            rep("Fixed delay", length(est_dates)))
)

# Convert to factors with desired order
df_plot <- df_plot %>%
  mutate(
    Method = factor(Method, 
                    levels = c("External estimate",
                               "Convolutional ratio",
                               "Lagged ratio")),
    Delay = factor(Delay,
                   levels = c("Time-varying", "Fixed delay"))
  )

# Create plot
p <- ggplot(data = df_plot, aes(x = Date, y = HFR, color = Method, linetype = Delay)) +
  geom_line(linewidth = 0.8) +
  ggtitle("Time-varying delay distribution") +
  xlab("") +
  ylab("HFR") +
  scale_x_date(breaks = seq(as.Date("2021-01-01"), max(df_plot$Date), by = "3 months"),
               date_labels = "%b %Y") +
  scale_y_continuous(limits = c(8, 40), expand = c(0, 0),  oob = scales::squish) +
  scale_color_manual(values = c("External estimate" = "red",
                                "Convolutional ratio" = "green3",
                                "Lagged ratio" = "dodgerblue")) +
  scale_linetype_manual(values = c("Time-varying" = "dashed",
                                   "Fixed delay" = "dotted"),
                        breaks = c("Time-varying", "Fixed delay"),
                        na.value = "solid") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -10),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 17),
    # axis.text.x = element_text(angle = 45, hjust = 1)
  )
p

# Save

ggsave(file.path(figure_folder, "time_varying_comparison.pdf"),
       plot = p, device = "pdf", width = 10, height = 5, dpi=300)

