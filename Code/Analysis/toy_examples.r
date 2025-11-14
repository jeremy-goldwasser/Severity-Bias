###############################################
## Unified Simulation Figures Script
## Produces 4 figures consistently formatted
###############################################

library(tidyverse)
library(here)
library(patchwork)

code_folder   <- here::here("Code", "Analysis")
data_folder   <- here::here("Data", "Real_data")
figure_folder <- here::here("Figures", "Simulated")
dir.create(figure_folder, showWarnings = FALSE, recursive = TRUE)

source(here::here(code_folder, "helper_functions.R"))

## ---------------------------------------------------------------------
## Load data
## ---------------------------------------------------------------------
hhs_hosps  <- readRDS(here::here(data_folder, "HHS_finalized.RData"))
jhu_deaths <- readRDS(here::here(data_folder, "JHU_finalized.RData"))

hosp_dates <- as.Date(names(hhs_hosps))

###############################################
## COMMON THEME
###############################################


theme_consistent <- theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    
    # bring legend closer
    legend.spacing.y = unit(-4, "pt"),
    legend.margin = margin(t = -2, b = 2),
    
    plot.title = element_text(size = 18),
    
    # KEEP x-axis ticks, REMOVE only the label
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 12),
    axis.ticks.x = element_line(),
    
    axis.title.y = element_text(size = 16),
    
    panel.grid.minor = element_blank()
  )

###############################################
## FIGURE 1 — sim_onehot
###############################################

hfrs_nhcs <- readRDS(here::here(data_folder, "HFRs_NHCS_rescaled.RData"))
nhcs_dates <- as.Date(names(hfrs_nhcs))

short_lag <- 14
long_lag  <- 28

estDates <- nhcs_dates[(long_lag + 1):length(nhcs_dates)]
gt_hfr   <- hfrs_nhcs[nhcs_dates %in% estDates]

est_hfr_short <- hfrs_nhcs[nhcs_dates %in% (estDates - short_lag)]
est_hfr_long  <- hfrs_nhcs[nhcs_dates %in% (estDates - long_lag)]

df_onehot <- tibble(
  Date   = rep(estDates, 3),
  HFR    = c(gt_hfr, est_hfr_short, est_hfr_long) * 100,
  Method = factor(rep(c("True HFR","Short delay","Long delay"), each = length(estDates)),
                  levels = c("True HFR","Short delay","Long delay"))
)

p_onehot <- ggplot(df_onehot,
                   aes(Date, HFR, color = Method, linetype = Method)) +
  geom_line() +
  labs(title = "One-hot delay: Bias from improper alignment",
       y = "HFR (%)", x = "") +
  scale_colour_manual(values = c(
    "True HFR"        = "red",
    "Short delay"       = "darkorange",
    "Long delay"     = "brown"
  )) +
  theme_consistent
p_onehot
ggsave(file.path(figure_folder, "sim_onehot.pdf"), p_onehot, width = 7, height = 5)

###############################################
## FIGURE 2 — sim_chging_primary
###############################################

n_dates <- 300
sim_hfrs <- seq(0.5, 0, length.out=n_dates)

sigmoid <- function(x) 1/(1+exp(-x))
xs <- seq(-9, 7, length.out=n_dates/2)

hosps <- (c(sigmoid(xs), rev(sigmoid(xs))) * 9000 + 1000)

# two-day delay
days <- c(0, 10)
two_day_d <- max(days)
two_day_delay <- rep(0, two_day_d+1)
two_day_delay[days+1] <- 0.5

est_indices <- (two_day_d+1):n_dates
n_est <- length(est_indices)

exp_deaths <- sapply(est_indices, function(t) {
  idx <- (t-two_day_d):t
  sum(hosps[idx] * sim_hfrs[idx] * rev(two_day_delay))
})

conv_hfrs <- sapply(est_indices, function(t) {
  idx <- (t-two_day_d):t
  sum(hosps[idx] * sim_hfrs[idx] * rev(two_day_delay)) /
    sum(hosps[idx] * rev(two_day_delay))
})

lag <- mean(days)
lagged_hfrs <- exp_deaths / hosps[est_indices - lag]

gt_hfrs <- sim_hfrs[est_indices]

conv_bias <- conv_hfrs - gt_hfrs
lagged_bias <- lagged_hfrs - gt_hfrs

bias_scale <- max(abs(c(conv_bias, lagged_bias)))
conv_bias_scaled <- conv_bias / bias_scale * 0.2
lagged_bias_scaled <- lagged_bias / bias_scale * 0.2

hosp_scaled <- hosps[est_indices]
hosp_scaled <- (hosp_scaled - min(hosp_scaled)) /
  (max(hosp_scaled) - min(hosp_scaled))
hosp_scaled <- hosp_scaled * (max(gt_hfrs) - min(gt_hfrs)) + min(gt_hfrs)


df <- data.frame(
  Date = rep(1:n_est, 4),
  Value = c(
    gt_hfrs * 100,
    hosp_scaled * 100,
    (lagged_bias_scaled + min(gt_hfrs)) * 100,
    (conv_bias_scaled   + min(gt_hfrs)) * 100
  ),
  Line = factor(
    rep(c("True HFR","Hospitalizations","Lagged Bias","Conv Bias"), each=n_est),
    levels=c("True HFR","Hospitalizations","Lagged Bias","Conv Bias")
  )
)

pChg <- ggplot(df, aes(x = Date, y = Value, colour = Line)) +
  geom_line() +
  labs(
    # title = "Changing Primary Incidence. True HFR and Estimation Bias.",
    title = "Linear HFR: Changing primary and estimation bias",
    x = "Date", 
    y = "HFR (%)"
  ) +
  scale_colour_manual(values = c(
    "True HFR"        = "red",
    "Hospitalizations" = "black",
    "Lagged Bias"     = "dodgerblue",
    "Conv Bias"       = "green3"
  )) +
  theme_consistent

pChg
ggsave(file.path(figure_folder, "sim_chging_primary.pdf"), pChg, width=7, height=5)

###############################################
## FIGURE 3 — Sinusoidal HFR
###############################################

oracle_lag <- compute_optimal_lag(hhs_hosps, jhu_deaths)

d <- 75
delay_distr <- make_delay_distr(
  Mean = oracle_lag,
  Sd   = round(0.9 * oracle_lag),
  d    = d
)

cfr_dates <- hosp_dates
n_dates   <- length(cfr_dates)
t_idx     <- seq_len(n_dates)

gt_hfrs <- 0.10 + 0.04 * sin(2*pi*t_idx / 180)
names(gt_hfrs) <- as.character(cfr_dates)

death_dates_for_hfr <- cfr_dates[(d+1):n_dates]
sim_deaths_const <- compute_noiseless_deaths(
  hhs_hosps, gt_hfrs, delay_distr, death_dates_for_hfr
)

est_dates <- seq(death_dates_for_hfr[1],
                 tail(death_dates_for_hfr, 1),
                 by = "7 days")

gt_hfrs_est <- gt_hfrs[cfr_dates %in% est_dates] * 100
hfr_conv <- compute_conv_hfrs(hhs_hosps, sim_deaths_const, est_dates, delay_distr, 1) * 100
hfr_lag  <- compute_lagged_hfrs(hhs_hosps, sim_deaths_const, l=oracle_lag,
                                w=1, dates=est_dates, real_time=TRUE) * 100

df_sine <- tibble(
  Date   = rep(est_dates, 3),
  HFR    = c(gt_hfrs_est, hfr_conv, hfr_lag),
  Method = factor(rep(c("True HFR","Convolutional ratio","Lagged ratio"), each = length(est_dates)),
                  levels = c("True HFR","Convolutional ratio","Lagged ratio"))
)

p_sine <- ggplot(df_sine,
                 aes(Date, HFR, color = Method, linetype = Method)) +
  geom_line() +
  labs(title = "Sinusoidal HFR: Convolutional vs lagged ratio",
       y = "HFR (%)", x = "") +
  theme_consistent

p_sine
ggsave(file.path(figure_folder, "sinusoidal_hfr_lag_vs_conv.pdf"),
       p_sine, width = 7, height = 5)

###############################################
## FIGURE 4 — Step-function HFR
###############################################

set.seed(123)
n_steps  <- 6
step_len <- floor(n_dates / n_steps)
step_vals <- runif(n_steps, 0.06, 0.14)

gt_step <- rep(step_vals, each = step_len)[seq_len(n_dates)]
names(gt_step) <- as.character(cfr_dates)

sim_deaths_step <- compute_noiseless_deaths(
  hhs_hosps, gt_step, delay_distr, death_dates_for_hfr
)

gt_step_est <- gt_step[cfr_dates %in% est_dates] * 100

hfr_conv_step <- compute_conv_hfrs(
  hhs_hosps, sim_deaths_step, est_dates, delay_distr, 1
) * 100

hfr_lag_step <- compute_lagged_hfrs(
  hhs_hosps, sim_deaths_step, l=oracle_lag, w=1, dates=est_dates, real_time=TRUE
) * 100

df_step <- tibble(
  Date   = rep(est_dates, 3),
  HFR    = c(gt_step_est, hfr_conv_step, hfr_lag_step),
  Method = factor(rep(c("True HFR","Convolutional ratio","Lagged ratio"),
                      each = length(est_dates)),
                  levels = c("True HFR","Convolutional ratio","Lagged ratio"))
)

p_step <- ggplot(df_step,
                 aes(Date, HFR, color = Method, linetype = Method)) +
  geom_line() +
  labs(title = "Step-function HFR: Convolutional vs lagged ratio",
       y = "HFR (%)", x = "") +
  theme_consistent

# save
p_step

ggsave(file.path(figure_folder, "step_hfr_lag_vs_conv.pdf"),
       p_step, width = 7, height = 5)

