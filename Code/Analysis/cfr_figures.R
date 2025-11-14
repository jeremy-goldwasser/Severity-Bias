#setwd("/Users/jeremygoldwasser/Desktop/HFR/Severity-Bias")
library(epidatr)
library(tidyverse)
library(patchwork)

code_folder <- here::here("Code", "Analysis")
data_folder <- here::here("Data", "Real_data")
figure_folder <- here::here("Figures")
source(here::here(code_folder, "helper_functions.R"))

geo_type <- "nation"
geo_value <- "*"
# Pull cases from JHU
# On their own, they have significant DOW effects; weekly smoothing is the easiest way to remove them
jhu_df <- data.frame(pub_covidcast(source="jhu-csse", signals="confirmed_7dav_incidence_num",
                                   geo_type = geo_type, geo_values=geo_value, 
                                   time_type="day"))
cases7 <- jhu_df$value; names(cases7) <- jhu_df$time_value
plot(cases7)

jhu_deaths <- readRDS(here::here(data_folder, "JHU_finalized.RData"))
plot(jhu_deaths)
hfrs_nhcs <- readRDS(here::here(data_folder, "HFRs_NHCS_rescaled.RData"))
hfr_dates <- as.Date(names(hfrs_nhcs))

n_dates <- length(hfr_dates)
oracle_lag_realtime <- compute_optimal_lag(cases7, jhu_deaths, verbose=TRUE)
total_cases <- sum(cases7[names(cases7) %in% hfr_dates[1:(n_dates-oracle_lag_realtime-1)]])
total_deaths <- sum(jhu_deaths[names(jhu_deaths) %in% hfr_dates[(oracle_lag_realtime+1):n_dates]])
total_cfr <- total_deaths/total_cases
total_cfr

cfrs <- hfrs_nhcs * total_cfr/mean(hfrs_nhcs)
mean(cfrs) # Good
cfr_dates <- hfr_dates
plot(cfrs)

# d <- 75
# delay_distr_jhu <- make_delay_distr(oracle_lag_realtime, oracle_lag_realtime*0.9, d)
# w <- 7
first_est_idx <- d+w
# 
cfr_dates_weekly <- cfr_dates[seq(first_est_idx, length(cfr_dates), by=7)]
# cfrs_weekly <- cfrs[seq(first_est_idx, length(cfr_dates), by=7)]



######################### Figure 1 — CFR analogues #########################

# Shared setup
d <- 75
delay_distr_med <- make_delay_distr(Mean = 20, Sd = round(20 * 0.9), d)
n_pts <- 200
n_ests <- n_pts - d
t0 <- as.Date("2022-02-01")

case_dates <- as.Date(names(cases7))
case_start_idx <- which(case_dates == t0)

# CFR sequence (replace HFR)
cfr_window <- cfrs[cfr_dates %in% as.Date((t0 - d):(t0 - d + n_pts - 1))]
flatter_cfrs <- (cfr_window - min(cfr_window)) * 0.25 + median(cfr_window) * 0.85
est_idx <- (d + 1):n_pts
cases_window <- cases7[(case_start_idx - d):(case_start_idx + n_ests - 1)]

calc_nish_cfrs <- function(cases, cfrs, delay_distr) {
  sapply(1:n_ests, function(t_idx) {
    trailing_cases <- cases[(t_idx):(d + t_idx)]
    trailing_cfrs  <- cfrs[(t_idx):(d + t_idx)]
    exp_deaths <- sum(rev(delay_distr) * trailing_cases * trailing_cfrs)
    denom      <- sum(rev(delay_distr) * trailing_cases)
    exp_deaths / denom
  })
}

################## Example A: Changing CFR ##################
est_flat   <- calc_nish_cfrs(cases_window, flatter_cfrs, delay_distr_med)
est_chging <- calc_nish_cfrs(cases_window, cfr_window, delay_distr_med)

dfA <- data.frame(
  t = t0 + rep(1:n_ests, 4) - 1,
  CFRs = c(cfr_window[est_idx], flatter_cfrs[est_idx], est_chging, est_flat),
  CFR  = factor(rep(c("True", "Estimated"), each = n_ests * 2)),
  Variation = factor(rep(rep(c("Original", "Flatter"), each = n_ests), 2))
)

plotA <- ggplot(dfA, aes(x = t, y = CFRs * 100, color = Variation, linetype = CFR)) +
  geom_line() +
  theme_bw() +
  ggtitle("Effect of changing CFR") +
  ylab("CFR (%)") + xlab("") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  scale_color_manual(values = c("Original" = "black", "Flatter" = "blue")) +
  theme(legend.position = "bottom")
ggsave(file.path(figure_folder, "Simulated", "toy_chging_cfr.pdf"),
       plot = plotA, device = "pdf", width = 6, height = 6)


################## Example B: Delay distributions ##################
short_mean <- 12; long_mean <- 28
delay_short <- make_delay_distr(Mean = short_mean, Sd = round(short_mean * 0.9), d)
delay_long  <- make_delay_distr(Mean = long_mean,  Sd = round(long_mean  * 0.9), d)

est_short <- calc_nish_cfrs(cases_window, cfr_window, delay_short)
est_med   <- calc_nish_cfrs(cases_window, cfr_window, delay_distr_med)
est_long  <- calc_nish_cfrs(cases_window, cfr_window, delay_long)

types <- c("12 d", "20 d", "28 d")
dfB <- data.frame(
  t = t0 + rep(1:n_ests, 4) - 1,
  CFRs = c(cfr_window[est_idx], est_short, est_med, est_long),
  CFR  = factor(c(rep("True", n_ests), rep("Estimated", n_ests * 3))),
  Type = factor(c(rep("True", n_ests), rep(types, each = n_ests)), levels = c("True", types))
)

# Main plot
plotB_main <- ggplot(dfB, aes(x = t, y = CFRs * 100, linetype = CFR, color = Type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Effect of true delay distribution") +
  ylab("CFR (%)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  scale_color_manual(values = c("True" = "black", "12 d" = "green3", "20 d" = "red2", "28 d" = "blue")) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 18),
    legend.text  = element_text(size = 14),
    axis.title.y = element_text(size = 15)
  )

# Delay distribution inset
dfDistr <- data.frame(
  t = rep(0:d, 3),
  Delay = c(delay_short, delay_distr_med, delay_long),
  Type = rep(types, each = d + 1)
)

delay_inset <- ggplot(dfDistr, aes(x = t, y = Delay, color = Type)) +
  geom_line(size = 0.9) +
  theme_bw() +
  scale_color_manual(
    breaks = types,
    values = c("12 d" = "green3", "20 d" = "red2", "28 d" = "blue")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin = unit(rep(1, 4), "mm")
  ) +
  xlab("Days after case onset") +
  ylab("") +
  ggtitle("Delay distribution")

# Combine with inset (upper-left corner)
plotB <- plotB_main +
  inset_element(
    delay_inset,
    left = 0.01,   # a bit more to the right
    right = 0.47,
    bottom = 0.03, # push it DOWN
    top = 0.73,    # lower top boundary too
    align_to = "panel"
  )
plotB
ggsave(file.path(figure_folder, "Simulated", "toy_delay_distr_cfr.pdf"),
       plot = plotB, device = "pdf", width = 6, height = 6)


################## Example C: Primary incidence ##################
max_obs <- max(cases_window[est_idx]); min_obs <- min(cases_window[est_idx])
sigmoid <- function(x) 1 / (1 + exp(-x))
alt_cases1 <- sigmoid(seq(-10, 5, length.out = n_pts)) * (max_obs - min_obs) + min_obs
alt_cases2 <- (1 - sigmoid(seq(-10, 5, length.out = n_pts))) * (max_obs - min_obs) + min_obs
names(alt_cases1) <- names(cases_window); names(alt_cases2) <- names(cases_window)

est_orig <- calc_nish_cfrs(cases_window, cfr_window, delay_distr_med)
est_rise <- calc_nish_cfrs(alt_cases1,   cfr_window, delay_distr_med)
est_fall <- calc_nish_cfrs(alt_cases2,   cfr_window, delay_distr_med)

dfC <- data.frame(
  t = t0 + rep(1:n_ests, 4) - 1,
  CFRs = c(cfr_window[est_idx], est_orig, est_rise, est_fall),
  CFR  = factor(c(rep("True", n_ests), rep("Estimated", n_ests * 3))),
  Type = factor(c(rep("True", n_ests),
                  rep(c("Original", "Rising", "Falling"), each = n_ests)),
                levels = c("True", "Original", "Rising", "Falling"))
)

# Main CFR plot
plotC_main <- ggplot(dfC, aes(x = t, y = CFRs * 100, linetype = CFR, color = Type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Effect of primary incidence") +
  ylab("CFR (%)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  scale_color_manual(values = c("Original" = "red2", "Rising" = "green3", "Falling" = "blue")) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 18),
    legend.text  = element_text(size = 14),
    axis.title.y = element_text(size = 15)
  )

# Inset: primary-incidence (cases) curves
dfCases <- data.frame(
  t = rep(as.Date(names(cases_window)[est_idx]), 3),
  Cases = c(cases_window[est_idx], alt_cases1[est_idx], alt_cases2[est_idx]),
  Type  = factor(rep(c("Original", "Rising", "Falling"),
                     each = length(est_idx)),
                 levels = c("Original", "Rising", "Falling"))
)

inset_cases <- ggplot(dfCases, aes(x = t, y = Cases, color = Type, group = Type)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = c("Original" = "red2", "Rising" = "green3", "Falling" = "blue")) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.margin = unit(rep(1, 4), "mm")
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cases")

# Combine with inset — same placement as your working Example B
plotC <- plotC_main +
  inset_element(
    inset_cases,
    left = 0.01,
    right = 0.47,
    bottom = 0.01,
    top = 0.73,
    align_to = "panel"
  )
plotC
ggsave(file.path(figure_folder, "Simulated", "toy_chging_primary_cfr.pdf"),
       plot = plotC, device = "pdf", width = 6, height = 6)


######################### Figure 2 — Misspecification (CFR) #########################

original_mean <- 20
shorter_mean  <- 16
longer_mean   <- 24

delay_distr_shorter  <- make_delay_distr(Mean = shorter_mean, Sd = round(shorter_mean * 0.9), d)
delay_distr_original <- make_delay_distr(Mean = original_mean, Sd = round(original_mean * 0.9), d)
delay_distr_longer   <- make_delay_distr(Mean = longer_mean, Sd = round(longer_mean * 0.9), d)

case_dates <- as.Date(names(cases7))
first_death_date <- max(cfr_dates[1], case_dates[1]) + d
last_death_date  <- min(cfr_dates[length(cfr_dates)], case_dates[length(cases7)])
death_dates      <- seq(first_death_date, last_death_date, by = "day")
est_dates        <- seq(as.Date("2021-11-01"), as.Date("2022-07-01"), by = "day")

# deaths simulated from true CFRs and "original" delay distribution
noiseless_deaths <- compute_noiseless_deaths(cases7, cfrs, delay_distr_original, est_dates)
n_ests <- length(est_dates)

gt_cfrs <- cfrs[cfr_dates %in% est_dates]
optimal_lag <- compute_optimal_lag(cases7, noiseless_deaths)
lagged_cfrs <- compute_lagged_hfrs(cases7, noiseless_deaths, l = optimal_lag, dates = est_dates)

conv_cfrs_oracle  <- compute_conv_hfrs(cases7, noiseless_deaths, est_dates, delay_distr_original)
conv_cfrs_shorter <- compute_conv_hfrs(cases7, noiseless_deaths, est_dates, delay_distr_shorter)
conv_cfrs_longer  <- compute_conv_hfrs(cases7, noiseless_deaths, est_dates, delay_distr_longer)

R_gamma_shorter  <- compute_R_gammas(cases7, delay_distr_original, delay_distr_shorter, est_dates)
R_gamma_original <- compute_R_gammas(cases7, delay_distr_original, delay_distr_original, est_dates)
R_gamma_longer   <- compute_R_gammas(cases7, delay_distr_original, delay_distr_longer, est_dates)

lagged_distr <- rep(0, d + 1)
lagged_distr[optimal_lag + 1] <- 1
R_gamma_lagged <- compute_R_gammas(cases7, delay_distr_original, lagged_distr, est_dates)

sources <- c("True", "Estimated")
types   <- c("Original", "Light-tailed", "Heavy-tailed", "Lagged")

df <- data.frame(
  t    = rep(est_dates, 5),
  CFRs = c(gt_cfrs, conv_cfrs_oracle, conv_cfrs_shorter, conv_cfrs_longer, lagged_cfrs),
  CFR  = factor(c(rep("True", n_ests), rep("Estimated", n_ests * 4)), levels = sources),
  Type = factor(c(rep("True", n_ests), rep(types, each = n_ests)), levels = c("True", types))
)

mispPlot <- ggplot(df, aes(x = t, y = CFRs*100, linetype = CFR, color = Type)) +
  geom_line(aes(group = interaction(CFR, Type))) +
  scale_color_manual(
    breaks = types,
    values = c("Light-tailed" = "green3", "Original" = "red2",
               "Heavy-tailed" = "blue", "Lagged" = "magenta")
  ) +
  ggtitle("Effect of misspecification") +
  labs(color = "Delay distribution") +
  theme_bw() +
  xlab("") + ylab("CFR (%)") +
  guides(linetype = guide_legend(order = 1), color = guide_legend(order = 2)) +
  scale_x_date(date_breaks = "1 month", labels = NULL) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(size = 17),
    legend.text = element_text(size = 13.5),
    legend.title = element_text(size = 13.5),
    axis.text.x = element_blank(),
    legend.spacing.y = unit(0.5, "cm"),
    legend.margin = margin(t = -10)
  )
mispPlot

delayTypes <- c("Light-tailed", "Original", "Heavy-tailed")
dfDistr <- data.frame(
  t     = rep(0:d, 3),
  Delay = c(delay_distr_shorter, delay_distr_original, delay_distr_longer),
  Type  = rep(delayTypes, each = d + 1)
)

delay_inset <- ggplot(dfDistr, aes(x = t, y = Delay, color = Type)) +
  geom_line() + theme_bw() +
  scale_color_manual(
    breaks = delayTypes,
    values = c("Light-tailed" = "green3", "Original" = "red2", "Heavy-tailed" = "blue")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("") + ylab("") + ggtitle("Delay distribution")

mispPlotInlaid <- mispPlot +
  annotation_custom(
    ggplotGrob(delay_inset),
    xmin = est_dates[1] - 15, xmax = est_dates[1] + 67,
    ymin = 0.35,
    ymax = 1.0
  )

df2 <- data.frame(
  t      = rep(est_dates, 4),
  Ratios = c(R_gamma_original, R_gamma_shorter, R_gamma_longer, R_gamma_lagged),
  Type   = factor(rep(types, each = n_ests), levels = types)
)

ratioPlot <- ggplot(df2, aes(x = t, y = Ratios, color = Type)) +
  geom_line() +
  scale_color_manual(
    breaks = types,
    values = c("Light-tailed" = "green3", "Original" = "red2",
               "Heavy-tailed" = "blue", "Lagged" = "magenta")
  ) +
  ggtitle("Misspecification factor") +
  theme_bw() +
  xlab("") + ylab("Factor") +
  scale_x_date(date_breaks = "1 month", labels = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15),
    axis.text.x = element_blank()
  )

df3 <- data.frame(t = est_dates, Cases = cases7[case_dates %in% est_dates])
casePlot <- ggplot(df3, aes(x = t, y = Cases)) +
  geom_line() +
  ylab("Cases") + ggtitle("Cases") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme_bw() + xlab("") + theme(legend.position = "none")

full_misp_plot <- mispPlotInlaid / ratioPlot / casePlot + plot_layout(heights = c(4.5, 2, 2))
full_misp_plot
ggsave(file.path(figure_folder, "Simulated", "toy_misp_cfr.pdf"),
       plot = full_misp_plot, device = "pdf", width = 9, height = 11)



############################## Figure 4 ##############################
########## Convolve NCHS-based CFR curve against cases to simulate deaths ##########

death_dates_for_cfr <- hfr_dates[(d+1):n_dates]
short_mean <- 12
long_mean <- 24
delay_distr_short <- make_delay_distr(Mean=short_mean, Sd=round(short_mean*.9), d)
delay_distr_long <- make_delay_distr(Mean=long_mean, Sd=round(long_mean*.9), d)
sim_deaths_short <- compute_noiseless_deaths(cases7, cfrs, delay_distr_short, death_dates_for_cfr)
sim_deaths_long <- compute_noiseless_deaths(cases7, cfrs, delay_distr_long, death_dates_for_cfr)


lag_short <- compute_optimal_lag(cases7, sim_deaths_short, verbose=TRUE)
est_dates <- cfr_dates_weekly
n_est <- length(est_dates)


########## Simulate estimated vs true CFRs ##########

# Short delay distribution
cfrs_lagged_short <- compute_lagged_hfrs(cases7, sim_deaths_short, lag_short, w=1, dates=est_dates)
cfrs_conv_short <- compute_conv_hfrs(cases7, sim_deaths_short, est_dates, delay_distr_short, w=1)

gt_cfrs <- cfrs[names(cfrs) %in% est_dates]
CFRtypes <- c("Ground truth", "Oracle convolutional ratio", "Lagged ratio")

dfCFRshort <- data.frame(
  Date = rep(est_dates, length(CFRtypes)),
  CFR = c(gt_cfrs, cfrs_conv_short, cfrs_lagged_short),
  Method = factor(rep(CFRtypes, each = n_est), levels = CFRtypes)
)


########## Long delay distribution ##########

lag_long <- compute_optimal_lag(cases7, sim_deaths_long)
cfrs_lagged_long <- compute_lagged_hfrs(cases7, sim_deaths_long, lag_long, w=1, dates=est_dates)
cfrs_conv_long <- compute_conv_hfrs(cases7, sim_deaths_long, est_dates, delay_distr_long, w=1)

dfCFRlong <- data.frame(
  Date = rep(est_dates, length(CFRtypes)),
  CFR = c(gt_cfrs, cfrs_conv_long, cfrs_lagged_long),
  Method = factor(rep(CFRtypes, each = n_est), levels = CFRtypes)
)

########## Inverted CFR ##########

inv_cfrs <- (1 / cfrs) * min(cfrs, na.rm = TRUE) * max(cfrs, na.rm = TRUE)

sim_deaths_short_inv <- compute_noiseless_deaths(cases7, inv_cfrs, delay_distr_short, death_dates_for_cfr)
lag_short <- compute_optimal_lag(cases7, sim_deaths_short_inv)

cfrs_lagged_short_inv <- compute_lagged_hfrs(cases7, sim_deaths_short_inv, lag_short, w = 1, dates = est_dates)
cfrs_conv_short_inv <- compute_conv_hfrs(cases7, sim_deaths_short_inv, est_dates, delay_distr_short, w = 1)
gt_cfrs_inv <- inv_cfrs[names(inv_cfrs) %in% est_dates]

dfCFRshort_inv <- data.frame(
  Date = rep(est_dates, length(CFRtypes)),
  CFR = c(gt_cfrs_inv, cfrs_conv_short_inv, cfrs_lagged_short_inv),
  Method = factor(rep(CFRtypes, each = n_est), levels = CFRtypes)
)


########## Inverted CFR – Long delay ##########

sim_deaths_long_inv <- compute_noiseless_deaths(cases7, inv_cfrs, delay_distr_long, death_dates_for_cfr)
lag_long <- compute_optimal_lag(cases7, sim_deaths_long_inv)

cfrs_lagged_long_inv <- compute_lagged_hfrs(cases7, sim_deaths_long_inv, lag_long, w = 1, dates = est_dates)
cfrs_conv_long_inv <- compute_conv_hfrs(cases7, sim_deaths_long_inv, est_dates, delay_distr_long, w = 1)

dfCFRlong_inv <- data.frame(
  Date = rep(est_dates, length(CFRtypes)),
  CFR = c(gt_cfrs_inv, cfrs_conv_long_inv, cfrs_lagged_long_inv),
  Method = factor(rep(CFRtypes, each = n_est), levels = CFRtypes)
)
########## Flat CFR ##########

flat_cfr <- 0.01
flat_cfrs <- rep(flat_cfr, length(cfrs))
names(flat_cfrs) <- names(cfrs)

sim_deaths_short_flat <- compute_noiseless_deaths(cases7, flat_cfrs, delay_distr_short, death_dates_for_cfr)
lag_short <- compute_optimal_lag(cases7, sim_deaths_short_flat)

cfrs_lagged_short_flat <- compute_lagged_hfrs(cases7, sim_deaths_short_flat, lag_short, w = 1, dates = est_dates)
cfrs_conv_short_flat <- compute_conv_hfrs(cases7, sim_deaths_short_flat, est_dates, delay_distr_short, w = 1)
gt_cfrs_flat <- flat_cfrs[names(flat_cfrs) %in% est_dates]

dfCFRshort_flat <- data.frame(
  Date = rep(est_dates, length(CFRtypes)),
  CFR = c(gt_cfrs_flat, cfrs_conv_short_flat, cfrs_lagged_short_flat),
  Method = factor(rep(CFRtypes, each = n_est), levels = CFRtypes)
)


########## Flat CFR – Long delay ##########

sim_deaths_long_flat <- compute_noiseless_deaths(cases7, flat_cfrs, delay_distr_long, death_dates_for_cfr)
lag_long <- compute_optimal_lag(cases7, sim_deaths_long_flat)

cfrs_lagged_long_flat <- compute_lagged_hfrs(cases7, sim_deaths_long_flat, lag_long, w = 1, dates = est_dates)
cfrs_conv_long_flat <- compute_conv_hfrs(cases7, sim_deaths_long_flat, est_dates, delay_distr_long, w = 1)

dfCFRlong_flat <- data.frame(
  Date = rep(est_dates, length(CFRtypes)),
  CFR = c(gt_cfrs_flat, cfrs_conv_long_flat, cfrs_lagged_long_flat),
  Method = factor(rep(CFRtypes, each = n_est), levels = CFRtypes)
)


########## Make grid figure ##########

p1 <- ggplot(data=dfCFRshort, aes(x = Date, y = CFR*100, color = Method, linetype=Method)) + 
  geom_line() + labs(color = "CFR", linetype = "CFR") + 
  ggtitle("Realistic CFR (%)", subtitle="Short delay distribution") +
  xlab("") + ylab("CFR") + theme_bw()

p2 <- ggplot(data=dfCFRlong, aes(x = Date, y = CFR*100, color = Method, linetype=Method)) + 
  geom_line() + labs(color = "CFR", linetype = "CFR") + 
  ggtitle("", subtitle="Long delay distribution") +
  xlab("") + ylab("CFR") + theme_bw()

p3 <- ggplot(data=dfCFRshort_inv, aes(x = Date, y = CFR*100, color = Method, linetype=Method)) + 
  geom_line() + labs(color = "CFR", linetype = "CFR") + 
  ggtitle("Inverted-Realistic CFR (%)", subtitle="Short delay distribution") +
  xlab("") + ylab("CFR") + theme_bw()

p4 <- ggplot(data=dfCFRlong_inv, aes(x = Date, y = CFR*100, color = Method, linetype=Method)) + 
  geom_line() + labs(color = "CFR", linetype = "CFR") + 
  ggtitle("", subtitle="Long delay distribution") +
  xlab("") + ylab("CFR") + theme_bw() 

p5 <- ggplot(data=dfCFRshort_flat, aes(x = Date, y = CFR*100, color = Method, linetype=Method)) + 
  geom_line() + labs(color = "CFR", linetype = "CFR") + 
  ggtitle(expression("Stationary CFR, ratios (%) and " * A[t]^gamma), subtitle = "Short delay distribution") +
  scale_y_continuous(sec.axis = sec_axis(~ . / (flat_cfr*100), name = expression(A[t]^gamma))) +
  xlab("") + ylab("CFR") + theme_bw()

p6 <- ggplot(data=dfCFRlong_flat, aes(x = Date, y = CFR*100, color = Method, linetype=Method)) + 
  geom_line() + labs(color = "CFR", linetype = "CFR") + 
  ggtitle("", subtitle="Long delay distribution") +
  scale_y_continuous(sec.axis = sec_axis(~ . / (flat_cfr*100), name = expression(A[t]^gamma))) +
  xlab("") + ylab("CFR") + theme_bw() 

Plot <- ((p1|p2) + plot_layout(tag_level = 'new')) /
  ((p3|p4) + plot_layout(tag_level = 'new')) /
  ((p5|p6) + plot_layout(tag_level = 'new')) +
  plot_annotation(title = 'Ratio estimates under various simulation models',
                  theme = theme(plot.title = element_text(size = 17))) +
  plot_layout(guides="collect") &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        plot.subtitle=element_text(size=13))

Plot
ggsave(file.path(figure_folder, "Simulated", "simulated_results_cfr.pdf"), plot = Plot, device = "pdf", width = 9, height = 8.5)






