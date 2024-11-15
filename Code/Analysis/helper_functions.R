library(lubridate)
library(Matrix)
library(extraDistr)

################ LAGGED HFR ################
compute_lagged_hfr <- function(hosps, deaths, t, l=19, w=1,
                               real_time=TRUE, centered=NA) {
  # real_time is same as forward-looking (previously)
  # also removed weekly, which made sure rounding to nearest multiple of 7
  w <- ifelse(w==0, 1, w)
  if (is.Date(date) == FALSE) {
    t <- as.Date(t)
  }
  
  if (is.na(centered)) {
    centered <- ifelse(real_time, FALSE, TRUE)
  }
  n_into_past <- ifelse(centered, floor((w-1)/2), w-1)
  n_into_future <- ifelse(centered, max(0, floor(w/2)), 0)
  if (real_time==FALSE) { # Y_{t+l}/X_{t}
    hosps_t <- sum(hosps[names(hosps)>=t-n_into_past & names(hosps)<=t+n_into_future])
    deaths_in_l <- sum(deaths[names(deaths)>=t+l-n_into_past & names(deaths)<=t+l+n_into_future])
    return(deaths_in_l / hosps_t)
  }
  # Y_{t}/X_{t-l}; different strategies for window
  hosps_l_ago <- sum(hosps[names(hosps)>=t-n_into_past-l & names(hosps)<=t+n_into_future-l])
  # print(t-n_into_past-l)
  # print(hosps[names(hosps)==(t-n_into_past-l)])
  deaths_t <- sum(deaths[names(deaths)>=t-n_into_past & names(deaths)<=t+n_into_future])
  return(deaths_t / hosps_l_ago)
}

compute_lagged_hfrs <- function(hosps, deaths, l, w=1, 
                                real_time=TRUE, dates=NULL, centered=NA) {
  if (is.null(dates)) {
    dates <- names(deaths)
  }
  lagged_hfrs <- sapply(dates, compute_lagged_hfr, l=l, w=w, real_time=real_time, 
                        hosps=hosps,deaths=deaths, centered=centered)
  names(lagged_hfrs) <- dates
  return(lagged_hfrs)
}

compute_optimal_lag <- function(hosps, deaths, verbose=FALSE) {
  hosp_dates <- names(hosps); death_dates <- names(deaths)
  intersect_dates <- as.Date(intersect(hosp_dates, death_dates))
  deaths_both <- deaths[death_dates %in% intersect_dates]
  hhs_hosps_both <- hosps[hosp_dates %in% intersect_dates]
  cc <- ccf(hhs_hosps_both, deaths_both, lag.max=40, plot=FALSE)
  max_correlation_lag <- which.max(abs(cc$acf))
  oracle_lag <- abs(cc$lag[max_correlation_lag])
  if (verbose) 
    print(paste0("Correlation-maximizing lag is ", oracle_lag, " days."))
  return(oracle_lag)
}


###### NISHIURA HFR ######
compute_conv_hfr <- function(hosps, deaths, t, delay_shape, w=1) {
  # W is length of trailing window ending at t
  d <- length(delay_shape) - 1
  w <- ifelse(w==0, 1, w)
  contributing_hosps <- sum(sapply(0:(w-1), function(i) {
    hosps_in_trailing_window <- hosps[names(hosps) %in% seq(t-d-i, t-i, by="day")]
    sum(rev(hosps_in_trailing_window) * delay_shape)
  }))
  deaths_in_window <- sum(deaths[(names(deaths) > t-w) & (names(deaths) <= t)])
  deaths_in_window/contributing_hosps
}

compute_conv_hfrs <- function(hosps, deaths, dates, delay_shape, w=1) {
  # Only a real-time estimator. Not implemented since the beginning.
  ahfrs <- sapply(dates, compute_conv_hfr, hosps=hosps, deaths=deaths, w=w, delay_shape=delay_shape)
  names(ahfrs) <- dates
  return(ahfrs)
}

compute_ahfrs_og <- function(hosps, deaths, dates, delay_shape, start_date) {
  d <- length(delay_shape)-1
  ref_date <- start_date
  for (i in 1:length(dates)) {
    t <- dates[i]
    window_dates <- seq(ref_date, t, by="day")
    deaths_in_window <- sum(deaths[names(deaths) %in% window_dates])
    contributing_hosps_in_window <- sum(sapply(window_dates, function(date) {
      hosps_in_trailing_window <- hosps[names(hosps) %in% seq(date-d, date, by="day")]
      sum(rev(hosps_in_trailing_window) * delay_shape)
    }))
    if (i==1) {
      nums <- c(deaths_in_window); denoms <- c(contributing_hosps_in_window)
      hfr <- deaths_in_window/contributing_hosps_in_window
      ahfrs <- c(hfr)
    } else {
      deaths_so_far <- nums[i-1] + deaths_in_window
      contr_hosps_so_far <- denoms[i-1] + contributing_hosps_in_window
      nums <- c(nums, deaths_so_far); denoms <- c(denoms, contr_hosps_so_far)
      hfr <- deaths_so_far/contr_hosps_so_far
      ahfrs <- c(ahfrs, hfr)
    }
    ref_date <- ref_date + 1
  }
  names(ahfrs) <- dates
  return(ahfrs)
}

###### Misc functions for preprocessing, estimation, and evaluation ######
mae <- function(vec1, vec2) {
  mean(abs(as.numeric(vec1) - as.numeric(vec2)), na.rm=T)
}

seven_day_smoothing <- function(data, centered=TRUE) {
  n <- length(data)
  if (centered==FALSE) {
    smoothed <- sapply(1:n, function(i) { # Will be flat for the first 7 days
      first <- ifelse(i-6 < 1, 1, i-6)
      last <- ifelse(i-6 < 1, 7, i)
      mean(data[first:last])
    })
  } else {
    smoothed <- sapply(1:n, function(i) { # Later, should account for boundary effects
      first <- ifelse(i-3 < 1, 1, i-3)
      last <- ifelse(i+3 > n, n, i+3)
      mean(data[first:last])
    })
  }
  
  if (!is.null(names(data))) {
    names(smoothed) <- names(data)
  }
  return(smoothed)
}

intersect_dates <- function(dates1, dates2) {
  return(as.Date(intersect(dates1, dates2), origin = "1970-01-01"))
}

get_gamma_params <- function(Mean, Sd) {
  Var <- Sd**2
  shape <- (Mean**2)/Var; rate <- Mean/Var
  return(c(shape, rate))
}  

make_delay_distr <- function(Mean, Sd, d) {
  params <- get_gamma_params(Mean=Mean, Sd=Sd)
  shape <- params[1]; rate <- params[2]
  DelayShape <- ddgamma(0:d, shape, rate)
  return(DelayShape)
}


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

compute_R_gammas <- function(hosps, true_delay_distr, plugin_delay_distr, est_dates) {
  d <- length(true_delay_distr) - 1
  R_gammas <- sapply(est_dates, function(t) {
    trailing_dates <- seq(t-d,t, by="day")
    trailing_hosps <- hosps[names(hosps) %in% trailing_dates]
    num <- sum(rev(true_delay_distr)*trailing_hosps)
    denom <- sum(rev(plugin_delay_distr)*trailing_hosps)
    num/denom
  })
  names(R_gammas) <- est_dates
  return(R_gammas)
}