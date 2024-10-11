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

#### Tune Lagged
tune_lagged_hfrs_retro <- function(hosps, deaths, lags, ws, dates) {
  lagged_cv_maes <- sapply(ws, function(w) {
    n_into_past <- floor(w/2)
    n_into_future <- floor((w-1)/2)
    maes_by_lag <- sapply(lags, function(lag) {
      lagged_hfrs_cv <- sapply(dates, function(t) {
        window <- seq(t-n_into_past, t+n_into_future, by="day")
        window_no_t <- window[!(window==t)]
        hosps_window <- hosps[as.Date(names(hosps)) %in% window_no_t]
        deaths_window <- deaths[as.Date(names(deaths)) %in% (window_no_t+lag)]
        total_hosps <- sum(hosps_window)
        
        # Check if total_hosps is 0, if so, return NA to avoid division by 0
        if (total_hosps == 0) {
          return(NA)
        }
        
        hfr_t <- sum(deaths_window)/total_hosps
      })
      
      # Remove NA values before calculating mae
      valid_values <- !is.na(lagged_hfrs_cv)
      Y_tl_ests <- hosps[as.Date(names(hosps)) %in% dates[valid_values]] * lagged_hfrs_cv[valid_values]
      Y_tl <- deaths[as.Date(names(deaths)) %in% (dates[valid_values] + lag)]
      mae(Y_tl_ests, Y_tl)
    })
    maes_by_lag
  })
  
  # Reshape in case only 1 lag is given
  lagged_cv_maes <- matrix(lagged_cv_maes, nrow=length(lags), ncol=length(ws))
  rownames(lagged_cv_maes) <- lags
  colnames(lagged_cv_maes) <- ws
  return(lagged_cv_maes)
}


###### NISHIURA HFR ######
compute_ahfr_dynamic <- function(hosps, deaths, t, delay_shape, w=1) {
  # W is length of trailing window ending at t
  d <- length(delay_shape) - 1
  contributing_hosps <- sum(sapply(0:(w-1), function(i) {
    hosps_in_trailing_window <- hosps[names(hosps) %in% seq(t-d-i, t-i, by="day")]
    sum(rev(hosps_in_trailing_window) * delay_shape)
  }))
  deaths_in_window <- sum(deaths[(names(deaths) > t-w) & (names(deaths) <= t)])
  deaths_in_window/contributing_hosps
}

compute_ahfrs_dynamic <- function(hosps, deaths, dates, delay_shape, w=1) {
  # Only a real-time estimator. Not implemented since the beginning.
  ahfrs <- sapply(dates, compute_ahfr_dynamic, hosps=hosps, deaths=deaths, w=w, delay_shape=delay_shape)
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
