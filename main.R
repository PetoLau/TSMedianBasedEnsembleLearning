# Header ----
# gc()

library(data.table)
library(TSrepr) # devtools::install_github("PetoLau/TSrepr")
library(forecast)
library(matrixStats)
library(kernlab)
library(wmtsa)
library(scales)
library(car)
library(lubridate)

# Read data ----
data <- fread("toronto.csv")
# subset last two years of data
data[, year := year(date)]
data <- data[year %in% c(2014, 2015)]
setkey(data, "dateTime")
# remove weekends and subset only working days without Monday - as in the our paper in ACTA Journal
data <- data[dayOfWeek %in% 2:5]

# Long term forecasts ----
# subset train curve for moving median and average computation
train_data <- data[year == 2014, load]
test_data <- data[year == 2015, load]
# plot it
plot(ts(train_data))
lines(ts(test_data), col = "red")
# helpers
freq <- 24
n_days <- length(train_data)/freq
  
# moving average
curve_mean <- rep(0, freq*n_days)
for(i in 0:(n_days-1)) {
  for(j in 1:freq) {
    curve_mean[(i*freq)+j] <- mean(train_data[((freq*i)+1):(freq*(i+3))])
  }
}

plot(ts(train_data, start = 0, freq = freq))
lines(ts(curve_mean, start = 1, freq = freq), col = "red", lwd = 2)

# daily enhancement
train_ts <- ts(train_data, start = 0, freq = freq)
mov_ave <- c(t(tapply(train_ts, cycle(train_ts), mean)))
train_ts_mean <- mean(mov_ave)
train_ts_diff <- c(mov_ave - train_ts_mean)

# compute final year ahead prediction
curve_mean_final <- rep(0, freq*n_days)
for(i in 0:(n_days-1)) {
  for(j in 1:freq) {
    curve_mean_final[(i*freq)+j] <- curve_mean[(i*freq)+j] + train_ts_diff[j]
  }
}

plot(ts(train_data, start = 0, freq = freq))
lines(ts(curve_mean_final, start = 1, freq = freq), col = "orange", lwd = 2)
lines(ts(curve_mean, start = 1, freq = freq), col = "red", lwd = 2)

# moving median
curve_median <- rep(0, freq*n_days)
for(i in 0:(n_days-1)) {
  for(j in 1:freq) {
    curve_median[(i*freq)+j] <- median(train_data[((freq*i)+1):(freq*(i+3))])
  }
}

plot(ts(train_data, start = 0, freq = freq))
lines(ts(curve_median, start = 1, freq = freq), col = "red", lwd = 2)

# daily enhancement
# train_ts <- ts(train_data, start = 0, freq = freq)
mov_med <- c(t(tapply(train_ts, cycle(train_ts), median)))
train_ts_med <- mean(mov_med)
train_ts_diff <- c(mov_med - train_ts_med)

# compute final year ahead prediction
curve_median_final <- rep(0, freq*n_days)
for(i in 0:(n_days-1)) {
  for(j in 1:freq) {
    curve_median_final[(i*freq)+j] <- curve_median[(i*freq)+j] + train_ts_diff[j]
  }
}

plot(ts(train_data, start = 0, freq = freq))
lines(ts(curve_median_final, start = 1, freq = freq), col = "orange", lwd = 2)
lines(ts(curve_median, start = 1, freq = freq), col = "red", lwd = 2)

# TBATS model
# split train_set to 3 parts for better forecasting accuracy
day_thirds <- n_days/3

tb_1 <- tbats(train_ts[1:(freq*day_thirds)], seasonal.periods = c(freq, 5*freq), use.box.cox = F,
              use.trend = T, use.damped.trend = T, use.arma.errors = T,
              use.parallel = T, num.cores = 4)
tb_1_f <- forecast(tb_1, freq*day_thirds)

tb_2 <- tbats(train_ts[((freq*day_thirds)+1):(freq*day_thirds*2)], seasonal.periods = c(freq, 5*freq),
              use.box.cox = F,
              use.trend = T, use.damped.trend = T, use.arma.errors = T,
              use.parallel = T, num.cores = 4)
tb_2_f <- forecast(tb_2, freq*day_thirds)

tb_3 <- tbats(train_ts[((freq*day_thirds*2)+1):(freq*day_thirds*3)], seasonal.periods = c(freq, 5*freq),
              use.box.cox = F,
              use.trend = T, use.damped.trend = T, use.arma.errors = T,
              use.parallel = T, num.cores = 4)
tb_3_f <- forecast(tb_3, freq*day_thirds)

tbats_for <- c(tb_1_f$mean, tb_2_f$mean, tb_3_f$mean)

# Main cycle - training and forecasting ----
source("forecasting_methods.R")
source("ensembleWeighting.R")

# First iteration - sliding window approach
# lengths of training sets
# for regression methods
n <- 4
k <- 5
# for time series analysis methods
n2 <- 6
k2 <- 7

# subset test set - to train and test part again
test_data <- c(tail(train_data, freq*n2), test_data)
# 1. time series analysis part
train_ts <- ts(test_data[1:(n2*freq)], start = 0, freq = freq)
test_ts <- ts(test_data[((n2*freq)+1):(k2*freq)], start = n2, freq = freq)

# wavelet transform for time series analysis methods
wavet <- wavShrink(train_ts, wavelet = "s8", shrink.fun = "soft", thresh.fun = "universal",
                   threshold = NULL, thresh.scale = 1, xform = "dwt", noise.variance = 0.0, reflect = TRUE)

train_wave <- ts(wavet, start = 0, freq = freq)
remain <- train_ts - train_wave

# 2. regression part
train_reg <- test_data[(((n2-n)*freq)+1):(n2*freq)]
train_reg_ts <- ts(train_reg, start = 0, freq = freq)

# wavelet transform for regression methods
wavet_reg <- wavShrink(train_reg_ts, wavelet = "s8", shrink.fun = "soft", thresh.fun = "universal",
                    threshold = NULL, thresh.scale = 1, xform = "dwt", noise.variance = 0.0, reflect = TRUE)
remain_reg <- train_reg - wavet_reg

# subset long-term forecasts
for_tbats <- tbats_for[1:freq]
for_ave <- curve_mean_final[1:freq]
for_med <- curve_median_final[1:freq]

# forecasting error tracking variables (MAPE)
# base_methods <- data.frame(rep(0, 11))
# ens_method <- data.frame(rep(0, 1))

# forecasting
ts_for <- ensemble_ts(elek = train_ts, elekwave = train_wave, eleknoise = remain,
                      for_tbats, for_ave, for_med, freq = freq)
reg_for <- ensemble_reg(train_reg, wavet_reg, remain_reg, n, k, freq = freq)

methods_for <- cbind(ts_for$methods, reg_for$methods)

# Ensemble learning
ens_wei <- ensembleWeighting(test_ts, methods_for, rep(1,11), rep(1,11))

# save data to variables
whole_forecast <- rep(0, freq*(n_days-2))
whole_forecast[1:freq] <- ens_wei$forecast
base_methods <- sapply(1:ncol(methods_for), function(x) mape(test_ts, methods_for[,x]))
ens_err <- mape(test_ts, ens_wei$forecast)

train_ts <- train_ts[-(1:freq)]
train_reg <- train_reg[-(1:freq)]

# continue cycle
for(i in 1:(n_days-3)) {

  train_ts <- ts(c(train_ts, test_ts), start=i, freq = freq)
  train_reg <- c(train_reg, test_ts)

  for_tbats <- tbats_for[((freq*i)+1):(freq*(i+1))]
  for_ave <- curve_mean_final[((freq*i)+1):(freq*(i+1))]
  for_med <- curve_median_final[((freq*i)+1):(freq*(i+1))]
  
  wavet_ts <- wavShrink(train_ts, wavelet="s8", shrink.fun="soft", thresh.fun="universal", threshold=NULL, thresh.scale=1, xform="dwt", noise.variance=0.0, reflect=TRUE)
  wavet_reg <- wavShrink(ts(train_reg, freq = freq), wavelet="s8", shrink.fun="soft", thresh.fun="universal", threshold=NULL, thresh.scale=1, xform="dwt", noise.variance=0.0, reflect=TRUE)
  
  train_wave <- ts(wavet_ts, start = i, freq = freq)
  remain <- train_ts - train_wave
  remain_reg <- train_reg - wavet_reg
  
  ts_for <- ensemble_ts(elek = train_ts, elekwave = train_wave, eleknoise = remain,
                        for_tbats, for_ave, for_med, freq = freq)
  reg_for <- ensemble_reg(train_reg, wavet_reg, remain_reg, n, k, freq = freq)
  
  methods_for <- cbind(ts_for$methods, reg_for$methods)
  
  test_ts <- ts(test_data[((n2*freq + (freq*i))+1):(n2*freq + (freq*(i+1)))], start = n2+i, freq = freq)
 
  ens_wei <- ensembleWeighting(test_ts, methods_for, ens_wei$weights, ens_wei$weights_scaled)
  
  whole_forecast[((freq*i)+1):(freq*(i+1))] <- ens_wei$forecast
  base_methods <- rbind(base_methods, sapply(1:ncol(methods_for), function(x) mape(test_ts, methods_for[,x])))
  ens_err <- c(ens_err, mape(test_ts, ens_wei$forecast))
  
  plot(train_ts, xlim = c(i, n2+i+1), main = paste(i))
  lines(test_ts, col = "blue", lwd = 2)
  lines(ts(ens_wei$forecast, start = n2+i, freq = freq), col = "red", lwd = 2)
  
  train_ts <- train_ts[-(1:freq)]
  train_reg <- train_reg[-(1:freq)]
}

colMeans(base_methods)
mean(ens_err)
apply(base_methods, 2, median)
median(ens_err)

wilcox.test(ens_err, as.vector(base_methods[,4]), paired = T, exact = T, conf.int = F, alternative = "less")
