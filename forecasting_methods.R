ensemble_ts <- function(elek, elekwave, eleknoise, for_long_1, for_long_2, for_long_3, freq = 24) {

  ar_w <- ar(elekwave, FALSE, 2*freq)
  ar_w_f <- forecast(ar_w, freq)
  ar_n <- ar(eleknoise, FALSE, 2*freq)
  ar_n_f <- forecast(ar_n, freq)
  
  HW_w <- HoltWinters(elekwave, beta = F, alpha = 0.15, gamma = 0.95, seasonal = "additive")
  HW_w_f <- forecast(HW_w, freq)
  HW_n <- HoltWinters(eleknoise, gamma = F, beta = F)
  HW_n_f <- forecast(HW_n, freq)
  
  stl.decom <- stl(elek, s.window = "periodic", robust = T)
  stl_exp_f <-forecast(stl.decom, h = freq, method = "ets", etsmodel = "ZZN")
  stl_arima_f <- forecast(stl.decom, h = freq, method = "arima")
  
  ann <- nnetar(elek)
  ann_f <- forecast(ann, freq)
  
  naive_f <-snaive(elek, freq)

  forecasts <- matrix(rep(0, 9*freq), nrow = freq, byrow = T)
  forecasts[,1] <- ar_w_f$mean + ar_n_f$mean
  forecasts[,2] <- HW_w_f$mean + HW_n_f$mean
  forecasts[,3] <- stl_exp_f$mean
  forecasts[,4] <- stl_arima_f$mean
  forecasts[,5] <- ann_f$mean
  forecasts[,6] <- naive_f$mean
  forecasts[,7] <- for_long_1
  forecasts[,8] <- for_long_2
  forecasts[,9] <- for_long_3
  
  list(methods = forecasts)
}

ensemble_reg <- function(train, wavet, noise, n, k, freq = 24) {
  
  # OLS regression
  seas_1 <- factor(rep(1:freq, n))
  data_train <- data.frame(Y = train, seas = seas_1)
  mat_train <- data.frame(Y = train, model.matrix(as.formula("Y ~ 0 +."), data = data_train))
  
  lm_m <- lm(Y ~ 0 + ., data = mat_train)
  
  new_data <- data.frame(mat_train[1:freq, -1])
  lm_f <- predict(lm_m, new_data)
  
  # SVR - wavelet + noise
  data_train <- data.frame(Y = wavet, seas = seas_1)
  mat_train <- data.frame(Y = wavet, model.matrix(as.formula("Y ~ 0 +."), data = data_train))
  
  svr_w <- ksvm(Y ~ ., data = mat_train, type = "eps-svr", kernel = "rbfdot", eps = 0.08)
  
  svr_w_f <- predict(svr_w, new_data)
  
  data_train <- data.frame(Y = noise, seas = seas_1)
  mat_train <- data.frame(Y = noise, model.matrix(as.formula("Y ~ 0 +."), data = data_train))
  
  svr_n <- ksvm(Y ~ ., data = mat_train, type = "eps-svr", kernel = "rbfdot", eps = 0.05)
  
  svr_n_f <- predict(svr_n, new_data)
  
  forecasts <- matrix(rep(0, 2*freq), nrow = freq, byrow = T)
  forecasts[,1] <- lm_f
  forecasts[,2] <- as.vector(svr_w_f) + as.vector(svr_n_f) 
  
  list(methods = forecasts)
}