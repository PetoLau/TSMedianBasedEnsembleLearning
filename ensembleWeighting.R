ensembleWeighting <- function(test, forecasts, weights, weights_scaled) {
  
  preds_ens <- rowWeightedMeans(forecasts, weights_scaled)
  
  mdae_for <- sapply(1:ncol(forecasts), function(x) mdae(test, forecasts[,x]))
  
  med <- median(mdae_for)
  multipliers <- vector(length = ncol(forecasts))
  
  for(j in 1:ncol(forecasts)) {
    multipliers[j] <- med / mdae_for[j]
  }
  
  new_weights <- weights * multipliers
  weights_scaled <- rescale(new_weights, to = c(1,10))
  
  list(weights_scaled = weights_scaled, weights = new_weights, forecast = preds_ens)
}
