TimeFunction <- function (n, dt = 1) {
  ## n: Number of observations. Positive integer
  ## dt: Time step between observations in seconds
  Time <- (1:n)/(3600/dt) ## Time in hours since tagging
}

AddTime <- function (data, dt = 1) {
  ### Time provides time since tagging
  data$Time <- NULL
  for(ind in unique(data$Ind)) {
    indicator <- (data$Ind == ind)
    data$Time[indicator] <- TimeFunction(n = sum(indicator), dt = dt)
  }
  return(data)
}

BeforeExposure <- function (data) {
  ## Find times of first buzz and first exposure time and
  ## restrict data to be between these for each whale
  for (k in unique(data$Ind)) {
    mintime <- min(data$Time[(data$Ind == k) & (data$Buzz == 1)])
    maxtime <- min(data$Time[(data$Ind == k) & (data$LOS == 1)])
    data <- subset(data, !(Ind == k & Time > maxtime))
    data <- subset(data, !(Ind == k & Time < mintime))
  }
  return(data)
}
