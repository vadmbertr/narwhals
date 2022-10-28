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

AfterExposure <- function (data, no_stress = TRUE) {
  if (no_stress) {
    ## We cut the first 24 hours of all individuals to avoid
    ## the effect of tagging
    data <- data[data$Time > 24,]
  } else {
    data <- data[data$Time > 0,]
  }
  return(data)
}

ExposureFunction <- function(X){
  ## X: Distance to ship in meters measured each second. Supposed to be positive

  ## Exposure variable
  Xtilde <- 1000/X ## Inverse of distance to ship in kilometers
  Xtilde[is.na(Xtilde)] <- 0 ## If not in line of sight, exposure is zero
  Xtilde
}

AddExposure <- function (data) {
  data$Xtilde <- ExposureFunction(X = data$Dist_Ship) ## Exposure variable based on 1/distance in km
  data$X <- data$Xtilde*(data$Seismik == 1)
  data$P <- data$Xtilde*(data$Seismik == 0)
  return(data)
}

OnlyAirgun <- function (data) {
  return(subset(data, P == 0))
}

RemoveNA <- function (data, cols) {
  return(data[complete.cases(dataAll[, cols])])
}