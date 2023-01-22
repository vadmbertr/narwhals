#--------------------------------------------------------------------------------
## Objective : fit glmer model using fixed AR and depth coefs
#---------------------------------------------------------------------------------

library(splines)
library(lme4)

glmer_buzz_ARDepth_expo <- function (dataArg, dataAR, ARcoefs, Depthcoefs) {
  data <- dataArg
  data$ARDepth <- as.matrix(dataAR) %*% ARcoefs
  data$ARDepth <- data$ARDepth + as.matrix(ns(data$Depth, knots = c(-323, -158, -54))) %*% Depthcoefs
  glmerAllBuzzDepth <- glmer(Buzz ~ offset(ARDepth) + ns(X, knots = quantile(data$X[data$X > 0], 1:2 / 3)) + (1 | Ind),
                             data = data,
                             nAGQ = 0,
                             weights = n,
                             family = poisson)
  return(glmerAllBuzzDepth)
}