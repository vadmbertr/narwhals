#--------------------------------------------------------------------------------
## Objective : fit glmer model using fixed AR and depth coefs
#---------------------------------------------------------------------------------

library(splines)
library(lme4)

glmer_buzz_AR_expo <- function (dataArg, dataAR, ARcoefs) {
  data <- dataArg
  data$AR <- as.matrix(dataAR) %*% ARcoefs
  glmerAllBuzz <- glmer(Buzz ~ offset(AR) + ns(X, knots = quantile(data$X[data$X > 0], 1:2 / 3)) + (1 | Ind),
                        data = data,
                        nAGQ = 0,
                        weights = n,
                        family = poisson)
  return(glmerAllBuzz)
}