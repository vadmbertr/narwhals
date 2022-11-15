#--------------------------------------------------------------------------------
## Objective : fit a bi exponential regression on the AR coefs
#---------------------------------------------------------------------------------

library(broom)

maxlag.bic <- readRDS("../data/glmER_buzz_depth_maxlag/maxlag.bic.rds")
ARcoef.best <- readRDS("../data/glmER_buzz_depth_maxlag/ARcoef.best.rds")

maxlag.opt <- as.integer(maxlag.bic[which.min(maxlag.bic[, 2]), 1])
ARcoef.est <- ARcoef.best$estimate[1:(maxlag.opt)+5]
ARcoef.se <- ARcoef.best$std.error[1:(maxlag.opt)+5]

biexp.reg <- nls(ARcoef.est ~ SSbiexp(1:maxlag.opt, A1, lrc1, A2, lrc2), weights = 1/ARcoef.se)
saveRDS(tidy(biexp.reg), file = "../data/glmER_buzz_depth_maxlag/ARcoef.RegBiExp.rds")

A1 <- summary(biexp.reg)$coefficients[1,1]
lrc1 <- exp(summary(biexp.reg)$coefficients[2,1])
A2 <- summary(biexp.reg)$coefficients[3,1]
lrc2 <- exp(summary(biexp.reg)$coefficients[4,1])
ARcoef.fit <- data.frame(maxlag = as.factor(rep(maxlag.opt, maxlag.opt)), lag = 1:maxlag.opt, truth = ARcoef.est)
ARcoef.fit$estimate <- A1 * exp(-lrc1 * ARcoef.fit$lag) + A2 * exp(-lrc2 * ARcoef.fit$lag)
saveRDS(ARcoef.fit, file = "../data/glmER_buzz_depth_maxlag/ARcoef.fit.rds")
