library(splines)
library(lme4)
library(ggplot2)
library(data.table)  ## For fast and easy reading of large data sets


#--------------------------------------------------------------------------------
## Objective : construct a model including memory and depth 
## estimated from a glm restricted to data before exposure
## and use them as offset in glmer restricted to data after exposure
#---------------------------------------------------------------------------------


## The two functions used in the analysis ##
TimeFunction <- function(n, dt = 1){
  ## n: Number of observations. Positive integer
  ## dt: Time step between observations in seconds
  Time <- (1:n)/(3600/dt) ## Time in hours since tagging
}

ExposureFunction <- function(X){
  ## X: Distance to ship in meters measured each second. Supposed to be positive
  
  ## Exposure variable
  Xtilde <- 1000/X ## Inverse of distance to ship in kilometers
  Xtilde[is.na(Xtilde)] <- 0 ## If not in line of sight, exposure is zero
  Xtilde
}

#---------------------------------------------------------------------------------
# Data importation, total dataset
#---------------------------------------------------------------------------------
#dataAll <- fread("~/Dropbox/Effect study EGRL 2018/Database/db_narwhal_2017_2018.txt")
#data <- dataAll
#
### Time provides time since tagging
#data$Time <- numeric()
#for(ind in unique(data$Ind)){
#  indicator <- (data$Ind == ind)
#  data$Time[indicator] <- TimeFunction(n = sum(indicator))
#}
### keep the 2018 data
#data = data[data$Year=="2018",]
#
#
#datasub = data[seq(1,dim(data)[1], by = 60),]
#fwrite(datasub, file="~/Dropbox/Adeline/Projet/Susanne Ditlevsen/Narval-Projet Tutore/subdb_narwhal_2018.txt")
#---------------------------------------------------------------------------------
# Data importation; sub-sample seq(1,dim(data)[1], by = 60), only 2018
#---------------------------------------------------------------------------------
dataAll <- fread("../../../data/1_buzzing/subdb_narwhal_2018.txt")
data <- dataAll


#---------------------------------------------------------------------------------
# 1. glm model before exposure to estimate the memory and influence of depth
#---------------------------------------------------------------------------------
## Define the first 100 time lags
temp <- as.data.frame(shift(data$Buzz, n = 1:100, give.names = TRUE))
data <- cbind(temp, data)

#---------------------------------------------------------------------------------
# We restrict the data to be before exposure
## Find times of first buzz and first exposure time and 
## restrict data to be between these for each whale
for(k in unique(data$Ind)){
  mintime <- min(data$Time[(data$Ind == k) & (data$Buzz == 1)])
  maxtime <- min(data$Time[(data$Ind == k) & (data$LOS == 1)])
  data <- subset(data, !(Ind == k & Time > maxtime))
  data <- subset(data, !(Ind == k & Time < mintime))
}


#---------------------------------------------------------------------------------
# Estimation of the model
# We fit a memory component to buzzes with lags from 1 to $k$, where we varied $k$ between 1 and 100, adjusted to Depth, entered as a nonlinear spline function. For each model, we computed the BIC and chose the model with the lowest BIC. Below, only a subset of k are chosen. This can be changed as desired in the definition of lagvector.
LagVariables <- names(data[,1:100])

glmresultsDepth <- list(NULL)

lagvector <- c(50,seq(55,65,1)) # TODO: 1- élargir la plage de recherche

n <- length(lagvector)
AICvector <- numeric(n) 
BICvector <- numeric(n)
# TODO: ? paraléliser ?
for(i in 1:n){
  form <- paste("Buzz ~ Ind + ns(Depth, knots = c(-323, -158, -54)) + ",
                paste(LagVariables[1:lagvector[i]], collapse = " + "))

  # TODO: 2- ajout effet aléatoire sur les individus
  glmAllBuzzDepth <- glm(form,
                         data = data, 
                         family = poisson)
  
  AIC <- AIC(glmAllBuzzDepth)
  BIC <- BIC(glmAllBuzzDepth)
  
  AICvector[i]    <- AIC
  BICvector[i]    <- BIC
  coeftable <- coefficients(summary(glmAllBuzzDepth))
  glmresultsDepth[[i]] <- list(lag = lagvector[i], coeftable = coeftable, AIC = AIC, BIC = BIC)
}
dataBICDepth <- data.frame(lag = lagvector, BIC = BICvector)


#---------------------------------------------------------------------------------
## Estimation of the AR coefficients.

k <- length(glmresultsDepth)
kvector <- 1:12#c(2:6,16:19)
temp    <- glmresultsDepth[[1]]
ARcoef  <- temp$coeftable[-(1:10), "Estimate"]
lag     <- temp$lag

fitresult <- data.frame(maxlag = as.factor(rep(lag, lag)), lag = 1:lag,
                        estimate = ARcoef)

curvefit  <- data.frame(maxlag = as.factor(rep(lag, lag)), lag = 1:lag,
                        b1 = rep(NA, lag), b2 = rep(NA, lag),
                        b3 = rep(NA, lag), b4 = rep(NA, lag))

for(i in kvector){
  temp    <- glmresultsDepth[[i]]
  lag     <- temp$lag
  # TODO: 2- différent pour GLMER
  ARcoef  <- temp$coeftable[-(1:10), "Estimate"]
  std.err <- temp$coeftable[-(1:10),"Std. Error"]
  w       <- 1/std.err
  n       <- length(ARcoef)

  # fm <- nls(ARcoef ~ SSbiexp(1:n, A1, lrc1, A2, lrc2), weights = w)
  # b1 <- summary(fm)$coefficients[1,1]
  # b2 <- exp(summary(fm)$coefficients[2,1])
  # b3 <- summary(fm)$coefficients[3,1]
  # b4 <- exp(summary(fm)$coefficients[4,1])
  fitresulti <- data.frame(maxlag = as.factor(rep(lag, lag)), lag = 1:lag,
                           estimate = ARcoef)
  # curvefiti  <- data.frame(maxlag = as.factor(rep(lag, lag)), lag = 1:lag,
  #                          b1 = rep(b1, lag), b2 = rep(b2, lag),
  #                          b3 = rep(b3, lag), b4 = rep(b4, lag))
  fitresult <- rbind(fitresult, fitresulti)
  # curvefit  <- rbind(curvefit, curvefiti)
}

# curvefit$estimate <- curvefit$b1 * exp(-curvefit$b2 * curvefit$lag) +
#   curvefit$b3 * exp(-curvefit$b4 * curvefit$lag)

ggplot(fitresult, aes(x = lag, y = estimate, color = maxlag)) +
  geom_point(alpha = 0.3) +
  # geom_line(aes(x = lag, y = estimate, color = maxlag), data = curvefit) +
  # geom_hline(yintercept = 0) +
  theme(legend.position="top")


#---------------------------------------------------------------------------------
## glmer model for the estimation of the exposure effect
## including the memory and depth as offset
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# importation of data and restriction to data after exposure
data <- dataAll

## We cut the first 24 hours of all individuals to avoid
## the effect of tagging
data <- data[data$Time > 24,]

#---------------------------------------------------------------------------------
## Defining exposures
data$Xtilde <- ExposureFunction(X = data$Dist_Ship) ## Exposure variable based on 1/distance in km
data$X = data$Xtilde*(data$Seismik == 1)
data$P = data$Xtilde*(data$Seismik == 0)

## Only data from airgun or no exposure at all (removing ship):
data <- subset(data, P == 0)


plotdist <- seq(2,50,0.1) # distances for plotting in km

#---------------------------------------------------------------------------------
## Estimation of the glmer model
# The generalized linear mixed model with a Poisson response distribution with a log-link to model the effect of exposure on buzzing rate (buzzes/min) with an autoregressive memory component. Exposure is defined as 1/distance (km).Exposure is entered non-linearly as an explanatory variable using natural cubic splines with 3 degrees of freedom (ns, package splines) with internal knots located at the 33th and 66th percentiles of the non-zero exposure values. Individual is included as a random effect allowing each animal to have a unique baseline (intercept) in their sound production rate. To obtain convergence the optimization is done by Adaptive Gauss-Hermite Quadrature, which is obtained by the option nAGQ = 0 in the glmer-call. The default is nAGQ = 1, the Laplace approximation, which does not reach convergence.

k=60; i=1
temp    <- glmresultsDepth[[i]]
lag     <- temp$lag
ARcoef  <- temp$coeftable[-(1:10), "Estimate"]
std.err <- temp$coeftable[-(1:10),"Std. Error"]
w       <- 1/std.err
n       <- length(ARcoef)

fm <- nls(ARcoef ~ SSbiexp(1:n, A1, lrc1, A2, lrc2), weights = w)
b1 <- summary(fm)$coefficients[1,1]
b2 <- exp(summary(fm)$coefficients[2,1])
b3 <- summary(fm)$coefficients[3,1]
b4 <- exp(summary(fm)$coefficients[4,1])

ARvec <- b1*exp(-b2*(1:60))+b3*exp(-b4*(1:60))
LagVariables <- names(data[,1:60])
dataAR  <- data[,LagVariables]
## Autoregressive component for offset in later analyses
data$ARDepth <- as.matrix(dataAR) %*% ARvec
## Depth coefficients for offset
Depthcoeff <- temp$coeftable[7:10,"Estimate"]
data$ARDepth <- data$ARDepth + as.matrix(ns(data$Depth, knots = c(-323, -158, -54)))%*% Depthcoeff

## Weights for the glmer analysis
data$n <- rep(0, length(data$Ind))
for(k in unique(data$Ind)){
  data$n[data$Ind == k] <- length(data$Ind[data$Ind == k])
}

glmerAllBuzzDepth <- glmer(Buzz ~ offset(ARDepth) + #ns(Depth, knots = c(-323, -158, -54)) +
                             ns(X, knots = quantile(data$X[data$X > 0],c(1:2)/3)) +
                             (1 | Ind),
                           data = data,
                           nAGQ=0,
                           weights = n,
                           family = poisson)

summary(glmerAllBuzzDepth)
 
# #---------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------
# ## Validation of the model with plots
# #---------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------
# #### Percentage of normal behavior
# 
# 
# # percentage of normal behavior
# ChangePop = predFramePop
# ChangePop$change = exp(ChangePop$predBuzzPop)
# ChangePop$change[1:length(plotdist)] = ChangePop$change[1:length(plotdist)]/exp(predFramePop0$predBuzzPop0[1])*100
# ChangePop$change[(length(plotdist)+1):(2*length(plotdist))] = ChangePop$change[(length(plotdist)+1):(2*length(plotdist))]/exp(predFramePop0$predBuzzPop0[2])*100
# 
# 
# ggplot(data = ChangePop,
#        aes(x = change, y = 1/X, color = model)) +
#   geom_line() +
#   ylim(0,40) + ylab("Distance to ship") +
#   xlab("Percentage of normal buzzing rate") +
#   #  facet_grid(~ model)  + 
#   theme(legend.position = "top", legend.title = element_blank())
# ggsave("PercentageNormalBehavior.pdf")
#  
# #### Model control for model with offset including AR +Depth   
# 
# data$predictDepth = predict(glmerAllBuzzDepth, type = "response")
# 
# ## Below we calculate the uniform residuals
# 
# Zall = list(NULL)
# m = 1
# #for(k in unique(data$Ind)){
# for (k in unique(data$Ind)){
#   datak = data[data$Ind == k,]
#   n     = length(datak$Buzz)
#   Z     = NULL  ## Uniform residuals
#   Buzzindices = (1:n)[datak$Buzz == 1]
#   nB          = length(Buzzindices)
#   dataki      = datak[1:Buzzindices[1],]
#   Z    = c(Z,(exp(-sum(dataki$predictDepth))))
#   for(i in 2:nB){
#     dataki = datak[(Buzzindices[i-1]+1):(Buzzindices[i]),]
#     Z    = c(Z, (1-exp(-sum(dataki$predictDepth))))
#   }
#   dataki = datak[Buzzindices[nB]:n,]
#   Z      = c(Z,(exp(-sum(dataki$predictDepth))))
#   Zall[[m]] = list(Z = Z, Ind = k)
#   m = m + 1
# }
# 
# 
# 
# ## Collect the uniform residuals in a data.frame
# Zdata = data.frame(Z = Zall[[1]]$Z, Ind = Zall[[1]]$Ind, 
#                    n = length(Zall[[1]]$Z))
# for(i in 2:6){
#   Zdata = rbind(Zdata, 
#                 data.frame(Z = Zall[[i]]$Z, Ind = Zall[[i]]$Ind, 
#                            n = length(Zall[[i]]$Z)))
# }
# 
# nn = length(Zdata$Z)
# 
# temp = acf(Zdata$Z, plot = FALSE)
# 
# acfplotdata = data.frame(acf = temp$acf, lag = temp$lag)
# 
# ggplot(data = acfplotdata, aes(x = lag, y = acf)) +
#   geom_col(size = 0.2) + coord_fixed(36) +
#   ylab("Autocorrelation of uniform residuals") +
#   xlab("Lag")
# 
# #ggsave("Zacf.png", width = 10, height = 8, units = "cm")
# 
# zplotdata = data.frame(Zlow = Zdata$Z[1:(nn-1)], Zupp = Zdata$Z[2:nn])
# 
# ggplot(data = zplotdata, aes(x = Zlow, y = Zupp)) +
#   geom_point(size = 0.5) + coord_fixed(1) +
#   xlab(expression(Z[i-1])) +
#   ylab(expression(Z[i])) +
#   geom_abline(slope = 1, intercept = 0, size = 1)
# #ggsave("Zplot.png", width = 10, height = 8, units = "cm")
# 
# 
# # QQplot for a given depth
# 
# Zall = list(NULL)
# m = 1
# for(k in unique(data$Ind)){
#   datak = data[data$Ind == k,]
#   n     = length(datak$Buzz)
#   Z     = NULL  ## Uniform residuals
#   Depthk = NULL
#   Xk = NULL
#   Pk = NULL
#   Buzzindices = (1:n)[datak$Buzz == 1]
#   nB          = length(Buzzindices)
#   for(i in 2:nB){
#     dataki = datak[(Buzzindices[i-1]+1):(Buzzindices[i]),]
#     Z    = c(Z, (1-exp(-sum(dataki$predictDepth))))
#     Depthk = c(Depthk, dataki$Depth[1])
#     Xk     = c(Xk, dataki$X[1])
#     Pk     = c(Pk, dataki$P[1])
#   }
#   Zall[[m]] = list(Z = Z, Ind = k, X = Xk, P = Pk , Depth = Depthk)
#   m = m + 1
# }
# 
# 
# ## Collect the uniform residuals in a data.frame
# Zdata = data.frame(Z = Zall[[1]]$Z, 
#                    Ind = Zall[[1]]$Ind, 
#                    n = length(Zall[[1]]$Z), 
#                    Depth = Zall[[1]]$Depth,
#                    X = Zall[[1]]$X, P = Zall[[1]]$P)
# for(i in 2:6){
#   Zdata = rbind(Zdata, 
#                 data.frame(Z = Zall[[i]]$Z, 
#                            Ind = Zall[[i]]$Ind, 
#                            n = length(Zall[[i]]$Z), 
#                            Depth = Zall[[i]]$Depth,
#                            X = Zall[[i]]$X, P = Zall[[i]]$P))
# }
# 
# ZdataDepthbt = Zdata[-Zdata$Depth>400 & Zdata$X==0 &Zdata$P==0, ]
# #ZdataDepthbt = Zdata[ Zdata$X==0 &Zdata$P==0, ]
# nresid = length(ZdataDepthbt$Z)
# Zorder = ZdataDepthbt$Z[order(ZdataDepthbt$Z)]
# 
# df = data.frame(qunif = (1:nresid)/nresid,
#                 qZ = Zorder)
# 
# ggplot(df, aes(x=qunif, y=Zorder))+
#   geom_point()+
#   xlab("Uniform theoretical quantiles")+
#   ylab("Empirical quantiles")
# #ggsave("QuantilesZ.png", width = 10, height = 8, units = "cm")
# 
# 
