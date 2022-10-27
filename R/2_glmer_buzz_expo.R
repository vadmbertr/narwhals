#--------------------------------------------------------------------------------
## Objective : find the optimal (BIC wise) memory lag
#---------------------------------------------------------------------------------

library(data.table)
library(splines)
library(lme4)
source("0_data.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read db path from command line
if (length(args) != 5) {
  print("Usage du script : Rscript 2_glmer_buzz_expo.R arg1 arg2 arg3 arg4 arg5")
  print("arg1 : le chemin vers la base de données")
  print("arg2 : la borne inférieure de recherche du lag maximum")
  print("arg3 : la borne supérieure de recherche du lag maximum")
  print("arg4 : le nombre de lag à considérer entre les deux bornes")
  print("arg5 : le chemin vers le dossier de sauvegarde des objets R")
  stop("Des arguments doivent être donnés au script.", call. = FALSE)
}

#---------------------------------------------------------------------------------
# Data preparation
dataAll <- fread(args[1])
data <- dataAll
## We restrict the data to be after exposure and 24 hours after tagging
data <- AfterExposure(data, no_stress = TRUE)
## We add exposure
data <- AddExposure(data)
### We restrict to airgun expositions
data <- OnlyAirgun(data)

#---------------------------------------------------------------------------------
# Estimation of the glmer model
# The generalized linear mixed model with a Poisson response distribution with a log-link to model the effect of exposure on buzzing rate (buzzes/min) with an autoregressive memory component.
# Exposure is defined as 1/distance (km).
# Exposure is entered non-linearly as an explanatory variable using natural cubic splines with 3 degrees of freedom (ns, package splines) with internal knots located at the 33th and 66th percentiles of the non-zero exposure values.
# Individual is included as a random effect allowing each animal to have a unique baseline (intercept) in their sound production rate.
# To obtain convergence the optimization is done by Adaptive Gauss-Hermite Quadrature, which is obtained by the option nAGQ = 0 in the glmer-call.
# The default is nAGQ = 1, the Laplace approximation, which does not reach convergence.

maxlag.bic <- readRDS("../data/glm_buzz_depth_maxlag/maxlag.bic.rds")
ARcoef.RegBiExp <- readRDS("../data/glm_buzz_depth_maxlag/ARcoef.RegBiExp.rds")

maxlag.opt <- as.integer(maxlag.bic[which.min(maxlag.bic[, 2]), 1])
ARvec <- ARcoef.RegBiExp$b1 * exp(-ARcoef.RegBiExp$b2 * (1:maxlag.opt)) +
  ARcoef.RegBiExp$b3 * exp(-ARcoef.RegBiExp$b4 * (1:maxlag.opt))
LagVariables <- names(data[, 1:maxlag.opt])
dataAR <- data[, LagVariables]

## Autoregressive component for offset in later analyses
data$ARDepth <- as.matrix(dataAR) %*% ARvec

## Depth coefficients for offset
Depthcoeff <- temp$coeftable[7:10, "Estimate"]
data$ARDepth <- data$ARDepth + as.matrix(ns(data$Depth, knots = c(-323, -158, -54))) %*% Depthcoeff

## Weights for the glmer analysis
data$n <- rep(0, length(data$Ind))
for (k in unique(data$Ind)) {
  data$n[data$Ind == k] <- length(data$Ind[data$Ind == k])
}

glmerAllBuzzDepth <- glmer(Buzz ~ offset(ARDepth) +
  ns(X, knots = quantile(data$X[data$X > 0], 1:2 / 3)) +
  (1 | Ind),
                           data = data,
                           nAGQ = 0,
                           weights = n,
                           family = poisson)

#---------------------------------------------------------------------------------
# For visual model validation

## Model control for model with offset including AR +Depth
predictDepth <- predict(glmerAllBuzzDepth, type = "response")
### Uniform residuals
Zall <- list(NULL)
m <- 1
for (k in unique(data$Ind)) {
  datak <- data[data$Ind == k,]
  n <- length(datak$Buzz)
  Z <- NULL  ## Uniform residuals
  Buzzindices <- (1:n)[datak$Buzz == 1]
  nB <- length(Buzzindices)
  dataki <- datak[1:Buzzindices[1],]
  Z <- c(Z, (exp(-sum(dataki$predictDepth))))
  for (i in 2:nB) {
    dataki <- datak[(Buzzindices[i - 1] + 1):(Buzzindices[i]),]
    Z <- c(Z, (1 - exp(-sum(dataki$predictDepth))))
  }
  dataki <- datak[Buzzindices[nB]:n,]
  Z <- c(Z, (exp(-sum(dataki$predictDepth))))
  Zall[[m]] <- list(Z = Z, Ind = k)
  m <- m + 1
}
Zdata <- data.frame(Z = Zall[[1]]$Z, Ind = Zall[[1]]$Ind,
                    n = length(Zall[[1]]$Z))
for (i in 2:6) {
  Zdata <- rbind(Zdata,
                 data.frame(Z = Zall[[i]]$Z, Ind = Zall[[i]]$Ind,
                            n = length(Zall[[i]]$Z)))
}
temp <- acf(Zdata$Z, plot = FALSE)
acfplotdata <- data.frame(acf = temp$acf, lag = temp$lag) # to save
zplotdata <- data.frame(Zlow = Zdata$Z[1:(nn - 1)], Zupp = Zdata$Z[2:nn]) # to save

# QQplot for a given depth
Zall <- list(NULL)
m <- 1
for (k in unique(data$Ind)) {
  datak <- data[data$Ind == k,]
  n <- length(datak$Buzz)
  Z <- NULL  ## Uniform residuals
  Depthk <- NULL
  Xk <- NULL
  Pk <- NULL
  Buzzindices <- (1:n)[datak$Buzz == 1]
  nB <- length(Buzzindices)
  for (i in 2:nB) {
    dataki <- datak[(Buzzindices[i - 1] + 1):(Buzzindices[i]),]
    Z <- c(Z, (1 - exp(-sum(dataki$predictDepth))))
    Depthk <- c(Depthk, dataki$Depth[1])
    Xk <- c(Xk, dataki$X[1])
    Pk <- c(Pk, dataki$P[1])
  }
  Zall[[m]] <- list(Z = Z, Ind = k, X = Xk, P = Pk, Depth = Depthk)
  m <- m + 1
}
Zdata <- data.frame(Z = Zall[[1]]$Z,
                    Ind = Zall[[1]]$Ind,
                    n = length(Zall[[1]]$Z),
                    Depth = Zall[[1]]$Depth,
                    X = Zall[[1]]$X, P = Zall[[1]]$P)
for (i in 2:6) {
  Zdata <- rbind(Zdata,
                 data.frame(Z = Zall[[i]]$Z,
                            Ind = Zall[[i]]$Ind,
                            n = length(Zall[[i]]$Z),
                            Depth = Zall[[i]]$Depth,
                            X = Zall[[i]]$X, P = Zall[[i]]$P))
}
ZdataDepthbt <- Zdata[-Zdata$Depth > 400 & Zdata$X == 0 & Zdata$P == 0,]
nresid <- length(ZdataDepthbt$Z)
Zorder <- ZdataDepthbt$Z[order(ZdataDepthbt$Z)]
df <- data.frame(qunif = (1:nresid) / nresid, qZ = Zorder) # to save

#---------------------------------------------------------------------------------
# Save R objects
# TODO
saveRDS(dataBICDepth, paste0(args[5], "/glmer_buzz_depth_maxlag_bic.rds"))
saveRDS(coeftable, paste0(args[5], "/glmer_buzz_depth_bestmaxlag_coefs.rds"))
