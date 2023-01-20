#--------------------------------------------------------------------------------
## Objective : fit glmer model using fixed AR coefs
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(splines)
library(lme4)
library(RhpcBLASctl)
blas_set_num_threads(1)
source("utils/0_data.R")
source("utils/1_glmer_buzz_AR_expo.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 2) {
  print("Usage du script : Rscript 1_glmer_buzz_AR_expo.R arg1 arg2")
  print("arg1 : le chemin vers la base de données")
  print("arg2 : le chemin vers le dossier de sauvegarde des objets R")
  stop("Des arguments doivent être donnés au script.", call. = FALSE)
}

#---------------------------------------------------------------------------------
# Data preparation
dataAll <- fread(args[1])
data <- dataAll
## keep the 2018 data
data <- data[data$Year == "2018", ]
## Add time since tagging
data <- AddTime(data)
## We restrict the data to be after exposure and 24 hours after tagging
data <- AfterExposure(data, no_stress = TRUE)
## We add exposure
data <- AddExposure(data)
### We restrict to airgun expositions
data <- OnlyAirgun(data)
## Weights for the glmer analysis
data$n <- rep(0, length(data$Ind))
for (k in unique(data$Ind)) {
  data$n[data$Ind == k] <- length(data$Ind[data$Ind == k])
}

#---------------------------------------------------------------------------------
# Estimation of the glmer model
# The generalized linear mixed model with a Poisson response distribution with a log-link to model the effect of exposure on buzzing rate (buzzes/min) with an autoregressive memory component.
# Exposure is defined as 1/distance (km).
# Exposure is entered non-linearly as an explanatory variable using natural cubic splines with 3 degrees of freedom (ns, package splines) with internal knots located at the 33th and 66th percentiles of the non-zero exposure values.
# Individual is included as a random effect allowing each animal to have a unique baseline (intercept) in their sound production rate.
# To obtain convergence the optimization is done by Adaptive Gauss-Hermite Quadrature, which is obtained by the option nAGQ = 0 in the glmer-call.
# The default is nAGQ = 1, the Laplace approximation, which does not reach convergence.

maxlag.bic <- readRDS("../data/glmER_buzz_depth_maxlag/maxlag.bic.rds")
ARcoef.best <- readRDS("../data/glmER_buzz_depth_maxlag/ARcoef.best.rds")
ARcoef.RegBiExp <- readRDS("../data/glmER_buzz_depth_maxlag/ARcoef.RegBiExp.rds")

maxlag.opt <- as.integer(maxlag.bic[which.min(maxlag.bic[, 2]), 1])

## Set ARcoef using optimal max lag
### Define the first maxlag.opt lags
temp <- as.data.frame(shift(data$Buzz, n = 1:maxlag.opt, give.names = TRUE))
data <- cbind(temp, data)
LagVariables <- names(data[, 1:maxlag.opt])
dataAR <- data[, LagVariables]

### Autoregressive component for offset
ARvec <- ARcoef.RegBiExp$estimate[1] * exp(-exp(ARcoef.RegBiExp$estimate[2]) * (1:maxlag.opt)) +
  ARcoef.RegBiExp$estimate[3] * exp(-exp(ARcoef.RegBiExp$estimate[4]) * (1:maxlag.opt))

glmerAllBuzz <- glmer_buzz_ARDepth_expo(data, dataAR, ARvec, Depthcoeff)

glmerAllBuzz.tidy <- tidy(glmerAllBuzz) # to save
glmerAllBuzz.glance <- glance(glmerAllBuzz) # to save

summary(glmerAllBuzz)

#---------------------------------------------------------------------------------
# Save R objects
saveRDS(glmerAllBuzz, paste0(args[2], "/glmerAllBuzz.rds"))
saveRDS(glmerAllBuzz.tidy, paste0(args[2], "/glmerAllBuzz.tidy.rds"))
saveRDS(glmerAllBuzz.glance, paste0(args[2], "/glmerAllBuzz.glance.rds"))
