#--------------------------------------------------------------------------------
## Objective : retrieve covariance matrix of the Depth coefficients
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(splines)
library(lme4)
source("utils/data.R")
source("utils/biexp.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read db path from command line
if (length(args)!=2) {
  print("Usage du script : RScript glmER_buzz_depth_maxlag.R arg1 arg2")
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
## We restrict the data to be before exposure
data <- BeforeExposure(data)

#---------------------------------------------------------------------------------
# Estimation of the model

maxlag.bic <- readRDS("../data/1_buzzing/glmER_buzz_depth_maxlag/maxlag.bic.rds")
ARcoef.RegBiExp <- readRDS("../data/1_buzzing/glmER_buzz_depth_maxlag/ARcoef.RegBiExp.rds")

maxlag.opt <- as.integer(maxlag.bic[which.min(maxlag.bic[, 2]), 1])

## Set ARcoef using optimal max lag
### Define the first maxlag.opt lags
temp <- as.data.frame(shift(data$Buzz, n = 1:maxlag.opt, give.names = TRUE))
data <- cbind(temp, data)
LagVariables <- names(data[, 1:maxlag.opt])
dataAR <- data[, LagVariables]

ARvec <- BiExp(ARcoef.RegBiExp$estimate[1], ARcoef.RegBiExp$estimate[2],
               ARcoef.RegBiExp$estimate[3], ARcoef.RegBiExp$estimate[4], maxlag = maxlag.opt)
data$AR <- as.matrix(dataAR) %*% ARvec

## Fit a glmer
glmERBuzzARDepth <- glmer(Buzz ~ offset(AR) + ns(Depth, knots = c(-323, -158, -54)) + (1 | Ind),
                           data = data,
                           family = poisson)

summary(glmERBuzzARDepth)

#---------------------------------------------------------------------------------
# Save R objects
saveRDS(tidy(glmERBuzzARDepth), paste0(args[2], "/glmERBuzzARDepth.tidy.rds"))
saveRDS(glance(glmERBuzzARDepth), paste0(args[2], "/glmERBuzzARDepth.glance.rds"))
saveRDS(vcov(glmERBuzzARDepth), paste0(args[2], "/glmERBuzzARDepth.vcov.rds"))
