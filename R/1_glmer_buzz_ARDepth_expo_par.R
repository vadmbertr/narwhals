#--------------------------------------------------------------------------------
## Objective : fit several glmer while drawing fixed coef and save expo coef
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(lme4)
library(mvtnorm)
library(splines)
library(parallel)
library(RhpcBLASctl)
source("0_data.R")
source("0_biexp.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 4) {
  print("Usage du script : Rscript 1_glmer_buzz_ARDepth_expo.R arg1 arg2")
  print("arg1 : le chemin vers la base de données")
  print("arg2 : le chemin vers le dossier de sauvegarde des objets R")
  print("arg3 : nombre de coefficients estimés souhaités")
  print("arg4 : le nombre de coeurs alloués") # cannot be retrieve from R
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
### Remove NA
# data <- RemoveNA(data, c("Ind", "Buzz", "Depth", "X"))

#---------------------------------------------------------------------------------
# Estimation of the glmer model
# The generalized linear mixed model with a Poisson response distribution with a log-link to model the effect of exposure on buzzing rate (buzzes/min) with an autoregressive memory component.
# Exposure is defined as 1/distance (km).
# Exposure is entered non-linearly as an explanatory variable using natural cubic splines with 3 degrees of freedom (ns, package splines) with internal knots located at the 33th and 66th percentiles of the non-zero exposure values.
# Individual is included as a random effect allowing each animal to have a unique baseline (intercept) in their sound production rate.
# To obtain convergence the optimization is done by Adaptive Gauss-Hermite Quadrature, which is obtained by the option nAGQ = 0 in the glmer-call.
# The default is nAGQ = 1, the Laplace approximation, which does not reach convergence.

maxlag.bic <- readRDS("../data/glmER_buzz_depth_maxlag/maxlag.bic.rds")
glmer.coefs <- readRDS("../data/glmER_buzz_depth_maxlag/ARcoef.best.rds")
biexp.coef.estimate <- readRDS("../data/glmer_biexp_AR_mc/biexp.coef.estimate.rds")
biexp.coef.cov <- readRDS("../data/glmer_biexp_AR_mc/biexp.coef.cov.rds")
Depth.vcov <- as.matrix(readRDS("../data/glmer_buzz_ARDepth/glmERBuzzARDepth.vcov.rds"))[2:5, 2:5]

maxlag.opt <- as.integer(maxlag.bic[which.min(maxlag.bic[, 2]), 1])
biexp.coef.estimate <- apply(biexp.coef.estimate, 2, mean)

## Set ARcoef using optimal max lag
### Define the first maxlag.opt lags
temp <- as.data.frame(shift(data$Buzz, n = 1:maxlag.opt, give.names = TRUE))
data <- cbind(temp, data)
LagVariables <- names(data[, 1:maxlag.opt])
dataAR <- data[, LagVariables]

fit.glmer <- function (i) {
  ### Autoregressive component for offset
  ARcoefs <- as.numeric(rmvnorm(1, mean = biexp.coef.estimate, sigma = biexp.coef.cov,
                                checkSymmetry = FALSE))
  ARvec <- BiExp(ARcoefs[[1]], ARcoefs[[3]], ARcoefs[[2]], ARcoefs[[4]], , maxlag = maxlag.opt)
  data$ARDepth <- as.matrix(dataAR) %*% ARvec

  ### Depth coefficients for offset
  Depthcoefs <- as.numeric(rmvnorm(1, mean = glmer.coefs$estimate[2:5], sigma = Depth.vcov,
                                   checkSymmetry = FALSE))
  data$ARDepth <- data$ARDepth + as.matrix(ns(data$Depth, knots = c(-323, -158, -54))) %*% Depthcoefs

  ## Weights for the glmer analysis
  data$n <- rep(0, length(data$Ind))
  for (k in unique(data$Ind)) {
    data$n[data$Ind == k] <- length(data$Ind[data$Ind == k])
  }

  glmerAllBuzzDepth <- glmer(Buzz ~ offset(ARDepth) + ns(X, knots = quantile(data$X[data$X > 0], 1:2 / 3)) + (1 | Ind),
                             data = data,
                             nAGQ = 0,
                             weights = n,
                             family = poisson)
  coefs <- tidy(glmerAllBuzzDepth)
  if (!any(grepl("Error", coefs$term, fixed = T))) {
    return(coefs[, c("term", "estimate", "std.error")])
  }
}

n.done <- function (df, nc) {
  if (is.na(nc)) return(0)
  else return(round(nrow(df) / nc))
}

expo.coef.path <- paste0(args[2], "/expo.coef.rds")
if (file.exists(expo.coef.path)) {
  expo.coef.all <- readRDS(expo.coef.path)
  n.coefs <- length(unique(expo.coef.all$term))
} else {
  expo.coef.all <- data.frame()
  n.coefs <- NA
}

## Parallelism
blas_set_num_threads(1)
### allocated RAM = RAM total * allocated cores / total cores
### n glmer // = min(allocated cores, allocated RAM / RAM per glmer)
ram.total <- 192
ram.per.job <- 7.5
n.cores <- detectCores()
n.cores.alloc <- as.numeric(args[4])
ram.alloc <- ram.total * n.cores.alloc / n.cores
n.jobs <- min(n.cores.alloc,
              floor(ram.alloc / ram.per.job),
              as.numeric(args[3]) - n.done(expo.coef.all, n.coefs))

while (as.numeric(args[3]) - n.done(expo.coef.all, n.coefs) > 0) {
  print(n.done(expo.coef.all, n.coefs))
  n.jobs <- min(n.jobs, as.numeric(args[3]) - n.done(expo.coef.all, n.coefs))
  expo.coef <- do.call(rbind, mclapply(1:n.jobs, fit.glmer, mc.cores = n.jobs))
  if (is.na(n.coefs)) {
    n.coefs <- as.integer(nrow(expo.coef) / n.jobs)
  }
  expo.coef.all <- rbind(expo.coef.all, expo.coef)
  saveRDS(expo.coef.all, expo.coef.path)
}
