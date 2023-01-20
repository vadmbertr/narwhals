#--------------------------------------------------------------------------------
## Objective : generate and save obj for plotting interval bands using MC procedure
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(splines)
library(lme4)
library(RhpcBLASctl)
source("utils/0_data.R")
source("utils/1_glmer_buzz_ARDepth_expo.R")
source("utils/1_glmer_buzz_AR_expo.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 4) {
  print("Usage du script : Rscript 1_glmer_buzz_ARDepth_expo_pred.R arg1 arg2 arg3 arg4")
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
data <- data[data$Year == "2018",]
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
# Estimation of the prediction bands

maxlag.bic <- readRDS("../data/glmER_buzz_depth_maxlag/maxlag.bic.rds")
coefs.estimate <- readRDS("../data/glmER_buzz_depth_maxlag/ARcoef.best.rds")
coefs.vcov <- as.matrix(readRDS("../data/glmer_buzz_ARDepth/glmERBuzzARDepth.vcov.rds"))

maxlag.opt <- as.integer(maxlag.bic[which.min(maxlag.bic[, 2]), 1])
coefs.idx <- 1:(4 + maxlag.opt) + 1
coefs.estimate <- coefs.estimate$estimate[coefs.idx]
coefs.vcov <- coefs.vcov[coefs.idx, coefs.idx]

## Set ARcoef using optimal max lag
### Define the first maxlag.opt lags
temp <- as.data.frame(shift(data$Buzz, n = 1:maxlag.opt, give.names = TRUE))
data <- cbind(temp, data)
LagVariables <- names(data[, 1:maxlag.opt])
dataAR <- data[, LagVariables]

## prediction matrices
plotdist <- seq(1, 70, 0.1) ## Plotting distances in km
Depthlevel <- -400
### ind - with expo - no Depth
predFrameMany <- NULL
for (k in unique(data$Ind)) {
  temp <- range(data$X[data$X > 0 & data$Ind == k])
  temp1 <- expand.grid(X = 1 / seq(1 / temp[2], 1 / temp[1], 0.1),
                       P = 0,
                       Depth = Depthlevel,
                       AR = 0,
                       Ind = k)
  predFrameMany <- rbind(predFrameMany, temp1)
}
### ind - with expo - with Depth
predFrameManyDepth <- NULL
Depthlevel <- -400  # TODO: ask Ad on this: Depth is not used directly but rather AR
for (k in unique(data$Ind)) {
  temp <- range(data$X[data$X > 0 & data$Ind == k])
  temp1 <- expand.grid(X = 1 / seq(1 / temp[2], 1 / temp[1], 0.1),
                       P = 0,
                       Depth = Depthlevel,
                       ARDepth = 0,
                       Ind = k)
  predFrameManyDepth <- rbind(predFrameManyDepth, temp1)
}
### pop - with expo - no Depth
predFramePop_temp <- expand.grid(X = 1 / plotdist,
                                 P = 0,
                                 Depth = Depthlevel,
                                 AR = 0,
                                 Ind = "Population")
### pop - with expo - no Depth
predFramePopDepth_temp <- expand.grid(X = 1 / plotdist,
                                      P = 0,
                                      Depth = Depthlevel,
                                      ARDepth = 0,
                                      Ind = "Population")
### ind - no expo - no Depth
predFrame0_temp <- expand.grid(X = 0,
                               P = 0,
                               Depth = Depthlevel,
                               AR = 0,
                               Ind = unique(data$Ind))
### ind - no expo - with Depth
predFrame0Depth_temp <- expand.grid(X = 0,
                                    P = 0,
                                    Depth = Depthlevel,
                                    ARDepth = 0,
                                    Ind = unique(data$Ind))
### pop - no expo - no Depth
predFramePop0_temp <- expand.grid(X = 0,
                                  P = 0,
                                  Depth = Depthlevel,
                                  AR = 0,
                                  Ind = "Population")
### pop - no expo - with Depth
predFramePop0Depth_temp <- expand.grid(X = 0,
                                  P = 0,
                                  Depth = Depthlevel,
                                  ARDepth = 0,
                                  Ind = "Population")

compute.pred.band <- function (i) {
  ### Components for offset
  coefs <- as.numeric(rmvnorm(1, mean = coefs.estimate, sigma = coefs.vcov, checkSymmetry = FALSE))
  ### Autoregressive component for offset
  ARcoefs <- coefs[1:maxlag.opt + 4]
  ### Depth coefficients for offset
  Depthcoefs <- coefs[1:4]

  ## fitting
  ### no depth
  glmerAllBuzz <- glmer_buzz_AR_expo(data, dataAR, ARcoefs)
  ### with depth
  glmerAllBuzzDepth <- glmer_buzz_ARDepth_expo(data, dataAR, ARcoefs, Depthcoefs)

  ## predictions
  ### individual prediction wrt exposure X for glmerAllBuzz
  predBuzz <- predict(glmerAllBuzz,
                      newdata = predFrameMany)
  predFrame <- cbind(predFrameMany, as.data.frame(predBuzz))
  predFrame$exposure <- "Trial"
  predFrame$model <- "Without Depth"
  ### individual prediction wrt exposure X for glmerAllBuzzDepth
  predBuzz <- predict(glmerAllBuzzDepth,
                      newdata = predFrameManyDepth)
  temp <- cbind(predFrameManyDepth, as.data.frame(predBuzz))
  temp$exposure <- "Trial"
  temp$model <- "With Depth"
  names(temp) <- c("X", "P", "Depth", "AR", "Ind", "predBuzz", "exposure", "model")
  predFrame <- rbind(predFrame, temp) # to save

  # population prediction wrt exposure X for glmerAllBuzz
  predBuzzPop <- predict(glmerAllBuzz,
                         newdata = predFramePop_temp,
                         re.form = NA)
  predFramePop <- cbind(predFramePop_temp, as.data.frame(predBuzzPop))
  predFramePop$model <- "Without Depth"
  # population prediction wrt exposure X for glmerAllBuzzDepth
  predBuzzPop <- predict(glmerAllBuzzDepth,
                         newdata = predFramePopDepth_temp,
                         re.form = NA)
  temp <- cbind(predFramePop_temp, as.data.frame(predBuzzPop))
  temp$model <- "With Depth"
  names(temp) <- c("X", "P", "Depth", "AR", "Ind", "predBuzzPop", "model")
  predFramePop <- rbind(predFramePop, temp) # to save

  # individual prediction with no exposure X for glmerAllBuzz
  predBuzz0 <- predict(glmerAllBuzz,
                       newdata = predFrame0_temp)
  predFrame0 <- cbind(predFrame0_temp, as.data.frame(predBuzz0))
  predFrame0$model <- "Without Depth"
  # individual prediction with no exposure X for glmerAllBuzzDepth
  predBuzz0 <- predict(glmerAllBuzzDepth,
                       newdata = predFrame0Depth_temp)
  temp <- cbind(predFrame0Depth_temp, as.data.frame(predBuzz0))
  temp$model <- "With Depth"
  names(temp) <- c("X", "P", "Depth", "AR", "Ind", "predBuzz0", "model")
  predFrame0 <- rbind(predFrame0, temp) # to save

  # population prediction with no exposure X for glmerAllBuzz
  predBuzzPop0 <- predict(glmerAllBuzz,
                          newdata = predFramePop0_temp,
                          re.form = NA)
  predFramePop0 <- cbind(predFramePop0_temp, as.data.frame(predBuzzPop0))
  predFramePop0$model <- "Without Depth"
  # population prediction with no exposure X for glmerAllBuzzDepth
  predBuzzPop0 <- predict(glmerAllBuzzDepth,
                          newdata = predFramePop0Depth_temp,
                          re.form = NA)
  temp <- cbind(predFramePop0Depth_temp, as.data.frame(predBuzzPop0))
  temp$model <- "With Depth"
  names(temp) <- c("X", "P", "Depth", "AR", "Ind", "predBuzzPop0", "model")
  predFramePop0 <- rbind(predFramePop0, temp) # to save

  # percentage of normal behavior
  ChangePop <- predFramePop # to save
  ChangePop$change <- exp(ChangePop$predBuzzPop)
  ChangePop$change[seq_along(plotdist)] <-
    ChangePop$change[seq_along(plotdist)] / exp(predFramePop0$predBuzzPop0[1]) * 100
  ChangePop$change[(length(plotdist) + 1):(2 * length(plotdist))] <-
    ChangePop$change[(length(plotdist) + 1):(2 * length(plotdist))] / exp(predFramePop0$predBuzzPop0[2]) * 100

  return(list(predFrame=predFrame, predFramePop=predFramePop, predFrame0=predFrame0, predFramePop0=predFramePop0,
              ChangePop=ChangePop))
}

n.done <- function (ml, nm) {
  if (is.na(nm)) return(0)
  else return(length(ml))
}

if (file.exists(paste0(args[2], "/predFrame.all.rds"))) {
  predFrame.all <- readRDS(paste0(args[2], "/predFrame.all.rds"))
  predFramePop.all <- readRDS(paste0(args[2], "/predFramePop.all.rds"))
  predFrame0.all <- readRDS(paste0(args[2], "/predFrame0.all.rds"))
  predFramePop0.all <- readRDS(paste0(args[2], "/predFramePop0.all.rds"))
  ChangePop.all <- readRDS(paste0(args[2], "/ChangePop.all.rds"))
  n.preds <- length(mat.list.all)
} else {
  predFrame.all <- list()
  predFramePop.all <- list()
  predFrame0.all <- list()
  predFramePop0.all <- list()
  ChangePop.all <- list()
  n.preds <- NA
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
              as.numeric(args[3]) - n.done(mat.list.all, n.preds))

while (as.numeric(args[3]) - n.done(mat.list.all, n.preds) > 0) {
  print(n.done(mat.list.all, n.preds))
  n.jobs <- min(n.jobs, as.numeric(args[3]) - n.done(mat.list.all, n.preds))
  res.list <- mclapply(1:n.jobs, compute.pred.band, mc.cores = n.jobs)
  if (is.na(n.preds)) {
    n.preds <- length(res.list)
  }

  for (i in seq_along(res.list)) {
    predFrame.all <- c(predFrame.all, res.list$predFrame)
    predFramePop.all <- c(predFramePop.all, res.list$predFramePop)
    predFrame0.all <- c(predFrame0.all, res.list$predFrame0)
    predFramePop0.all <- c(predFramePop0.all, res.list$predFramePop0)
    ChangePop.all <- c(ChangePop.all, res.list$ChangePop)
  }

  # Save objects
  saveRDS(predFrame.all, paste0(args[2], "/predFrame.all.rds"))
  saveRDS(predFramePop.all, paste0(args[2], "/predFramePop.all.rds"))
  saveRDS(predFrame0.all, paste0(args[2], "/predFrame0.all.rds"))
  saveRDS(predFramePop0.all, paste0(args[2], "/predFramePop0.all.rds"))
  saveRDS(ChangePop.all, paste0(args[2], "/ChangePop.all.rds"))
}

#---------------------------------------------------------------------------------
# Save objects
predFrame <- simplify2array(predFrame.all)
predFramePop <- simplify2array(predFramePop.all)
predFrame0 <- simplify2array(predFrame0.all)
predFramePop0 <- simplify2array(predFramePop0.all)
ChangePop <- simplify2array(ChangePop.all)
saveRDS(predFrame, paste0(args[2], "/predFrame.rds"))
saveRDS(predFramePop, paste0(args[2], "/predFramePop.rds"))
saveRDS(predFrame0, paste0(args[2], "/predFrame0.rds"))
saveRDS(predFramePop0, paste0(args[2], "/predFramePop0.rds"))
saveRDS(ChangePop, paste0(args[2], "/ChangePop.rds"))
