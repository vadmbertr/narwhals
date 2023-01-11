#--------------------------------------------------------------------------------
## Objective : generate and save obj for plots
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(splines)
library(lme4)
library(RhpcBLASctl)
blas_set_num_threads(1)
source("0_data.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 2) {
  print("Usage du script : Rscript 1_glmer_buzz_ARDepth_expo_plots.R arg1 arg2")
  print("arg1 : le chemin vers la base de données")
  print("arg2 : le chemin vers le dossier de sauvegarde des objets R")
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
### Remove NA
# data <- RemoveNA(data, c("Ind", "Buzz", "Depth", "X"))

#---------------------------------------------------------------------------------
# Load models
glmerAllBuzz <- readRDS(paste0(dirname(args[2]), "/glmer_buzz_AR_expo/glmerAllBuzz.rds"))
glmerAllBuzzDepth <- readRDS(paste0(dirname(args[2]), "/glmer_buzz_ARDepth_expo/glmerAllBuzzDepth.rds"))

#---------------------------------------------------------------------------------
# Obj for ploting
plotdist <- seq(1, 70, 0.1) ## Plotting distances in km
# individual prediction wrt exposure X for glmerAllBuzz
predFrameMany <- NULL
Depthlevel <- -400
for (k in unique(data$Ind)) {
  temp <- range(data$X[data$X > 0 & data$Ind == k])
  temp1 <- expand.grid(X = 1 / seq(1 / temp[2], 1 / temp[1], 0.1),
                       P = 0,
                       Depth = Depthlevel,
                       AR = 0,
                       Ind = k)
  predFrameMany <- rbind(predFrameMany, temp1)
}
predBuzz <- predict(glmerAllBuzz,
                    newdata = predFrameMany)
predFrame <- cbind(predFrameMany, as.data.frame(predBuzz))
predFrame$exposure <- "Trial"
predFrame$model <- "Without Depth"

# individual prediction wrt exposure X for glmerAllBuzzDepth
predFrameMany <- NULL
for (k in unique(data$Ind)) {
  temp <- range(data$X[data$X > 0 & data$Ind == k])
  temp1 <- expand.grid(X = 1 / seq(1 / temp[2], 1 / temp[1], 0.1),
                       P = 0,
                       Depth = Depthlevel,
                       ARDepth = 0,
                       Ind = k)
  predFrameMany <- rbind(predFrameMany, temp1)
}
predBuzz <- predict(glmerAllBuzzDepth,
                    newdata = predFrameMany)
temp <- cbind(predFrameMany, as.data.frame(predBuzz))
temp$exposure <- "Trial"
temp$model <- "With Depth"
names(temp) <- c("X", "P", "Depth", "AR", "Ind", "predBuzz", "exposure", "model")
predFrame <- rbind(predFrame, temp) # to save

# population prediction wrt exposure X for glmerAllBuzz
predFramePop_temp <- expand.grid(X = 1 / plotdist,
                                 P = 0,
                                 Depth = Depthlevel,
                                 AR = 0,
                                 Ind = "Population")
predBuzzPop <- predict(glmerAllBuzz,
                       newdata = predFramePop_temp,
                       re.form = NA)
predFramePop <- cbind(predFramePop_temp, as.data.frame(predBuzzPop))
predFramePop$model <- "Without Depth"
# population prediction wrt exposure X for glmerAllBuzzDepth
predFramePop_temp <- expand.grid(X = 1 / plotdist,
                                 P = 0,
                                 Depth = Depthlevel,
                                 ARDepth = 0,
                                 Ind = "Population")
predBuzzPop <- predict(glmerAllBuzzDepth,
                       newdata = predFramePop_temp,
                       re.form = NA)
temp <- cbind(predFramePop_temp, as.data.frame(predBuzzPop))
temp$model <- "With Depth"
names(temp) <- c("X", "P", "Depth", "AR", "Ind", "predBuzzPop", "model")
predFramePop <- rbind(predFramePop, temp) # to save

# individual prediction with no exposure X for glmerAllBuzz
predFrame0_temp <- expand.grid(X = 0,
                               P = 0,
                               Depth = Depthlevel,
                               AR = 0,
                               Ind = unique(data$Ind))
predBuzz0 <- predict(glmerAllBuzz,
                     newdata = predFrame0_temp)
predFrame0 <- cbind(predFrame0_temp, as.data.frame(predBuzz0))
predFrame0$model <- "Without Depth"
# individual prediction with no exposure X for glmerAllBuzzDepth
predFrame0_temp <- expand.grid(X = 0,
                               P = 0,
                               Depth = Depthlevel,
                               ARDepth = 0,
                               Ind = unique(data$Ind))
predBuzz0 <- predict(glmerAllBuzzDepth,
                     newdata = predFrame0_temp)
temp <- cbind(predFrame0_temp, as.data.frame(predBuzz0))
temp$model <- "With Depth"
names(temp) <- c("X", "P", "Depth", "AR", "Ind", "predBuzz0", "model")
predFrame0 <- rbind(predFrame0, temp) # to save

# population prediction with no exposure X for glmerAllBuzz
predFramePop0_temp <- expand.grid(X = 0,
                                  P = 0,
                                  Depth = Depthlevel,
                                  AR = 0,
                                  Ind = "Population")
predBuzzPop0 <- predict(glmerAllBuzz,
                        newdata = predFramePop0_temp,
                        re.form = NA)
predFramePop0 <- cbind(predFramePop0_temp, as.data.frame(predBuzzPop0))
predFramePop0$model <- "Without Depth"
# population prediction with no exposure X for glmerAllBuzzDepth
predFramePop0_temp <- expand.grid(X = 0,
                                  P = 0,
                                  Depth = Depthlevel,
                                  ARDepth = 0,
                                  Ind = "Population")
predBuzzPop0 <- predict(glmerAllBuzzDepth,
                        newdata = predFramePop0_temp,
                        re.form = NA)
temp <- cbind(predFramePop0_temp, as.data.frame(predBuzzPop0))
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

#---------------------------------------------------------------------------------
# Save objects
saveRDS(predFrame, paste0(args[2], "/predFrame.rds"))
saveRDS(predFramePop, paste0(args[2], "/predFramePop.rds"))
saveRDS(predFrame0, paste0(args[2], "/predFrame0.rds"))
saveRDS(predFramePop0, paste0(args[2], "/predFramePop0.rds"))
saveRDS(ChangePop, paste0(args[2], "/ChangePop.rds"))
