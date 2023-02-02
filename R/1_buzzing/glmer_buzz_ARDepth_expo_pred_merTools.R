#--------------------------------------------------------------------------------
## Objective : generate and save obj for plotting prediction band using merTools
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(lme4)
library(merTools)
library(splines)
source("utils/data.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 2) {
  print("Usage du script : RScript glmer_buzz_ARDepth_expo_pred.R arg1 arg2")
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

#---------------------------------------------------------------------------------
# Load models
glmerAllBuzzDepth <- readRDS("../../data/1_buzzing/glmer_buzz_ARDepth_expo/glmerAllBuzzDepth.rds")

#---------------------------------------------------------------------------------
# Obj for ploting
plotdist <- seq(1, 70, 0.1) ## Plotting distances in km
# individual prediction wrt exposure X for glmerAllBuzzDepth
predFrame <- NULL
for (k in unique(data$Ind)) {
  temp <- range(data$X[data$X > 0 & data$Ind == k])
  temp1 <- expand.grid(X = 1 / seq(1 / temp[2], 1 / temp[1], 0.1),
                       ARDepth = 0,
                       Ind = k)
  predFrame <- rbind(predFrame, temp1)
}
predBuzz <- predictInterval(merMod = glmerAllBuzzDepth, newdata = predFrame,
                            which = "fixed", level = 0.95, include.resid.var = FALSE)
predFrame <- cbind(predFrame, as.data.frame(predBuzz)) # to save

# population prediction wrt exposure X for glmerAllBuzzDepth
predFramePop <- expand.grid(X = 1 / plotdist,
                            ARDepth = 0,
                            Ind = "Population")
predBuzzPop <- predictInterval(merMod = glmerAllBuzzDepth, newdata = predFramePop,
                            which = "fixed", level = 0.95, include.resid.var = FALSE)
predFramePop <- cbind(predFramePop, as.data.frame(predBuzzPop)) # to save

# individual prediction with no exposure X for glmerAllBuzzDepth
predFrame0 <- expand.grid(X = 0,
                          ARDepth = 0,
                          Ind = unique(data$Ind))
predBuzz0 <- predictInterval(merMod = glmerAllBuzzDepth, newdata = predFrame0,
                            which = "fixed", level = 0.95, include.resid.var = FALSE)
predFrame0 <- cbind(predFrame0, as.data.frame(predBuzz0)) # to save

# population prediction with no exposure X for glmerAllBuzzDepth
predFramePop0 <- expand.grid(X = 0,
                             ARDepth = 0,
                             Ind = "Population")
predBuzzPop0 <- predictInterval(merMod = glmerAllBuzzDepth, newdata = predFramePop0,
                                which = "fixed", level = 0.95, include.resid.var = FALSE)
predFramePop0 <- cbind(predFramePop0, as.data.frame(predBuzzPop0)) # to save

# percentage of normal behavior
ChangePop <- predFramePop # to save
ChangePop$fit <- exp(ChangePop$fit) / exp(predFramePop0$fit) * 100
## CI
ChangePop$upr <- exp(ChangePop$upr) / exp(predFramePop0$lwr) * 100
ChangePop$lwr <- exp(ChangePop$lwr) / exp(predFramePop0$upr) * 100

#---------------------------------------------------------------------------------
# Save objects
saveRDS(predFrame, paste0(args[2], "/predFrame.rds"))
saveRDS(predFramePop, paste0(args[2], "/predFramePop.rds"))
saveRDS(predFrame0, paste0(args[2], "/predFrame0.rds"))
saveRDS(predFramePop0, paste0(args[2], "/predFramePop0.rds"))
saveRDS(ChangePop, paste0(args[2], "/ChangePop.rds"))
