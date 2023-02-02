#---------------------------------------------------------------------------------
## Objective : generate and save obj for plotting prediction band using delta meth
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(lme4)
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
glmerAllBuzzDepth <- readRDS(paste0(args[2], "/../glmer_buzz_ARDepth_expo/glmerAllBuzzDepth.rds"))

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
predBuzz <- predict(glmerAllBuzzDepth,
                    newdata = predFrame)
predFrame <- cbind(predFrame, as.data.frame(predBuzz)) # to save

# population prediction wrt exposure X for glmerAllBuzzDepth
predFramePop <- expand.grid(X = 1 / plotdist,
                            ARDepth = 0,
                            Ind = "Population")
predBuzzPop <- predict(glmerAllBuzzDepth,
                       newdata = predFramePop,
                       re.form = NA)
predFramePop <- cbind(predFramePop, as.data.frame(predBuzzPop)) # to save

# individual prediction with no exposure X for glmerAllBuzzDepth
predFrame0 <- expand.grid(X = 0,
                          ARDepth = 0,
                          Ind = unique(data$Ind))
predBuzz0 <- predict(glmerAllBuzzDepth,
                     newdata = predFrame0)
predFrame0 <- cbind(predFrame0, as.data.frame(predBuzz0)) # to save

# population prediction with no exposure X for glmerAllBuzzDepth
predFramePop0 <- expand.grid(X = 0,
                             ARDepth = 0,
                             Ind = "Population")
predBuzzPop0 <- predict(glmerAllBuzzDepth,
                        newdata = predFramePop0,
                        re.form = NA)
predFramePop0 <- cbind(predFramePop0, as.data.frame(predBuzzPop0)) # to save

# percentage of normal behavior
ChangePop <- predFramePop # to save
ChangePop$change <- exp(ChangePop$predBuzzPop)
ChangePop$change <- ChangePop$change / exp(predFramePop0$predBuzzPop0) * 100
## CI
alpha <- .05
### mean, var estimates
expo.coef <- readRDS("../../data/1_buzzing/glmer_buzz_ARDepth_expo_par/expo.coef.mvnorm.mc.rds")
#### no intercept as we are looking at the percentage of normal behaviour
# expo.coef <- expo.coef[!expo.coef$term == "(Intercept)",]
expo.coef.estimate <- expo.coef[, c("term", "estimate")]
expo.coef.estimate$seq <- with(expo.coef.estimate, ave(estimate, term, FUN = seq_along))
# expo.coef.estimate <- dcast(expo.coef.estimate, seq ~ term, value.var = "estimate")[, 2:4]
expo.coef.estimate <- dcast(expo.coef.estimate, seq ~ term, value.var = "estimate")[, 1:4]
Beta.hat <- apply(expo.coef.estimate, 2, mean)
Sigma.hat <- cov(expo.coef.estimate)
### f.hat
expo <- 1 / plotdist
X <- ns(expo, knots = quantile(data$X[data$X > 0], 1:2 / 3))
X <- cbind(rep(1, nrow(X)), X)
# f2.hat <- (exp(X %*% Beta.hat) * 100)^2
f2.hat <- (exp(X %*% Beta.hat - predFramePop0$predBuzzPop0) * 100)^2
### var
sigma2.hat <- f2.hat * diag(X %*% Sigma.hat %*% t(X))
### CI
ChangePop$CI <- qnorm(1 - alpha / 2) * sqrt(sigma2.hat / length(Beta.hat))

#---------------------------------------------------------------------------------
# Save objects
saveRDS(predFrame, paste0(args[2], "/predFrame.rds"))
saveRDS(predFramePop, paste0(args[2], "/predFramePop.rds"))
saveRDS(predFrame0, paste0(args[2], "/predFrame0.rds"))
saveRDS(predFramePop0, paste0(args[2], "/predFramePop0.rds"))
saveRDS(ChangePop, paste0(args[2], "/ChangePop.rds"))
