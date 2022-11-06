#--------------------------------------------------------------------------------
## Objective : find the optimal (BIC wise) memory lag
#---------------------------------------------------------------------------------

library(broom.mixed)
library(data.table)
library(splines)
library(lme4)
source("0_data.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read db path from command line
if (length(args)!=5) {
  print("Usage du script : Rscript 1_glmER_buzz_depth_maxlag.R arg1 arg2 arg3 arg4 arg5")
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
## keep the 2018 data
data <- data[data$Year == "2018", ]
## Add time since tagging
data <- AddTime(data)
## We restrict the data to be before exposure
data <- BeforeExposure(data)

#---------------------------------------------------------------------------------
# Estimation of the model
# We fit a memory component to buzzes with lags from 1 to $k$, where we varied $k$ between 1 and maxlag.to,
# adjusted to Depth, entered as a nonlinear spline function.
# For each model, we computed the BIC and chose the model with the lowest BIC.

maxlag.from <- as.numeric(args[2])
maxlag.to <- as.numeric(args[3])
maxlag.n <- as.numeric(args[4])

maxlag.bic.path <- paste0(args[5], "/maxlag.bic.rds")
if (file.exists(maxlag.bic.path)) {
  maxlag.bic <- readRDS(maxlag.bic.path)
  BICvector <- maxlag.bic$BIC
  names(BICvector) <- maxlag.bic$maxlag
} else {
  BICvector <- NULL
}

## Define the first maxlag.n lags
temp <- as.data.frame(shift(data$Buzz, n = 1:maxlag.to, give.names = TRUE))
data <- cbind(temp, data)
LagVariables <- names(data[, 1:maxlag.to])

splineDepth <- ns(data$Depth, knots = c(-323, -158, -54))

## Search
if (maxlag.to == 0) { # allow to bypass search
  maxlag.best <- maxlag.from
} else {
  maxlag.best <- Inf
}
BIC.best.idx <- NULL
### We stop when the interval range indicates a sufficiently small search step
while (maxlag.to - maxlag.from > 2 | is.infinite(maxlag.best)) {
  ### maxlag.n integers between maxlag.from and maxlag.to
  lagvector <- unique(round(seq(from = maxlag.from, to = maxlag.to, length.out = maxlag.n)))
  lagvector <- setdiff(lagvector, as.numeric(names(BICvector)))
  print(lagvector)

  ### fit a glm for every maxlag in lagvector and compute its BIC
  temp <- NULL
  for (maxlag in lagvector) {
    form <- paste("Buzz ~ (1 | Ind) + splineDepth + ",
                  paste(LagVariables[1:maxlag], collapse = " + "))

    glmERAllBuzzDepth <- glmer(form,
                               data = data,
                               family = poisson)
    bic <- BIC(glmERAllBuzzDepth)
    print(paste0(maxlag, ": ", bic))

    #### online save, not really nice
    temp <- c(temp, bic)
    names(temp) <- lagvector[1:length(temp)]
    BICvector <- c(BICvector, temp)
    BICvector <- BICvector[as.character(sort(as.numeric(unique(names(BICvector)))))]
    saveRDS(data.frame(maxlag = as.numeric(names(BICvector)), BIC = as.numeric(BICvector)),
            paste0(args[5], "/maxlag.bic.rds"))
  }

  ### concatenate, ensure uniqueness and sort on the maxlags values
  BICvector <- c(BICvector, temp)
  BICvector <- BICvector[as.character(sort(as.numeric(unique(names(BICvector)))))]

  ### smallest BIC and corresponding maxlag
  BIC.best.idx <- which.min(BICvector)
  maxlag.best <- as.numeric(names(BICvector)[[BIC.best.idx]])

  ### update the search interval
  if (BIC.best.idx == 1) {
    maxlag.from <- as.numeric(names(BICvector)[[BIC.best.idx]])
    maxlag.to <- as.numeric(names(BICvector)[[BIC.best.idx + 1]])
  } else if (BIC.best.idx == length(BICvector)) {
    maxlag.from <- as.numeric(names(BICvector)[[BIC.best.idx - 1]])
    maxlag.to <- as.numeric(names(BICvector)[[BIC.best.idx]])
  } else {
    maxlag.from <- as.numeric(names(BICvector)[[BIC.best.idx - 1]])
    maxlag.to <- as.numeric(names(BICvector)[[BIC.best.idx + 1]])
  }
}

print(paste0("maxlag.best: ", maxlag.best))

#---------------------------------------------------------------------------------
# Build results objects
dataBICDepth <- data.frame(maxlag = as.numeric(names(BICvector)), BIC = as.numeric(BICvector))

## Fit a glm with the best maxlag
form <- paste("Buzz ~ (1 | Ind) + splineDepth + ",
              paste(LagVariables[1:maxlag.best], collapse = " + "))
glmERAllBuzzDepth <- glmer(form,
                           data = data,
                           family = poisson)

glmERAllBuzzDepth.tidy <- tidy(glmERAllBuzzDepth) # to save
glmERAllBuzzDepth.glance <- glance(glmERAllBuzzDepth) # to save

summary(glmERAllBuzzDepth)

#---------------------------------------------------------------------------------
# Save R objects
saveRDS(glmERAllBuzzDepth.tidy, paste0(args[5], "/glmERAllBuzzDepth.tidy.rds"))
saveRDS(glmERAllBuzzDepth.glance, paste0(args[5], "/glmERAllBuzzDepth.glance.rds"))
saveRDS(dataBICDepth, paste0(args[5], "/maxlag.bic.rds"))
