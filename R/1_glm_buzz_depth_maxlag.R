#--------------------------------------------------------------------------------
## Objective : find the optimal (BIC wise) memory lag
#---------------------------------------------------------------------------------

library(data.table)
library(splines)
source("0_data.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read db path from command line
if (length(args)!=4) {
  print("Usage du script : Rscript 1_glm_buzz_depth_maxlag.R arg1 arg2 arg3 arg4")
  print("arg1 : le chemin vers la base de données")
  print("arg2 : la borne inférieure de recherche du lag maximum")
  print("arg3 : la borne supérieure de recherche du lag maximum")
  print("arg4 : le nombre de lag à considérer entre les deux bornes")
  stop("Des arguments doivent être donnés au script.", call. = FALSE)
}

#---------------------------------------------------------------------------------
# Data preparation
dataAll <- fread(args[1])
data <- dataAll
## keep the 2018 data TODO: why?
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

## Define the first maxlag.n lags
temp <- as.data.frame(shift(data$Buzz, n = maxlag.from:maxlag.to, give.names = TRUE))
data <- cbind(temp, data)
LagVariables <- names(data[, 1:maxlag.to])

splineDepth <- ns(data$Depth, knots = c(-323, -158, -54))

## Search
BICvector <- NULL
maxlag.best <- Inf
BIC.best.idx <- NULL
### We stop when the interval range indicates a sufficiently small search step
while (maxlag.to - maxlag.from > 2 | is.infinite(maxlag.best)) {
  ### maxlag.n integers between maxlag.from and maxlag.to
  lagvector <- unique(round(seq(from = maxlag.from, to = maxlag.to, length.out = maxlag.n)))
  lagvector <- setdiff(lagvector, as.numeric(names(BICvector)))
  print(lagvector)

  ### fit a glm for every maxlag in lagvector and compute its BIC
  temp <- sapply(lagvector, function (maxlag) {
    form <- paste("Buzz ~ Ind + splineDepth + ",
                  paste(LagVariables[1:maxlag], collapse = " + "))

    # TODO: 2- ajout effet aléatoire sur les individus
    glmAllBuzzDepth <- glm(form,
                           data = data,
                           family = poisson)
    BIC(glmAllBuzzDepth)
  })
  names(temp) <- lagvector

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

#---------------------------------------------------------------------------------
# Build results objects
dataBICDepth <- data.frame(maxlag = as.numeric(names(BICvector)), BIC = as.numeric(BICvector))

## Fit a glm with the best maxlag
form <- paste("Buzz ~ Ind + splineDepth + ",
              paste(LagVariables[1:maxlag.best], collapse = " + "))
# TODO: 2- ajout effet aléatoire sur les individus
glmAllBuzzDepth <- glm(form,
                       data = data,
                       family = poisson)
coeftable <- coefficients(summary(glmAllBuzzDepth))

#---------------------------------------------------------------------------------
# Save R objects
saveRDS(dataBICDepth, "../data/glm_buzz_depth_maxlag_bic.rds")
saveRDS(coeftable, "../data/glm_buzz_depth_bestmaxlag_coefs.rds")
