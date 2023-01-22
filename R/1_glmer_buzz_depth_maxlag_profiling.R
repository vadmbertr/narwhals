#--------------------------------------------------------------------------------
## Objective : keep the time and memory information for GLMER with 1:90 lag (by=10)
#---------------------------------------------------------------------------------

library(bench)
library(data.table)
library(splines)
library(lme4)
source("utils/0_data.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read db path from command line
if (length(args) != 2) {
  print("Usage du script : Rscript 1_glmer_buzz_depth_memory_time.R arg1 arg2")
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

#-------------------------------------------------------------------------------------
# Estimation of the model
# We fit a memory component to buzzes with lags from 1 to $k$, where we varied $k$ between 1 and maxlag.to,
# adjusted to Depth, entered as a nonlinear spline function.
# For each model, we computed the BIC and chose the model with the lowest BIC.

maxlag.to <- 90

## Define the first maxlag.n lags
temp <- as.data.frame(shift(data$Buzz, n = 1:maxlag.to, give.names = TRUE))
data <- cbind(temp, data)
LagVariables <- names(data[, 1:maxlag.to])

splineDepth <- ns(data$Depth, knots = c(-323, -158, -54))

lagvector <- c(1, seq(10, 90, by = 10))

# retrieve previous runs
maxlag.memory.time.path <- paste0(args[2], "/profiling.rds")
if (file.exists(maxlag.memory.time.path)) {
  glmer.maxlag.memory.time <- readRDS(maxlag.memory.time.path)
  lagvector <- setdiff(lagvector,
                       as.numeric(rownames(glmer.maxlag.memory.time)))
} else {
  glmer.maxlag.memory.time <- NULL
}

print(lagvector)

for (maxlag in lagvector) {
  print(maxlag)
  form <- paste("Buzz ~ (1 | Ind) + splineDepth + ",
                  paste(LagVariables[1:maxlag], collapse = " + "))
  bench_res <- mark(glmer = glmer(form, data = data, family = poisson),
                    iterations = 1, time_unit = "s")
  temp <- as.matrix(bench_res[, c("total_time", "mem_alloc")])
  rownames(temp) <- maxlag
  glmer.maxlag.memory.time <- rbind(glmer.maxlag.memory.time, temp)
  # online save
  saveRDS(glmer.maxlag.memory.time, maxlag.memory.time.path)
}
