#-------------------------------------------------------------------------------------------------
## Objective : MC approach to estimate AR bi-exponential regression coefficients mean and var/cov
#--------------------------------------------------------------------------------------------------

library(broom)
library(mvtnorm)
library(parallel)
library(RhpcBLASctl)
source("0_biexp.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 3) {
  print("Usage du script : Rscript 1_biexp_AR_mc.R arg1 arg2 arg3")
  print("arg1 : le chemin vers le dossier de sauvegarde des objets R")
  print("arg2 : nombre de coefficients estimés souhaités")
  print("arg3 : le nombre de coeurs alloués") # cannot be retrieve from R
  stop("Des arguments doivent être donnés au script.", call. = FALSE)
}

#---------------------------------------------------------------------------------
# MC

maxlag.bic <- readRDS("../data/glmER_buzz_depth_maxlag/maxlag.bic.rds")
glmer.coefs <- readRDS("../data/glmer_buzz_ARDepth/glmERBuzzARDepth.tidy.rds")
glmer.coefs.vcov <- as.matrix(readRDS("../data/glmer_buzz_ARDepth/glmERBuzzARDepth.vcov.rds"))

maxlag.opt <- as.integer(maxlag.bic[which.min(maxlag.bic[, 2]), 1])
AR.est <- glmer.coefs$estimate[1:(maxlag.opt)+5]
AR.vcov <- glmer.coefs.vcov[1:(maxlag.opt)+5, 1:(maxlag.opt)+5]
AR.se <- diag(AR.vcov)

fit.nls <- function (i) {
  ARcoefs <- as.numeric(rmvnorm(1, mean = AR.est, sigma = AR.vcov,
                                checkSymmetry = FALSE))
  biexp.reg <- nls(ARcoefs ~ SSbiexp(1:maxlag.opt, A1, lrc1, A2, lrc2), weights = 1/AR.se)
  return(tidy(biexp.reg)[, c("term", "estimate", "std.error")])
}

n.done <- function (df, nc) {
  if (is.na(nc)) return(0)
  else return(round(nrow(df) / nc))
}

biexp.coef.path <- paste0(args[2], "/biexp.AR.est.rds")
if (file.exists(biexp.coef.path)) {
  biexp.coef.all <- readRDS(biexp.coef.path)
  n.coefs <- length(unique(biexp.coef.all$term))
} else {
  biexp.coef.all <- data.frame()
  n.coefs <- NA
}

## Parallelism
blas_set_num_threads(1) # otherwise some kind of deadlock can occured
n.jobs <- as.numeric(args[3]) - 1

while (as.numeric(args[2]) - n.done(biexp.coef.all, n.coefs) > 0) {
  print(n.done(biexp.coef.all, n.coefs))
  n.jobs <- min(n.jobs, as.numeric(args[2]) - n.done(biexp.coef.all, n.coefs))
  biexp.coef <- do.call(rbind, mclapply(1:n.jobs, fit.nls, mc.cores = n.jobs))
  if (is.na(n.coefs)) {
    n.coefs <- as.integer(nrow(biexp.coef) / n.jobs)
  }
  biexp.coef.all <- rbind(expo.coef.all, expo.coef)
  saveRDS(biexp.coef.all, expo.coef.path)
}
