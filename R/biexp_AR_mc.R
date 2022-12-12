#-------------------------------------------------------------------------------------------------
## Objective : MC approach to estimate AR bi-exponential regression coefficients mean and var/cov
#--------------------------------------------------------------------------------------------------

library(broom)
library(mvtnorm)
library(reshape2)
source("0_biexp.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 2) {
  print("Usage du script : Rscript 1_biexp_AR_mc.R arg1 arg2")
  print("arg1 : le chemin vers le dossier de sauvegarde des objets R")
  print("arg2 : nombre de coefficients estimés souhaités")
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
AR.se <- glmer.coefs$std.error[1:(maxlag.opt)+5]

biexp.coef <- replicate(args[2], {
  ARcoefs <- as.numeric(rmvnorm(1, mean = AR.est, sigma = AR.vcov,
                                checkSymmetry = FALSE))
  biexp.reg <- tryCatch({
    nls(ARcoefs ~ SSbiexp(1:maxlag.opt, A1, lrc1, A2, lrc2),
                   weights = 1/AR.se)
  }, error = function (c) return(NULL))

  if (is.null(biexp.reg)) return(NULL)
  else return(as.matrix(tidy(biexp.reg)[, c("term", "estimate", "std.error")]))
}, simplify = FALSE)
biexp.coef <- do.call(rbind, biexp.coef)

biexp.coef.estimate <- as.data.frame(biexp.coef)
biexp.coef.estimate$seq <- with(biexp.coef.estimate,
                                ave(std.error, term, FUN = seq_along))
biexp.coef.estimate <- apply(dcast(biexp.coef.estimate, seq ~ term, value.var = "estimate")[2:5],
                             2, as.numeric)
biexp.coef.cov <- cov(biexp.coef.estimate)

saveRDS(biexp.coef, paste0(args[1], "/biexp.coef.rds"))
saveRDS(biexp.coef.estimate, paste0(args[1], "/biexp.coef.estimate.rds"))
saveRDS(biexp.coef.cov, paste0(args[1], "/biexp.coef.cov.rds"))
