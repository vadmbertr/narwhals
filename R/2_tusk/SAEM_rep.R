#--------------------------------------------------------------------------------
## Objective : repeat SAEM algo several times and produce error metrics
#---------------------------------------------------------------------------------

library(parallel)
library(RhpcBLASctl)
source("utils/simulation.R")
source("utils/MCMC.R")
source("utils/SAEM.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 3) {
  print("arg1 : le nombre de répétitions de l'expérience")
  print("arg2 : le nombre de coeurs alloués")
  print("arg3 : le chemin vers le dossier de sauvegarde des objets R")
  stop("Des arguments doivent être donnés au script.", call. = FALSE)
}

#---------------------------------------------------------------------------------
set.seed(15)

run.saem <- function (i) {
  xi <- rxi()
  Y <- f(x, xi) + rnorm(n, 0, omega)
  saem.obj <- saem.alg(Y)
  return(data.frame(omega = saem.obj$omega.c, psi = saem.obj$psi.c, gamma = saem.obj$gamma.c,
                    A = saem.obj$A.c, B = saem.obj$B.c, a = saem.obj$a.c, b = saem.obj$b.c))
}

n.rep <- as.numeric(args[[1]])
blas_set_num_threads(1)
n.jobs <- as.numeric(args[[2]])
parameters.est <- do.call(rbind, mclapply(1:n.rep, run.saem, mc.cores = n.jobs))
saveRDS(parameters.est, paste0(args[[3]], "/parameters.est.smc.rds"))