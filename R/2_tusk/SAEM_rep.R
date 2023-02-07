#--------------------------------------------------------------------------------
## Objective : repeat SAEM algo several times and produce error metrics
#---------------------------------------------------------------------------------

source("utils/SAEM.R")

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 2) {
  print("arg1 : le chemin vers le dossier de sauvegarde des objets R")
  print("arg2 : le nombre de répétitions de l'expérience")
  stop("Des arguments doivent être donnés au script.", call. = FALSE)
}

#---------------------------------------------------------------------------------
set.seed(15)

# calculation of parameter vectors
R <- args(args[2])
omega.vec <- rep(0, R)
gamma.vec <- rep(0, R)
psi.vec <- rep(0, R)

for (i in 1:R) {
  xi <- rxi(rep(0, n))
  Y <- f(xi) + rnorm(n, 0, omega)
  saem.obj <- saem.alg(Y)
  omega.vec[[r]] <- saem.obj$omega.c
  psi.vec[[r]] <- saem.obj$psi.c
  gamma.vec[[r]] <- saem.obj$gamma.c
}

#---------------------------------------------------------------------------------
# Save R objects
saveRDS(omega.vec, paste0(args[1], "/omega.vec.rds"))
saveRDS(psi.vec, paste0(args[1], "/psi.Vec.rds"))
saveRDS(gamma.vec, paste0(args[1], "/gamma.Vec.rds"))