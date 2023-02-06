

#---------------------------------------------------------------------------------
# Read script arguments
args <- commandArgs(trailingOnly = TRUE) # read args from command line
if (length(args) != 2) {
  print("arg1 : le chemin vers le dossier de sauvegarde des objets R")
  print("arg2 : le nombre d'itération pour le l'erreur")
  stop("Des arguments doivent Ãªtre donnÃ©s au script.", call. = FALSE)
}

#---------------------------------------------------------------------------------
# parameters
A <- 1 / 2
B <- -1 / 4
b <- 1
a <- 0.1
beta <- 0.05
sigma <- 0.1
omega <- 0.01
delta <- 1
psi <- exp(-delta * beta)
gamma <- sigma / sqrt(2 * beta) * sqrt(1 - psi^2)

n <- 500
x <- 1:n

#---------------------------------------------------------------------------------
# Functions
g <- function (xi, a.arg = a) {
  a.arg * x + xi
}
f <- function (xi, A.arg = A, a.arg = a, B.arg = B, b.arg = b) {
  A.arg * sin(g(xi, a.arg) + b.arg) +
    B.arg * sin(2 * g(xi, a.arg) + 2 * b.arg + pi / 2)
}
rxi <- function (xi=rep(0,n), psi.arg = psi, gamma.arg = gamma) {
  for (i in 2:length(xi)) {
    xi[[i]] <- xi[[i - 1]] * psi.arg + rnorm(1, 0, gamma.arg)
  }
  return(xi)
}

xi <- rxi(rep(0, n))
Y <- f(xi) + rnorm(n, 0, omega)


acceptance.target <- .23
delta.ar <- 0.1

mcmc.log.lik <- function (xi.arg, omega.arg, psi.arg, gamma.arg) {
  return(sum(dnorm(Y, f(xi.arg), omega.arg, log = TRUE)) +
           sum(dnorm(xi.arg[2:length(xi.arg)], xi.arg[1:(length(xi.arg) - 1)] * psi.arg,
                     gamma.arg, log = TRUE)))
}

mcmc.init <- function (n.rep) {
  return(list(xi = matrix(nrow = n.rep, ncol = n),
              xi.p = matrix(nrow = n.rep, ncol = n),
              delta = matrix(nrow = n.rep, ncol = n),
              acceptance.rate = matrix(nrow = n.rep, ncol = n),
              xi.c = rep(0, n),
              delta.c = rep(.05, n),
              accepted.n = rep(0, n),
              n.it = 0))
}

mcmc.alg <- function (n.rep, mcmc.arg = NULL, omega.arg = omega, psi.arg = psi, gamma.arg = gamma) {
  if (is.null(mcmc.arg)) {
    mcmc.arg <- mcmc.init(n.rep)
  }
  
  for(j in 1:n.rep) {
    mcmc.arg$n.it <- mcmc.arg$n.it + 1
    for (i in 1:n) {
      xi.p <- mcmc.arg$xi.c
      xi.p[[i]] <- rnorm(1, xi.p[[i]], mcmc.arg$delta.c[[i]])
      mcmc.arg$xi.p[mcmc.arg$n.it, i] <- xi.p[[i]]
      log.prob <- min(1, mcmc.log.lik(xi.p, omega.arg, psi.arg, gamma.arg) -
                        mcmc.log.lik(mcmc.arg$xi.c, omega.arg, psi.arg, gamma.arg))
      if (log(runif(1)) < log.prob) {
        mcmc.arg$xi.c <- xi.p
        mcmc.arg$accepted.n[[i]] <- mcmc.arg$accepted.n[[i]] + 1
      }
      
      mcmc.arg$acceptance.rate[mcmc.arg$n.it, i] <- mcmc.arg$accepted.n[[i]] / mcmc.arg$n.it
      if (mcmc.arg$acceptance.rate[mcmc.arg$n.it, i] < acceptance.target * (1 - .1)) {
        mcmc.arg$delta.c[[i]] <- mcmc.arg$delta.c[[i]] * (1 - delta.ar)
      } else if (mcmc.arg$acceptance.rate[mcmc.arg$n.it, i] > acceptance.target * (1 + .1)) {
        mcmc.arg$delta.c[[i]] <- mcmc.arg$delta.c[[i]] * (1 + delta.ar)
      }
    }
    mcmc.arg$xi[mcmc.arg$n.it, ] <- mcmc.arg$xi.c
    mcmc.arg$delta[mcmc.arg$n.it, ] <- mcmc.arg$delta.c
  }
  
  return(mcmc.arg)
}



saem.alg <- function (n.rep = 500, n.mcmc = 20, n.alpha = 90) {
  s1.c <- 0
  s1 <- rep(s1.c, n.rep)
  s2.c <- 0
  s2 <- rep(s2.c, n.rep)
  s3.c <- 0
  s3 <- rep(s3.c, n.rep)
  s4.c <- 0
  s4 <- rep(s4.c, n.rep)
  
  omega.c <- .5
  omega <- rep(omega.c, n.rep)
  psi.c <- .5
  psi <- rep(psi.c, n.rep)
  gamma.c <- .5
  gamma <- rep(gamma.c, n.rep)
  
  mcmc.rep <- rep(5, n.rep)
  mcmc.rep[n.mcmc:n.rep] <- 1
  
  alpha <- rep(1, n.rep)
  alpha[n.alpha:n.rep] <- 1 / ((1:(n.rep - n.alpha + 1))^0.8)
  
  mcmc.obj <- mcmc.init(sum(mcmc.rep))
  
  for (k in 1:n.rep) {
    # MCMC
    mcmc.obj <- mcmc.alg(mcmc.rep[[k]], mcmc.obj, omega.c, psi.c, gamma.c)
    
    # E
    S1 <- mean((Y - f(mcmc.obj$xi.c))^2)  # TODO: add A, B, a, b
    S2 <- sum(mcmc.obj$xi.c[1:(n - 1)] * mcmc.obj$xi.c[2:n])
    S3 <- sum(mcmc.obj$xi.c[1:(n - 1)]^2)
    S4 <- sum(mcmc.obj$xi.c[2:n]^2)
    
    # SA
    s1.c <- s1.c + alpha[[k]] * (S1 - s1.c)
    s2.c <- s2.c + alpha[[k]] * (S2 - s2.c)
    s3.c <- s3.c + alpha[[k]] * (S3 - s3.c)
    s4.c <- s4.c + alpha[[k]] * (S4 - s4.c)
    s1[[k]] <- s1.c
    s2[[k]] <- s2.c
    s3[[k]] <- s3.c
    s4[[k]] <- s4.c
    
    # M
    omega.c <- sqrt(s1.c)
    psi.c <- s2.c / s3.c
    gamma.c <- sqrt((psi.c^2 * s3.c - 2 * psi.c * s2.c + s4.c) / n)
    omega[[k]] <- omega.c
    psi[[k]] <- psi.c
    gamma[[k]] <- gamma.c
  }
  
  return(list(mcmc = mcmc.obj,
              s1.c = s1.c, s2.c = s2.c, s3.c = s3.c, s4.c = s4.c,
              s1 = s1, s2 = s2, s3 = s3, s4 = s4,
              omega.c = omega.c, psi.c = psi.c, gamma.c = gamma.c,
              omega = omega, psi = psi, gamma = gamma))
}

#---------------------------------------------------------------------------------
# calculation of parameter vectors
r<-args(args[2])
omega.vec<-rep(0,r)
gamma.vec<-rep(0,r)
psi.vec<-rep(0,r)


for (w in 1:r){
  saem.obj <- saem.alg()
  omega.vec[w]<-saem.obj$omega.c
  psi.vec[w]<-saem.obj$psi.c
  gamma.vec[w]<-saem.obj$gamma.c
}

#---------------------------------------------------------------------------------
# Save R objects
saveRDS(omega.vec, paste0(args[1], "/omega.vec.rds"))
saveRDS(psi.vec, paste0(args[1], "/psi.Vec.rds"))
saveRDS(gamma.vec, paste0(args[1], "/gamma.Vec.rds"))


