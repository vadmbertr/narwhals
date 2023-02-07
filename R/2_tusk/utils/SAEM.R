source("MCMC.R")

saem.alg <- function(Y.obs, n.rep = 500, n.mcmc = 20, n.alpha = 90) {
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

  n.time <- length(Y.obs)
  mcmc.obj <- mcmc.init(n.time, sum(mcmc.rep))

  for (k in 1:n.rep) {
    # MCMC
    mcmc.obj <- mcmc.alg(Y.obs, mcmc.rep[[k]], mcmc.obj, omega.c, psi.c, gamma.c)

    # E
    S1 <- mean((Y.obs - f(mcmc.obj$xi.c))^2)  # TODO: add A, B, a, b
    S2 <- sum(mcmc.obj$xi.c[1:(n.time - 1)] * mcmc.obj$xi.c[2:n.time])
    S3 <- sum(mcmc.obj$xi.c[1:(n.time - 1)]^2)
    S4 <- sum(mcmc.obj$xi.c[2:n.time]^2)

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
    gamma.c <- sqrt((psi.c^2 * s3.c - 2 * psi.c * s2.c + s4.c) / n.time)
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