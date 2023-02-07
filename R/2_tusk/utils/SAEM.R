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
  A.B.max <- max(Y) + omega.c
  A.c <- runif(1, -A.B.max, A.B.max)
  A <- rep(A.c, n.rep)
  B.c <- runif(1, -A.B.max, A.B.max)
  B <- rep(B.c, n.rep)
  b.max <- 2
  b.c <- runif(1, 0, b.max)
  b <- rep(b.c, n.rep)
  # we initialize "a" using the frequency of "Y" with the highest spectral density. Otherwise "nls" is not converging
  ssp <- spectrum(Y, plot = FALSE)
  a.c <- 2 * pi * ssp$freq[[which.max(ssp$spec)]]
  a <- rep(a.c, n.rep)

  mcmc.rep <- rep(5, n.rep)
  mcmc.rep[n.mcmc:n.rep] <- 1

  alpha <- rep(1, n.rep)
  alpha[n.alpha:n.rep] <- 1 / ((1:(n.rep - n.alpha + 1))^0.8)

  n.time <- length(Y.obs)
  mcmc.obj <- mcmc.init(n.time, sum(mcmc.rep))

  for (k in 1:n.rep) {
    # MCMC
    mcmc.obj <- mcmc.alg(Y.obs, mcmc.rep[[k]], mcmc.obj, omega.c, psi.c, gamma.c, A.c, B.c, a.c, b.c)

    # E
    S1 <- mean((Y.obs - f(mcmc.obj$xi.c, A.arg = A.c, B.arg = B.c, a.arg = a.c, b.arg = b.c))^2)
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
    # function to minimize
    f.min <- function(A, B, a, b) {
      return(f(mcmc.obj$xi.c, A.arg = A, B.arg = B, a.arg = a, b.arg = b))
    }
    phi.nls <- nls(Y ~ f.min(A, B, a, b),
                   start = list(A = A.c, B = B.c, a = a.c, b = b.c),
                   control = nls.control(warnOnly = TRUE))
    A.c <- coef(phi.nls)[[1]]
    B.c <- coef(phi.nls)[[2]]
    a.c <- coef(phi.nls)[[3]]
    b.c <- coef(phi.nls)[[4]]

    omega[[k]] <- omega.c
    psi[[k]] <- psi.c
    gamma[[k]] <- gamma.c
    A[[k]] <- A.c
    B[[k]] <- B.c
    a[[k]] <- a.c
    b[[k]] <- b.c
  }

  return(list(mcmc = mcmc.obj,
              s1.c = s1.c, s2.c = s2.c, s3.c = s3.c, s4.c = s4.c,
              s1 = s1, s2 = s2, s3 = s3, s4 = s4,
              omega.c = omega.c, psi.c = psi.c, gamma.c = gamma.c, A.c = A.c, B.c = B.c, a.c = a.c, b.c = b.c,
              omega = omega, psi = psi, gamma = gamma, A = A, B = B, a = a, b = b))
}