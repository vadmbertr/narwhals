# "a" is the pulsation
init.a <- function (Y.arg) {
  ssp <- spectrum(Y.arg, plot = FALSE)
  return(2 * pi * ssp$freq[[which.max(ssp$spec)]])
}

find.x0 <- function (Y.arg) {
  sign.p <- sign(Y.arg[[1]])
  for (i in seq_len(length(Y.arg))) {
    sign.c <- sign(Y.arg[[i]])
    if (sign.c == 0) {
      return(x[[i]])
    } else if (sign.p != sign.c) {
      w.c <- abs(1 / Y.arg[[i]])
      w.p <- abs(1 / Y.arg[[i - 1]])
      return((x[[i]] * w.c + x[[i - 1]] * w.p) / (w.c + w.p))
    }
  }
}

# "b" is the phase shift
# if "b" = 0 then: Y(x) = 0 <=> x = k * (7 * pi / (8 * a)) := x0 with k = 1
# hence (if no noise): b = (7 * pi / 8) - x0 * a
init.b <- function (Y.arg, a.arg) {
  return((7 * pi / 8) - find.x0(Y.arg) * a.arg)
}

saem.alg <- function(Y.obs, n.rep = 500, n.part = 500, smc = FALSE, n.mcmc = 20, n.alpha = 90) {
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
  A.B.max <- max(abs(Y.obs)) + omega.c
  A.c <- runif(1, -A.B.max, A.B.max)
  A <- rep(A.c, n.rep)
  B.c <- runif(1, -A.B.max, A.B.max)
  B <- rep(B.c, n.rep)
  # we initialize "a" using the frequency of "Y" with the highest spectral density. Otherwise "nls" is not converging
  a.est <- init.a(Y.obs)
  a.c <- a.est
  a <- rep(a.c, n.rep)
  # we can also derive an approximate starting value for "b". ?? Is it robust to noise ??
  b.est <- init.b(Y.obs, a.est)
  b.c <- b.est
  b <- rep(b.c, n.rep)

  mcmc.rep <- rep(5, n.rep)
  mcmc.rep[n.mcmc:n.rep] <- 1

  alpha <- rep(1, n.rep)
  alpha[n.alpha:n.rep] <- 1 / ((1:(n.rep - n.alpha + 1))^0.8)

  n.time <- length(Y.obs)
  xi.obj <- mcmc.init(n.time, sum(mcmc.rep))

  for (k in 1:n.rep) {
    # XI
    if (smc) {
      xi.obj <- smc.alg(Y.obs, n.part, omega.c, psi.c, gamma.c, A.c, B.c, a.c, b.c)
    } else {
      xi.obj <- mcmc.alg(Y.obs, mcmc.rep[[k]], xi.obj, omega.c, psi.c, gamma.c, A.c, B.c, a.c, b.c)
    }

    # E
    S1 <- mean((Y.obs - f(x, xi.obj$xi.c, A.arg = A.c, B.arg = B.c, a.arg = a.c, b.arg = b.c))^2)
    S2 <- sum(xi.obj$xi.c[1:(n.time - 1)] * xi.obj$xi.c[2:n.time])
    S3 <- sum(xi.obj$xi.c[1:(n.time - 1)]^2)
    S4 <- sum(xi.obj$xi.c[2:n.time]^2)

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
      return(f(x, xi.obj$xi.c, A.arg = A, B.arg = B, a.arg = a, b.arg = b))
    }
    phi.nls <- nls(Y.obs ~ f.min(A, B, a, b),
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

  return(list(mcmc = xi.obj,
              s1.c = s1.c, s2.c = s2.c, s3.c = s3.c, s4.c = s4.c,
              s1 = s1, s2 = s2, s3 = s3, s4 = s4,
              omega.c = omega.c, psi.c = psi.c, gamma.c = gamma.c, A.c = A.c, B.c = B.c, a.c = a.c, b.c = b.c,
              omega = omega, psi = psi, gamma = gamma, A = A, B = B, a = a, b = b))
}