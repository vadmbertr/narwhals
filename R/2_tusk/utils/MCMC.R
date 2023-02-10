acceptance.target <- .23
delta.ar <- .1

mcmc.log.lik <- function(Y.obs, xi.arg, omega.arg, psi.arg, gamma.arg, A.arg, B.arg, a.arg, b.arg) {
  return(sum(dnorm(Y.obs, f(x, xi.arg, A.arg = A.arg, B.arg = B.arg, a.arg = a.arg, b.arg = b.arg), omega.arg, log = TRUE)) +
           sum(dnorm(xi.arg[2:length(xi.arg)], xi.arg[1:(length(xi.arg) - 1)] * psi.arg,
                     gamma.arg, log = TRUE)))
}

mcmc.init <- function(n.time, n.rep) {
  return(list(xi = matrix(nrow = n.rep, ncol = n.time),
              xi.p = matrix(nrow = n.rep, ncol = n.time),
              delta = matrix(nrow = n.rep, ncol = n.time),
              acceptance.rate = matrix(nrow = n.rep, ncol = n.time),
              xi.c = rep(0, n.time),
              delta.c = rep(.05, n.time),
              accepted.n = rep(0, n.time),
              n.it = 0))
}

mcmc.alg <- function(Y.obs, n.rep, mcmc.arg = NULL, omega.arg = omega, psi.arg = psi, gamma.arg = gamma,
                     A.arg = A, B.arg = B, a.arg = a, b.arg = b) {
  if (is.null(mcmc.arg)) {
    mcmc.arg <- mcmc.init(length(Y.obs), n.rep)
  }

  for (j in 1:n.rep) {
    mcmc.arg$n.it <- mcmc.arg$n.it + 1
    for (i in seq_len(length(Y.obs))) {
      xi.p <- mcmc.arg$xi.c
      xi.p[[i]] <- rnorm(1, xi.p[[i]], mcmc.arg$delta.c[[i]])
      mcmc.arg$xi.p[mcmc.arg$n.it, i] <- xi.p[[i]]
      log.prob <- min(1, mcmc.log.lik(Y.obs, xi.p, omega.arg, psi.arg, gamma.arg, A.arg, B.arg, a.arg, b.arg) -
        mcmc.log.lik(Y.obs, mcmc.arg$xi.c, omega.arg, psi.arg, gamma.arg, A.arg, B.arg, a.arg, b.arg))
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
    mcmc.arg$xi[mcmc.arg$n.it,] <- mcmc.arg$xi.c
    mcmc.arg$delta[mcmc.arg$n.it,] <- mcmc.arg$delta.c
  }

  return(mcmc.arg)
}

pmcmc.init <- function(n.time, n.rep) {
  return(list(xi = matrix(nrow = n.rep, ncol = n.time),
              xi.p = matrix(nrow = n.rep, ncol = n.time),
              acceptance.rate = matrix(nrow = n.rep, ncol = n.time),
              xi.c = rep(0, n.time),
              accepted.n = rep(0, n.time),
              n.it = 0))
}

pmcmc.alg <- function(Y.obs, n.rep, n.part = 100, mcmc.arg = NULL, omega.arg = omega, psi.arg = psi, gamma.arg = gamma,
                     A.arg = A, B.arg = B, a.arg = a, b.arg = b) {
  if (is.null(mcmc.arg)) {
    mcmc.arg <- pmcmc.init(length(Y.obs), n.rep)
  }

  for (j in 1:n.rep) {
    mcmc.arg$n.it <- mcmc.arg$n.it + 1
    apf.obj <- apf.alg(Y.obs, n.part)
    for (i in seq_len(length(Y.obs))) {
      xi.p <- mcmc.arg$xi.c
      xi.p[[i]] <- apf.obj$xi.c[[i]]
      mcmc.arg$xi.p[mcmc.arg$n.it, i] <- xi.p[[i]]
      log.prob <- min(1, mcmc.log.lik(Y.obs, xi.p, omega.arg, psi.arg, gamma.arg, A.arg, B.arg, a.arg, b.arg) -
        mcmc.log.lik(Y.obs, mcmc.arg$xi.c, omega.arg, psi.arg, gamma.arg, A.arg, B.arg, a.arg, b.arg))
      if (log(runif(1)) < log.prob) {
        mcmc.arg$xi.c <- xi.p
        mcmc.arg$accepted.n[[i]] <- mcmc.arg$accepted.n[[i]] + 1
      }

      mcmc.arg$acceptance.rate[mcmc.arg$n.it, i] <- mcmc.arg$accepted.n[[i]] / mcmc.arg$n.it
    }
    mcmc.arg$xi[mcmc.arg$n.it,] <- mcmc.arg$xi.c
  }

  return(mcmc.arg)
}