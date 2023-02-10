smc.init <- function(n.time, n.part) {
  return(list(xi.part = matrix(nrow = n.part, ncol = n.time),
              w.part = matrix(nrow = n.part, ncol = n.time),
              xi.c = rep(0, n.time)))
}

smc.alg <- function (Y.obs, n.part, omega.arg = omega, psi.arg = psi, gamma.arg = gamma,
                     A.arg = A, B.arg = B, a.arg = a, b.arg = b) {
  n.time <- length(Y.obs)
  smc.obj <- smc.init(n.time, n.part)

  smc.obj$xi.part[, 1] <- rep(0, n.part)
  smc.obj$w.part[, 1] <- rep(1 / n.part, n.part)

  for (i in 2:n.time) {
    xi.part <- f.xi(smc.obj$xi.part[, i - 1], psi.arg, gamma.arg)  # draw \hat{Xi(t)} cond to Xi(t-1)
    W1 <- dnorm(Y.obs[[i]], f(i, xi.part, A.arg = A.arg, B.arg = B.arg, a.arg = a.arg, b.arg = b.arg), omega.arg)  # prob of Y(t) cond to \hat{Xi(t)}

    # exclude almost zero probalities
    plausible.idx <- which(W1 >= 1e-8)
    xi.part <- xi.part[plausible.idx]
    W1 <- W1[plausible.idx]

    w1 <- W1 / sum(W1)
    idx <- sample(seq_along(plausible.idx), n.part, replace = TRUE, prob = w1)  # particles sampling
    smc.obj$xi.part[, i] <- xi.part[idx]
    smc.obj$w.part[, i] <- w1[idx]
  }

  # trajectory sampling
  smc.obj$xi.c <- smc.obj$xi.part[sample(seq_len(n.part), 1, prob = smc.obj$w.part[, n.time]), ]
  # smc.obj$xi.c <- apply(smc.obj$xi.part * smc.obj$w.part, 2, sum) / apply(smc.obj$w.part, 2, sum)
  return(smc.obj)
}