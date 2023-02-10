apf.init <- function(n.time, n.part) {
  return(list(xi.part = matrix(nrow = n.part, ncol = n.time),
              w.part = matrix(nrow = n.part, ncol = n.time),
              xi.c = rep(0, n.time)))
}

apf.alg <- function (Y.obs, n.part, omega.arg = omega, psi.arg = psi, gamma.arg = gamma,
                     A.arg = A, B.arg = B, a.arg = a, b.arg = b) {
  n.time <- length(Y.obs)
  apf.obj <- apf.init(n.time, n.part)

  apf.obj$xi.part[, 1] <- rep(0, n.part)
  apf.obj$w.part[, 1] <- rep(1 / n.part, n.part)

  for (i in 2:n.time) {
    # 1st stage
    xi.part <- f.xi(apf.obj$xi.part[, i - 1], psi.arg, gamma.arg)
    p1 <- dnorm(Y.obs[[i]], f(i, xi.part, A.arg = A.arg, B.arg = B.arg, a.arg = a.arg, b.arg = b.arg), omega.arg)
    p1 <- sapply(p1, function (p) { max(p, 1e-8) })
    w1 <- p1 * apf.obj$w.part[, i - 1]
    w1 <- w1 / sum(w1)
    # 2nd stage
    idx <- sample(1:n.part, replace = TRUE, prob = w1)
    xi.part <- apf.obj$xi.part[idx, i - 1]
    apf.obj$xi.part[, i] <- f.xi(xi.part, psi.arg, gamma.arg)
    p2 <- dnorm(Y.obs[[i]], f(i, apf.obj$xi.part[, i], A.arg = A.arg, B.arg = B.arg, a.arg = a.arg, b.arg = b.arg), omega.arg)
    p2 <- sapply(p2, function (p) { max(p, 1e-8) })
    w2 <- p2 / p1
    apf.obj$w.part[, i] <- w2 / sum(w2)
  }

  apf.obj$xi.c <- apf.obj$xi.part[sample(1:n.part, 1, prob = apf.obj$w.part[, n.time]), ]
  return(apf.obj)
}