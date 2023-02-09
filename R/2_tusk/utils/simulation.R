# Parameters
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

# Functions
g <- function(x.arg, xi.arg, a.arg = a) {
  a.arg * x.arg + xi.arg
}
f <- function(x.arg, xi.arg, A.arg = A, B.arg = B, a.arg = a, b.arg = b) {
  A.arg * sin(g(x.arg, xi.arg, a.arg) + b.arg) +
    B.arg * sin(2 * (g(x.arg, xi.arg, a.arg) + b.arg) + pi / 2)
}
f.xi <- function (xi.p, psi.arg = psi, gamma.arg = gamma) {
  return(xi.p * psi.arg + rnorm(length(xi.p), 0, gamma.arg))
}
rxi <- function(xi.arg = rep(0, n), psi.arg = psi, gamma.arg = gamma) {
  for (i in 2:length(xi.arg)) {
    xi.arg[[i]] <- f.xi(xi.arg[[i - 1]], psi.arg, gamma.arg)
  }
  return(xi.arg)
}