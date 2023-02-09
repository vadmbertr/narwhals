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
g <- function(xi.arg, a.arg = a) {
  a.arg * x + xi.arg
}
f <- function(xi.arg, A.arg = A, B.arg = B, a.arg = a, b.arg = b) {
  A.arg * sin(g(xi.arg, a.arg) + b.arg) +
    B.arg * sin(2 * (g(xi.arg, a.arg) + b.arg) + pi / 2)
}
rxi <- function(xi.arg = rep(0, n), psi.arg = psi, gamma.arg = gamma) {
  for (i in 2:length(xi.arg)) {
    xi.arg[[i]] <- xi.arg[[i - 1]] * psi.arg + rnorm(1, 0, gamma.arg)
  }
  return(xi.arg)
}