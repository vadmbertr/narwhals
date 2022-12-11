BiExp.nls <- function (y, x, w) {
  biexp.reg <- nls(y ~ SSbiexp(x, A1, lrc1, A2, lrc2), weights = w)
  coefs <- summary(biexp.reg)$coefficients[1:4, 1]
  names(coefs) <- c("A1", "lrc1", "A2", "lrc2")
  return(coefs)
}

BiExp <- function (A1, lrc1, A2, lrc2, maxlag = NULL, lag = NULL) {
  ## maxlag: Number of maximum memory lag. Positive integer
  ## A1, lrc1; A2, lrc2: Coefficients
  if (is.null(lag)) {
    lag <- 1:maxlag
  }
  A1 * exp(-exp(lrc1) * lag) + A2 * exp(-exp(lrc2) * lag)
}