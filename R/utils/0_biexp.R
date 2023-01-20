BiExp <- function (A1, lrc1, A2, lrc2, maxlag = NULL, lag = NULL) {
  ## maxlag: Number of maximum memory lag. Positive integer
  ## A1, lrc1; A2, lrc2: Coefficients
  if (is.null(lag)) {
    lag <- 1:maxlag
  }
  A1 * exp(-exp(lrc1) * lag) + A2 * exp(-exp(lrc2) * lag)
}