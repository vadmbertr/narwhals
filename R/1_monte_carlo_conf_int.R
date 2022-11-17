percentile_conf_int <- function (dist, level = .05) {
  conf.int <- quantile(dist, c(level/2, 1-level/2))
  names(conf.int) <- c("lower", "upper")
  return(conf.int)
}
