percentile_conf_int <- function (dist, level = .05) {
  conf.int <- quantile(as.numeric(dist), c(level/2, 1-level/2), na.rm = TRUE)
  names(conf.int) <- c("lower", "upper")
  return(conf.int)
}
