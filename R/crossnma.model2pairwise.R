crossnma.model2pairwise <- function(x) {
  sm <- x$sm
  ##
  dat <- x$all.data.ad
  trt <- dat$trt
  n <- dat$n
  outcome <- dat$outcome
  study <- dat$study
  se <- dat$se
  ##
  if (sm %in% c("OR", "RR")) {
    suppressWarnings(
      res <-
        pairwise(treat = trt, event = outcome, n = n,
                 studlab = study, sm = sm, warn = FALSE)
    )
  }
  else if (sm %in% c("MD", "SMD")) {
    suppressWarnings(
      res <-
        pairwise(treat = trt,
                 mean = outcome, n = n, sd = se,
                 studlab = study, sm = sm, warn = FALSE)
    )
  }
  ##
  res
}
