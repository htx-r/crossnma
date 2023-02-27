## Calculate treatment effect of each observed comparison using MLE
## estimator
##
## Produce a single value to be used in constructing default priors


max_TE <- function(data, sm) {
  if (sm %in% c("OR", "RR"))
    deltas <- pairwise(treat = data$trt, event = data$outcome, n = data$n,
                       studlab = data$study, sm = sm,
                       incr = 0.5, addincr = TRUE, RR.Cochrane = TRUE)$TE
  ##
  else if (sm %in% c("MD", "SMD")) {
    deltas <- pairwise(treat = data$trt, mean = data$outcome, n = data$n,
                       sd = rep(1, length(data$trt)),
                       studlab = data$study, sm = sm)$TE
  }
  ##
  if (sm == "SMD") {
    ## this is not true, the MD should be divided by s.pooled
    s.pooled <- rep(1, length(deltas))
    deltas <- deltas / s.pooled
  }
  
  ## Return maximum delta for priors
  max(abs(deltas))
}
