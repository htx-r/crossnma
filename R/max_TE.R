# This produces a single value to be used in constructing default priors
max_TE <- function(data,sm){
  #Bind variables to function
  adj_r.e <- adj_r.c <-theta.e <- theta.c <-
    delta <- time.e <- time.c <- NULL
  # create a table with the treatment comparisons for which we have data
  table <- by.comparison(data, sm = sm)

  # calculate treatment effect of each observed comparison using MLE estimator
  if (sm == "OR"){
    deltas <- table %>%
      mutate(adj_r.e = (outcome.e+0.5)/(n.e + 1), # add 0.5 to ensure ratio is non-zero
             adj_r.c = (outcome.c+0.5)/(n.c + 1)) %>%
      mutate(theta.e = log(adj_r.e/(1-adj_r.e)),
             theta.c = log(adj_r.c/(1-adj_r.c))) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    # make numeric
    deltas <- pull(deltas, delta)

  } else if (sm =="RR"){
    deltas <- table %>%
      mutate(adj_r.e = (outcome.e+0.5)/(n.e + 1),
             adj_r.c = (outcome.c+0.5)/(n.c + 1)) %>%
      mutate(theta.e = log(adj_r.e),
             theta.c = log(adj_r.c)) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    # make numeric
    deltas <- pull(deltas, delta)

  } else if (sm == "MD"){
    deltas <- table %>%
      mutate(delta = as.numeric(outcome.e)-as.numeric(outcome.c)) %>%
      select(delta)
    # make numeric
    deltas <- pull(deltas, delta)

  } else if (sm == "SMD"){ # this is not true, the MD should be divided by s.pooled

    deltas <- table %>%
      mutate(delta = (as.numeric(outcome.e)-as.numeric(outcome.c))/as.numeric(s.pooled)) %>%
      select(delta)
    # make numeric
    deltas <- pull(deltas, delta)
  }

  # return maximum delta for priors
  return(max(abs(deltas)))

}
