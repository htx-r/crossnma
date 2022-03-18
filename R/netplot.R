#' Produce a network plot
#' @description This function creates a plot of evidence network using an object produced by \code{crossnma.model}.
#'
#' @param model A \code{crossnmaModel} object produced by \code{crossnma.model}.
#' @param ... \dots Additional arguments from \code{\link{netgraph.netmeta}}
#'
#' @examples
#' # Two datasets
#' data(prt.data) #  participant-level data
#' data(std.data) # study-level data
#'  #=========================#
#'   # Create a jags model  #
#'  #=========================#
#'  # We conduct a network meta-analysis assuming a random-effects model.
#'  # The data comes from randomized-controlled trials and non-randomized studies (combined naively)
#' mod <- crossnma.model(prt.data=prt.data,
#'                       std.data=std.data,
#'                       trt="trt",
#'                       study="study",
#'                       outcome="outcome",
#'                       n="n",
#'                       design="design",
#'                       reference="A",
#'                       trt.effect="random",
#'                       covariate = NULL,
#'                       method.bias="naive"
#'                       )
#'  #=========================#
#'     # Network plot  #
#'  #=========================#
#' netplot(mod)
#'
#'@export
netplot <- function(model,
                    ...){
  # Bind variables to function
  trt <- NULL
  r <- NULL
  study <- NULL
data.pair = suppressWarnings(pairwise(treat = trt, event = r, n = n,
                          data = model$all.data.ad, studlab = study, sm = "OR",warn = FALSE))
net1 = suppressWarnings(netmeta(data.pair,warn = FALSE))

## ------------------------------------------------------------------------
return(netgraph(net1,...))
}

