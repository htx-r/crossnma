#' Produce a network plot
#' @description This function creates a plot of evidence network using an object produced by \code{crossnma.model}.
#'
#' @param x A \code{crossnmaModel} object produced by \code{crossnma.model}.
#' @param ... \dots Additional arguments from \code{\link{netgraph.netmeta}}
#'
#' @examples
#' # We conduct a network meta-analysis assuming a random-effects
#' # model.
#' # The data comes from randomized-controlled trials and
#' # non-randomized studies (combined naively)
#' head(ipddata) # participant-level data
#' head(stddata) # study-level data
#' 
#' #=========================#
#' # Create a jags model     #
#' #=========================#
#' mod <- crossnma.model(treat, id, relapse, n, design,
#'   prt.data = ipddata, std.data = stddata,
#'   reference = "A", trt.effect = "random", method.bias = "naive")
#' 
#' #=========================#
#' # Fit jags model          #
#' #=========================#
#' fit <-
#'   crossnma.run(model = mod, n.adapt = 20,
#'     n.iter = 50, thin = 1, n.chains = 3)
#'
#' #=========================#
#' # Network plot            #
#' #=========================#
#' netgraph(mod)
#'
#'@method netgraph crossnmaModel
#'@export

netgraph.crossnmaModel <- function(x, ...) {
  ## Bind variables to function
  trt <- NULL
  r <- NULL
  study <- NULL
  data.pair <-
    suppressWarnings(pairwise(treat = trt, event = r, n = n,
                              data = x$all.data.ad, studlab = study,
                              sm = "OR", warn = FALSE))
  net1 <- suppressWarnings(netmeta(data.pair,
                                   warn = FALSE))
  return(netgraph(net1,...))
}
