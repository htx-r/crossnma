#' @description This function creates a plot of evidence network using an object produced by \code{crossnma.model}.
#'
#' @param model A \code{crossnmaModel} object produced by \code{crossnma.model}.
#' @param ... \dots Additional arguments from netgraph.netmeta
#'
#' See \code{\link{netgraph.netmeta}} for more details.
#'
#' @return \code{crossnma.run} returns an object of class \code{crossrun} which is a list containing the following components:
#' @return \code{samples}  The MCMC samples produced by running the JAGS model.
#' @return \code{model}  The \code{crossnmaModel} object obtained from \code{crossnma.model} which was used to run \code{jags}.
#' @return \code{trt.key}  A table of the treatment names and their correspondence to integers used in the JAGS model.
#' @examples
#' # An example from participant-level data and study-level data.
#' # data
#' data(prt.data)
#' data(std.data)
#'  #=========================#
#'   # Create a jags model  #
#'  #=========================#
#' mod <- crossnma.model(prt.data=prt.data,
#'                   std.data=std.data,
#'                   trt='trt',
#'                   study='study',
#'                   outcome='outcome',
#'                   n='n',
#'                   design='design',
#'                   reference='A',
#'                   trt.effect='random',
#'                   covariate = NULL,
#'                   method.bias='naive'
#'                    )
#'  #=========================#
#'     # Network plot  #
#'  #=========================#
#' netplot(mod)
#'
#'@export
netplot <- function(model,
                    ...){
data.pair = suppressWarnings(pairwise(treat = trt, event = r, n = n,
                          data = model$all.data.ad, studlab = study, sm = "OR",warn = FALSE))
net1 = suppressWarnings(netmeta(data.pair,warn = FALSE))

## ------------------------------------------------------------------------
return(netgraph(net1,...))
}

