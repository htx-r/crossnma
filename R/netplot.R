#' @description This function creates a plot of evidence network using an object produced by \code{crossnma.model}.
#'
#' @param model A \code{crossnmaModel} object produced by \code{crossnma.model}.
#' @param ... \dots Additional arguments from netgraph.netmeta
#'
#' See \code{\link{netgraph.netmeta}} for more details.
#'
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

