#' Print call used to create JAGS model for cross-design & -format
#' network meta-analysis or regression
#' 
#' @description
#' Print call used to create JAGS model for cross-design & -format
#' network meta-analysis or regression
#' 
#' @param x An object of class \code{crossnma.model}.
#' @param \dots Additional arguments (ignored).
#'
#' @return
#' No return value (print function).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{crossnma.model}}
#' 
#' @keywords print
#' 
#' @examples
#' # We conduct a network meta-analysis assuming a random-effects
#' # model.
#' # The data comes from randomized-controlled trials and
#' # non-randomized studies (combined naively)
#' head(ipddata) # participant-level data
#' stddata # study-level data
#' 
#' # Create a JAGS model
#' mod <- crossnma.model(treat, id, relapse, n, design,
#'   prt.data = ipddata, std.data = stddata,
#'   reference = "A", trt.effect = "random", method.bias = "naive")
#' mod
#' 
#' @method print crossnma.model
#' @export


print.crossnma.model <- function(x, ...) {
  
  print(x$call)
  
  invisible(NULL)
}
