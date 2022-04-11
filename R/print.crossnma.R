#' Print call used to create JAGS model for cross-design & -format
#' network meta-analysis or regression
#' 
#' @description
#' Print call used to create JAGS model for cross-design & -format
#' network meta-analysis or regression
#' 
#' @param x An object of class \code{crossnma}.
#' @param \dots Additional arguments (ignored).
#'
#' @return
#' No return value (print function).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#'
#' @seealso \code{\link{crossnma}}
#' 
#' @keywords print
#' 
#' @examples
#' # We conduct a network meta-analysis assuming a random-effects
#' # model.
#' # The data comes from randomized-controlled trials and
#' # non-randomized studies (combined naively)
#' head(ipddata) # participant-level data
#' head(stddata) # study-level data
#' 
#' # Create a JAGS model
#' mod <- crossnma.model(treat, id, relapse, n, design,
#'   prt.data = ipddata, std.data = stddata,
#'   reference = "A", trt.effect = "random", method.bias = "naive")
#' 
#' # Fit JAGS model
#' # (suppress warning 'Adaptation incomplete' due to n.adapt = 20)
#' fit <-
#'   suppressWarnings(crossnma(mod, n.adapt = 20,
#'     n.iter = 50, thin = 1, n.chains = 3))
#' fit
#'
#' @method print crossnma
#' @export


print.crossnma <- function(x, ...) {

  cat("Model:\n")
  print(x$model)

  cat("\nTreatment coding:\n")
  print(x$trt.key)
  
  invisible(NULL)
}
