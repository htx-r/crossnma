#' Print summary of cross-design & -format network meta-analysis or
#' regression
#' 
#' @description
#' Print summary results for cross-design and cross-format network
#' meta-analysis or meta-regression
#' 
#' @param x An object of class \code{crossnma}.
#' @param digits The number of significant digits printed. The default
#'   value is 3.
#' @param \dots Additional arguments.
#'
#' @return
#' No return value (print function).
#' 
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}
#'
#' @seealso \code{\link{summary.crossnma}}
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
#'
#' # Display the output (with 5 digits)
#' print(summary(fit), digits = 5)
#' 
#' @method print summary.crossnma
#' @export


print.summary.crossnma <- function(x, digits = 3, ...) {
  
  if (attr(x, "exp"))
    message("Mean and quantiles are exponentiated")
  
  prmatrix(round(x, digits = digits), quote = FALSE, right = TRUE, ...)
  
  invisible(NULL)
}
