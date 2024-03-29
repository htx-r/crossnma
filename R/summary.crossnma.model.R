#' Summary function for crossnma.model object
#'
#' @description
#' Summary function for crossnma.model object
#'
#' @param object An object generated by the
#'   \code{\link{crossnma.model}}.
#' @param \dots Additional arguments (ignored)
#'
#' @author Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{print.summary.crossnma.model}}
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
#'
#' summary(mod)
#'
#' @method summary crossnma.model
#' @export


summary.crossnma.model <- function(object, ...) {

  chkclass(object, "crossnma.model")
  ##
  res <- object
  class(res) <- c("summary.crossnma.model", class(res))
  ##
  res
}
