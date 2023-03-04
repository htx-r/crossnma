#' Print results of cross-design & -format network meta-analysis or
#' regression
#' 
#' @description
#' Print call used to create JAGS model for cross-design & -format
#' network meta-analysis or regression
#' 
#' @param x An object of class \code{crossnma}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed. If \code{backtransf = TRUE}, results for
#'   \code{sm = "OR"} are presented as odds ratios rather than log
#'   odds ratios, for example.
#' @param digits The number of significant digits printed.
#' @param \dots Additional arguments.
#'
#' @return
#' No return value (print function).
#' 
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}, Guido
#'   Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
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
#' stddata # study-level data
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
#'     n.iter = 50, n.thin = 1, n.chains = 3))
#' fit
#'
#' @method print crossnma
#' @export


print.crossnma <- function(x,
                           backtransf = x$model$backtransf,
                           digits = gs("digits"),
                           ...) {
  
  chkclass(x, "crossnma")
  ##
  chklogical(backtransf)
  chknumeric(digits, min = 0)
  
  sx <- summary(x, quantiles = x$model$quantiles, backtransf = backtransf)
  ##
  n.eff <- formatN(sx[, "n.eff"], digits = 0, text.NA = ".")
  ##
  mat <- formatN(sx, digits = digits, text.NA = ".")
  mat[, "n.eff"] <- n.eff
  ##
  ##sel.d <- startsWith(rownames(mat), "d.")
  ##mat <- rbind(mat[sel.d, , drop = FALSE],
  ##             rep("", ncol(mat)),
  ##             mat[!sel.d, , drop = FALSE])
  ##
  if (attr(sx, "exp")) {
    sel.tau <- startsWith(row.names(mat), "tau")
    row.names(mat)[!sel.tau] <-
      paste0("exp(", row.names(mat)[!sel.tau], ")")
  }
  ##
  prmatrix(mat, quote = FALSE, right = TRUE, ...)
  ##
  if (attr(sx, "exp"))
    message("Mean and quantiles are exponentiated")
  
  invisible(NULL)
}
