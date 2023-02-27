#' Print summary of cross-design & -format network meta-analysis or
#' regression
#'
#' @description
#' Print results of cross-design and cross-format network
#' meta-analysis or meta-regression. In addition, the call used to
#' create the JAGS model is printed.
#'
#' @param x An object of class \code{crossnma}.
#' @param digits The number of significant digits printed. The default
#'   value is 3.
#' @param \dots Additional arguments.
#'
#' @return
#' No return value (print function).
#'
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}, Guido
#'   Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
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
#' stddata # study-level data
#'
#' # Create a JAGS model
#' mod <- crossnma.model(treat, id, relapse, n, design,
#'   prt.data = ipddata, std.data = stddata,
#'   reference = "A", trt.effect = "random", method.bias = "naive")
#'
#' # Fit JAGS model
#' set.seed(1909)
#' fit <- crossnma(mod, n.burnin = 10, n.iter = 50,
#'   n.thin = 1, n.chains = 3)
#'
#' # Display the output (with 5 digits)
#' print(summary(fit), digits = 5)
#'
#' @method print summary.crossnma
#' @export


print.summary.crossnma <- function(x, digits = gs("digits"), ...) {
  
  chkclass(x, "summary.crossnma")
  ##
  object <- attr(x, "object")
  
  chkclass(x, "summary.crossnma")
  cat("Model:\n")
  print(object$model)
  
  cat("\nTreatment coding:\n")
  print(object$model$trt.key)
  cat("\n")
  
  pw <- crossnma.model2pairwise(object$model)
  net <- suppressWarnings(netmeta(pw, warn = FALSE))
  ##
  cat(paste("Number of studies: k = ", net$k, "\n", sep = ""))
  cat(paste0("Number of pairwise comparisons: m = ", net$m, "\n"))
  if (!is.null(net$n.trts))
    cat(paste0("Number of observations: o = ",
               round(sum(net$n.trts, na.rm = TRUE), 1),
               "\n"))
  cat(paste0("Number of treatments: n = ", net$n, "\n"))
  cat(paste0("Number of designs: d = ", net$d, "\n"))
  cat("\n")
  
  n.eff <- formatN(x[, "n.eff"], digits = 0, text.NA = ".")
  ##
  mat <- formatN(x, digits = digits, text.NA = ".")
  mat[, "n.eff"] <- n.eff
  ##
  ##sel.d <- startsWith(rownames(mat), "d.")
  ##mat <- rbind(mat[sel.d, , drop = FALSE],
  ##             rep("", ncol(mat)),
  ##             mat[!sel.d, , drop = FALSE])
  ##
  if (attr(x, "exp")) {
    sel.tau <- startsWith(row.names(mat), "tau")
    row.names(mat)[!sel.tau] <-
      paste0("exp(", row.names(mat)[!sel.tau], ")")
  }
  ##
  prmatrix(mat, quote = FALSE, right = TRUE, ...)
  ##
  if (attr(x, "exp"))
    message("Mean and quantiles are exponentiated")
  
  invisible(NULL)
}
