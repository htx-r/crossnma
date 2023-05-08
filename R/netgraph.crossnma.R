#' Produce a network plot
#'
#' @description
#' Create a network plot of the cross network meta-analysis or
#' meta-regression
#'
#' @param x An object produced by \code{\link{crossnma}}.
#' @param ... \dots Additional arguments (passed on to
#'   \code{\link{netgraph.netmeta}})
#' @return A data frame containing the following columns:
#' \item{labels}{Treatment labels.}
#' \item{seq}{Sequence of treatment labels.}
#' \item{xpos}{Position of treatment / edge on x-axis.}
#' \item{ypos}{Position of treatment / edge on y-axis.}
#' \item{zpos}{Position of treatment / edge on z-axis (for 3-D
#'   plots).}
#' \item{xpos.labels}{Position of treatment labels on x-axis (for 2-D
#'   plots).}
#' \item{ypos.labels}{Position of treatment labels on y-axis (for 2-D
#'   plots).}
#' \item{adj.x}{Adjustment for treatment label on x-axis.}
#' \item{adj.y}{Adjustment for treatment label on y-axis.}
#' \item{adj.z}{Adjustment for treatment label on z-axis (for 3-D
#'   plots).}
#'
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}
#'
#' @seealso \code{\link[netmeta]{netgraph.netmeta}}
#'
#' @examples
#' \dontrun{
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
#' fit <- crossnma(mod)
#'
#' # Create network plot
#' netgraph(fit, plastic = FALSE, cex.points = 7, adj = 0.5)
#' }
#'
#'@method netgraph crossnma
#'@export


netgraph.crossnma <- function(x, ...) {
  
  chkclass(x, "crossnma")
  ##
  pw <- crossnma.model2pairwise(x$model)
  ##
  netgraph(suppressWarnings(netmeta(pw, warn = FALSE)), ...)
}
