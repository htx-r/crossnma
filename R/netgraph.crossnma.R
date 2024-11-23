#' Produce a network plot
#'
#' @description
#' Create a network plot of the cross network meta-analysis or
#' meta-regression
#'
#' @param x An object produced by \code{\link{crossnma}}.
#' @param labels An optional vector with treatment labels.
#' @param adj One, two, or three values in [0, 1] (or a vector /
#'   matrix with length / number of rows equal to the number of
#'   treatments) specifying the x (and optionally y and z) adjustment
#'   for treatment labels.
#' @param offset Distance between edges (i.e. treatments) in graph and
#'   treatment labels for 2-D plots (value of 0.0175 corresponds to a
#'   difference of 1.75\% of the range on x- and y-axis).
#' @param points A logical indicating whether points should be printed
#'   at nodes (i.e. treatments) of the network graph.
#' @param cex.points Corresponding size for points. Can be a vector
#'   with length equal to the number of treatments.
#' @param ... \dots Additional arguments (passed on to
#'   \code{\link{netgraph.netmeta}})
#' 
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


netgraph.crossnma <- function(x,
                              ##
                              labels,
                              adj = NULL,
                              offset =
                                if (!is.null(adj) &&
                                    all(unique(adj) == 0.5))
                                  0
                                else
                                  0.0175,
                              ##
                              points = !missing(cex.points),
                              cex.points = 1,
                              ##
                              ...) {
  
  chkclass(x, "crossnma")
  ##
  pw <- crossnma.model2pairwise(x$model)
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  net <- suppressWarnings(netmeta(pw, warn = FALSE))
  ##
  if (!missing(labels)) {
    ##
    labels <- catch("labels", mc, x, sfsp)
    if (is.null(labels))
      labels <- catch("labels", mc, net, sfsp)
    ##
    if (is.null(labels))
      stop("Argument 'labels' must be not NULL.")
  }
  else
    labels <- net$trts
  ##
  if (!missing(adj)) {
    adj <- catch("adj", mc, net, sfsp)
    if (is.data.frame(adj))
      adj <- as.matrix(adj)
    if (is.logical(adj))
      adj <- 1L * adj
  }
  ##
  if (!missing(offset))
    offset <- catch("offset", mc, net, sfsp)
  ##
  cex.points <- replaceNULL(catch("cex.points", mc, net, sfsp), 1)
  ##
  netgraph(net,
           labels = labels, adj = adj, offset = offset,
           points = points, cex.points = cex.points,
           ...)
}
