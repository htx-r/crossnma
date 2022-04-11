#' Produce a network plot
#' 
#' @description
#' Create a network plot of the cross network meta-analysis or
#' meta-regression
#'
#' @param x An object produced by \code{\link{crossnma.model}}.
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
#' # Create network plot
#' netgraph(mod)
#'
#'@method netgraph crossnma.model
#'@export


netgraph.crossnma.model <- function(x, ...) {
  
  chkclass(x, "crossnma.model")
  ##
  ## Bind variables to function
  trt <- r <- study <- NULL
  ##
  dat <-
    suppressWarnings(pairwise(treat = trt, event = r, n = n, studlab = study,
                              data = x$all.data.ad,
                              sm = "OR", warn = FALSE))
  ##
  netgraph(suppressWarnings(netmeta(dat, warn = FALSE)), ...)
}
