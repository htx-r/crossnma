#' Run JAGS to fit cross NMA and NMR
#'
#' @description
#' This function takes the JAGS model from an object produced by
#' \code{\link{crossnma.model}} and runs it using \code{jags.model} in
#' rjags package.
#'
#' @param x An object produced by \code{\link{crossnma.model}}.
#' @param inits A list of lists with \code{n.chains} elements; each
#'   element contains initial values for each model parameter or a
#'   function that generates starting values. Default is different
#'   numbers in \code{.RNG.seed} and \code{.RNG.name =
#'   "base::Mersenne-Twister"}.
#' @param n.adapt Number of adaptations for the MCMC chains.
#' @param n.burnin Number of burnin iterations for the MCMC
#'   chains. Default is \code{n.iter / 2} which discards the first
#'   half of the iterations.
#' @param n.iter Number of iterations to run each MCMC chain.
#' @param n.thin Number of thinning for the MCMC chains. Default is
#'   max(1, floor((n.iter - n.burnin) / 1000)), that is only thinning
#'   if there are more than 2000 iterations.
#' @param n.chains Number of MCMC chains.
#' @param monitor A character vector of the names of the parameters to
#'   be monitored. Basic parameters (depends on the analysis) will be
#'   automatically monitored and only additional parameters need to be
#'   specified.
#' @param level.ma The level used to calculate credible intervals for
#'   network estimates.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param quiet A logical passed on to \code{\link{jags.model}}.
#'
#' @return
#' An object of class \code{crossnma} which is a list containing the
#' following components:
#' \item{jagsfit}{An "rjags" object produced when rjags package used
#'   to run the JAGS model.}
#' \item{model}{The \code{crossnma.model} object obtained from
#'   \code{\link{crossnma.model}} which was used to run JAGS.}
#' \item{trt.key}{A table of treatment names and their correspondence
#'   to integers used in the JAGS model.}
#' \item{inits, n.adapt, n.burnin, n.iter}{As defined above.}
#' \item{n.thin, n.chains}{As defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{crossnma} used to create
#'   object.}
#'
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}, Guido
#'   Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{crossnma.model}},
#'   \code{\link[rjags]{jags.model}}
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
#' # Display the output
#' summary(fit)
#' plot(fit)
#' }
#'
#' @export


crossnma <- function(x,
                     inits = NULL,
                     n.adapt = 1000,
                     n.burnin = floor(n.iter / 2),
                     n.iter = 2000,
                     n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
                     n.chains = 2,
                     monitor = NULL,
                     level.ma = x$level.ma,
                     backtransf = x$backtransf,
                     quiet = TRUE
                     ) {

  chkclass(x, "crossnma.model")
  ##
  chknumeric(n.adapt, min = 1, length = 1)
  chknumeric(n.burnin, min = 1, length = 1)
  chknumeric(n.iter, min = 1, length = 1)
  chknumeric(n.thin, min = 1, length = 1)
  chknumeric(n.chains, min = 1, length = 1)
  chklevel(level.ma)
  chklogical(backtransf)
  chklogical(quiet)


  if (!missing(level.ma)) {
    x$level.ma <- level.ma
    x$quantiles =
      c((1 - level.ma) / 2, 0.5, 1 - (1 - level.ma) / 2)
  }


  if (is.null(inits)) {
    seeds <- sample(.Machine$integer.max, n.chains)
    inits <- vector("list", n.chains)
    for (i in seq_len(n.chains))
      inits[[i]] <- list(.RNG.seed = seeds[i],
                         .RNG.name = "base::Mersenne-Twister")
  }


  suppressWarnings(
    jagsfit <-
      jags.model(textConnection(x$model),
                 x$data, inits = inits,
                 n.chains = n.chains, n.adapt = n.adapt,
                 quiet = quiet))
  ##
  if (n.burnin != 0)
    update(jagsfit, n.burnin)


  ##
  ## Monitor (basics)
  ##
  monitor <- c("d", monitor)
  ##
  if (x$trt.effect == "random")
    monitor <- c(monitor, "tau")
  ##
  ## Monitor (meta-regression)
  ##
  if (!is.null(x$covariate)) {
    n.covs <- length(x$covariate) # [[1]]
    id.cov <- seq_len(n.covs)
    ##
    monitor.reg <- character(0)
    ##
    if (x$split.regcoef) {
      if (x$regb.effect == 'independent')
        monitor.reg <- c(monitor.reg, paste0("betab.t_", id.cov))
      ##
      if (x$regw.effect == 'independent')
        monitor.reg <- c(monitor.reg, paste0("betaw.t_", id.cov))
      ##
      monitor.reg <-
        c(monitor.reg, paste0("bb_", id.cov), paste0("bw_", id.cov))
      ##
      if (x$regb.effect == 'random')
        monitor.reg <- c(monitor.reg, paste0("tau.bb_", id.cov))
      ##
      if (x$regw.effect == 'random')
        monitor.reg <- c(monitor.reg, paste0("tau.bw_", id.cov))
    }
    else {
      if (x$regb.effect == 'independent' && x$regw.effect == 'independent')
        monitor.reg <- c(monitor.reg, paste0("beta.t_", id.cov))
      ##
      monitor.reg <- c(monitor.reg, paste0("b_", id.cov))
      ##
      if (x$regb.effect == 'random' && x$regw.effect == 'random')
        monitor.reg <- c(monitor.reg, paste0("tau.b_", id.cov))
    }
    ##
    monitor <- c(monitor, monitor.reg)
  }
  ##
  ## Monitor (bias)
  ##
  if (!is.null(x$method.bias)) {
    monitor.bias <- character(length = 0)
    ##
    if (x$method.bias == 'adjust1') {
      if (x$bias.type == 'both') {
        monitor.bias <- c("g1", "g2")
        ##
        if (length(x$jagsdata$std.act.yes) != 0)
          monitor.bias <- c(monitor.bias, "g.act1", "g.act2")
        ##
        if (x$bias.effect == 'random')
          monitor.bias <- c(monitor.bias, "tau.gamma1", "tau.gamma2")
      }
      else if (x$bias.type %in% c('add', 'mult')) {
        monitor.bias <- c("g")
        ##
        if (length(x$jagsdata$std.act.yes) != 0)
          monitor.bias <- c(monitor.bias, "g.act")
        ##
        if (x$bias.effect == 'random')
          monitor.bias <- c(monitor.bias, "tau.gamma")
      }
    }
    else if (x$method.bias == 'adjust2') {
      monitor.bias <- "g"
      ##
      if (length(x$jagsdata$std.act.yes) != 0)
        monitor.bias <- c(monitor.bias, "g.act")
      ##
      if (x$bias.effect == 'random')
        monitor.bias <- c(monitor.bias, "tau.gamma")

    }
    monitor <- c(monitor, monitor.bias)
  }


  res <- list(
    samples = coda.samples(jagsfit,
                           variable.names = monitor,
                           n.iter = n.iter, n.thin = n.thin),
    jagsfit = jagsfit,
    model = x,
    trt.key = x$trt.key,
    ##
    inits = inits,
    n.adapt = n.adapt,
    n.burnin = n.burnin,
    n.iter = n.iter,
    n.thin = n.thin,
    n.chains = n.chains,
    ##
    call = match.call(),
    version = packageDescription("crossnma")$Version
  )
  ##
  class(res) <- "crossnma"
  ##
  res
}
