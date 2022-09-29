#' Run JAGS to fit cross NMA and NMR
#'
#' @description
#' This function takes the JAGS model from an object produced by
#' \code{\link{crossnma.model}} and runs it using \code{jags.model} in
#' rjags package.
#'
#' @param x An object produced by \code{\link{crossnma.model}}.
#' @param n.adapt Number of adaptations for the MCMC chains. Default
#'   is 1000.
#' @param n.burnin Number of burnin iterations for the MCMC chains.
#' @param n.iter Number of iterations for the MCMC chains.
#' @param thin Number of thinning for the MCMC chains. Default is 1.
#' @param n.chains Number of MCMC chains. Default is 2.
#' @param quiet A logical passed on to
#'   \code{\link[rjags]{jags.model}}.
#' @param monitor A vector of additional parameters to
#'   monitor. Default is NULL.
#'
#' @return
#' An object of class \code{crossnma} which is a list containing the
#' following components:
#' \item{samples}{The MCMC samples produced by running the JAGS
#'   model.}
#' \item{model}{The \code{crossnma.model} object obtained from
#'   \code{\link{crossnma.model}} which was used to run JAGS.}
#' \item{trt.key}{A table of treatment names and their correspondence
#'   to integers used in the JAGS model.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{crossnma} used to create
#'   object.}
#'
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#'
#' @seealso \code{\link{crossnma.model}},
#'   \code{\link[rjags]{jags.model}}
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
#'   prt.data = ipddata, std.data = stddata,sm="OR",
#'   reference = "A", trt.effect = "random", method.bias = "naive")
#'
#' # Fit JAGS model
#' # (suppress warning 'Adaptation incomplete' due to n.adapt = 20)
#' fit <-
#'   suppressWarnings(crossnma(mod, n.adapt = 20,
#'     n.iter = 50, thin = 1, n.chains = 3))
#'
#' # Display the output
#' summary(fit)
#' plot(fit)
#'
#' @export


crossnma <- function(x,
                     n.adapt = 1000,
                     n.burnin = floor(n.iter / 2),
                     n.iter = 10000,
                     thin = 1,
                     n.chains = 2,
                     quiet = TRUE,
                     monitor = NULL
                     ) {

  chkclass(x, "crossnma.model")
  ##
  chknumeric(n.adapt, min = 1, length = 1)
  chknumeric(n.burnin, min = 1, length = 1)
  chknumeric(n.iter, min = 1, length = 1)
  chknumeric(thin, min = 1, length = 1)
  chknumeric(n.chains, min = 1, length = 1)
  chklogical(quiet)


  seeds <- sample(.Machine$integer.max, n.chains, replace = FALSE)
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(.RNG.seed = seeds[i],
                       .RNG.name = "base::Mersenne-Twister")


  suppressWarnings(jagsfit <-
                     jags.model(textConnection(x$model),
                                x$data,
                                n.chains = n.chains,
                                n.adapt = n.adapt,
                                inits = inits,
                                quiet = quiet))


  ## runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
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
    n.covs <- length(x$covariate[[1]])
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


  res <- list(samples = coda.samples(jagsfit,
                                     variable.names = monitor,
                                     n.iter = n.iter, thin = thin),
              model = x,
              trt.key = x$trt.key,
              call = match.call(),
              version = packageDescription("crossnma")$Version)
  ##
  class(res) <- "crossnma"
  ##
  res
}
