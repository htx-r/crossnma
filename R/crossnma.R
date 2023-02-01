#' Run JAGS to fit cross NMA and NMR
#'
#' @description
#' This function takes the JAGS model from an object produced by
#' \code{\link{crossnma.model}} and runs it using \code{jags.parallel} in
#' R2jags package.
#'
#' @param x An object produced by \code{\link{crossnma.model}}.
#' @param inits A list with n.chains elements; each element is itself a list of initial values for each model parameter,
#' or a function that generates starting values. Default is NULL where jags.parallel() function
#' will internally generate starting values for each parameter.
#' @param n.burnin Number of burnin iterations for the MCMC chains. Default is n.iter/2 which discards the first half of the iterations.
#' @param n.iter Number of iterations to run each MCMC chain. Default is 2000.
#' @param n.thin Number of thinning for the MCMC chains. Default is max(1, floor((n.iter - n.burnin) / 1000)),
#' that is only thinning if there are more than 2000 iterations.
#' @param n.chains Number of MCMC chains. Default is 2.
#' @param monitor A character vector of the names of the parameters to
#'   be monitored. Basic parameters (depends on the analysis) will be automatically monitored
#'   and only additional parameters need to be specified.
#'   Default is NULL.
#'
#' @return
#' An object of class \code{crossnma} which is a list containing the
#' following components:
#' \item{jagsfit}{ An rjags-class object produced when R2jags package used to run the JAGS
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
#'   \code{\link[R2jags]{jags.parallel}}
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
#' fit <-crossnma(mod,
#' n.burnin =10,n.iter = 50, n.thin = 1, n.chains = 3)
#'
#' # Display the output
#' summary(fit)
#' plot(fit)
#'
#' @export


crossnma <- function(x,
                     inits=NULL,
                     #n.adapt = 1000,
                     n.burnin = floor(n.iter / 2),
                     n.iter = 2000,
                     n.thin = max(1, floor((n.iter - n.burnin)/1000)),
                     n.chains = 2,
                     #quiet = TRUE,
                     monitor = NULL
                     ) {

  chkclass(x, "crossnma.model")
  ##
  #chknumeric(n.adapt, min = 1, length = 1)
  chknumeric(n.burnin, min = 1, length = 1)
  chknumeric(n.iter, min = 1, length = 1)
  chknumeric(n.thin, min = 1, length = 1)
  chknumeric(n.chains, min = 1, length = 1)
  #chklogical(quiet)


  # inits <- list()
  # for (i in 1:n.chains) inits[[i]] <- seeds[i]
    # inits[[i]] <- list(.RNG.seed = seeds[i],
    #                    .RNG.name = "base::Mersenne-Twister")


  # suppressWarnings(jagsfit <-
  #                    jags.model(textConnection(x$model),
  #                               x$data,
  #                               n.chains = n.chains,
  #                               n.adapt = n.adapt,
  #                               inits = inits,
  #                               quiet = quiet))


  ## runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  ##
  # if (n.burnin != 0)
  #   update(jagsfit, n.burnin)

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

  # random seed for JAGS model (to produce identical results)
  seeds <- sample(.Machine$integer.max, n.chains, replace = FALSE)
  # Run JAGS model
  jmodel <- x$model
  jagsfit <- jags.parallel(data=x$data,
                inits = inits,
                parameters.to.save = monitor,
                model.file=jmodel,
                n.chains = n.chains,
                n.iter = n.iter,
                n.burnin = n.burnin,
                n.thin = n.thin,
                jags.seed=seeds,
                DIC=FALSE)
  res <- list(# samples = coda.samples(jagsfit,
    #                                  variable.names = monitor,
    #                                  n.iter = n.iter, thin = thin),
              jagsfit=jagsfit,
              model = x,
              trt.key = x$trt.key,
              call = match.call(),
              version = packageDescription("crossnma")$Version)
  ##
  class(res) <- "crossnma"
  ##
  res
}
