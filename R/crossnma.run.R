#' Run JAGS to synthesize cross-design evidence and cross-format data in NMA and NMR for dichotomous outcomes
#' @description Takes the JAGS model from an object produced by \code{crossnma.model} and runs it using \code{jags} package.
#'
#' @param model A \code{crossnmaModel} object produced by \code{crossnma.model}.
#' @param n.adapt Number of adaptations for the MCMC chains.
#' @param n.burnin Number of burnin iterations for the MCMC chains.
#' @param n.iter Number of iterations for the MCMC chains.
#' @param thin Number of thinning for the MCMC chains. Default is 1.
#' @param n.chains Number of MCMC chains. Default is 2.
#' @param quiet A logical. If TRUE, the warning message will not be displayed
#' @param monitor A vector of additional parameters to monitor. Default is NULL
#' See \code{\link{jags.model}} for more info.
#'
#' @return \code{crossnma.run} returns an object of class \code{crossrun} which is a list containing the following components:
#' @return \code{samples}  The MCMC samples produced by running the JAGS model.
#' @return \code{model}  The \code{crossnmaModel} object obtained from \code{crossnma.model} which was used to run \code{jags}.
#' @return \code{trt.key}  A table of the treatment names and their correspondence to integers used in the JAGS model.
#' @examples
#' # An example from participant-level data and study-level data.
#' # data
#' data(prt.data)
#' data(std.data)
#'  #=========================#
#'   # Create a jags model  #
#'  #=========================#
#'  # We conduct a network meta-analysis assuming a random effect model.
#'  # The data comes from randomised-controlled trials and non-randomised studies. They will be combined naively.
#'  # The data has 2 different formats: individual participant data (prt.data) and study-level data (std.data).
#' mod <- crossnma.model(prt.data=prt.data,
#'                   std.data=std.data,
#'                   trt="trt",
#'                   study="study",
#'                   outcome="outcome",
#'                   n="n",
#'                   design="design",
#'                   reference="A",
#'                   trt.effect="random",
#'                   covariate = NULL,
#'                   method.bias="naive"
#'                    )
#'  #=========================#
#'     # Fit jags model  #
#'  #=========================#
#' fit <- crossnma.run(model=mod,
#'                 n.adapt = 20,
#'                 n.iter=50,
#'                 thin=1,
#'                 n.chains=3)
#'
#'  #=========================#
#'    # Display the output   #
#'  #=========================#
#' summary(fit)
#' plot(fit)
#'
#'
#' @seealso \code{\link{crossnma.model}},\code{\link{jags.model}}
#' @export

crossnma.run <- function(model,
                         n.adapt = 1000,
                         n.burnin = floor(n.iter / 2),
                         n.iter,
                         thin=1,
                         n.chains=2,
                         quiet=TRUE,
                         monitor=NULL
                         ){

  if(class(model) != "crossnmaModel")
    stop("\'model\' must be a valid crossnmaModel object created using the crossnma.model function.")


  seeds <- sample(.Machine$integer.max, n.chains, replace = FALSE)
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(.RNG.seed = seeds[i], .RNG.name = "base::Mersenne-Twister")


  jagsfit <- jags.model(textConnection(model$jags),        #Create a connection so JAGS can access the variables
                        model$data,
                        n.chains=n.chains,
                        n.adapt=n.adapt,
                        inits = inits,
                        quiet=quiet)

  # runjags.options(silent.jags=TRUE, silent.runjags=TRUE)
  if(n.burnin!=0) update(jagsfit, n.burnin)

  # ---- monitor
  # basics
  make.monitor <- c("d",monitor)
  if(model$trt.effect=="random") make.monitor <- c(make.monitor, "tau")
  # meta-regression
  if(!is.null(model$covariate)){
    if(model$split.regcoef){
      make.monitor.reg <- c()
      if(model$regb.effect=='independent'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg0 <- paste0("betab.t_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg0)
        }
      }
      if(model$regw.effect=='independent'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg0 <- paste0("betaw.t_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg0)
        }
      }
      for (i in 1:length(model$covariate[[1]])) {
        make.monitor.reg0 <- c(paste0("bb_",i),paste0("bw_",i))
        make.monitor.reg <- c(make.monitor.reg,make.monitor.reg0)
      }
      if(model$regb.effect=='random'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg1 <- paste0("tau.bb_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg1)
        }
      }
      if(model$regw.effect=='random'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg2 <- paste0("tau.bw_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg2)
        }
      }

    }else{
      make.monitor.reg <- c()
      if(model$regb.effect=='independent' && model$regw.effect=='independent'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg0 <- paste0("beta.t_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg0)
        }
      }
      for (i in 1:length(model$covariate[[1]])) {
        make.monitor.reg0 <- paste0("b_",i)
        make.monitor.reg <- c(make.monitor.reg,make.monitor.reg0)
      }
      if(model$regb.effect=='random'&&model$regw.effect=='random'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg1 <- paste0("tau.b_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg1)
        }
      }

    }
    make.monitor <- c(make.monitor, make.monitor.reg)
  }

  if(!is.null(model$method.bias)){
    make.monitor.bias <- c()
    if(model$method.bias=='adjust1'){
      if(model$bias.type=='both'){
        make.monitor.bias <- c("g1","g2")
        if(length(model$jagsdata$std.act.yes)!=0) make.monitor.bias <- c(make.monitor.bias,c("g.act1","g.act2"))
        if(model$bias.effect=='random') make.monitor.bias <- c(make.monitor.bias,"tau.gamma1","tau.gamma2")
      }else if(model$bias.type%in%c('add','mult')){
        make.monitor.bias <- c("g")
        if(length(model$jagsdata$std.act.yes)!=0) make.monitor.bias <- c(make.monitor.bias,"g.act")
        if(model$bias.effect=='random') make.monitor.bias <- c(make.monitor.bias,"tau.gamma")
      }
    }
    if(model$method.bias=='adjust2'){
      make.monitor.bias <- c("g")
      if(length(model$jagsdata$std.act.yes)!=0) make.monitor.bias <- c(make.monitor.bias,"g.act")
      if(model$bias.effect=='random') make.monitor.bias <- c(make.monitor.bias,"tau.gamma")

    }
    make.monitor <- c(make.monitor, make.monitor.bias)
  }


  jagssamples <- coda.samples(jagsfit,
                              variable.names=make.monitor,
                              n.iter=n.iter,
                              thin=thin)

  crossrun <- structure(list(samples=jagssamples,
                             model=model,
                             "trt.key"=model$trt.key),
                        class = "crossnma")
  return(crossrun)
}


