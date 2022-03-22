#' Creates a JAGS model and the needed data to perform cross NMA and
#' NMR (dichotomous outcomes)
#' @description This function creates a JAGS model and the needed
#'   data. The JAGS code is created from the internal function
#'   \code{crossnma.code}.
#' @param trt Treatment variable in prt.data and std.data.
#' @param study Study variable in prt.data and std.data.
#' @param outcome Outcome variable in prt.data and std.data.
#' @param n Number of participants in std.data.
#' @param design Design variable in prt.data and std.data.
#' @param prt.data An object of class data.frame containing the
#'   individual participant dataset. Each row contains the data of a
#'   single participant.  The dataset needs to have the following
#'   columns: treatment, study identification, outcome (event and
#'   non-event), design. Additional columns might be required for
#'   certain analyses.
#' @param std.data An object of class data.frame containing the
#'   study-level dataset. Each row represents the information of study
#'   arm.  The dataset needs to have the following columns: treatment,
#'   study identification, outcome (number of events), sample size and
#'   design. Additional columns might be required for certain
#'   analyses.
#' @param reference A character indicating the name of the reference
#'   treatment. When the reference is not specified, the first
#'   alphabetic treatment will be used as a reference in the analysis.
#' @param trt.effect A character defining the model for the
#'   study-specific treatment effects. Options are 'random' (default)
#'   or 'common'.
#' @param cov1 Optional first covariate in prt.data and std.data to
#'   conduct network meta-regression (see Details).
#' @param cov2 Optional second covariate in prt.data and std.data to
#'   conduct network meta-regression (see Details).
#' @param cov3 Optional third covariate in prt.data and std.data to
#'   conduct network meta-regression (see Details).
#' @param cov1.ref An optional value to center the first covariate
#'   which is only useful for a continuous covariate. Dichotomous
#'   covariates should be given NA value. The default is the overall
#'   minimum covariate value from all studies.
#' @param cov2.ref An optional value to center the second covariate
#'   which is only useful for a continuous covariate. Dichotomous
#'   covariates should be given NA value. The default is the overall
#'   minimum covariate value from all studies.
#' @param cov3.ref An optional value to center the third covariate
#'   which is only useful for a continuous covariate. Dichotomous
#'   covariates should be given NA value. The default is the overall
#'   minimum covariate value from all studies.
#' @param reg0.effect An optional character (needed when
#'   \code{covariate} is not NULL) indicating the relationship across
#'   studies for the prognostic effects expressed by the regression
#'   coefficient, (\eqn{\beta_0}), in a study \eqn{j}.  Options are
#'   'independent' or 'random'. We recommend using 'independent'
#'   (default).
#' @param regb.effect An optional character (needed when
#'   \code{covariate} is not NULL) indicating the relationship across
#'   treatments for the between-study regression coefficient
#'   (\eqn{\beta^B}). This parameter quantifies the treatment-mean
#'   covariate interaction.  Options are 'independent', 'random' or
#'   'common'. Default is 'random'.
#' @param regw.effect An optional character (needed when
#'   \code{covariate} is not NULL) indicating the relationship across
#'   treatments for the within-study regression coefficient
#'   (\eqn{\beta^W}). This parameter quantifies the
#'   treatment-covariate interaction effect at the individual level.
#'   Options are 'independent', 'random' and 'common'. Default is
#'   'random'.
#' @param split.regcoef A logical value (needed when \code{covariate}
#'   is not NULL). If TRUE (default) the within- and between-study
#'   coefficients will be splitted in the analysis of prt.data.  When
#'   the split.regcoef = FALSE, only a single regression coefficient
#'   will be estimated to represent both the between-studies and
#'   within-studies covariate effects.  In this case, both arguments
#'   \code{regb.effect} and \code{regw.effect} need to be given the
#'   same option to model the single regression effect.
#' @param method.bias A character for defining the method to combine
#'   randomized clinical trials (RCT) and non-randomized studies
#'   (NRS).  Options are 'naive' for naive or unadjusted synthesize,
#'   'prior' for using NRS evidence to construct priors for the
#'   relative treatment effects in RCTs analysis, or 'adjust1' and
#'   'adjust2' to allow a bias adjustment. When only one design is
#'   available (either rct or nrs), this argument needs also to be
#'   specified to indicate whether unadjusted (naive) or bias-adjusted
#'   analysis (adjust1 or adjust2) should be applied.
#' @param bias A character indicating the name of the variable
#'   (required when method.bias= 'adjust1' or 'adjust2') that includes
#'   the risk of bias of the study in prt.data and std.data.  The
#'   values of this variable should be a character with entries that
#'   need to be spelled as such 'low', 'high' or 'unclear'. These
#'   values need to be repeated for the participants that belong to
#'   the same study.
#' @param bias.group An optional character for defining the name of
#'   the variable in prt.data and std.data that indicates the bias
#'   effect in each study (can be provided when method.bias= 'adjust1'
#'   or 'adjust2').  The entries of these variables should be either 1
#'   (study has inactive treatment and its estimate should be adjusted
#'   for bias effect), 2 (study has only active treatments and its
#'   estimate should be adjusted for bias effect (different from
#'   inactive bias effect) or 0 (study does not need any bias
#'   adjustment). The values need to be repeated for the participants
#'   assigned to the same treatment. Default is 1.
#' @param unfav A character for defining the name of the variable
#'   which indicates the unfavored treatment in each study (should be
#'   provided when method.bias= 'adjust1' or 'adjust2') in prt.data
#'   and std.data.  The entries of this variable are either 0
#'   (unfavored treatment) or 1 (favorable treatment or
#'   treatments). Each study should include only one 0 entry. The
#'   values need to be repeated for participants who take the same
#'   treatment.
#' @param bias.type An optional character defining of bias on the
#'   treatment effect (required when method.bias='adjust1').  Three
#'   options are possible: 'add' to add the additive bias
#'   effect,'mult' for multiplicative bias effect and 'both' includes
#'   both an additive and a multiplicative terms.
#' @param bias.covariate A character of the variable name that will be
#'   used in estimating the probability of bias (can be provided when
#'   method.bias='adjust1' or 'adjust2')
#' @param bias.effect An optional character indicating the
#'   relationship for the bias coefficients across studies.  Options
#'   are 'random' or 'common' (default). It is required when
#'   method.bias='adjust1' or 'adjust2'.
#' @param down.wgt An optional numeric indicating the percent to which
#'   studies at high risk of bias will be downweighed on average. The
#'   value ranges between 0 and 1. It can be provided when
#'   method.bias='adjust1' or 'adjust2'.
#' @param prior An optional list to control the prior for various
#'   parameters in JAGS model. When effects are set as 'random', we
#'   can set the heterogeneity parameters for: tau.trt for the
#'   treatment effects, tau.reg0 for the effect of prognostic
#'   covariates, tau.regb and tau.regw for within- and between-study
#'   covariate effect, respectively.  and tau.gamma for bias
#'   effect. The default of all heterogeneity parameters is
#'   'dunif(0,2)'. Currently only the uniform distribution is
#'   supported.  When the method.bias= 'adjust1' or 'adjust2', the
#'   user may provide priors to control the bias probability.  For the
#'   bias probabilities, beta distributions are assumed with the
#'   following default values: RCT with low
#'   (pi.low.rct='dbeta(1,10)'), high (pi.high.rct='dbeta(10,1)')
#'   bias, NRS with low (pi.low.rct='dbeta(1,30)') / high
#'   (pi.high.rct='dbeta(30,1)') bias (pi.low.nrs, pi.high.nrs).
#' @param run.nrs An optional list is needed when the NRS used as a
#'   prior (method.bias='prior').  The list consists of the following:
#'   (\code{var.infl}) controls the common inflation of the variance
#'   of NRS estimates (\eqn{w}) and its values range between 0 (NRS
#'   does not contribute at all and the prior is vague) and 1 (the NRS
#'   evidence is used at face value, default approach).  The parameter
#'   (\code{mean.shift}) is the bias shift (\eqn{\zeta}) to be
#'   added/subtracted from the estimated mean treatment effects (on
#'   the log-scale) from NRS network (0 is the default).
#'   \code{trt.effect} is a character indicates how to combine
#'   treatment effects across NRS studies .Options are 'random' or
#'   'common' (default).  Here you can also specify the arguments to
#'   control the MCMC chains with default value is in the parentheses:
#'   the number of adaptions n.adapt (500), number of iterations
#'   n.iter(10000), number of burn in n.burnin (4000), number of
#'   thinning thin (1) and number of chains n.chains
#'   (2). \code{\link{jags.model}} from rjags package describes these
#'   arguments.
#'
#' @details
#'  Covariates provided in arguments \code{cov1}, \code{cov2} and
#'  \code{cov3} can be either numeric or dichotomous (should be
#'  provided as factor or character) variables. By default, no
#'  covariate adjustment is applied (network meta-analysis).
#'
#' 
#' @return \code{crossnma.model} returns an object of class
#'   \code{crossnmaModel} which is a list containing the following
#'   components (the JAGS model to be ran under these settings):
#' @return \code{jagsmodel} A long character string containing JAGS
#'   code that will be run in \code{\link{jags.model}}.
#' @return \code{data} The data to be used to run JAGS model.
#' @return \code{trt.key} A table of the treatments and its mapped
#'   integer number (as used in JAGS model).
#' @return \code{study.key} A table of the studies and its mapped
#'   integer number (as used in JAGS model).
#' @return \code{trt.effect} A character defining the model for the
#'   study-specific treatment effects.
#' @return \code{method.bias} A character for defining the method to
#'   combine randomized clinical trials (RCT) and non-randomized
#'   studies (NRS).
#' @return \code{covariate} A list of the the names of the covariates
#'   in prt.data and std.data used in network meta-regression.
#' @return \code{cov1.ref, cov2.ref, cov3.ref} Values to center
#'   continuous covariates. Dichotomous covariates take NA.
#' @return \code{cov1.labels, cov2.labels, cov3.labels} A matrix with the levels of each
#'   dichotomous covariate and the corresponding assigned 0/1 values.
#' @return \code{split.regcoef} A logical value. If FALSE the within-
#'   and between-study regression coefficients will be considered
#'   equal.
#' @return \code{regb.effect} A character indicating the model for the
#'   between-study regression coefficients across studies.
#' @return \code{regw.effect} A character indicating the model for the
#'   within-study regression coefficients across studies.
#' @return \code{bias.effect} A character indicating the model for the
#'   bias coefficients across studies.
#' @return \code{bias.type} A character indicating the effect of bias
#'   on the treatment effect; additive ('add') or multiplicative
#'   ('mult') or both ('both').
#' @return \code{all.data.ad} A data.frame object with the prt.data
#'   (after it is aggregated) and std.data in a single dataset.
#' 
#' @examples
#' # We conduct a network meta-analysis assuming a random-effects
#' # model.
#' # The data comes from randomized-controlled trials and
#' # non-randomized studies (combined naively)
#' head(ipddata) # participant-level data
#' head(stddata) # study-level data
#' 
#' #=========================#
#' # Create a jags model     #
#' #=========================#
#' mod <- crossnma.model(treat, id, relapse, n, design,
#'   prt.data = ipddata, std.data = stddata,
#'   reference = "A", trt.effect = "random", method.bias = "naive")
#' 
#' #=========================#
#' # Fit jags model          #
#' #=========================#
#' fit <-
#'   crossnma.run(model = mod, n.adapt = 20,
#'     n.iter = 50, thin = 1, n.chains = 3)
#'
#' #=========================#
#' # Display the output      #
#' #=========================#
#' summary(fit)
#' plot(fit)
#'
#' @seealso \code{\link{crossnma.run}}, \code{\link{jags.model}}
#' @export

crossnma.model <- function(trt,
                           study,
                           outcome,
                           n,
                           design,
                           ##
                           cov1 = NULL,
                           cov2 = NULL,
                           cov3 = NULL,
                           ##
                           bias = NULL,
                           unfav = NULL,
                           bias.covariate = NULL,
                           bias.group = NULL,
                           ##
                           prt.data = NULL,
                           std.data = NULL,
                           ##
                           reference=NULL,
                           trt.effect = 'random',
                           #---------- meta regression ----------
                           cov1.ref = NULL,
                           cov2.ref = NULL,
                           cov3.ref = NULL,
                           reg0.effect = 'independent',
                           regb.effect = 'random',
                           regw.effect = 'random',
                           split.regcoef = TRUE,
                           #---------- bias adjustment ----------
                           method.bias = NULL,
                           bias.type=NULL,
                           bias.effect = 'common',
                           down.wgt=NULL,
                           # ---------- prior ----------
                           prior=list(tau.trt=NULL,
                                      tau.reg0=NULL,
                                      tau.regb=NULL,
                                      tau.regw=NULL,
                                      tau.gamma=NULL,
                                      pi.high.rct=NULL,
                                      pi.low.rct=NULL,
                                      pi.high.nrs=NULL,
                                      pi.low.nrs=NULL
                           ),
                           # ---------- when method.bias='prior' ----------
                           run.nrs=list(var.infl=1,
                                        mean.shift=0,
                                        trt.effect="common",
                                        n.adapt = 2000,
                                        n.iter=10000,
                                        n.burnin = 4000,
                                        thin=1,
                                        n.chains=2)) {
  
  ## Check and set variables
  ##
  if (missing(trt))
    stop("Mandatory argument 'trt' missing.")
  if (missing(study))
    stop("Mandatory argument 'study' missing.")
  if (missing(outcome))
    stop("Mandatory argument 'outcome' missing.")
  if (missing(design))
    stop("Mandatory argument 'design' missing.")
  ##
  trt.effect <- setchar(trt.effect, c("common", "random"))
  ##
  reg0.effect <- setchar(reg0.effect, c("independent", "random"))
  regb.effect <- setchar(regb.effect, c("common", "independent", "random"))
  regw.effect <- setchar(regw.effect, c("common", "independent", "random"))
  ##
  if (!is.null(split.regcoef) && !is.logical(split.regcoef))
    stop("Argument 'split.regcoef' must be a logical value.")
  ##
  if (!is.null(method.bias))
    method.bias <-
      setchar(method.bias, c("naive", "prior", "adjust1", "adjust2"))
  ##
  if (!is.null(bias.type))
    bias.type <-
      setchar(bias.type, c("add", "mult", "both"))
  else if (!is.null(method.bias) && method.bias == "adjust1")
    stop("Argument 'bias.type' must be provided if method.bias = \"adjust1\".")
  ##
  bias.effect <- setchar(bias.effect, c("common", "random"))
  ##
  if (!is.null(down.wgt)) {
    if (!length(down.wgt == 1) && !(down.wgt > 0 & down.wgt < 1))
      stop("Values of argument 'down.wgt' should be between 0 and 1")
    if (!method.bias %in% c('adjust1','adjust2'))
      stop("'down.wgt' should be specified only if method.bias = 'adjust1' or 'adjust2'")
  }
  ##
  if (!is.null(prior)) {
    if (trt.effect == 'common' & !is.null(prior$tau.trt))
      warning("The prior of the heterogeneity between relative treatments parameters is ignored")
    if (reg0.effect == 'common' & !is.null(prior$tau.reg0))
      warning("The prior of the heterogeneity between progonostic parameters is ignored")
    if (regw.effect == 'common' & !is.null(prior$tau.regw))
      warning("The prior of the heterogeneity between within-study interaction parameters is ignored")
    if (regb.effect == 'common' & !is.null(prior$tau.regb))
      warning("The prior of the heterogeneity between between-study interaction parameters is ignored")
    if (bias.effect == 'common' & !is.null(prior$tau.gamma))
      warning("The prior of the heterogeneity between bias effect parameters is ignored")
    if (!is.null(down.wgt) & !is.null(prior$tau.gamma))
      message("The assigned prior for the heterogeneity parameters of bias effect is ignored when 'down.wgt' is provided")
  }
  
  
  options(warn=-1)
  ## Bind variables to function
  ##
  study.jags <- trt.ini <- trt.jags <-
    arm <- value <- variable <- bias_index <-
      x.bias <- x1 <- x1f <- x2 <- x2f <- x3 <- x3f <-
        ref.trt.std <- n.arms <- NULL
  
  ## Extract names of covariates
  ##
  covariates <- NULL
  if (!missing(cov1)) {
    covariates <- deparse(substitute(cov1))
    if (!missing(cov2)) {
      covariates <- c(covariates, deparse(substitute(cov2)))
      if (!missing(cov3)) {
        covariates <- c(covariates, deparse(substitute(cov3)))
      }
    }
  }
  
  
  ## Prepare IPD dataset
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  excl1 <- FALSE
  ##
  if (!is.null(prt.data)) {

    trt <- catch("trt", mc, prt.data, sfsp)
    ##
    study <- catch("study", mc, prt.data, sfsp)
    if (is.factor(study))
      study <- as.character(study)
    ##
    outcome <- catch("outcome", mc, prt.data, sfsp)
    if (!is.numeric(outcome) | any(!outcome %in% 0:1))
      stop("Outcome values must be either 0 or 1 (IPD dataset).")
    ##
    design <- catch("design", mc, prt.data, sfsp)
    design <- as.character(design)
    design <- setchar(design, c("nrs", "rct"),
                      text = paste("must contain values",
                                   '"nrs" or "rct"',
                                   "(IPD dataset)"))
    ##
    bias <- catch("bias", mc, prt.data, sfsp)
    if (is.null(bias) && !is.null(method.bias) &&
        method.bias %in% c("adjust1", "adjust2"))
      stop("Argument 'bias' must be provided if ",
           "method.bias = \"adjust1\" (IPD dataset).")
    ##
    bias.covariate <- catch("bias.covariate", mc, prt.data, sfsp)
    bias.group <- catch("bias.group", mc, prt.data, sfsp)
    unfav <- catch("unfav", mc, prt.data, sfsp)
    ##
    data11 <- data.frame(trt = trt, study = study,
                         r = outcome,
                         design = design)
    ##
    if (!is.null(bias)) {
      bias <- as.character(bias)
      bias <- setchar(bias, c("low", "high", "unclear"),
                      text = paste("must contain values",
                                   '"low", "high", or "unclear"',
                                   "(IPD dataset)"))
      data11$bias <- bias
    }
    if (!is.null(bias.covariate))
      data11$x.bias <- bias.covariate
    if (!is.null(bias.group))
      data11$bias.group <- bias.group
    ##
    if (!is.null(unfav)) {
      if (!is.numeric(unfav) | any(!unfav %in% 0:1))
        stop("Values of argument 'unfav' must be either 0 or 1 (IPD dataset).")
      data11$unfav <- unfav
    }
    ##
    cov1 <- catch("cov1", mc, prt.data, sfsp)
    cov2 <- catch("cov2", mc, prt.data, sfsp)
    cov3 <- catch("cov3", mc, prt.data, sfsp)
    ##
    if (!is.null(cov1))
      data11$x1 <- cov1
    if (!is.null(cov2))
      data11$x2 <- cov2
    if (!is.null(cov3))
      data11$x3 <- cov3
    ##
    data11$study <- paste0(data11$study,".ipd")
    ##
    ## Delete NAs
    ##
    nam <- c("trt", "study", "outcome", "design", "bias", "unfav", "bias.group")
    nam <- nam[nam %in% names(data11)]
    excl1 <- apply(is.na(data11[, nam]), 1, anyNA)
    if (any(excl1))
      data11 <- data11[!excl1, ]
    ##
    ## check unfav: unique value 0 per study (repeated for the same treatment)
    ##
    if (!is.null(data11$unfav)){
      chk.unfav1 <- data11 %>%
        group_by(study) %>%
        group_map(~length(unique(subset(.x, unfav == 0, select = c(trt)))) != 1) %>%
        unlist()
      if (any(chk.unfav1)) stop("For 'unfav' variable in prt.data, each study could be provided by only a unique 0 for a specific treatment")
    }
    ##
    ## check unique bias per study
    ##
    if (!is.null(data11$bias)){
      chk.bias1 <- data11 %>%
        group_by(study) %>%
        group_map(~length(unique(.x$bias)) !=1 ) %>%
        unlist()
      ##
      if (any(chk.bias1)) stop("The 'bias' should be a vector of length 2 where the first element is the name of the variable in prt.data and the second for the std.data")
    }
  }
  else {
    data11 <- NULL
  }
  
  
  ## Prepare AD dataset
  ##
  excl2 <- FALSE
  ##
  if (!is.null(std.data)){

    if (missing(outcome))
      stop("Mandatory argument 'outcome' missing.")
    
    trt <- catch("trt", mc, std.data, sfsp)
    ##
    study <- catch("study", mc, std.data, sfsp)
    if (is.factor(study))
      study <- as.character(study)
    ##
    outcome <- catch("outcome", mc, std.data, sfsp)
    if (!is.numeric(outcome) | any(outcome < 0) | any(!outcome %% 1 == 0))
      stop("Outcome values must be integers greater than or equal to 0 ",
           "(study-level dataset).")
    ##
    n <- catch("n", mc, std.data, sfsp)
    if (!is.numeric(n) | any(n <= 0) | any(!n %% 1 == 0))
      stop("Sample sizes must be integers greater than 0 ",
           "(study-level dataset).")
    if (any(outcome > n))
      stop("Sample sizes must be larger than number of events ",
           "(study-level dataset).")
    ##
    design <- catch("design", mc, std.data, sfsp)
    design <- as.character(design)
    design <- setchar(design, c("nrs", "rct"),
                      text = paste("must contain values",
                                   '"nrs" or "rct"',
                                   "(study-level dataset)"))
    ##
    bias <- catch("bias", mc, std.data, sfsp)
    if (is.null(bias) && !is.null(method.bias) &&
        method.bias %in% c("adjust1", "adjust2"))
      stop("Argument 'bias' must be provided if ",
           "method.bias = \"adjust1\" (study-level dataset).")
    ##
    ##
    bias.covariate <- catch("bias.covariate", mc, std.data, sfsp)
    bias.group <- catch("bias.group", mc, std.data, sfsp)
    unfav <- catch("unfav", mc, std.data, sfsp)
    ##
    data22 <- data.frame(trt = trt, study = study,
                         r = outcome, n = n,
                         design = design)
    ##
    if (!is.null(bias)) {
      bias <- as.character(bias)
      bias <- setchar(bias, c("low", "high", "unclear"),
                      text = paste("must contain values",
                                   '"low", "high", or "unclear"',
                                   "(study-level dataset)"))
      data22$bias <- bias
    }
    if (!is.null(bias.covariate))
      data22$x.bias <- bias.covariate
    if (!is.null(bias.group))
      data22$bias.group <- bias.group
    ##
    if (!is.null(unfav)) {
      if (!is.numeric(unfav) | any(!unfav %in% 0:1))
        stop("Values of argument 'unfav' must be either 0 or 1 ",
             "(study-level dataset).")
      data22$unfav <- unfav
    }
    ##
    cov1 <- catch("cov1", mc, std.data, sfsp)
    cov2 <- catch("cov2", mc, std.data, sfsp)
    cov3 <- catch("cov3", mc, std.data, sfsp)
    ##
    if (!is.null(cov1))
      data22$x1 <- cov1
    if (!is.null(cov2))
      data22$x2 <- cov2
    if (!is.null(cov3))
      data22$x3 <- cov3
    ##
    data22$study <- paste0(data22$study, ".ad")
    ##
    ## Delete NAs
    ##
    nam <- c("trt", "study", "outcome", "n", "design", "bias", "unfav", "bias.group")
    nam <- nam[nam %in% names(data22)]
    excl2 <- apply(is.na(data22[, nam]), 1, anyNA)
    if (any(excl2))
      data22 <- data22[!excl2, ]
    ##
    ## check unfav: unique value 0 per study (repeated for the same treatment)
    ##
    if (!is.null(data22$unfav)){
      chk.unfav2 <- data22 %>%
        group_by(study) %>%
        group_map(~length(unique(subset(.x, unfav == 0, select = c(trt)))) != 1) %>%
        unlist()
      ##
      if (any(chk.unfav2))
        stop("For 'unfav' variable in std.data, each study could be provided by only a unique 0 for a specific treatment")
    }
    ##
    ## check unique bias per study
    ##
    if (!is.null(data22$bias)){
      chk.bias2 <- data22 %>%
        group_by(study) %>%
        group_map(~length(unique(.x$bias)) !=1 ) %>%
        unlist()
      ##
      if (any(chk.bias2))
        stop("The 'bias' should be a vector of length 2 where the first element is the name of the variable in prt.data and the second for the std.data")
    }
  }
  else{
    data22 <- NULL
  }
  
  
  ## Messages for missing values
  ##
  if (any(excl1))
    message('Participants with missing data in one of these variables: ',
            'outcome, trt, design, study, bias, unfav or bias.group are ',
            'discarded from the analysis')
  ##
  if (any(excl1) | any(excl2))
    message('Arms with missing data in these variables: ',
            'outcome, n, bias, unfav or bias.group are ',
            'discarded from the analysis')
  
  
  ## jagsdata for IPD
  ##
  
  ## pull relevant fields from the data and apply naming convention
  
  ## include / exclude NRS
  if (is.null(method.bias)) {
    if (any(data11$design == 'nrs') | any(data22$design == 'nrs'))
      stop('You should specify the method to combine RCT and NRS.')
    ##
    data1 <- data11
    data2 <- data22
  }
  else {
    if (method.bias %in% c('naive','adjust1','adjust2')) {
      data1 <- data11
      data2 <- data22
    }
    else if (method.bias == 'prior') {
      data1 <- data11[data11$design != 'nrs', ]
      data2 <- data22[data22$design != 'nrs', ]
      data1.nrs <- data11[data11$design == 'nrs', ]
      data2.nrs <- data22[data22$design == 'nrs', ]
    }
  }
  
  
  cov.ref <- NULL
  ##
  ## Set reference covariate values if missing
  if (isCol(data1, "x1")) {
    if (missing(cov1.ref)) {
      if (is.numeric(data1$x1) & is.numeric(data2$x1))
        cov1.ref <- min(c(min(data1$x1, na.rm = TRUE),
                          min(data2$x1, na.rm = TRUE)))
      else
        cov1.ref <- NA
    }
    else {
      if (length(cov1.ref) != 1)
        stop("Argument 'cov1.ref' must be of length 1.")
      ##
      if (!(is.numeric(data1$x1) & is.numeric(data2$x1)) & !is.na(cov1.ref)) {
        warning("Argument 'cov1.ref' set to NA as first covariate ",
                "is not continuous.")
        cov1.ref <- NA
      }
    }
    cov.ref <- cov1.ref
    ##
    if (isCol(data1, "x2")) {
      if (missing(cov2.ref)) {
        if (is.numeric(data1$x2) & is.numeric(data2$x2))
          cov2.ref <- min(c(min(data1$x2, na.rm = TRUE),
                            min(data2$x2, na.rm = TRUE)))
        else
          cov2.ref <- NA
      }
      else {
        if (length(cov2.ref) != 1)
          stop("Argument 'cov2.ref' must be of length 1.")
        ##
        if (!(is.numeric(data1$x2) & is.numeric(data2$x2)) & !is.na(cov2.ref)) {
          warning("Argument 'cov2.ref' set to NA as first covariate ",
                  "is not continuous.")
          cov2.ref <- NA
        }
      }
      cov.ref <- c(cov.ref, cov2.ref)
      ##
      if (isCol(data1, "x3")) {
        if (missing(cov3.ref)) {
          if (is.numeric(data1$x3) & is.numeric(data2$x3))
            cov3.ref <- min(c(min(data1$x3, na.rm = TRUE),
                              min(data2$x3, na.rm = TRUE)))
          else
            cov3.ref <- NA
        }
        else {
          if (length(cov3.ref) != 1)
            stop("Argument 'cov3.ref' must be of length 1.")
          ##
          if (!(is.numeric(data1$x3) & is.numeric(data2$x3)) &
              !is.na(cov3.ref)) {
            warning("Argument 'cov3.ref' set to NA as first covariate ",
                    "is not continuous.")
            cov3.ref <- NA
          }
        }
      } 
      cov.ref <- c(cov.ref, cov3.ref)
    }
  }
  ##
  cov1.labels <- NULL
  cov2.labels <- NULL
  cov3.labels <- NULL
  
  
  ## set a trt key from the two datasets
  ##
  trts <- sort(as.character(unique(c(data1$trt, data2$trt))))
  if (is.null(reference))
    reference <- trts[1]
  else
    reference <- setref(reference, trts)
  ##
  trt.key <- data.frame(trt.ini = c(reference, trts[trts != reference]),
                        trt.jags = seq_along(trts))
  
  
  ## set a study key from the two datasets
  ##
  study.key <- data.frame(std.id = c(unique(data1$study), unique(data2$study)))
  study.key$study.jags <- seq_len(nrow(study.key))
  
  
  if (!is.null(prt.data)){
    #Trt mapping
    #add treatment mapping to data
    data1 %<>% mutate(trt.jags=mapvalues(as.character(trt),
                                         from=trt.key$trt.ini,
                                         to=trt.key$trt.jags,
                                         warn_missing = FALSE) %>% as.integer)
    # Study mapping
    #add study mapping to data
    data1 %<>% mutate(study.jags=mapvalues(study,
                                           from=study.key$std.id,
                                           to=study.key$study.jags,
                                           warn_missing = FALSE) %>%  as.integer)

    # create bias_index or x.bias based on RoB and study type RCT or NRS when  method.bias= 'adjust1' or 'adjust2'
    if (!is.null(bias)){
      data1%<>%mutate(bias_index=case_when(
        design=='rct'&bias=='high'~ 1,
        design=='rct'&bias=='low'~ 2,
        design=='nrs'&bias=='high'~ 3,
        design=='nrs'&bias=='low'~ 4,
        bias=='unclear'~ 5
      ))
      bias_index.ipd<- data1 %>% group_by(study,bias_index) %>% group_keys() %>% select('bias_index')
      if (!is.null(bias.covariate)){
        # continuous
        if (is.numeric(data1$x.bias)) {
          # mean covariate if continuous
          xbias.ipd <- data1 %>%
            group_by(study.jags) %>%
            group_map(~mean(.x$x.bias, na.rm = TRUE)) %>% unlist()
          # Factor with 2 levels
        } else if (is.factor(data1$x.bias) || is.character(data1$x.bias)) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data1$x.bias)) > 2)
            stop("crossnma does not currently support bias-regression with categorical variables that have more than two levels.")
          if (length(unique(data1$x.bias)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data1$x.bias))
            data1$x.bias <- as.factor(data1$x.bias)
          data1$x.bias <- as.numeric(data1$x.bias != levels(data1$x.bias)[1])
          data1 <- data1 %>%
            group_by(study.jags) %>%
            dplyr::mutate(x.bias=mean(x.bias, na.rm = TRUE))
        } else {stop("Invalid datatype for bias covariate.")}

        # bias_index.ipd <- NULL
      }else{
        xbias.ipd <- NULL
      }

    }else{
      bias_index.ipd <- NULL
      xbias.ipd <- NULL
    }
    # pre-process the covariate if specified
    if (!is.null(cov1)) {
      if (is.numeric(data1$x1)) {
        ## mean covariate
        data1 <- data1 %>%
          group_by(study.jags) %>%
          dplyr::mutate(xm1.ipd = mean(x1, na.rm = TRUE))
        ## Center the covariate and the mean covariate
        data1$x1 <- data1$x1 - cov1.ref
        data1$xm1.ipd <- data1$xm1.ipd - cov1.ref
      }
      else if (is.factor(data1$x1) || is.character(data1$x1)) {
        ## check that covariate has fewer than 3 levels and convert
        ## strings and factors to binary covariates
        if (length(unique(data1$x1)) > 2)
          stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
        if (length(unique(data1$x1)) == 1)
          stop("Covariate should have more than one unique value.")
        if (is.character(data1$x1))
        data1$x1f <- as.factor(data1$x1) # represent the covariate as a factor
        data1$x1 <- as.numeric(data1$x1f != levels(data1$x1f)[1]) # tranfer it to numeric to be used in JAGS
        data1 <- data1 %>%
          group_by(study.jags) %>%
          dplyr::mutate(xm1.ipd=mean(x1, na.rm = TRUE))
        cov1.labels <- data1 %>% group_by(x1f,x1) %>% group_keys()
        data1$x1f <- NULL # no need for the factor version of x1
      }
      else
        stop("Invalid datatype for first covariate.")
      ## covariate2
      if (!is.null(data1$x2)){
        if (is.numeric(data1$x2)) {
          ## mean covariate if continuous
          data1 <- data1 %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm2.ipd=mean(x2, na.rm = TRUE))
          ## Center the covariate and the mean covariate
          data1$x2 <- data1$x2 - cov2.ref
          data1$xm2.ipd <- data1$xm2.ipd - cov2.ref
        }
        else if (is.factor(data1$x2) || is.character(data1$x2)) {
          ## check that covariate has fewer than 3 levels and convert
          ## strings and factors to binary covariates
          if (length(unique(data1$x2)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if (length(unique(data1$x2)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data1$x2))
            data1$x2f <- as.factor(data1$x2) # represent the covariate as a factor
          data1$x2 <- as.numeric(data1$x2f != levels(data1$x2f)[1]) # tranfer it to numeric to be used in JAGS
          data1 <- data1 %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm2.ipd=mean(x2, na.rm = TRUE))
          cov2.labels <- data1 %>% group_by(x2f,x2) %>% group_keys()
          data1$x2f <- NULL # no need for the factor version of x2
        }
        else
          stop("Invalid datatype for second covariate.")
      }
      else {
        xm2.ipd <- NULL
      }
      # covariate3
      if (!is.null(data1$x3)){
        if (is.numeric(data1$x3)) {
          data1 <- data1 %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm3.ipd=mean(x3, na.rm = TRUE))
          ## Center the covariate and the mean covariate
          data1$x3 <- data1$x3 - cov3.ref
          data1$xm3.ipd <- data1$xm3.ipd - cov3.ref
        }
        else if (is.factor(data1$x3) || is.character(data1$x3)) {
          ## check that covariate has fewer than 3 levels and convert
          ## strings and factors to binary covariates
          if (length(unique(data1$x3)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if (length(unique(data1$x3)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data1$x3))
            data1$x3f <- as.factor(data1$x3) # represent the covariate as a factor
          data1$x3 <- as.numeric(data1$x3f != levels(data1$x3f)[1]) # tranfer it to numeric to be used in JAGS
          data1 <- data1 %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm3.ipd=mean(x3, na.rm = TRUE))
          cov3.labels <- data1 %>% group_by(x3f,x3) %>% group_keys()
          data1$x3f <- NULL # no need for the factor version of x3
        }
        else
          stop("Invalid datatype for third covariate.")
      }
      else {
        xm3.ipd <- NULL
      }
    }
    else{
      xm1.ipd <- NULL
      xm2.ipd <- NULL
      xm3.ipd <- NULL
    }
    
    
    # create a matrix of treatment per study row
    jagsdata1 <- list()

    # create the matrix of trt index following the values of unfav column (adjust 1&2)
    if (method.bias%in%c("adjust1","adjust2")) {
      if (is.null(bias.group)){data1$bias.group <- 1} # Default, make bias adjustment when bias.group is not provided

      # From the unfav column create new ref treatment per study
      dd0 <- data1 %>%
        group_by(study.jags) %>%
        dplyr::mutate(ref.trt.std=.data[["trt"]][unfav==0][1])
      # For each study, arrange treatments by the new ref
      ns <- length(unique(dd0$study.jags))
      dd1 <-sapply(1:ns,
                   function(i){
                     dstd0 <- dd0[dd0$study.jags==unique(dd0$study.jags)[i],]
                     dstd  <- dstd0 %>% arrange(match(trt,ref.trt.std))
                   }
                   ,simplify = FALSE)
      dd2 <- do.call(rbind,dd1)
      # create a matrix with the treatment index
      jagsdata1$t.ipd <- dd2 %>%
        select(trt.jags,study.jags) %>% unique() %>%
        group_by(study.jags) %>%
        dplyr::mutate(arm = row_number()) %>% ungroup() %>%
        spread(arm, trt.jags) %>%
        select(-study.jags) %>%
        as.matrix()

      # generate JAGS data object
      jagstemp <- data1 %>% select(-c(study,trt,design,bias.group,unfav,bias_index,bias))
      for (v in names(jagstemp)){
        jagsdata1[[v]] <- jagstemp %>% pull(v)
      }
    }else{
      jagsdata1$t.ipd <- data1 %>% group_by(study.jags,trt.jags) %>% group_keys() %>%
        group_by(study.jags) %>%
        dplyr::mutate(arm = row_number()) %>% ungroup() %>%
        spread(arm, trt.jags) %>%
        select(-study.jags) %>%
        as.matrix()
      # generate JAGS data object
      jagstemp <- data1 %>% select(-c(study,trt,design))
      for (v in names(jagstemp)){
        jagsdata1[[v]] <- jagstemp %>% pull(v)
      }
    }
    # add number of treatments, studies, and arms to JAGS data object
    jagsdata1$nt <- trt.key %>% nrow()
    jagsdata1$ns.ipd <- ifelse(!is.null(data1),data1$study %>% unique() %>% length(),0)
    jagsdata1$na.ipd <- data1 %>% arrange(study.jags) %>% group_by(study.jags) %>% group_map(~length(unique(.x$trt))) %>%
      unlist()
    jagsdata1$np <- data1 %>% nrow()
    jagsdata1$x.bias <- NULL

    # modify names in JAGS object
    names(jagsdata1)[names(jagsdata1) == "r"]  <- "y"
    names(jagsdata1)[names(jagsdata1) == "trt.jags"] <- "trt"
    names(jagsdata1)[names(jagsdata1) == "study.jags"] <- "study"
  } else{
    jagsdata1 <- list(ns.ipd=0, nt = trt.key %>% nrow())
    xbias.ipd <- NULL
    bias_index.ipd <- NULL
  }
  #=======================
  # AD
  if (!is.null(std.data)){
    if (!is.numeric(data2$n)) stop("Sample size must be an integer greater than 0.")
    ifelse(floor(data2$n) != data2$n | data2$n<1, stop("Sample size must be an integer greater than 0."), 1)

    if (!is.numeric(data2$r)) stop("Outcome must be numeric.")

    #add treatment mapping to data
    data2 %<>% mutate(trt.jags=mapvalues(trt,
                                         from=trt.key$trt.ini,
                                         to=trt.key$trt.jags,
                                         warn_missing = FALSE) %>% as.integer)
    #add study mapping to data
    data2 %<>% mutate(study.jags=mapvalues(study,
                                           from=study.key$std.id,
                                           to=study.key$study.jags,
                                           warn_missing = FALSE) %>% as.integer)

    # add bias_index based on RoB and study design RCT or NRS. when method.bias ='adjust1' or 'adjust2'
    if (!is.null(bias)){
      data2%<>%mutate(bias_index=case_when(
        design=='rct'&bias=='high'~ 1,
        design=='rct'&bias=='low'~ 2,
        design=='nrs'&bias=='high'~ 3,
        design=='nrs'&bias=='low'~ 4,
        bias=='unclear'~ 5
      ))
      bias_index.ad<- data2 %>% arrange(study.jags) %>% group_by(study.jags,bias_index) %>% group_keys() %>% select('bias_index')
      if (!is.null(bias.covariate)){
        # continuous
        if (is.numeric(data2$x.bias)) {
          # mean covariate if continuous
          xbias.ad <- data2 %>%
            arrange(study.jags) %>%
            group_by(study.jags) %>%
            group_map(~mean(.x$x.bias, na.rm = TRUE)) %>% unlist()
          # Factor with 2 levels
        } else if (is.factor(data2$x.bias) || is.character(data2$x.bias)) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data2$x.bias)) > 2)
            stop("crossnma does not currently support bias-regression with categorical variables that have more than two levels.")
          if (length(unique(data2$x.bias)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data2$x.bias))
            data2$x.bias <- as.factor(data2$x.bias)
          data2$x.bias <- as.numeric(data2$x.bias != levels(data2$x.bias)[1])
          data2 <- data2 %>%
            arrange(study.jags) %>%
            group_by(study.jags) %>%
            dplyr::mutate(x.bias=mean(x.bias, na.rm = TRUE))
        } else {stop("Invalid datatype for bias covariate.")}

      }else{
        xbias.ad <- NULL
      }

    }
    else{
      bias_index.ad <- NULL
      xbias.ad <- NULL
    }
    
    # pre-process the covariate if specified
    if (!is.null(cov1)) {
      ## continuous
      if (is.numeric(data2$x1)) {
        ## mean covariate
        data2 <- data2 %>%
          arrange(study.jags) %>%
          group_by(study.jags) %>%
          dplyr::mutate(xm1.ad=mean(x1, na.rm = TRUE))
        ## Center the mean covariate
        ## cov1.ref if specified
        data2$xm1.ad <- data2$xm1.ad - cov1.ref
      }
      else if (is.factor(data2$x1) || is.character(data2$x1)) {
        ## check that covariate has fewer than 3 levels and convert
        ## strings and factors to binary covariates
        if (length(unique(data2$x1)) > 2)
          stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
        if (length(unique(data2$x1)) == 1)
          stop("Covariate should have more than one unique value.")
        if (is.character(data2$x1))
          data2$x1 <- as.factor(data2$x1)
        data2$x1 <- as.numeric(data2$x1 != levels(data2$x1)[1])
        data2 <- data2 %>%
          arrange(study.jags) %>%
          group_by(study.jags) %>%
          dplyr::mutate(xm1.ad=mean(x1, na.rm = TRUE))
      }
      else
        stop("Invalid datatype for covariate.")

      # covariate2
      if (!is.null(data2$x2)){
        # continuous
        if (is.numeric(data2$x2)) {
          # mean covariate
          data2 <- data2 %>%
            arrange(study.jags) %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm2.ad=mean(x2, na.rm = TRUE))
          ## Center the covariate
          data2$xm2.ad <- data2$xm2.ad - cov2.ref
        }
        else if (is.factor(data2$x2) || is.character(data2$x2)) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data2$x2)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if (length(unique(data2$x2)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data2$x2))
            data2$x2 <- as.factor(data2$x2)
          data2$x2 <- as.numeric(data2$x2 != levels(data2$x2)[1])
          data2 <- data2 %>%
            arrange(study.jags) %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm2.ad=mean(x2, na.rm = TRUE))
        }
        else
          stop("Invalid datatype for second covariate.")
      }
      else {
        xm2.ad <- NULL
      }
      # covariate3
      if (!is.null(data2$x3)){
        # continuous
        if (is.numeric(data2$x3)) {
          # mean covariate
          data2 <- data2 %>%
            arrange(study.jags) %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm3.ad=mean(x3, na.rm = TRUE))
          ## Center the mean covariate
          data2$xm3.ad <- data2$xm3.ad - cov3.ref
        }
        else if (is.factor(data2$x3) || is.character(data2$x3)) {
          ## check that covariate has fewer than 3 levels and convert
          ## strings and factors to binary covariates
          if (length(unique(data2$x3)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if (length(unique(data2$x3)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data2$x3))
            data2$x3 <- as.factor(data2$x3)
          data2$x3 <- as.numeric(data2$x3 != levels(data2$x3)[1])
          data2 <- data2 %>%
            arrange(study.jags) %>%
            group_by(study.jags) %>%
            dplyr::mutate(xm3.ad=mean(x3, na.rm = TRUE))
        }
        else
          stop("Invalid datatype for third covariate.")
      }
      else {
        xm3.ad <- NULL
      }
    }
    else{
      xm1.ad <- NULL
      xm2.ad <- NULL
      xm3.ad <- NULL
    }


    # generate JAGS data object

    # create the matrix of trt index following the values of unfav column (adjust 1&2)
    if (method.bias%in%c("adjust1","adjust2")) {
      if (is.null(bias.group)){data2$bias.group <- 1} # Default, make bias adjustment when bias.group is no provided

      # From the unfav column create new ref treatment per study
      dd0 <- data2 %>%
        arrange(study.jags) %>%
        group_by(study.jags) %>%
        dplyr::mutate(ref.trt.std=.data[["trt"]][unfav==0])
      # For each study, arrange treatments by the new ref
      ns <- length(unique(dd0$study.jags))
      dd1 <-sapply(1:ns,
                   function(i){
                     dstd0 <- dd0[dd0$study.jags==unique(dd0$study.jags)[i],]
                     dstd<- dstd0 %>% arrange(match(trt,ref.trt.std))
                   }
                   ,simplify = FALSE)
      dd2 <- do.call(rbind,dd1)


      # create a matrix with the treatment index
      jagstemp2 <- dd2 %>% arrange(study.jags) %>% group_by(study.jags) %>% dplyr::mutate(arm = row_number()) %>%
        ungroup() %>% dplyr::select(-c(trt,design,bias,ref.trt.std,unfav,bias.group,bias_index,study))  %>% gather("variable", "value",-study.jags, -arm) %>% spread(arm, value)
      jagsdata2 <- list()
      for (v in unique(jagstemp2$variable)){
        jagsdata2[[v]] <- as.matrix(jagstemp2 %>% filter(variable == v) %>% select(-study.jags, -variable))
      }
    }else{
      jagstemp2 <- data2 %>% arrange(study.jags,trt.jags) %>% group_by(study.jags) %>% dplyr::mutate(arm = row_number()) %>%
        ungroup() %>% select(-c(trt,design,bias,study))  %>% gather("variable", "value",-study.jags, -arm) %>% spread(arm, value)

      jagsdata2 <- list()
      for (v in unique(jagstemp2$variable)){
        jagsdata2[[v]] <- as.matrix(jagstemp2 %>% filter(variable == v) %>% select(-study.jags, -variable))
      }
    }

    # add number of treatments, studies, and arms to JAGS data object
    jagsdata2$ns.ad <- ifelse(!is.null(data2),data2$study.jags %>% unique() %>% length(),0)
    jagsdata2$na.ad <- data2 %>% arrange(study.jags) %>% group_by(study.jags) %>% dplyr::summarize(n.arms = n()) %>%
      ungroup() %>% select(n.arms) %>% t() %>% as.vector
    # add covariate
    jagsdata2$x1 <- jagsdata2$x2 <- jagsdata2$x3 <- NULL
    jagsdata2$xm1.ad <- unique(data2$xm1.ad)
    jagsdata2$xm2.ad <- unique(data2$xm2.ad)
    jagsdata2$xm3.ad <- unique(data2$xm3.ad)
    jagsdata2$x.bias <- NULL

    # change names
    names(jagsdata2)[names(jagsdata2) == "trt.jags"] <- "t.ad"
  }else{
    jagsdata2 <- list(ns.ad=0)
    xbias.ad <- NULL
    bias_index.ad <- NULL
  }
  # combine jagsdata of IPD and AD
  jagsdata <- c(jagsdata1,jagsdata2)

  # combine bias_index and bias covariate from IPD and AD
  jagsdata$bias_index <- c(bias_index.ipd$bias_index,bias_index.ad$bias_index)
  jagsdata$xbias <- c(xbias.ipd,xbias.ad)

  # when method.bias is adjust1 or adjust 2: add studies index:
  # 1. studies need bias adjustment and has inactive treatment (bias.group=1)
  # 2. studies need bias adjustment but has only active treatment (bias.group=2)
  # 3. studies don't need any bias adjustment
  bmat <- rbind(ifelse(!is.null(data1),list(data1 %>%
                                              arrange(study.jags) %>%
                                              group_by(study.jags) %>%
                                              select(bias.group) %>% unique() %>% select(bias.group)), list(NULL))[[1]],
                ifelse(!is.null(data2),list(data2 %>%
                                              arrange(study.jags) %>%
                                              group_by(study.jags) %>% select(bias.group) %>% unique()), list(NULL))[[1]]
  )
  jagsdata$std.in <-bmat$study.jags[bmat$bias.group==1]
  jagsdata$std.act.no <-bmat$study.jags[bmat$bias.group==0]
  jagsdata$std.act.yes <-bmat$study.jags[bmat$bias.group==2]


  #====================================
  # use NRS as prior, jags need to be run for only NRS

  if (method.bias == 'prior') {
    ## data NRS
    trts.nrs <- sort(as.character(unique(c(data1.nrs$trt, data2.nrs$trt))))
    if (is.null(reference))
      reference <- trts.nrs[1]
    else
      reference <-
        setref(reference, trts.nrs, varname = "reference",
               error.text =
                 paste("Reference treatment should be present in the list of",
                       "treatments in NRS."))
    ##
    trt.key.nrs <-
      data.frame(trt.ini = c(reference, trts.nrs[trts.nrs != reference]),
                 trt.jags = seq_along(trts.nrs))
    
    
    ## set a study key from the two datasets
    study.key.nrs <-
      data.frame(std.id = c(unique(data1.nrs$study), unique(data2.nrs$study)))
    study.key.nrs$study.jags <- seq_len(nrow(study.key.nrs))
    
    
    #====================================
    # 1. IPD-NRS
    if (!is.null(data1.nrs)){
      if (!(reference %in% data1.nrs$trt)) stop("Reference treatment is not present in the list of treatments.")

      #Trt mapping
      data1.nrs %<>% mutate(trt.jags=mapvalues(trt,
                                               from=trt.key.nrs$trt.ini,
                                               to=trt.key.nrs$trt.jags,
                                               warn_missing = FALSE) %>% as.integer) # check as.character

      #add study mapping to data
      data1.nrs %<>% mutate(study.jags=mapvalues(study,
                                                 from=study.key.nrs$std.id,
                                                 to=study.key.nrs$study.jags,
                                                 warn_missing = FALSE) %>% as.integer) # check as.character

      # add a matrix of treatment per study row
      jagsdata1.nrs <- list()

      jagsdata1.nrs$t.ipd <- data1.nrs %>%
        arrange(study.jags) %>%
        group_by(study.jags,trt.jags) %>% group_keys() %>%
        group_by(study.jags) %>%
        dplyr::mutate(arm = row_number()) %>% ungroup() %>%
        spread(arm, trt.jags) %>%
        select(-study.jags) %>%
        as.matrix()

      #generate JAGS data object

      jagstemp.nrs1 <- data1.nrs %>% select(-c(study,trt,design))
      for (v in names(jagstemp.nrs1)){
        jagsdata1.nrs[[v]] <- jagstemp.nrs1 %>% pull(v) # %>% as.vector() # %>% select(-trial, -variable)
      }
      #modify BUGS object for the various family/link combinations
      names(jagsdata1.nrs)[names(jagsdata1.nrs) == "r"]  <- "y"
      names(jagsdata1.nrs)[names(jagsdata1.nrs) == "trt.jags"] <- "trt"
      names(jagsdata1.nrs)[names(jagsdata1.nrs) == "study.jags"] <- "study"

      #add number of treatments, studies, and arms to BUGS data object
      jagsdata1.nrs$nt <- trt.key.nrs %>% nrow()
      jagsdata1.nrs$ns.ipd <- data1.nrs$study %>% unique() %>% length()
      jagsdata1.nrs$na.ipd <- data1.nrs %>% arrange(study.jags) %>% group_by(study.jags) %>% group_map(~length(unique(.x$trt))) %>%
        unlist()
      jagsdata1.nrs$np <- data1.nrs %>% nrow()
    } else {
      jagsdata1.nrs <- list(ns.ipd=0, nt = trt.key.nrs %>% nrow())
    }
    #====================================
    # AD - NRS
    if (!is.null(data2.nrs)){
      jagsdata2.nrs <- list()
      #add treatment mapping to data
      data2.nrs %<>% mutate(trt.jags=mapvalues(trt,
                                               from=trt.key.nrs$trt.ini,
                                               to=trt.key.nrs$trt.jags,
                                               warn_missing = FALSE) %>% as.integer)
      #add study mapping to data
      data2.nrs %<>% mutate(study.jags=mapvalues(study,
                                                 from=study.key.nrs$std.id,
                                                 to=study.key.nrs$study.jags,
                                                 warn_missing = FALSE) %>% as.integer)


      jagstemp.nrs2 <- data2.nrs %>% arrange(study.jags) %>% group_by(study.jags) %>% dplyr::mutate(arm = row_number()) %>%
        ungroup() %>% select(-c(trt,design,study))  %>% gather("variable", "value", -study.jags, -arm) %>% spread(arm, value)

      for (v in unique(jagstemp.nrs2$variable)){
        jagsdata2.nrs[[v]] <- as.matrix(jagstemp.nrs2 %>% filter(variable == v) %>% select(-study.jags, -variable))
      }
      names(jagsdata2.nrs)[names(jagsdata2.nrs) == "trt.jags"] <- "t.ad"

      #add number of treatments, studies, and arms to JAGS data object
      jagsdata2.nrs$ns.ad <- data2.nrs$study %>% unique() %>% length()
      jagsdata2.nrs$na.ad <- data2.nrs %>% arrange(study.jags) %>% group_by(study.jags) %>% dplyr::summarize(n.arms = n()) %>%
        ungroup() %>% select(n.arms) %>% t() %>% as.vector
    } else{
      jagsdata2.nrs <- list(ns.ad=0)
    }

    # combine jagsdata of IPD and AD
    jagsdata.nrs <- c(jagsdata1.nrs,jagsdata2.nrs)
    # set default values to the list in run.nrs
    var.infl.nrs <- ifelse(is.null(run.nrs[["var.infl"]]),1,run.nrs[["var.infl"]])
    mean.shift.nrs <- ifelse(is.null(run.nrs[["mean.shift"]]),0,run.nrs[["mean.shift"]])
    trt.effect.nrs <- ifelse(is.null(run.nrs[["trt.effect"]]),"common",run.nrs[["trt.effect"]])
    n.chains.nrs <- ifelse(is.null(run.nrs[["n.chains"]]),2,run.nrs[["n.chains"]])
    n.adapt.nrs <- ifelse(is.null(run.nrs[["n.adapt"]]),500,run.nrs[["n.adapt"]])
    n.iter.nrs <- ifelse(is.null(run.nrs[["n.iter"]]),10000,run.nrs[["n.iter"]])
    thin.nrs <- ifelse(is.null(run.nrs[["thin"]]),1,run.nrs[["thin"]])
    n.burnin.nrs <-ifelse(is.null(run.nrs[["n.burnin"]]),4000,run.nrs[["n.burnin"]])
    # jags code NRS
    model.nrs <- crossnma.code(ipd = ifelse(nrow(data1.nrs)==0,F,T),
                               ad = ifelse(nrow(data2.nrs)==0,F,T),
                               trt.effect=trt.effect.nrs,
                               covariate=NULL,
                               split.regcoef=F,
                               reg0.effect=NULL,
                               regb.effect=NULL,
                               regw.effect=NULL,
                               bias.effect=NULL,
                               bias.type=NULL,
                               prior.tau.trt=NULL,
                               prior.tau.reg0=NULL,
                               prior.tau.regb=NULL,
                               prior.tau.regw=NULL,
                               prior.tau.gamma=NULL,
                               prior.pi.high.rct=NULL,
                               prior.pi.low.rct=NULL,
                               prior.pi.high.nrs=NULL,
                               prior.pi.low.nrs=NULL,
                               method.bias = NULL,
                               d.prior.nrs=NULL)
    # jags run NRS

    seeds <- sample(.Machine$integer.max, n.chains.nrs , replace = FALSE)
    inits <- list()
    for (i in 1:n.chains.nrs)
      inits[[i]] <- list(.RNG.seed = seeds[i], .RNG.name = "base::Mersenne-Twister")

    jagsmodel.nrs <- jags.model(textConnection(model.nrs),        #Create a connection so JAGS can access the variables
                                jagsdata.nrs,
                                n.chains=n.chains.nrs,
                                n.adapt=n.adapt.nrs,
                                inits = inits)
    if (n.burnin.nrs!=0) update(jagsmodel.nrs, n.burnin.nrs)

    jagssamples.nrs <- coda.samples(jagsmodel.nrs,
                                    variable.names="d",
                                    n.iter=n.iter.nrs,
                                    thin=thin.nrs
    )

    # # # # # # # # # # # #
    # Output: prior for d's

    # map NRS trt to RCT trt
    trt.key2 <- trt.key %>% dplyr::mutate(trt.jags.nrs=mapvalues(trt.ini,
                                                                 from=trt.key.nrs$trt.ini,
                                                                 to=trt.key.nrs$trt.jags,
                                                                 warn_missing = FALSE) %>% as.integer)
    d.nrs <- summary(jagssamples.nrs)[[1]][,'Mean']+mean.shift.nrs
    prec.nrs <- ifelse(is.null(var.infl.nrs),1,var.infl.nrs)/(summary(jagssamples.nrs)[[1]][,'SD']^2)

    d.prior.nrs <- ""
    for (i in 2:nrow(trt.key2)) {
      d.prior0 <- paste0("d[",trt.key2$trt.jags[i],
                         "]~dnorm(",
                         ifelse(is.na(d.nrs[trt.key2$trt.jags.nrs[i]])|d.nrs[trt.key2$trt.jags.nrs[i]]==0,0,d.nrs[trt.key2$trt.jags.nrs[i]]),",",
                         ifelse(is.na(prec.nrs[trt.key2$trt.jags.nrs[i]])|prec.nrs[trt.key2$trt.jags.nrs[i]]==Inf,10^-2,prec.nrs[trt.key2$trt.jags.nrs[i]]),
                         ")
                    ")
      d.prior.nrs <- paste0(d.prior.nrs,d.prior0)
    }


  }else {d.prior.nrs <- NULL}

  #===================================================
  # data to be used in netgraph.crossnma() to plot the network

  # aggregate IPD dataset
  if (!is.null(data1)&!is.null(data2)){
    prt.data.ad0 <- sapply(1:length(unique(data1$study)),
                           function(i){
                             with(data1,
                                  data.frame(
                                    study=unique(data1$study)[i],
                                    trt=unique(data1[data1$study==unique(data1$study)[i],]$trt),
                                    r=sum(data1[data1$study==unique(data1$study)[i],]$r),
                                    n =nrow(data1[data1$study==unique(data1$study)[i],]),
                                    design=unique(data1[data1$study==unique(data1$study)[i],]$design)
                                  )
                             )
                           }, simplify = F)
    prt.data.ad <- do.call(rbind,prt.data.ad0)

    # combine IPD and AD in a single dataset
    all.data.ad <- rbind(data2[c('study','trt','r','n','design')],prt.data.ad)
  } else if (!is.null(data1)){
    prt.data.ad0 <- sapply(1:length(unique(data1$study)),
                           function(i){
                             with(data1,
                                  data.frame(
                                    study=unique(data1$study)[i],
                                    trt=unique(data1[data1$study==unique(data1$study)[i],]$trt),
                                    r=sum(data1[data1$study==unique(data1$study)[i],]$r),
                                    n =nrow(data1[data1$study==unique(data1$study)[i],]),
                                    design=unique(data1[data1$study==unique(data1$study)[i],]$design)
                                  )
                             )
                           }, simplify = F)
    prt.data.ad <- do.call(rbind,prt.data.ad0)
    all.data.ad <- prt.data.ad
  } else if (!is.null(data2)){
    all.data.ad <- data2[c('study','trt','r','n','design')]
  }
  #=======================
  # jags code
  model <- crossnma.code(ipd = !is.null(prt.data),
                         ad = !is.null(std.data),
                         trt.effect=trt.effect,
                         covariate=cov1,
                         split.regcoef=split.regcoef,
                         reg0.effect=reg0.effect,
                         regb.effect=regb.effect,
                         regw.effect=regw.effect,
                         bias.effect=bias.effect,
                         bias.type=bias.type,
                         bias.covariate=bias.covariate,
                         add.std.in=ifelse(length(jagsdata$std.in)==0,FALSE, TRUE),
                         add.std.act.no=ifelse(length(jagsdata$std.act.no)==0,FALSE, TRUE),
                         add.std.act.yes=ifelse(length(jagsdata$std.act.yes)==0,FALSE, TRUE),
                         prior.tau.trt=prior[['tau.trt']],
                         prior.tau.reg0=prior[['tau.reg0']],
                         prior.tau.regb=prior[['tau.regb']],
                         prior.tau.regw=prior[['tau.regw']],
                         prior.tau.gamma=prior[['tau.gamma']],
                         v=ifelse(is.null(down.wgt),
                                  list(NULL),
                                  list(down.wgt/(1-down.wgt)))[[1]],
                         prior.pi.high.rct=prior[['pi.high.rct']],
                         prior.pi.low.rct=prior[['pi.low.rct']],
                         prior.pi.high.nrs=prior[['pi.high.nrs']],
                         prior.pi.low.nrs=prior[['pi.low.nrs']],
                         method.bias = method.bias,
                         d.prior.nrs=d.prior.nrs
  )

  crossmodel <- structure(list(jagsmodel=model,
                               data=jagsdata,
                               trt.key=trt.key,
                               study.key=study.key,
                               trt.effect=trt.effect,
                               method.bias=method.bias,
                               ##
                               covariate = covariates,
                               cov.ref = cov.ref,
                               dich.cov.labels = cov1.labels,
                               ##
                               split.regcoef=split.regcoef,
                               regb.effect=regb.effect,
                               regw.effect=regw.effect,
                               bias.effect=bias.effect,
                               bias.type=bias.type,
                               all.data.ad=all.data.ad
                               ),
                          class = "crossnmaModel")
  
  return(crossmodel)
}
