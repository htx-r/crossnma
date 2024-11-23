#' Create JAGS model and data to perform cross network meta analysis
#' or meta-regression
#'
#' @description
#' This function creates a JAGS model and the needed data for
#' cross-design and cross-format network meta-analysis or
#' meta-regression for different types of outcome
#'
#' @param trt Treatment variable in \code{prt.data} and
#'   \code{std.data}.
#' @param study Study variable in \code{prt.data} and \code{std.data}.
#' @param outcome Outcome variable in \code{prt.data} and
#'   \code{std.data}.
#' @param n Number of participants in \code{std.data}.
#' @param design Design variable in \code{prt.data} and
#'   \code{std.data}.
#' @param se Standard error variable in \code{std.data} (required only
#'   for continuous outcome when \code{sm = "MD"} or \code{"SMD"}).
#' @param cov1 Optional first covariate in \code{prt.data} and
#'   \code{std.data} to conduct network meta-regression (see Details).
#' @param cov2 Optional second covariate in \code{prt.data} and
#'   \code{std.data} to conduct network meta-regression (see Details).
#' @param cov3 Optional third covariate in \code{prt.data} and
#'   \code{std.data} to conduct network meta-regression (see Details).
#' @param bias Optional variable with information on risk of bias in
#'   \code{prt.data} and \code{std.data}. Possible values for this
#'   variable are "low", "high" or "unclear" (can be
#'   abbreviated). These values must be identical for all participants
#'   from the same study. Either this variable or bias.covariate
#'   variable should be provided when method.bias = "adjust1" or
#'   "adjust2".
#' @param unfav An optional variable in \code{prt.data} and
#'   \code{std.data} indicating the unfavored treatment in each study
#'   (should be provided when method.bias = "adjust1" or
#'   "adjust2"). The entries of this variable are either 0 (unfavored
#'   treatment) or 1 (favorable treatment or treatments). Each study
#'   should include only one 0 entry. The values need to be repeated
#'   for participants who take the same treatment.
#' @param bias.covariate An optional variable in \code{prt.data} and
#'   \code{std.data} indicate the covariate used to estimate the
#'   probability of bias. Either this variable or bias variable should
#'   be provided when method.bias = "adjust1" or "adjust2".
#' @param bias.group An optional variable in \code{prt.data} and
#'   \code{std.data} that indicates the bias effect in each study (can
#'   be provided when method.bias = "adjust1" or "adjust2"). The
#'   entries of these variables should be either 1 (study has inactive
#'   treatment and its estimate should be adjusted for bias effect), 2
#'   (study has only active treatments and its estimate should be
#'   adjusted for bias effect (different from inactive bias effect) or
#'   0 (study does not need any bias adjustment). The values need to
#'   be repeated for the participants assigned to the same
#'   treatment. Default is 1.
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
#' @param sm A character indicating the underlying summary measure.
#'   Options are: Odds Ratio "OR" (default), Risk Ratio "RR", Mean
#'   Difference "MD" or Standardised Mean Difference "SMD".
#' @param reference A character indicating the name of the reference
#'   treatment. When the reference is not specified, the first
#'   alphabetic treatment will be used as a reference in the analysis.
#' @param trt.effect A character defining the model for the
#'   study-specific treatment effects. Options are "random" (default)
#'   or "common".
#' @param level.ma The level used to calculate credible intervals for
#'   network estimates.
#' @param sucra Logical. If TRUE SUCRA (Surface Under the Cumulative
#'   Ranking) values will be calculated within JAGS.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated. This
#'   argument is required when \code{sucra} is TRUE.
#' @param cov1.value The participant covariate value of \code{cov1}
#'   for which to report the results. Must be specified for network
#'   meta-regression, \code{sucra} is TRUE and when individual
#'   participant dataset is used in the analysis. For dichotomous
#'   covariates, a character of the level (used in the data) should be
#'   indicated.
#' @param cov2.value The participant covariate value of \code{cov2}
#'   for which to report the results. Must be specified for network
#'   meta-regression, \code{sucra} is TRUE and when individual
#'   participant dataset is used in the analysis. For dichotomous
#'   covariates, a character of the level (used in the data) should be
#'   indicated.
#' @param cov3.value The participant covariate value of \code{cov3}
#'   for which to report the results. Must be specified for network
#'   meta-regression, \code{sucra} is TRUE and when individual
#'   participant dataset is used in the analysis. For dichotomous
#'   covariates, a character of the level (used in the data) should be
#'   indicated.
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
#' @param reg0.effect An optional character (can by provided when at
#'   least \code{cov1} is not NULL) indicating the relationship across
#'   studies for the prognostic effects expressed by the regression
#'   coefficient, (\eqn{\beta_0}), in a study \eqn{j}. Options are
#'   "independent" or "random". We recommend using "independent"
#'   (default).
#' @param regb.effect An optional character (can by provided when at
#'   least \code{cov1} is not NULL) indicating the relationship across
#'   treatments for the between-study regression coefficient
#'   (\eqn{\beta^B}). This parameter quantifies the treatment-mean
#'   covariate interaction.  Options are "independent", "random" or
#'   "common". Default is "random".
#' @param regw.effect An optional character (can by provided when at
#'   least \code{cov1} is not NULL) indicating the relationship across
#'   treatments for the within-study regression coefficient
#'   (\eqn{\beta^W}). This parameter quantifies the
#'   treatment-covariate interaction effect at the individual level.
#'   Options are "independent", "random" and "common". Default is
#'   "random".
#' @param split.regcoef A logical value (needed when at least
#'   \code{cov1} is not NULL). If TRUE (default) the within- and
#'   between-study coefficients will be splitted in the analysis of
#'   \code{prt.data}.  When the split.regcoef = FALSE, only a single
#'   regression coefficient will be estimated to represent both the
#'   between-studies and within-studies covariate effects. In this
#'   case, both arguments \code{regb.effect} and \code{regw.effect}
#'   need to be given the same option to model the single regression
#'   effect.
#' @param method.bias A character for defining the method to combine
#'   randomized clinical trials (RCT) and non-randomized studies
#'   (NRS).  Options are "naive" for naive or unadjusted synthesize,
#'   "prior" for using NRS evidence to construct priors for the
#'   relative treatment effects in RCTs analysis, or "adjust1" and
#'   "adjust2" to allow a bias adjustment. When only one design is
#'   available (either rct or nrs), this argument needs also to be
#'   specified to indicate whether unadjusted (naive) or bias-adjusted
#'   analysis (adjust1 or adjust2) should be applied.
#' @param bias.type An optional character defining the relationship
#'   between the bias effect and the treatment effect (required when
#'   method.bias = "adjust1"). Three options are possible: "add" to
#'   add the additive bias effect, "mult" for multiplicative bias
#'   effect and "both" includes both an additive and a multiplicative
#'   terms.
#' @param bias.effect An optional character indicating the
#'   relationship for the bias coefficients across studies.  Options
#'   are "random" or "common" (default). It can be provided when
#'   method.bias = "adjust1" or "adjust2".
#' @param down.wgt An optional numeric indicating the percent to which
#'   studies at high risk of bias will be downweighted on average. The
#'   value ranges between 0 and 1. It can be provided when method.bias
#'   = "adjust1" or "adjust2".
#' @param prior.tau.trt Optional string to specify the prior for the
#'   between-study heterogeneity in treatment effects in JAGS model
#'   (when trt.effect="random").  The default prior is constructed
#'   from the data (see Details).
#' @param prior.tau.reg0 Optional string to specify the prior for the
#'   between-study heterogeneity in prognostic effects in JAGS model
#'   (when reg0.effect="random").  The default prior is constructed
#'   from the data (see Details).
#' @param prior.tau.regb Optional string to specify the prior for the
#'   between-study heterogeneity in between-study covariate effects in
#'   JAGS model (when regb.effect="random").  The default prior is
#'   constructed from the data (see Details).
#' @param prior.tau.regw Optional string to specify the prior for the
#'   between-study heterogeneity in within-study covariate effects in
#'   JAGS model (when regw.effect="random").  The default prior is
#'   constructed from the data (see Details).
#' @param prior.tau.bias Optional string to specify the prior for the
#'   between-study heterogeneity in bias effects in JAGS model (when
#'   bias.effect="random").
#' @param prior.pi.low.rct Optional string to provide the prior for
#'   the bias probability of randomised clinical trials (RCT) with low
#'   risk of bias in JAGS model (when the method.bias = "adjust1" or
#'   "adjust2" and the variable "bias" is provided).  The default is
#'   the beta distribution "dbeta(1,10)".
#' @param prior.pi.high.rct Optional string to provide the prior for
#'   the bias probability of randomised clinical trials (RCT) with
#'   high risk of bias in JAGS model (when the method.bias = "adjust1"
#'   or "adjust2" and the variable "bias" is provided).  The default
#'   is the beta distribution "dbeta(10,1)".
#' @param prior.pi.low.nrs Optional string to provide the prior for
#'   the bias probability of non-randomised studies (NRS) with low
#'   risk of bias in JAGS model (when the method.bias = "adjust1" or
#'   "adjust2" and the variable "bias" is provided).  The default is
#'   the beta distribution "dbeta(1,30)".
#' @param prior.pi.high.nrs Optional string to provide the prior for
#'   the bias probability of non-randomised studies (NRS) with high
#'   risk of bias in JAGS model (when the method.bias = "adjust1" or
#'   "adjust2" and the variable "bias" is provided).  The default is
#'   the beta distribution "dbeta(30,1)".
#' @param run.nrs.var.infl Optional numeric controls the common
#'   inflation of the variance of NRS estimates (\eqn{w}) and its
#'   values range between 0 (NRS does not contribute at all and the
#'   prior is vague) and 1 (the NRS evidence is used at face value,
#'   default approach). This argument can be provided when the NRS
#'   used as a prior (method.bias = "prior").
#' @param run.nrs.mean.shift Optional numeric controls the bias shift
#'   (\eqn{\zeta}) to be added / subtracted from the estimated mean
#'   treatment effects (on the log-scale when \code{sm = "OR"} or
#'   \code{"RR"}) from NRS network (0 is the default).  This argument
#'   can be provided when the NRS used as a prior (\code{method.bias =
#'   "prior"}).
#' @param run.nrs.n.adapt Optional numeric specifies the number of iterations for adaptation.
#' This determines how many steps the algorithm takes to adjust its parameters before starting
#' the main sampling process. Default is 1000. This argument can be provided when the NRS used as a prior
#'   (method.bias = "prior").
#' @param run.nrs.trt.effect Optional character indicates how to
#'   combine treatment effects across NRS studies. Options are
#'   "random" or "common" (default).  This argument can be provided
#'   when the NRS used as a prior (method.bias = "prior").
#' @param run.nrs.n.iter Optional numeric specifies the number of
#'   iterations to run MCMC chains for NRS network.  Default is
#'   10000. This argument can be provided when the NRS used as a prior
#'   (method.bias = "prior").
#' @param run.nrs.n.burnin Optional numeric specifies the number of
#'   burn-in to run MCMC chains for NRS network.  Default is
#'   4000. This argument can be provided when the NRS used as a prior
#'   (method.bias = "prior").
#' @param run.nrs.thin Optional numeric specifying thinning to run
#'   MCMC chains for NRS network. Default is 1. This argument can be
#'   provided when the NRS used as a prior (method.bias = "prior").
#' @param run.nrs.n.chains Optional numeric specifies the number of
#'   chains to run MCMC chains for NRS network.  Default is 2. This
#'   argument can be provided when the NRS used as a prior
#'   (method.bias = "prior").
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param run.nrs.n.thin Deprecated argument (replaced by
#'   \code{run.nrs.thin}).
#'
#' @details
#' This function creates a JAGS model and the needed data. The JAGS
#' code is created from the internal function \code{crossnma.code}.
#'
#' Covariates provided in arguments \code{cov1}, \code{cov2} and
#' \code{cov3} can be either numeric or dichotomous (should be
#' provided as factor or character) variables. By default, no
#' covariate adjustment is applied (network meta-analysis).
#'
#' The default prior for the between-study heterogeneity parameters
#' (prior.tau.trt, prior.tau.reg0, prior.tau.regb, prior.tau.regw and
#' prior.tau.bias) is a uniform distribution over the range 0 to ML,
#' where ML is the largest maximum likelihood estimates of all
#' relative treatment effects in all studies.
#'
#' @return
#' An object of class \code{crossnma.model} containing information on
#' the JAGS model, which is a list containing the following
#' components:
#'
#' \item{model}{A long character string containing JAGS code that
#'   will be run in \code{\link[R2jags]{jags.parallel}}.}
#' \item{data}{The data to be used to run JAGS model.}
#' \item{trt.key}{A table of the treatments and its mapped integer
#'   number (as used in JAGS model).}
#' \item{study.key}{A table of the studies and its mapped integer
#'   number (as used in JAGS model).}
#' \item{trt.effect}{A character defining the model for the
#'   study-specific treatment effects.}
#' \item{method.bias}{A character for defining the method to analyse combine
#'   randomized clinical trials (RCT) or \/ and non-randomized studies
#'   (NRS).}
#' \item{covariate}{A vector of the the names of the covariates
#'   (\code{cov1, cov2 and cov3}) in prt.data and std.data used in
#'   network meta-regression.}
#' \item{cov.ref}{A vector of values of \code{cov1.ref, cov2.ref,
#'   cov3.ref} to center continuous covariates. Dichotomous covariates
#'   take NA.}
#' \item{dich.cov.labels}{A matrix with the levels of each dichotomous
#'   covariate and the corresponding assigned 0 / 1 values.}
#' \item{split.regcoef}{A logical value. If FALSE the within- and
#'   between-study regression coefficients will be considered equal.}
#' \item{regb.effect}{A character indicating the model for the
#'   between-study regression coefficients across studies.}
#' \item{regw.effect}{A character indicating the model for the
#'   within-study regression coefficients across studies.}
#' \item{bias.effect}{A character indicating the model for the bias
#'   coefficients across studies.}
#' \item{bias.type}{A character indicating the effect of bias on the
#'   treatment effect; additive ("add") or multiplicative ("mult") or
#'   both ("both").}
#' \item{all.data.ad}{A data.frame object with the prt.data (after it
#'   is aggregated) and std.data in a single dataset.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{crossnma} used to create
#'   object.}
#'
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}, Guido
#'   Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{crossnma}}, \code{\link[R2jags]{jags.parallel}}
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
#' # Print call of JAGS model
#' mod
#'
#' # Print JAGS code
#' summary(mod)
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


crossnma.model <- function(trt,
                           study,
                           outcome,
                           n,
                           design,
                           se,
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
                           sm,
                           reference = NULL,
                           trt.effect = "random",
                           level.ma = gs("level.ma"),
                           ## ---------- SUCRA score ----------
                           sucra = FALSE,
                           small.values = NULL,
                           cov1.value = NULL,
                           cov2.value = NULL,
                           cov3.value = NULL,
                           ## ---------- meta regression ----------
                           cov1.ref = NULL,
                           cov2.ref = NULL,
                           cov3.ref = NULL,
                           reg0.effect = "independent",
                           regb.effect = "random",
                           regw.effect = "random",
                           split.regcoef = TRUE,
                           ## ---------- bias adjustment ----------
                           method.bias = NULL,
                           bias.type = NULL,
                           bias.effect = "common",
                           down.wgt = NULL,
                           ## ---------- prior ----------
                           prior.tau.trt = NULL,
                           prior.tau.reg0 = NULL,
                           prior.tau.regb = NULL,
                           prior.tau.regw = NULL,
                           prior.tau.bias = NULL,
                           prior.pi.high.rct = NULL,
                           prior.pi.low.rct = NULL,
                           prior.pi.high.nrs = NULL,
                           prior.pi.low.nrs = NULL,
                           ## ---------- when method.bias = "prior" ----------
                           run.nrs.var.infl = 1,
                           run.nrs.mean.shift = 0,
                           run.nrs.trt.effect = "common",
                           run.nrs.n.adapt = 1000,
                           run.nrs.n.iter = 10000,
                           run.nrs.n.burnin = 4000,
                           run.nrs.thin = 1,
                           run.nrs.n.chains = 2,
                           ##
                           backtransf = gs("backtransf"),
                           ##
                           run.nrs.n.thin = NULL
                           ) {

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
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  missing.se <- missing(se)
  ##
  if (missing(sm)) {
    if (missing.se)
      sm <- "OR"
    else
      sm <- "MD"
  }
  else
    sm <- setchar(sm, c("OR", "RR", "MD", "SMD"))
  ##
  if (missing.se & sm %in% c("MD", "SMD"))
    stop("Argument 'se' must be provided for continuous outcome (sm = \"",
         sm, ")\".")
  if (!missing.se & sm %in% c("OR", "RR"))
    warning("Argument 'se' ignored for binary outcome (sm = \"",
            sm, ")\".")
  ##
  trt.effect <- setchar(trt.effect, c("common", "random"))
  chklevel(level.ma)
  ##
  cov1.prt <- cov2.prt <- cov3.prt <-
    cov1.std <- cov2.std <- cov3.std <- NULL
  ##
  if (!is.null(prt.data)) {
    cov1.prt <- catch("cov1", mc, prt.data, sfsp)
    cov2.prt <- catch("cov2", mc, prt.data, sfsp)
    cov3.prt <- catch("cov3", mc, prt.data, sfsp)
  }
  ##
  if (!is.null(std.data)) {
    cov1.std <- catch("cov1", mc, std.data, sfsp)
    cov2.std <- catch("cov2", mc, std.data, sfsp)
    cov3.std <- catch("cov3", mc, std.data, sfsp)
  }
  ##
  avail.cov1 <- !is.null(cov1.prt) | !is.null(cov1.std)
  avail.cov2 <- !is.null(cov2.prt) | !is.null(cov2.std)
  avail.cov3 <- !is.null(cov3.prt) | !is.null(cov3.std)
  ##
  chklogical(sucra)
  ##
  if (sucra & !is.null(small.values))
    small.values <- setsv(small.values)
  else
    small.values <- NULL
  ##
  if (sucra & is.null(small.values))
    stop("Argument 'small.values' must be provided if sucra = TRUE.")
  ##
  if (sucra & avail.cov1 & is.null(cov1.value))
    stop("Argument 'cov1.value' must be provided if ",
         "sucra = TRUE and 'cov1' is specified.")
  ##
  if (sucra & avail.cov2 & is.null(cov1.value))
    stop("Argument 'cov2.value' must be provided if ",
         "sucra = TRUE and 'cov2' is specified.")
  ##
  if (sucra & avail.cov3 & is.null(cov1.value))
    stop("Argument 'cov3.value' must be provided if ",
         "sucra = TRUE and 'cov3' is specified.")
  ##
  reg0.effect <- setchar(reg0.effect, c("independent", "random"))
  regb.effect <- setchar(regb.effect, c("common", "independent", "random"))
  regw.effect <- setchar(regw.effect, c("common", "independent", "random"))
  ##
  chklogical(split.regcoef)
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
    if (!method.bias %in% c("adjust1", "adjust2"))
      stop("'down.wgt' should be specified only if ",
           "method.bias = 'adjust1' or 'adjust2'")
  }
  ##
  if (trt.effect == "common" & !is.null(prior.tau.trt))
    warning("The prior of the heterogeneity between ",
            "relative treatments parameters is ignored")
  if (reg0.effect == "common" & !is.null(prior.tau.reg0))
    warning("The prior of the heterogeneity between ",
            "progonostic parameters is ignored")
  if (regw.effect == "common" & !is.null(prior.tau.regw))
    warning("The prior of the heterogeneity between ",
            "within-study interaction parameters is ignored")
  if (regb.effect == "common" & !is.null(prior.tau.regb))
    warning("The prior of the heterogeneity between ",
            "between-study interaction parameters is ignored")
  if (bias.effect == "common" & !is.null(prior.tau.bias))
    warning("The prior of the heterogeneity between ",
            "bias effect parameters is ignored")
  if (!is.null(down.wgt) & !is.null(prior.tau.bias))
    message("The assigned prior for the heterogeneity parameters of ",
            "bias effect is ignored when 'down.wgt' is provided")
  ##
  chklogical(backtransf)

  ## Bind variables to function
  ##
  study.jags <- trt.ini <- trt.jags <-
    arm <- value <- variable <- bias_index <-
      x.bias <- x1 <- x1f <- x2 <- x2f <- x3 <- x3f <-
        ref.trt.std <- n.arms <-
          prec <- index <- num <- den <- na <- . <-
            events <- nonevents <- NULL

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
  excl1 <- FALSE
  ##
  if (!is.null(prt.data)) {
    trt <- catch("trt", mc, prt.data, sfsp)
    if (is.factor(trt))
      trt <- as.character(trt)
    ##
    study <- catch("study", mc, prt.data, sfsp)
    if (is.factor(study))
      study <- as.character(study)
    ##
    outcome <- catch("outcome", mc, prt.data, sfsp)

    if (sm %in% c("OR", "RR") &&
        (!is.numeric(outcome) | any(!outcome %in% 0:1)))
      stop("Binary outcome values must be either 0 or 1 (IPD dataset).")
    ##
    design <- catch("design", mc, prt.data, sfsp)
    design <- as.character(design)
    design <-
      setchar(design, c("nrs", "rct"),
              text = paste("must contain values \"nrs\" or \"rct\"",
                           "(IPD dataset)"))
    ##
    bias <- catch("bias", mc, prt.data, sfsp)
    bias.covariate <- catch("bias.covariate", mc, prt.data, sfsp)
    ##
    if (!is.null(method.bias) && method.bias %in% c("adjust1", "adjust2") &&
        is.null(bias) && is.null(bias.covariate))
      stop("Argument 'bias' or 'bias.covariate' must be provided if ",
           "method.bias = \"adjust1\" or \"adjust2\" (IPD dataset).")
    ##
    bias.group <- catch("bias.group", mc, prt.data, sfsp)
    ##
    unfav <- catch("unfav", mc, prt.data, sfsp)
    ##
    if (!is.null(method.bias) && method.bias %in% c("adjust1", "adjust2") &&
        is.null(unfav))
      stop("Argument 'unfav' must be provided if ",
           "method.bias = \"adjust1\" or \"adjust2\" (IPD dataset).")
    ##
    data11 <-
      data.frame(trt = trt, study = study, outcome = outcome, design = design,
                 stringsAsFactors = FALSE)
    ##
    if (!is.null(bias))
      data11$bias <-
        setchar(as.character(bias),
                c("low", "high", "unclear"),
                text = paste("must contain values \"low\", \"high\", or",
                             "\"unclear\" (IPD dataset)"))
    ##
    if (!is.null(bias.covariate))
      data11$x.bias <- bias.covariate
    ##
    if (!is.null(bias.group))
      data11$bias.group <- bias.group
    ##
    if (!is.null(unfav)) {
      if (!is.numeric(unfav) | any(!unfav %in% 0:1))
        stop("Values of argument 'unfav' must be either 0 or 1 (IPD dataset).")
      data11$unfav <- unfav
    }
    ##
    if (!is.null(cov1.prt))
      data11$x1 <- cov1.prt
    if (!is.null(cov2.prt))
      data11$x2 <- cov2.prt
    if (!is.null(cov3.prt))
      data11$x3 <- cov3.prt
    ##
    ## Exclude studies without or all events from network meta-analysis
    ##
    if (sm %in% c("RR", "OR")) {
      data11 %<>%
        group_by(study) %>%
        mutate(n = length(outcome),
               events = sum(outcome),
               nonevents = sum(n - outcome)) %>%
        as.data.frame()
      ##
      if (any(data11$events == 0)) {
        study00 <- data11 %>% filter(events == 0) %>%
          select(study) %>% unique() %>% as.character()
        ##
        if (length(study00) == 1)
          warning("Study '", study00,
                  "' without any events excluded from network meta-analysis.",
                  call. = FALSE)
        else if (length(study00) > 1)
          warning("Studies without any events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", study00, "'"),
                        collapse = " - "),
                  call. = FALSE)
        ##
        data11 %<>% filter(study != study00) %>% as.data.frame()
      }
      ##
      if (any(data11$nonevents == 0)) {
        study11 <- data11 %>% filter(nonevents == 0) %>%
          select(study) %>% unique() %>% as.character()
        ##
        if (length(study11) == 1)
          warning("Study '", study11,
                  "' with all events excluded from network meta-analysis.",
                  call. = FALSE)
        else if (length(study11) > 1)
          warning("Studies with all events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", study11, "'"),
                        collapse = " - "),
                  call. = FALSE)
        ##
        data11 %<>% filter(study != study11)
      }
      ##
      data11 %<>% select(-n) %>% select(-events) %>% select(-nonevents) %>%
        as.data.frame()
    }
    ##
    data11$study <- paste0(data11$study, ".ipd")
    ##
    ## Delete missing values
    ##
    nam <-
      c("trt", "study", "outcome", "design", "bias", "unfav", "bias.group")
    ##
    nam <- nam[nam %in% names(data11)]
    excl1 <- apply(data11[, nam], 1, anyNA)
    if (any(excl1))
      data11 <- data11[!excl1, ]
    ##
    ## Check unfav: unique value 0 per study (repeated for the same
    ## treatment)
    ##
    if (isCol(data11, "unfav")) {
      chk.unfav1 <- data11 %>%
        group_by(study) %>%
        group_map(~ length(unique(subset(.x, unfav == 0,
                                         select = c(trt)))) != 1) %>%
        unlist()
      if (any(chk.unfav1))
        stop("For 'unfav' variable in prt.data, each study could be ",
             "provided by only a unique 0 for a specific treatment")
    }
    ##
    ## Check unique bias per study
    ##
    if (isCol(data11, "bias")) {
      chk.bias1 <- data11 %>%
        group_by(study) %>%
        group_map(~length(unique(.x$bias)) !=1 ) %>%
        unlist()
      ##
      if (any(chk.bias1))
        stop("The 'bias' should be a vector of length 2 where the ",
             "first element is the name of the variable in prt.data and the ",
             "second for the std.data")
    }
  }
  else
    data11 <- NULL


  ## Prepare AD dataset
  ##
  excl2 <- FALSE
  ##
  if (!is.null(std.data)) {

    if (missing(outcome))
      stop("Mandatory argument 'outcome' missing.")

    trt <- catch("trt", mc, std.data, sfsp)
    if (is.factor(trt))
      trt <- as.character(trt)
    ##
    study <- catch("study", mc, std.data, sfsp)
    if (is.factor(study))
      study <- as.character(study)
    ##
    outcome <- catch("outcome", mc, std.data, sfsp)
    if (sm %in% c("OR", "RR") &&
        (!is.numeric(outcome) |
         any(outcome < 0) | any(!(outcome %% 1 == 0))))
      stop("Binary outcome values must be integers greater than or ",
           "equal to 0 (study-level dataset).")
    ##
    n <- catch("n", mc, std.data, sfsp)
    if (!is.numeric(n) | any(n <= 0) | any(!(n %% 1 == 0)))
      stop("Sample sizes must be integers greater than 0 ",
           "(study-level dataset).")
    if (sm %in% c("OR", "RR") & any(outcome > n))
      stop("Sample sizes must be larger equal than number of events ",
           "(study-level dataset).")
    ##
    design <- catch("design", mc, std.data, sfsp)
    design <- as.character(design)
    design <-
      setchar(design, c("nrs", "rct"),
              text = paste("must contain values \"nrs\" or \"rct\"",
                           "(study-level dataset)"))
    ##
    bias <- catch("bias", mc, std.data, sfsp)

    bias.covariate <- catch("bias.covariate", mc, std.data, sfsp)
    if (is.null(bias) && is.null(bias.covariate) && !is.null(method.bias) &&
        method.bias %in% c("adjust1", "adjust2"))
      stop("Argument 'bias' or 'bias.covariate' must be provided if ",
           "method.bias = \"adjust1\" or \"adjust2\" (study-level dataset).")
    ##
    bias.group <- catch("bias.group", mc, std.data, sfsp)
    unfav <- catch("unfav", mc, std.data, sfsp)
    if (is.null(unfav) && !is.null(method.bias) &&
        method.bias %in% c("adjust1", "adjust2"))
      stop("Argument 'unfav' must be provided if ",
           "method.bias = \"adjust1\" or \"adjust2\" (study-level dataset).")

    ##
    data22 <- data.frame(trt = trt, study = study,
                         outcome = outcome, n = n,
                         design = design,
                         stringsAsFactors = FALSE)


    ## ** se (standard error) needed for continuous outcome
    ##
    if (sm %in% c("MD", "SMD")) {
      se <- catch("se", mc, std.data, sfsp)
      data22$prec.delta.ad <- se^-2
    }

    ##
    if (!is.null(bias)) {
      data22$bias <-
        setchar(as.character(bias),
                c("low", "high", "unclear"),
                text = paste("must contain values \"low\", \"high\", or",
                             "\"unclear\" (study-level dataset)"))
    }
    ##
    if (!is.null(bias.covariate))
      data22$x.bias <- bias.covariate
    ##
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
    if (!is.null(cov1.std))
      data22$x1 <- cov1.std
    if (!is.null(cov2.std))
      data22$x2 <- cov2.std
    if (!is.null(cov3.std))
      data22$x3 <- cov3.std
    ##
    ## Exclude studies without or all events from network meta-analysis
    ##
    if (sm %in% c("RR", "OR")) {
      data22 %<>%
        group_by(study) %>%
        mutate(events = sum(outcome),
               nonevents = sum(n - outcome)) %>%
        as.data.frame()
      ##
      if (any(data22$events == 0)) {
        study00 <- data22 %>% filter(events == 0) %>%
          select(study) %>% unique() %>% as.character()
        ##
        if (length(study00) == 1)
          warning("Study '", study00,
                  "' without any events excluded from network meta-analysis.",
                  call. = FALSE)
        else if (length(study00) > 1)
          warning("Studies without any events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", study00, "'"),
                        collapse = " - "),
                  call. = FALSE)
        ##
        data22 %<>% filter(study != study00) %>% as.data.frame()
      }
      ##
      if (any(data22$nonevents == 0)) {
        study11 <- data22 %>% filter(nonevents == 0) %>%
          select(study) %>% unique() %>% as.character()
        ##
        if (length(study11) == 1)
          warning("Study '", study11,
                  "' with all events excluded from network meta-analysis.",
                  call. = FALSE)
        else if (length(study11) > 1)
          warning("Studies with all events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", study11, "'"),
                        collapse = " - "),
                  call. = FALSE)
        ##
        data22 %<>% filter(study != study11)
      }
      ##
      data22 %<>% select(-events) %>% select(-nonevents) %>%
        as.data.frame()
    }
    ##
    data22$study <- paste0(data22$study, ".ad")
    ##
    ## Delete missing values
    ##
    nam <- c("trt", "study", "outcome", "n", "design", "bias", "unfav",
             "bias.group", "prec.delta.ad")
    ##
    nam <- nam[nam %in% names(data22)]
    excl2 <- apply(data22[, nam], 1, anyNA)
    if (any(excl2))
      data22 <- data22[!excl2, ]
    ##
    ## Check unfav: unique value 0 per study (repeated for the same
    ## treatment)
    ##
    if (isCol(data22, "unfav")) {
      chk.unfav2 <- data22 %>%
        group_by(study) %>%
        group_map(~ length(unique(subset(.x, unfav == 0,
                                         select = c(trt)))) != 1) %>%
        unlist()
      ##
      if (any(chk.unfav2))
        stop("For 'unfav' variable in std.data, each study could be provided ",
             "by only a unique 0 for a specific treatment")
    }
    ##
    ## Check unique bias per study
    ##
    if (isCol(data22, "bias")) {
      chk.bias2 <- data22 %>%
        group_by(study) %>%
        group_map(~length(unique(.x$bias)) !=1 ) %>%
        unlist()
      ##
      if (any(chk.bias2))
        stop("The 'bias' should be a vector of length 2 where the ",
             "first element is the name of the variable in prt.data and the ",
             "second for the std.data")
    }
  }
  else
    data22 <- NULL


  ## Messages for missing values
  ##
  if (any(excl1))
    message("Participants with missing data in one of these variables: ",
            "outcome, trt, design, study, bias, unfav or bias.group are ",
            "discarded from the analysis")
  ##
  if (any(excl1) | any(excl2))
    message("Arms with missing data in these variables: ",
            "outcome, n, bias, unfav or bias.group are ",
            "discarded from the analysis")


  ## jagsdata for IPD
  ##

  ## Pull relevant fields from the data and apply naming convention

  ## Include / exclude NRS
  ##
  if (is.null(method.bias)) {
    if (any(data11$design == "nrs") | any(data22$design == "nrs"))
      stop("You should specify the method to combine RCT and NRS.")
    ##
    data1 <- data11
    data2 <- data22
  }
  else {
    if (method.bias %in% c("naive", "adjust1", "adjust2")) {
      data1 <- data11
      data2 <- data22
    }
    else if (method.bias == "prior") {
      data1 <- data11[data11$design != "nrs", ]
      data2 <- data22[data22$design != "nrs", ]
      data1.nrs <- data11[data11$design == "nrs", ]
      data2.nrs <- data22[data22$design == "nrs", ]
    }
  }

  cov.ref <- NULL
  ##
  ## Set reference covariate values if missing
  if (!is.null(prt.data) & !is.null(std.data)) {
    ## IPD and AD are provided
    if (isCol(data1, "x1") & isCol(data2, "x1")) {
      if (missing(cov1.ref)) {
        if (is.numeric(data1$x1) & is.numeric(data2$x1) &
            !(all(data2$x1 < 1) & all(data2$x1 > 0)))
          cov1.ref <- min(c(min(data1$x1, na.rm = TRUE),
                            min(data2$x1, na.rm = TRUE)))
        else
          cov1.ref <- NA
      }
      else {
        if (length(cov1.ref) != 1)
          stop("Argument 'cov1.ref' must be of length 1.")
        ##
        if (!(is.numeric(data1$x1) & is.numeric(data2$x1)) &
            !is.na(cov1.ref)) {
          warning("Argument 'cov1.ref' set to NA as first covariate ",
                  "is not continuous.")
          cov1.ref <- NA
        }
      }
      cov.ref <- cov1.ref
      ##
      if (isCol(data1, "x2")&isCol(data2, "x2")) {
        if (missing(cov2.ref)) {
          if (is.numeric(data1$x2) & is.numeric(data2$x2) &
              !(all(data2$x2 < 1) & all(data2$x2 > 0)))
            cov2.ref <- min(c(min(data1$x2, na.rm = TRUE),
                              min(data2$x2, na.rm = TRUE)))
          else
            cov2.ref <- NA
        }
        else {
          if (length(cov2.ref) != 1)
            stop("Argument 'cov2.ref' must be of length 1.")
          ##
          if (!(is.numeric(data1$x2) & is.numeric(data2$x2)) &
              !is.na(cov2.ref)) {
            warning("Argument 'cov2.ref' set to NA as first covariate ",
                    "is not continuous.")
            cov2.ref <- NA
          }
        }
        cov.ref <- c(cov.ref, cov2.ref)
        ##
        if (isCol(data1, "x3") & isCol(data2, "x3")) {
          if (missing(cov3.ref)) {
            if (is.numeric(data1$x3) & is.numeric(data2$x3) &
                !(all(data2$x3 < 1) &
                  all(data2$x3 > 0)))
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
  }
  else if (is.null(prt.data) | missing(prt.data) & !is.null(std.data)) {
    ## only ADs
    if (isCol(data2, "x1")) {
      if (missing(cov1.ref)) {
        if (is.numeric(data2$x1)&!(all(data2$x1 < 1)&all(data2$x1 > 0)))
          cov1.ref <- min(data2$x1, na.rm = TRUE)
        else
          cov1.ref <- NA
      }
      else {
        if (length(cov1.ref) != 1)
          stop("Argument 'cov1.ref' must be of length 1.")
        ##
        if (!is.numeric(data2$x1) |
            (all(data2$x1 < 1) & all(data2$x1 > 0)) &
            !is.na(cov1.ref)) {
          warning("Argument 'cov1.ref' set to NA as first covariate ",
                  "is not continuous.")
          cov1.ref <- NA
        }
      }
      cov.ref <- cov1.ref
      ##
      if (isCol(data2, "x2")) {
        if (missing(cov2.ref)) {
          if (is.numeric(data2$x2) &
              !(all(data2$x2 < 1) & all(data2$x2 > 0)))
            cov2.ref <- min(data2$x2, na.rm = TRUE)
          else
            cov2.ref <- NA
        }
        else {
          if (length(cov2.ref) != 1)
            stop("Argument 'cov2.ref' must be of length 1.")
          ##
          if (!is.numeric(data2$x2) |
              (all(data2$x2 < 1) & all(data2$x2 > 0)) &
              !is.na(cov2.ref)) {
            warning("Argument 'cov2.ref' set to NA as first covariate ",
                    "is not continuous.")
            cov2.ref <- NA
          }
        }
        cov.ref <- c(cov.ref, cov2.ref)
        ##
        if (isCol(data2, "x3")) {
          if (missing(cov3.ref)) {
            if (is.numeric(data2$x3) &
                !(all(data2$x3 < 1) & all(data2$x3 > 0)))
              cov3.ref <- min(data2$x3, na.rm = TRUE)
            else
              cov3.ref <- NA
          }
          else {
            if (length(cov3.ref) != 1)
              stop("Argument 'cov3.ref' must be of length 1.")
            ##
            if (!is.numeric(data2$x3) |
                (all(data2$x3 < 1) & all(data2$x3 > 0)) &
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
  }
  else if (!is.null(prt.data) & is.null(std.data) | missing(std.data)) {
    ## only IPDs
    if (isCol(data1, "x1")) {
      if (missing(cov1.ref)) {
        if (is.numeric(data1$x1))
          cov1.ref <- min(data1$x1, na.rm = TRUE)
        else
          cov1.ref <- NA
      }
      else {
        if (length(cov1.ref) != 1)
          stop("Argument 'cov1.ref' must be of length 1.")
        ##

        if (!is.numeric(data1$x1) & !is.na(cov1.ref)) {
          warning("Argument 'cov1.ref' set to NA as first covariate ",
                  "is not continuous.")
          cov1.ref <- NA
        }
      }
      cov.ref <- cov1.ref
      ##
      if (isCol(data1, "x2")) {
        if (missing(cov2.ref)) {
          if (is.numeric(data1$x2))
            cov2.ref <- min(data1$x2, na.rm = TRUE)
          else
            cov2.ref <- NA
        }
        else {
          if (length(cov2.ref) != 1)
            stop("Argument 'cov2.ref' must be of length 1.")
          ##
          if (!is.numeric(data1$x2) & !is.na(cov2.ref)) {
            warning("Argument 'cov2.ref' set to NA as first covariate ",
                    "is not continuous.")
            cov2.ref <- NA
          }
        }
        cov.ref <- c(cov.ref, cov2.ref)
        ##
        if (isCol(data1, "x3")) {
          if (missing(cov3.ref)) {
            if (is.numeric(data1$x3))
              cov3.ref <- min(data1$x3, na.rm = TRUE)
            else
              cov3.ref <- NA
          }
          else {
            if (length(cov3.ref) != 1)
              stop("Argument 'cov3.ref' must be of length 1.")
            ##
            if (!is.numeric(data1$x3) & !is.na(cov3.ref)) {
              warning("Argument 'cov3.ref' set to NA as first covariate ",
                      "is not continuous.")
              cov3.ref <- NA
            }
          }
        }
        cov.ref <- c(cov.ref, cov3.ref)
      }
    }
  }
  else{
    stop("Either the individual participant dataset (prt.data) or ",
         "the study-level dataset (std.data) should be provided.")
  }
  ##
  cov1.labels <- NULL
  cov2.labels <- NULL
  cov3.labels <- NULL


  ## Set a trt key from the two datasets
  ##
  trts <- sort(as.character(unique(c(as.character(data1$trt),
                                     as.character(data2$trt)))))
  if (is.null(reference))
    reference <- trts[1]
  else
    reference <- setref(reference, trts)
  ##
  trt.key <- data.frame(trt.ini = c(reference, trts[trts != reference]),
                        trt.jags = seq_along(trts),
                        stringsAsFactors = FALSE)


  ## Set a study key from the two datasets
  ##
  study.key <-
    data.frame(std.id = c(unique(data1$study), unique(data2$study)),
               stringsAsFactors = FALSE)
  study.key$study.jags <- seq_len(nrow(study.key))

  error.metacat <-
    paste("crossnma does not currently support meta-regression with",
          "categorical variables that have more than two levels.")
  ##
  error.biascat <-
    paste("crossnma does not currently support bias-regression with",
          "categorical variables that have more than two levels.")


  ##
  ## Individual participant data
  ##
  bias_index.ipd <- NULL
  xbias.ipd <- NULL
  ##
  xm1.ipd <- NULL
  xm2.ipd <- NULL
  xm3.ipd <- NULL
  ##
  if (!is.null(prt.data)) {
    ##
    ## Add treatment and study mapping to data
    ##
    data1 <- addmapvars(data1, trt.key, study.key)
    ##
    ## Add bias_index or x.bias based on RoB and study design RCT or
    ## NRS when method.bias = "adjust1" or "adjust2"
    ##
    data1 <- addbiasvars(data1, txt = error.biascat)
    xbias.ipd <- attr(data1, "x.bias")
    bias_index.ipd <- attr(data1, "bias_index")
    ##
    ## Pre-process first covariate if available
    ##
    if (isCol(data1, "x1")) {
      data1 <- addmeancov(data1, "x1", cov1.ref, txt = error.metacat)
      ##
      data1$x1 <- data1$mytempvar
      data1$mytempvar <- NULL
      ##
      data1$xm1.ipd <- xm1.ipd <- attr(data1, "cov.mean")
      ##
      cov1.labels <- attr(data1, "cov.labels")
      ##
      ## Second covariate
      ##
      if (isCol(data1, "x2")) {
        data1 <- addmeancov(data1, "x2", cov2.ref, txt = error.metacat)
        ##
        data1$x2 <- data1$mytempvar
        data1$mytempvar <- NULL
        ##
        data1$xm2.ipd <- xm2.ipd <- attr(data1, "cov.mean")
        ##
        cov2.labels <- attr(data1, "cov.labels")
      }
      ##
      ## Third covariate
      ##
      if (isCol(data1, "x3")) {
        data1 <- addmeancov(data1, "x3", cov3.ref, txt = error.metacat)
        ##
        data1$x3 <- data1$mytempvar
        data1$mytempvar <- NULL
        ##
        data1$xm3.ipd <- xm3.ipd <- attr(data1, "cov.mean")
        ##
        cov3.labels <- attr(data1, "cov.labels")
      }
      ##
      attr(data1, "cov.mean") <- NULL
      attr(data1, "cov.labels") <- NULL
    }
    
    
    ## Create a matrix of treatment per study row
    jagsdata1 <- list()

    ## Create the matrix of trt index following the values of unfav
    ## column (adjust 1 & 2)
    if (method.bias %in% c("adjust1", "adjust2")) {

      ## Default, make bias adjustment when bias.group is not provided
      if (is.null(bias.group))
        data1$bias.group <- 1

      ## From the unfav column create new ref treatment per study
      suppressMessages(
        dd0 <- data1 %>%
          group_by(study.jags) %>%
          mutate(ref.trt.std = .data[["trt"]][unfav == 0][1]))
      ## For each study, arrange treatments by the new ref
      ns <- length(unique(dd0$study.jags))
      dd1 <-
        sapply(1:ns,
               function(i) {
                 dstd0 <- dd0[dd0$study.jags == unique(dd0$study.jags)[i], ]
                 dstd <- dstd0 %>%
                   arrange(match(trt, ref.trt.std))
               },
               simplify = FALSE)
      dd2 <- do.call(rbind, dd1)
      ## Create a matrix with the treatment index
      suppressMessages(
        jagsdata1$t.ipd <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          select(trt.jags, study.jags) %>%
          unique() %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          spread(arm, trt.jags) %>%
          select(-study.jags) %>%
          as.matrix())

      ## Generate JAGS data object
      if (!is.null(bias)) {
        ## bias not needed and bias_index added later on study-level
        jagstemp <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          select(-c(study, trt, design, bias.group, unfav, bias_index, bias))

      }
      else{
        ## Added later on study-level (no need on participant-level)
        jagstemp <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          select(-c(study, trt, design, bias.group, unfav, x.bias))

      }
      for (v in names(jagstemp)) {
        jagsdata1[[v]] <- jagstemp %>%
          pull(v)
      }
      ## Additionally for continuous outcome, when sm = "MD" or "SMD",
      if (sm %in% c("MD", "SMD")) {
        ## 1. compute sd
        sd_jk <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          do(prec = 1 / var(.$outcome)) %>%
          unnest(cols = prec)

        ##  2. Represent "sd" column as a matrix with dim: study X
        ##  treatment arm
        jagsdata1$prec.delta.ipd <- sd_jk %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          select(-c(trt.jags)) %>%
          gather("variable", "value", -study.jags, -arm) %>%
          spread(arm, value) %>%
          select(-c(study.jags, variable)) %>%
          as.matrix()
        ## 3. Create a vector with treatment index (to be used for the
        ## precision matrix)
        jagsdata1$trt.index <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(index = cumsum(!duplicated(trt.jags))) %>%
          pull(index)
      }
      if (sm == "SMD") {
        s_n_jk <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          do(sd = sd(.$outcome),
             n = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(sd, n))
        ## For continuous outcome (doesn't matter the order of
        ## treatments)
        s.pool0.ipd <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          do(num = sum(.$sd^2 * (.$n - 1)),
             den = sum(.$n),
             na = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(num, den, na))
        ##
        jagsdata1$s.pool.ipd <- with(s.pool0.ipd, sqrt(num / (den - na)))
      }
    }
    else {
      suppressMessages(
        jagsdata1$t.ipd <- data1 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          group_keys() %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          spread(arm, trt.jags) %>%
          select(-study.jags) %>%
          as.matrix())
      ## Generate JAGS data object
      jagstemp <- data1 %>%
        arrange(study.jags, trt.jags) %>%
        select(-c(study, trt, design))
      for (v in names(jagstemp)) {
        jagsdata1[[v]] <- jagstemp %>%
          pull(v)
      }
      ##! Additionally for continuous outcome, when sm="MD" or "SMD",
      if (sm %in% c("MD", "SMD")) {
        ## 1. Compute SD
        sd_jk <- data1%>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          do(prec =1 / var(.$outcome)) %>%
          unnest(cols = prec)
        ##  Represent "sd" column as a matrix with dim: study X
        ##  treatment arm
        jagsdata1$prec.delta.ipd <- sd_jk %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          select(-c(trt.jags)) %>%
          gather("variable", "value", -study.jags, -arm) %>%
          spread(arm, value) %>%
          select(-c(study.jags, variable)) %>%
          as.matrix()
        ## 3. Create a vector with treatment index (to be used for the
        ## precision matrix)
        jagsdata1$trt.index <- data1 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(index = cumsum(!duplicated(trt.jags))) %>%
          pull(index)
      }
      if (sm == "SMD") {
        s_n_jk <- data1 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          do(sd = sd(.$outcome),
             n = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(sd, n))
        ## For continuous outcome (doesn't matter the order of
        ## treatments)
        s.pool0.ipd <- data1 %>% #!! data1 should be  replaced by s_n_jk?
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          do(num = sum(.$sd^2 * (.$n - 1)),
             den = sum(.$n),
             na = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(num, den, na))
        ##
        jagsdata1$s.pool.ipd <- with(s.pool0.ipd, sqrt(num / (den - na)))
      }
    }
    ## Add number of treatments, studies, and arms to JAGS data object
    jagsdata1$nt <- trt.key %>% nrow()
    jagsdata1$ns.ipd <-
      if (!is.null(data1))
        data1$study %>%
          unique() %>%
          length()
      else
        0
    suppressMessages(
      jagsdata1$na.ipd <-
        data1 %>%
        arrange(study.jags, trt.jags) %>%
        group_by(study.jags) %>%
        group_map(~length(unique(.x$trt))) %>%
        unlist())
    jagsdata1$np <- data1 %>%
      nrow()
    ##
    ## Add cov1.value, cov2.value, cov3.value, needed when sucra = TRUE
    ## and cov1, cov2, cov3 are provided
    if (sucra) {
      labs <- rbind(cov1.labels, cov2.labels, cov3.labels)
      jagsdata1$cov1.value <-
        if (is.numeric(cov1.value))
          cov1.value
        else
          as.numeric(labs[labs[, 1] == cov1.value, 2])
      ##
      jagsdata1$cov2.value <-
        if (is.numeric(cov2.value))
          cov2.value
        else
          as.numeric(labs[labs[, 1] == cov2.value, 2])
      ##
      jagsdata1$cov3.value <-
        if (is.numeric(cov3.value))
          cov3.value
        else
          as.numeric(labs[labs[, 1] == cov3.value, 2])
    }
    else {
      jagsdata1$cov1.value <- NULL
      jagsdata1$cov2.value <- NULL
      jagsdata1$cov3.value <- NULL
    }
    ## Modify names in JAGS object
    names(jagsdata1)[names(jagsdata1) == "outcome"]  <- "y"
    names(jagsdata1)[names(jagsdata1) == "trt.jags"] <- "trt"
    names(jagsdata1)[names(jagsdata1) == "study.jags"] <- "study"
    jagsdata1$x.bias <- NULL
  }
  else {
    jagsdata1 <- list(ns.ipd = 0, nt = trt.key %>% nrow())
    xbias.ipd <- NULL
    bias_index.ipd <- NULL
  }
  
  
  ##
  ## Aggregated data
  ##
  bias_index.ad <- NULL
  xbias.ad <- NULL
  ##
  xm1.ad <- NULL
  xm2.ad <- NULL
  xm3.ad <- NULL
  ##
  if (!is.null(std.data)) {
    if (!is.numeric(data2$n))
      stop("Sample size must be an integer greater than 0.")
    if (any(floor(data2$n) != data2$n | data2$n < 1))
      stop("Sample size must be an integer greater than 0.")
    if (!is.numeric(data2$outcome))
      stop("Outcome must be numeric.")
    if (sm %in% c("MD", "SMD")) {
      if (!is.numeric(data2$prec.delta.ad))
        stop("Standard error must be numeric.")
    }
    ##
    ## Add treatment and study mapping to data
    ##
    data2 <- addmapvars(data2, trt.key, study.key)
    ##
    ## Add bias_index based on RoB and study design RCT or NRS when
    ## method.bias = "adjust1" or "adjust2"
    ##
    data2 <- addbiasvars(data2, txt = error.biascat)
    xbias.ad <- attr(data2, "x.bias")
    bias_index.ad <- attr(data2, "bias_index")
    ##
    ## Pre-process first covariate if available
    ##
    if (isCol(data2, "x1")) {
      data2 <- addmeancov(data2, "x1", cov1.ref, ipd = FALSE)
      ##
      data2$x1 <- data2$mytempvar
      xm1.ad <- attr(data2, "cov.mean")
      ##
      data2$mytempvar <- data2$mytempvar.mean <- NULL
      ##
      ## Second covariate
      ##
      if (isCol(data2, "x2")) {
        data2 <- addmeancov(data2, "x2", cov2.ref, ipd = FALSE)
        ##
        data2$x2 <- data2$mytempvar
        xm2.ad <- attr(data2, "cov.mean")
        ##
        data2$mytempvar <- data2$mytempvar.mean <- NULL
      }
      ##
      ## Third covariate
      ##
      if (isCol(data2, "x3")) {
        data2 <- addmeancov(data2, "x3", cov3.ref, ipd = FALSE)
        ##
        data2$x3 <- data2$mytempvar
        xm3.ad <- attr(data2, "cov.mean")
        ##
        data2$mytempvar <- data2$mytempvar.mean <- NULL
      }
      ##
      attr(data2, "cov.mean") <- NULL
    }
    ## Generate JAGS data object
    ## Create the matrix of trt index following the values of unfav
    ## column (adjust 1 & 2)
    if (method.bias %in% c("adjust1", "adjust2")) {
      ## Default, make bias adjustment when bias.group is no provided
      if (is.null(bias.group))
        data2$bias.group <- 1

      ## From the unfav column create new ref treatment per study
      suppressMessages(
        dd0 <- data2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(ref.trt.std = .data[["trt"]][unfav == 0]))
      ## For each study, arrange treatments by the new ref
      ns <- length(unique(dd0$study.jags))
      dd1 <-
        sapply(1:ns,
               function(i) {
                 dstd0 <- dd0[dd0$study.jags == unique(dd0$study.jags)[i], ]
                 dstd <- dstd0 %>%
                   arrange(match(trt, ref.trt.std))
               },
               simplify = FALSE)
      dd2 <- do.call(rbind, dd1)

      ## Create a matrix with the treatment index
      if (!is.null(bias)) {
        ## bias not needed, bias_index added later as bias_index.ad as
        ## it needs to be combined with bias_index.ipd
        suppressMessages(
          jagstemp2 <- dd2 %>%
            arrange(study.jags, trt.jags) %>%
            group_by(study.jags) %>%
            mutate(arm = row_number()) %>%
            ungroup() %>%
            select(-c(trt, design, bias, ref.trt.std,
                             unfav, bias.group, bias_index, study)) %>%
            gather("variable", "value", -study.jags, -arm) %>%
            spread(arm, value))
      }
      else{
        ## x.bias added later as xbias.ad as it needs to be combined
        ## with xbias.ipd
        suppressMessages(
          jagstemp2 <- dd2 %>%
            arrange(study.jags, trt.jags) %>%
            group_by(study.jags) %>%
            mutate(arm = row_number()) %>%
            ungroup() %>%
            select(-c(trt, design, ref.trt.std, unfav,
                             bias.group, x.bias, study)) %>%
            gather("variable", "value", -study.jags, -arm) %>%
            spread(arm, value))
      }
      jagsdata2 <- list()
      ##
      for (v in unique(jagstemp2$variable)) {
        suppressMessages(
          jagsdata2[[v]] <-
            as.matrix(jagstemp2 %>%
                      filter(variable == v) %>%
                      select(-study.jags, -variable)))

      }
      ## For continuous outcome with SMD
      if (sm == "SMD") {
        s.pool0.ad <- dd2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          do(num = sum((1 / .$prec.delta.ad) *(.$n - 1)),
             ## num = sum(.$se^2 * .$n),
             den = sum(.$n),
             na = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(num, den, na))
        ##
        jagsdata2$s.pool.ad <- with(s.pool0.ad, sqrt(num / (den - na)))
      }

    }
    else {
      suppressMessages(
        jagstemp2 <- data2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          select(-c(trt, design, bias, study)) %>%
          gather("variable", "value", -study.jags, -arm) %>%
          spread(arm, value))
      ##
      jagsdata2 <- list()
      for (v in unique(jagstemp2$variable)) {
        suppressMessages(
          jagsdata2[[v]] <-
            as.matrix(jagstemp2 %>%
                      filter(variable == v) %>%
                      select(-study.jags, -variable)))
      }
      ## Calculate SMD for continuous outcome
      if (sm == "SMD") {
        s.pool0.ad <- data2 %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          do(num = sum((1 / .$prec.delta.ad) * (.$n - 1)),
             den = sum(.$n),
             na = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(num, den, na))
        ##
        jagsdata2$s.pool.ad <- with(s.pool0.ad, sqrt(num / (den - na)))
      }
    }

    ## Add number of treatments, studies, and arms to JAGS data object
    suppressMessages(
      jagsdata2$ns.ad <-
        if (!is.null(data2))
          data2$study.jags %>%
          unique() %>%
          length()
        else
          0)
    ##
    suppressMessages(
      jagsdata2$na.ad <-
        data2 %>%
        arrange(study.jags, trt.jags) %>%
        group_by(study.jags) %>%
        summarize(n.arms = n()) %>%
        ungroup() %>%
        select(n.arms) %>%
        t() %>%
        as.vector)
    ## Add covariate

    jagsdata2$x1 <- jagsdata2$x2 <- jagsdata2$x3 <- NULL
    if (isCol(data2, "x1"))
      jagsdata2$xm1.ad <- xm1.ad
    if (isCol(data2, "x2"))
      jagsdata2$xm2.ad <- xm2.ad
    if (isCol(data2, "x3"))
      jagsdata2$xm3.ad <- xm3.ad
    jagsdata2$x.bias <- NULL

    ## Change names
    if (sm %in% c("OR", "RR"))
      names(jagsdata2)[names(jagsdata2) == "outcome"] <- "r"
    names(jagsdata2)[names(jagsdata2) == "trt.jags"] <- "t.ad"
    if (sm %in% c("MD", "SMD"))
      names(jagsdata2)[names(jagsdata2) == "outcome"] <- "ybar"
  }
  else {
    jagsdata2 <- list(ns.ad = 0)
    xbias.ad <- NULL
    bias_index.ad <- NULL
  }

  ## Combine jagsdata of IPD and AD
  jagsdata <- c(jagsdata1, jagsdata2)

  ## Add mean study to calculate sucra (for betab * std.mean)
  if (sucra){
    jagsdata$stds.mean1 <-
      suppressWarnings(mean(c(xm1.ad, xm1.ipd), na.rm = TRUE))
    ##
    jagsdata$stds.mean2 <-
      suppressWarnings(mean(c(xm2.ad, xm2.ipd), na.rm = TRUE))
    ##
    jagsdata$stds.mean3 <-
      suppressWarnings(mean(c(xm3.ad, xm3.ipd), na.rm = TRUE))
    ##
    jagsdata$cov.ref <- cov.ref
  }
  else {
    jagsdata$stds.mean1 <- jagsdata$stds.mean2 <- jagsdata$stds.mean3 <- NULL
    jagsdata$cov.ref <- NULL
  }
  ## Combine bias_index and bias covariate from IPD and AD
  jagsdata$bias_index <-
    c(bias_index.ipd$bias_index, bias_index.ad$bias_index)
  jagsdata$xbias <- c(xbias.ipd, xbias.ad)
  
  
  ## when method.bias is adjust1 or adjust 2: add studies index:
  ## 1. studies need bias adjustment and has inactive treatment
  ##    (bias.group  =1)
  ## 2. studies need bias adjustment but has only active treatment
  ##    (bias.group = 2)
  ## 3. studies don't need any bias adjustment
  ##
  bmat <- rbind(
    if (!is.null(data1))
      suppressMessages(
        list(data1 %>%
             arrange(study.jags, trt.jags) %>%
             group_by(study.jags) %>%
             select(bias.group) %>%
             unique() %>%
             select(bias.group))[[1]])
    else
      NULL,
    if (!is.null(data2))
      suppressMessages(
        list(data2 %>%
             arrange(study.jags, trt.jags) %>%
             group_by(study.jags) %>%
             select(bias.group) %>%
             unique())[[1]])
    else
      NULL
  )
  ##
  if (method.bias %in% c("adjust1", "adjust2") &
      isCol(bmat, "study.jags") &
      isCol(bmat, "bias.group")) {
    jagsdata$std.in <- bmat$study.jags[bmat$bias.group == 1]
    jagsdata$std.act.no <- bmat$study.jags[bmat$bias.group == 0]
    jagsdata$std.act.yes <- bmat$study.jags[bmat$bias.group == 2]
  }

  ##
  ## Data to be used in netgraph.crossnma() to plot the network
  ##
  ## Aggregate IPD dataset
  if (!is.null(data1) & !is.null(data2)) {
    ## When IPD & AD provided

    ## prt.data.ad0 <- sapply(1:length(unique(data1$study)),
    ##                        function(i) {
    ##                          with(data1,
    ##                               data.frame(
    ##                                 study=unique(data1$study)[i],
    ##                                 trt=unique(data1[data1$study==unique(data1$study)[i],]$trt),
    ##                                 outcome = sum(data1[data1$study==unique(data1$study)[i],]$outcome),
    ##                                 n =nrow(data1[data1$study==unique(data1$study)[i],]),
    ##                                 design = unique(data1[data1$study==unique(data1$study)[i],]$design)
    ##                               )
    ##                          )
    ##                        }, simplify = F)
    ## prt.data.ad <- do.call(rbind,prt.data.ad0)


    if (sm %in% c("MD", "SMD")) {
      ## For continuous outcome
      prt.data.ad <- data1 %>%
        arrange(study, trt) %>%
        group_by(study, trt) %>%
        do(outcome = mean(.$outcome),
           n = nrow(.),
           se = sd(.$outcome),
           design = unique(.$design)) %>%
        unnest(cols = c(outcome, n, se, design)) %>%
        as.data.frame()
      ## Combine IPD and AD in a single dataset
      data2.ad <- data2
      data2.ad[, "se"] <- 1 / sqrt(data2.ad[, "prec.delta.ad"])
      all.data.ad <-
        rbind(data2.ad[c("study", "trt", "outcome", "n", "se", "design")],
              prt.data.ad)
      if (sm == "SMD") {
        ## Calculate s.pooled for SMD
        s.pooled <- all.data.ad %>%
          group_by(study) %>%
          do(s.pooled =
               sqrt(sum(.$se^2 * (.$n - 1)) /
                    (sum(.$n) - summarize(., n.arms = group_size(.)) %>%
                     pull(n.arms))) ) %>%
          unnest(cols = c(s.pooled))

        ## Add s.pooled per study to the dataset
        all.data.ad %<>%
          mutate(s.pooled =
                   mapvalues(study,
                             from = s.pooled$study,
                             to = s.pooled$s.pooled,
                             warn_missing = FALSE) )# %>% select(-se)
      }
    }
    else{
      ## IPD & AD with binary outcome
      prt.data.ad <- data1 %>%
        arrange(study, trt) %>%
        group_by(study, trt) %>%
        do(outcome = sum(.$outcome),
           n = nrow(.),
           design = unique(.$design)) %>%
        unnest(cols = c(outcome, n, design)) %>%
        as.data.frame()
      ## Combine IPD and AD in a single dataset
      all.data.ad <-
        rbind(data2[c("study", "trt", "outcome", "n", "design")],
              prt.data.ad)
    }

  }
  else if (!is.null(data1)) {
    ## Only IPD provided
    if (sm %in% c("MD", "SMD")) {
      ## Continuous outcome
      prt.data.ad <- data1 %>%
        arrange(study, trt) %>%
        group_by(study, trt) %>%
        do(outcome = mean(.$outcome),
           n = nrow(.),
           se = sd(.$outcome),
           design = unique(.$design)) %>%
        unnest(cols = c(outcome, n, se, design)) %>%
        as.data.frame()
      ##
      all.data.ad <- prt.data.ad
      ##
      if (sm == "SMD") {
        ## Calculate s.pooled for SMD
        s.pooled <- all.data.ad %>%
          group_by(study) %>%
          do(s.pooled = sqrt(sum(.$se^2 * (.$n - 1)) /
                             (sum(.$n) -
                              summarize(., n.arms = group_size(.)) %>%
                              pull(n.arms)))) %>%
          unnest(cols = c(s.pooled))
        ## Add s.pooled per study to the dataset
        all.data.ad %<>%
          mutate(s.pooled =
                   mapvalues(study,
                             from = s.pooled$study,
                             to = s.pooled$s.pooled,
                             warn_missing = FALSE)) # %>% select(-se)
      }
    }
    else{
      ## Binary outcome
      prt.data.ad <- data1 %>%
        arrange(study, trt) %>%
        group_by(study, trt) %>%
        do(outcome = sum(.$outcome),
           n = nrow(.),
           design = unique(.$design)) %>%
        unnest(cols = c(outcome, n, design)) %>%
        as.data.frame()
      ##
      all.data.ad <- prt.data.ad
    }

  }
  else if (!is.null(data2)) {
    ## Only AD provided
    if (sm %in% c("MD", "SMD")) {
      ## Continuous outcome
      data2.ad <- data2
      data2.ad[, "se"] <- 1 / sqrt(data2.ad[, "prec.delta.ad"])
      all.data.ad <-
        data2.ad[c("study", "trt", "outcome", "n", "se", "design")]
      ##
      if (sm == "SMD") {
        ## Calculate s.pooled for SMD
        s.pooled <- all.data.ad %>%
          group_by(study) %>%
          do(s.pooled = sqrt(sum(.$se^2 * (.$n - 1)) /
                             (sum(.$n) -
                              summarize(., n.arms = group_size(.)) %>%
                              pull(n.arms))) ) %>%
          unnest(cols = c(s.pooled))
        ## Add s.pooled per study to the dataset
        all.data.ad %<>%
          mutate(s.pooled =
                   mapvalues(study,
                             from = s.pooled$study,
                             to = s.pooled$s.pooled,
                             warn_missing = FALSE)) # %>% select(-se)
      }
    }
    else{
      ## Binary data
      all.data.ad <- data2[, c("study", "trt", "outcome", "n", "design")]
    }

  }
  ##
  ## Construct default priors
  ##
  max.d <- max_TE(all.data.ad, sm = sm)


  ##
  ## use NRS as prior, JAGS needs to be run for only NRS
  ##
  if (method.bias == "prior") {
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
                 trt.jags = seq_along(trts.nrs),
                 stringsAsFactors = FALSE)


    ## Set a study key from the two datasets
    study.key.nrs <-
      data.frame(std.id = c(unique(data1.nrs$study), unique(data2.nrs$study)),
                 stringsAsFactors = FALSE)
    study.key.nrs$study.jags <- seq_len(nrow(study.key.nrs))


    ##
    ## 1. IPD-NRS
    ##
    if (!is.null(data1.nrs)) {
      if (!(reference %in% data1.nrs$trt))
        stop("Reference treatment is not present in the list of treatments.")

      ## Trt mapping
      data1.nrs %<>%
        mutate(trt.jags =
                 mapvalues(trt,
                           from = trt.key.nrs$trt.ini,
                           to = trt.key.nrs$trt.jags,
                           warn_missing = FALSE) %>%
                 as.integer)
      ## Add study mapping to data
      suppressMessages(
        data1.nrs %<>%
        mutate(study.jags =
                 mapvalues(study,
                           from = study.key.nrs$std.id,
                           to = study.key.nrs$study.jags,
                           warn_missing = FALSE) %>%
                 as.integer))
      ## Add a matrix of treatment per study row
      jagsdata1.nrs <- list()

      suppressMessages(
        jagsdata1.nrs$t.ipd <- data1.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          group_keys() %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          spread(arm, trt.jags) %>%
          select(-study.jags) %>%
          as.matrix())

      ## Generate JAGS data object

      jagstemp.nrs1 <- data1.nrs %>%
        arrange(study.jags, trt.jags) %>%
        select(-c(study, trt, design))
      for (v in names(jagstemp.nrs1)) {
        jagsdata1.nrs[[v]] <- jagstemp.nrs1 %>%
          pull(v)
        ## %>% as.vector() # %>% select(-trial, -variable)
      }
      ## ! Additionally for continuous outcome
      if (sm %in% c("MD", "SMD")) {
        ## 1. compute sd
        sd_jk <- data1.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          select(outcome) %>%
          summarise_all(sd)
        ##  2. represent 'sd' column as a matrix with dim: study X treatment arm
        jagsdata1.nrs$sd <- sd_jk %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          select(-c(trt.jags)) %>%
          gather("variable", "value", -study.jags, -arm) %>%
          spread(arm, value) %>%
          select(-c(study.jags, variable)) %>%
          as.matrix()
      }
      if (sm == "SMD") {
        s_n_jk <- data1.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags, trt.jags) %>%
          do(sd = sd(.$outcome),
             n = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(sd, n))
        ## For continuous outcome (doesn't matter the order of
        ## treatments)
        s.pool0.ipd <- data1.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          do(num = sum(.$sd^2 * (.$n - 1)),
             den = sum(.$n),
             na = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(num, den, na))
        ##
        jagsdata1.nrs$s.pool.ipd <- with(s.pool0.ipd, sqrt(num / (den - na)))
      }
      ## Modify BUGS object for the various family / link combinations
      names(jagsdata1.nrs)[names(jagsdata1.nrs) == "outcome"]  <- "y"
      names(jagsdata1.nrs)[names(jagsdata1.nrs) == "trt.jags"] <- "trt"
      names(jagsdata1.nrs)[names(jagsdata1.nrs) == "study.jags"] <- "study"

      ## Add number of treatments, studies, and arms to BUGS data object
      jagsdata1.nrs$nt <- trt.key.nrs %>%
        nrow()
      jagsdata1.nrs$ns.ipd <- data1.nrs$study %>%
        unique() %>%
        length()
      suppressMessages(
        jagsdata1.nrs$na.ipd <- data1.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          group_map(~length(unique(.x$trt))) %>%
          unlist())
      jagsdata1.nrs$np <- data1.nrs %>%
        nrow()
    }
    else
      jagsdata1.nrs <- list(ns.ipd = 0, nt = trt.key.nrs %>% nrow())
    ##
    ## AD - NRS
    ##
    if (!is.null(data2.nrs)) {
      jagsdata2.nrs <- list()
      ##add treatment mapping to data
      data2.nrs %<>%
        mutate(trt.jags =
                 mapvalues(trt,
                           from = trt.key.nrs$trt.ini,
                           to = trt.key.nrs$trt.jags,
                           warn_missing = FALSE) %>%
                 as.integer)
      ##add study mapping to data
      suppressMessages(
        data2.nrs %<>%
        mutate(study.jags =
                 mapvalues(study,
                           from = study.key.nrs$std.id,
                           to = study.key.nrs$study.jags,
                           warn_missing = FALSE) %>%
                 as.integer))

      suppressMessages(
        jagstemp.nrs2 <- data2.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          mutate(arm = row_number()) %>%
          ungroup() %>%
          select(-c(trt, design, study)) %>%
          gather("variable", "value", -study.jags, -arm) %>%
          spread(arm, value))

      for (v in unique(jagstemp.nrs2$variable)) {
        suppressMessages(
          jagsdata2.nrs[[v]] <-
            as.matrix(jagstemp.nrs2 %>%
                      filter(variable == v) %>%
                      select(-study.jags, -variable)))
      }

      ## Calculate pooled standard deviation for continuous outcome
      if (sm == "SMD") {
        s.pool0.ad <- data2.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          do(num = sum((1 / .$prec.delta.ad) * (.$n - 1)),
             ## num = sum(.$se^2 * .$n),
             den = sum(.$n),
             na = summarize(., n.arms = group_size(.)) %>%
               pull(n.arms)) %>%
          unnest(cols = c(num, den, na))
        ##
        jagsdata2.nrs$s.pool.ad <- with(s.pool0.ad, sqrt(num / (den - na)))
      }
      if (sm %in% c("OR", "RR"))
        names(jagsdata2.nrs)[names(jagsdata2.nrs) == "outcome"]  <- "r"
      names(jagsdata2.nrs)[names(jagsdata2.nrs) == "trt.jags"] <- "t.ad"
      if (sm %in% c("MD", "SMD"))
        names(jagsdata2.nrs)[names(jagsdata2) == "outcome"] <- "ybar"

      ## Add number of treatments, studies, and arms to JAGS data
      ## object
      jagsdata2.nrs$ns.ad <- data2.nrs$study %>%
        unique() %>%
        length()
      suppressMessages(
        jagsdata2.nrs$na.ad <- data2.nrs %>%
          arrange(study.jags, trt.jags) %>%
          group_by(study.jags) %>%
          summarize(n.arms = n()) %>%
          ungroup() %>%
          select(n.arms) %>%
          t() %>%
          as.vector)
    }
    else
      jagsdata2.nrs <- list(ns.ad = 0)


    ## Combine jagsdata of IPD and AD
    jagsdata.nrs <- c(jagsdata1.nrs, jagsdata2.nrs)
    ## Set default values to the list in run.nrs
    var.infl.nrs <- replaceNULL(run.nrs.var.infl, 1)
    mean.shift.nrs <- replaceNULL(run.nrs.mean.shift, 0)
    trt.effect.nrs <- replaceNULL(run.nrs.trt.effect, "common")
    n.adapt.nrs <- replaceNULL(run.nrs.n.adapt, 1000)
    n.chains.nrs <- replaceNULL(run.nrs.n.chains, 2)
    ##
    n.iter.nrs <- replaceNULL(run.nrs.n.iter, 10000)
    if (!missing(run.nrs.n.thin)) {
      if (!missing(run.nrs.thin))
        warning("Deprecated argument 'run.nrs.n.thin' ignored as ",
                "argument 'run.nrs.thin' is provided.")
      else {
        warning("Argument 'run.nrs.n.thin' is deprecated; please use ",
                "argument 'run.nrs.thin' instead.")
        run.nrs.thin <- run.nrs.n.thin
      }
    }
    thin.nrs <- replaceNULL(run.nrs.n.thin, 1)
    ##
    n.burnin.nrs <- replaceNULL(run.nrs.n.burnin, 4000)


    ## JAGS code NRS
    ##
    model.nrs <- crossnma.code(ipd = !is.null(data1.nrs),
                               ad = !is.null(data2.nrs),
                               sm = sm,
                               max.d = max.d,
                               trt.effect = trt.effect.nrs,
                               ## ---------- SUCRA score ----------
                               sucra = FALSE,
                               small.values = NULL,
                               cov1.value = NULL,
                               cov2.value = NULL,
                               cov3.value = NULL,
                               cov.ref = NULL,
                               ## ---------- meta regression ----------
                               covariate = NULL,
                               split.regcoef = FALSE,
                               reg0.effect = NULL,
                               regb.effect = NULL,
                               regw.effect = NULL,
                               bias.effect = NULL,
                               bias.type = NULL,
                               prior.tau.trt = NULL,
                               prior.tau.reg0 = NULL,
                               prior.tau.regb = NULL,
                               prior.tau.regw = NULL,
                               prior.tau.gamma = NULL,
                               prior.pi.high.rct = NULL,
                               prior.pi.low.rct = NULL,
                               prior.pi.high.nrs = NULL,
                               prior.pi.low.nrs = NULL,
                               method.bias = NULL,
                               d.prior.nrs = NULL)
    ## JAGS run NRS
    ## seeds <- sample(.Machine$integer.max, n.chains.nrs)
    ## jmodel <- model.nrs
    ## jagsmodel.nrs <- jags.parallel(data = jagsdata.nrs,
    ##                          inits = NULL,
    ##                          parameters.to.save = "d",
    ##                          model.file = jmodel,
    ##                          n.chains = n.chains.nrs,
    ##                          n.iter = n.iter.nrs,
    ##                          n.burnin = n.burnin.nrs,
    ##                          thin = thin.nrs,
    ##                          jags.seed = seeds,
    ##                          DIC = FALSE)


    inits <- list()
    seeds <- sample(.Machine$integer.max, n.chains.nrs)
    for (i in 1:n.chains.nrs)
      inits[[i]] <- list(.RNG.seed = seeds[i],
                         .RNG.name = "base::Mersenne-Twister")
    ##
    jagsmodel.nrs <-
      suppressWarnings(jags.model(textConnection(model.nrs),
                                  jagsdata.nrs,
                                  n.chains = n.chains.nrs,
                                  n.adapt = n.adapt.nrs,
                                  inits = inits,
                                  quiet = TRUE))
    ##
    if (n.burnin.nrs != 0)
      suppressWarnings(update(jagsmodel.nrs, n.burnin.nrs))


    jagssamples.nrs <-
      coda.samples(jagsmodel.nrs, variable.names = "d",
                   n.iter = n.iter.nrs, thin = thin.nrs)

    ## Output: prior for d's
    ##
    ## map NRS trt to RCT trt
    ##
    trt.key2 <-
      trt.key %>%
      mutate(
        trt.jags.nrs =
          mapvalues(trt.ini,
                    from = trt.key.nrs$trt.ini,
                    to = trt.key.nrs$trt.jags,
                    warn_missing = FALSE) %>%
          as.integer)
    ## d.nrs <-
    ##   summary(as.mcmc(jagsmodel.nrs))[[1]][, "Mean"] + mean.shift.nrs
    d.nrs <- summary(jagssamples.nrs)[[1]][, "Mean"] + mean.shift.nrs

    prec.nrs <-
      if (!is.null(var.infl.nrs))
        var.infl.nrs
      else
        1
    ## prec.nrs <-
    ##   prec.nrs / (summary(as.mcmc(jagsmodel.nrs))[[1]][, "SD"]^2)
    prec.nrs <- prec.nrs / (summary(jagssamples.nrs)[[1]][, "SD"]^2)

    ##
    d.prior.nrs <- "\n"
    for (i in 2:nrow(trt.key2)) {
      d.prior0 <-
        paste0("d[", trt.key2$trt.jags[i], "] ~ dnorm(",
               ifelse(is.na(d.nrs[trt.key2$trt.jags.nrs[i]]),
                      0,
                      d.nrs[trt.key2$trt.jags.nrs[i]]),
               ", ",
               ifelse(is.na(prec.nrs[trt.key2$trt.jags.nrs[i]]) |
                      prec.nrs[trt.key2$trt.jags.nrs[i]] == Inf,
               (max.d * 15)^(-2),
               prec.nrs[trt.key2$trt.jags.nrs[i]]),
               ")",
               if (i < nrow(trt.key2))
                 "\n"
               else
                 "")
      ##
      d.prior.nrs <- paste0(d.prior.nrs, d.prior0)
    }
  }
  else
    d.prior.nrs <- NULL


  ##
  ## JAGS code
  ##
  model <- crossnma.code(ipd = !is.null(prt.data),
                         ad = !is.null(std.data),
                         sm = sm,
                         max.d = max.d,
                         trt.effect = trt.effect,
                         ## ---------- SUCRA score ----------
                         sucra = sucra,
                         small.values = small.values,
                         cov1.value = cov1.value,
                         cov2.value = cov2.value,
                         cov3.value = cov3.value,
                         cov.ref = cov.ref,
                         covariate = covariates,
                         split.regcoef = split.regcoef,
                         reg0.effect = reg0.effect,
                         regb.effect = regb.effect,
                         regw.effect = regw.effect,
                         bias.effect = bias.effect,
                         bias.type = bias.type,
                         bias.covariate = bias.covariate,
                         add.std.in =
                           method.bias %in% c("adjust1", "adjust2") &&
                           length(jagsdata$std.in) > 0,
                         add.std.act.no =
                           method.bias %in% c("adjust1", "adjust2") &&
                           length(jagsdata$std.act.no) > 0,
                         add.std.act.yes =
                           method.bias %in% c("adjust1", "adjust2") &&
                           length(jagsdata$std.act.yes) > 0,
                         prior.tau.trt = prior.tau.trt,
                         prior.tau.reg0 = prior.tau.reg0,
                         prior.tau.regb = prior.tau.regb,
                         prior.tau.regw = prior.tau.regw,
                         prior.tau.gamma = prior.tau.bias,
                         v = if (!is.null(down.wgt))
                               list(down.wgt / (1 - down.wgt))[[1]]
                             else
                               NULL,
                         prior.pi.high.rct = prior.pi.high.rct,
                         prior.pi.low.rct = prior.pi.low.rct,
                         prior.pi.high.nrs = prior.pi.high.nrs,
                         prior.pi.low.nrs = prior.pi.low.nrs,
                         method.bias = method.bias,
                         d.prior.nrs = d.prior.nrs
                         )


  res <- list(model = model,
              data = jagsdata,
              sm = sm,
              reference = reference,
              trt.key = trt.key,
              study.key = study.key,
              trt.effect = trt.effect,
              sucra=sucra,
              method.bias = method.bias,
              level.ma = level.ma,
              quantiles =
                c((1 - level.ma) / 2, 0.5, 1 - (1 - level.ma) / 2),
              backtransf = if (sm %in% c("MD", "SMD")) FALSE else backtransf,
              ##
              covariate = covariates,
              cov.ref = cov.ref,
              dich.cov.labels = rbind(cov1.labels, cov2.labels, cov3.labels),
              ##
              split.regcoef = split.regcoef,
              regb.effect = regb.effect,
              regw.effect = regw.effect,
              bias.effect = bias.effect,
              bias.type = bias.type,
              all.data.ad = all.data.ad,
              ##
              call = match.call(),
              version = packageDescription("crossnma")$Version)
  ##
  class(res) <- "crossnma.model"

  res
}
