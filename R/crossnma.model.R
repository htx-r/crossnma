# delete the message from study.jags once I run the adjust1 and adjust2 models
# delete NA from data is not applied! CHECK
#' Create JAGS model to synthesize cross-design evidence and cross-format data in NMA and NMR for dichotomous outcomes
#' @description This function creates a JAGS model and the needed data. The JAGS code is created from the internal function \code{crossnma.code}.
#'
#' @param prt.data An object of class data.frame containing the individual participant dataset. Each row contains the data of a single participant.
#' The data frame needs to have the following columns: treatment, study identification, outcome (event and non-event), design. Additional columns might be required for certain analyses.
#' @param std.data An object of class data.frame containing the study-level dataset. Each row represents the information of study arm.
#' The data frame needs to have the following columns: treatment, study identification, outcome (number of events), sample size and design. Additional columns might be required for certain analyses.
#' @param trt A charachter of the name of the treatment variable in prt.data and std.data.
#' @param study A charachter of the name of the study variable in prt.data and std.data.
#' @param outcome A charachter of the name of the outcome variable in prt.data and std.data.
#' @param n A character of the name of the number of participants variable in std.data.
#' @param design A charachter of the name of the design variable in prt.data and std.data.
#' @param reference A character indicating the name of the reference treatment. When the reference is not specified, the first alphabetic treatment will be used as a reference in the analysis.
#' @param trt.effect A character defining the model for the study-specific treatment effects. Options are 'random' (default) or 'common'.
#' @param covariate An optional vector indicating the name of the covariates in prt.data and std.data (to conduct a network meta-regression)
#' The covariates can be either numeric or dichotomous variables. The user can provide up to 3 covariates.
#' For example, we set `covariate=c(‘age’, ‘sex’)` to adjust for 2 covariates.
#' The default option is `covariate=NULL` where no covariate adjustment is applied (network meta-analysis).
#' @param reg0.effect An optional character (when \code{covariate} is not NULL) indicating the relationship across studies for the prognostic effects expressed by the regression coefficient, (\eqn{\beta_0}), in a study \eqn{j}.
#' Options are 'independent' or 'random'. We recommend using 'independent' (default).
#' @param regb.effect An optional character (when \code{covariate} is not NULL) indicating the relationship across studies for the between-study regression coefficient (\eqn{\beta^B}). This parameter quantifies the treatment-mean covariate interaction.
#'  Options are 'random' or 'common'. Default is 'random'.
#' @param regw.effect An optional character (when \code{covariate} is not NULL)  indicating the relationship across studies for the within-study regression coefficient (\eqn{\beta^W}). This parameter quantifies the treatment-covariate interaction effect at the individual level.
#' Options are 'random' and 'common'. Default is 'random'.
#' @param split.regcoef A logical value (when \code{covariate} is not NULL). If TRUE the within- and between-study coefficients will be splitted in the analysis of prt.data.
#' The default is TRUE. When the split.regcoef = FALSE, only a single regression coefficient will be estimated to represent both the between-studies and within-studies covariate effects.
#' @param method.bias A character for defining the method to combine randomised clinical trials (RCT) and non-randomised studies (NRS) (when \code{design} has both designs nrs and rct).
#' Options are 'naive' for naive synthesize, 'prior' for using NRS to inform priors for the relative treatment effects in RCTs.
#' or 'adjust1' and 'adjust2' to allow a bias adjustment. When only one design is available (either rct or nrs), this argument needs also to be specified to indicate whether unadjusted (naive) or bias-adjusted analysis (adjust1 or adjust2) should be applied.
#' @param bias A charachter indicating the name of the variable (required when method.bias='adjust1' or 'adjust2') that includes the risk of bias adjustment in prt.data and std.data.
#' The values of this variable should be a character with entries that need to be spelled as such 'low', 'high' or 'unclear'. These values need to be repeated for the participants that belong to the same study.
#' @param bias.type An optional character defining of bias on the treatment effect (required when method.bias='adjust1').
#' Three options are possible: 'add' to add the additive bias effect,'mult' for multiplicative bias effect and 'both' includes both an additive and a multiplicative terms.
#' @param bias.covariate A charachter of the variable name that will be used in estimating the probability of bias (can be provided when method.bias='adjust1' or 'adjust2')
#' @param bias.effect An optional character indicating the relationship for the bias coefficients across studies.
#' Options are 'random' or 'common' (default). It is required when method.bias='adjust1' or 'adjust2'.
#' @param unfav A charachter which defines the names of the variables in prt.data and std.data that include an indicator of the unfavoured treatment in each study (should be provided when method.bias='adjust1' or 'adjust2')
#' The entries of these variables should be either 0 (unfavoured treatment) or 1 (favourable treatment or treatments). Each study should include only one 0. The values need to be repeated for the participants that belong to the same study.
#' @param bias.group An optional charachter which defines the names of the variables in prt.data and std.data that indicates the bias effect in each study (can be provided when method.bias='adjust1' or 'adjust2')
#' The entries of these variables should be either 1 (study has inactive treatment and its estimate should be adjusted for bias effect), 2 (study has only active treatments and its estimate should be adjusted for bias effect (different from inactive bias effect)
#' or 0 (study doesn't need bias adjustment). The values need to be repeated for the participants that belong to the same study. Default is 1 which means applying bias adjustment to all study-specific estimates if the study is at high risk of bias.
#' @param prior An optional list to control the prior for various parameters in JAGS model. When effects are set as 'random', we can set the heterogeneity parameters for: tau.trt for the treatment effects,
#' tau.reg0 for the effect of prognostic covariates, tau.regb and tau.regw for within- and between-study covariate effect, respectively.
#' and tau.gamma for bias effect. The default of all heterogeneity parameters is 'dunif(0,2)'. Currently only the uniform distribution is supported.
#' When the method.bias='adjust1' or 'adjust2', the user may provide priors to control the bias probability.
#' For the bias probabilities, beta distributions are assumed with the following default values: RCT with low (pi.low.rct='dbeta(1,10)')/high (pi.high.rct='dbeta(10,1)') bias, NRS with low(pi.low.rct='dbeta(1,30)')/high (pi.high.rct='dbeta(30,1)') bias (pi.low.nrs, pi.high.nrs).
#' @param run.nrs An optional list is needed when the NRS used as a prior (method.bias='prior').
#' The list consists of the following: (\code{var.infl}) controls the common inflation of the varaince of NRS estimates (\eqn{w}) and its values range between 0 (NRS do not contribute at all and the prior is vague) and 1 (the NRS evidence is used at face value, default approach).
#' The parameter (\code{mean.shift}) is the bias shift (\eqn{\zeta}) to be added/subtracted from the estimated mean treatment effects (on the log-scale) from NRS network (0 is the default). Either (\code{var.infl}) or (\code{mean.shift}) should be provided but not both.
#' Here you can also specify the arguments to control the MCMC chains with default value is in the parentheses: the number of adaptions n.adapt (500), number of iterations n.iter(10000), number of burn in n.burnin (4000),
#' number of thinning thin (1) and number of chains n.chains (2), see \code{\link{jags.model}} arguments from rjags package.
#' @return \code{crossnma.model} returns an object of class \code{crossnmaModel} which is a list containing the following components:
#' @return \code{jagsmodel}  A long character string containing JAGS code that will be run in \code{\link{jags}}.
#' @return \code{data}  The data used in the JAGS code.
#' @return \code{trt.key}  A table of the treatments and its mapped integer number (used in JAGS model).
#' @return \code{study.key}  A table of the studies and its mapped integer number (used in JAGS model).
#' @return \code{trt.effect} A character defining the model for the study-specific treatment effects.
#' @return \code{method.bias}  A character for defining the method to combine randomised clinical trials (RCT) and non-randomised studies (NRS).
#' @return \code{covariate}  A list of the the names of the covariates in prt.data and std.data used in network meta-regression.
#' @return \code{split.regcoef} A logical value. If FALSE the within- and between-study regression coefficients will be considered equal.
#' @return \code{regb.effect} An optional character  indicating the model for the between-study regression coefficients across studies.
#' @return \code{regw.effect} An optional character indicating the model for the within-study regression coefficients across studies.
#' @return \code{bias.effect} An optional character indicating the model for the bias coefficients across studies.
#' @return \code{bias.type} A character indicating the effect of bias on the treatment effect; additive ('add') or multiplicative ('mult') or both ('both').
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
#'                   trt='trt',
#'                   study='study',
#'                   outcome=outcome',
#'                   n='n',
#'                   design='design',
#'                   reference='A',
#'                   trt.effect='random',
#'                   covariate = NULL,
#'                   method.bias='naive'
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
#' @export crossnma.model
#' @importFrom magrittr %>%
#' @importFrom magrittr "%<>%"
#' @importFrom plyr mapvalues
#' @importFrom rlang quo
#' @import dplyr
#' @import ggplot2
#' @import rjags
#' @import tidyr
#' @import coda
#' @import netmeta
#' @seealso \code{\link{crossnma.run}}, \code{\link{jags.model}}


crossnma.model <- function(prt.data,
                       std.data,
                       trt,
                       study,
                       outcome,
                       n,
                       design,
                       reference=NULL,
                       trt.effect='random',
                       #---------- meta regression ----------
                       covariate = NULL,#list(c('age','sex','EDSS'),c('age','sex','EDSS')), #default NULL
                       reg0.effect='independent',
                       regb.effect='random',
                       regw.effect='random',
                       split.regcoef=T,
                       #---------- bias adjustment ----------
                       method.bias = NULL,
                       bias=NULL, #optional, required for adjust1. It can be either: low, high,
                       bias.group = NULL,
                       unfav=NULL,
                       bias.type=NULL,#c('add','mult','both'),
                       bias.covariate=NULL,
                       bias.effect='common',
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
                                    n.adapt = 500,
                                    n.iter=10000,
                                    n.burnin = 4000,
                                    thin=1,
                                    n.chains=2)
){
  options(warn=-1)
  # Bind variables to function
  trt.ini <- NULL
  trt.jags <- NULL
  std.id <- NULL
  study.jags <- NULL
  arm <- NULL
  value <- NULL
  variable <- NULL

  #============================================
  # prepare IPD and AD
  if(!is.null(prt.data)){
    varlist1 <- c(trt = trt,
                  study = study,
                  r = outcome,
                  design=design,
                  bias=bias,
                  x.bias=bias.covariate,
                  bias.group=bias.group,
                  unfav=unfav)

    if(!is.null(covariate)){
      x11 <- covariate[1]
      x21<- ifelse(is.na(covariate[2]),list(NULL),covariate[2])[[1]]
      x31<- ifelse(is.na(covariate[3]),list(NULL),covariate[3])[[1]]
    } else {
      x11 <- NULL
      x21 <- NULL
      x31 <- NULL
    }
    covlist1 <- c(x1=x11,x2=x21,x3=x31)

    data11 <- prt.data[,c(varlist1,covlist1)]
    names(data11) <- c(names(varlist1), names(covlist1))
  } else{
    data11 <- NULL
  }
  if(!is.null(std.data)){
    varlist2 <- c(trt = trt,
                  study = study,
                  r = outcome,
                  n = n,
                  design=design,
                  bias=bias,
                  x.bias=bias.covariate,
                  bias.group=bias.group,
                  unfav=unfav)
    if(!is.null(covariate)){
      x12 <- covariate[1]
      x22<- ifelse(is.na(covariate[2]),list(NULL),covariate[2])[[1]]
      x32<- ifelse(is.na(covariate[3]),list(NULL),covariate[3])[[1]]
    }else {
      x12 <- NULL
      x22 <- NULL
      x32 <- NULL
    }
    covlist2 <- c(x1=x12,x2=x22,x3=x32)
    data22 <-std.data[, c(varlist2,covlist2)]
    names(data22) <- c(names(varlist2), names(covlist2))
  }else{
    data22 <- NULL
  }
  # if(!is.null(prt.data)){
  #   varlist1 <- c(trt = trt[1],
  #                 study = study[1],
  #                 r = outcome[1],
  #                 design=design[1],
  #                 bias=bias[1],
  #                 x.bias=bias.covariate[1],
  #                 bias.group=bias.group[1],
  #                 unfav=unfav[1])
  #
  #   if(!is.null(covariate)){
  #     x11 <- covariate[[1]][1]
  #     x21<- ifelse(is.na(covariate[[1]][2]),list(NULL),covariate[[1]][2])[[1]]
  #     x31<- ifelse(is.na(covariate[[1]][3]),list(NULL),covariate[[1]][3])[[1]]
  #   } else {
  #     x11 <- NULL
  #     x21 <- NULL
  #     x31 <- NULL
  #   }
  #   covlist1 <- c(x1=x11,x2=x21,x3=x31)
  #
  #   data11 <- prt.data[,c(varlist1,covlist1)]
  #   names(data11) <- c(names(varlist1), names(covlist1))
  # } else{
  #   data11 <- NULL
  # }
  # if(!is.null(std.data)){
  #   if(!is.null(prt.data)){
  #     varlist2 <- c(trt = trt[2],
  #                   study = study[2],
  #                   r = outcome[2],
  #                   n = n,
  #                   design=design[2],
  #                   bias=bias[2],
  #                   x.bias=bias.covariate[2],
  #                   bias.group=bias.group[2],
  #                   unfav=unfav[2])
  #     if(!is.null(covariate)){
  #       x12 <- covariate[[2]][1]
  #       x22<- ifelse(is.na(covariate[[2]][2]),list(NULL),covariate[[2]][2])[[1]]
  #       x32<- ifelse(is.na(covariate[[2]][3]),list(NULL),covariate[[2]][3])[[1]]
  #     } else {
  #       x12 <- NULL
  #       x22 <- NULL
  #       x32 <- NULL
  #     }
  #     covlist2 <- c(x1=x12,x2=x22,x3=x32)
  #     data22 <-std.data[, c(varlist2,covlist2)]
  #     names(data22) <- c(names(varlist2), names(covlist2))
  #   }else{
  #     varlist2 <- c(trt = trt[1],
  #                   study = study[1],
  #                   r = outcome[1],
  #                   n = n,
  #                   design=design[1],
  #                   bias=bias[1],
  #                   x.bias=bias.covariate[1],
  #                   bias.group=bias.group[1],
  #                   unfav=unfav[1])
  #     if(!is.null(covariate)){
  #       x12 <- covariate[[1]][1]
  #       x22<- ifelse(is.na(covariate[[1]][2]),list(NULL),covariate[[1]][2])[[1]]
  #       x32<- ifelse(is.na(covariate[[1]][3]),list(NULL),covariate[[1]][3])[[1]]
  #     } else {
  #       x12 <- NULL
  #       x22 <- NULL
  #       x32 <- NULL
  #     }
  #     covlist2 <- c(x1=x12,x2=x22,x3=x32)
  #     data22 <-std.data[, c(varlist2,covlist2)]
  #     names(data22) <- c(names(varlist2), names(covlist2))
  #   }
  # } else{
  #   data22 <- NULL
  # }

  # discard NA's
  excl1 <- is.na(data11$r) |is.na(data11$study) |is.na(data11$trt) |is.na(data11$design) | if(!is.null(data11$bias)){is.na(data11$bias)}else{FALSE} | if(!is.null(data11$unfav)){is.na(data11$unfav)}else{FALSE} | if(!is.null(data11$bias.group2)){is.na(data11$bias.group2)}else{FALSE}
  excl2 <- is.na(data22$r) | is.na(data22$n) |is.na(data22$study) |is.na(data22$trt) |is.na(data22$design) | if(!is.null(data22$bias)){is.na(data22$bias)}else{FALSE} | if(!is.null(data22$unfav)){is.na(data22$unfav)}else{FALSE} | if(!is.null(data22$bias.group2)){is.na(data22$bias.group2)}else{FALSE}
  if (sum(excl1)>0) message('Participants with missing data in one of these variables: outcome, trt, design, study, bias, unfav or bias.group are discarded from the analysis')
  if (sum(excl1)>0|sum(excl2)>0) message('Arms with missing data in these variables: outcome, n, bias, unfav or bias.group are discarded from the analysis')
  data11 <- data11[!excl1,]
  data22 <- data22[!excl2,]

  # checking ------------------------
  # 1. formatting
  # factor to character: trt & study
  if(is.factor(data11$trt)) data11$trt <- as.character(data11$trt)
  if(is.factor(data22$trt)) data22$trt <- as.character(data22$trt)

  if(is.factor(data11$study)) data11$study <- as.character(data11$study)
  if(is.factor(data22$study)) data22$study <- as.character(data22$study)

  if(!is.null(bias)) if(is.factor(data11$bias)) data11$bias <- as.character(data11$bias)
  if(!is.null(bias)) if(is.factor(data22$bias)) data22$bias <- as.character(data22$bias)

  if(is.factor(data11$design)) data11$design <- as.character(data11$design)
  if(is.factor(data22$design)) data22$design <- as.character(data22$design)


  # 3. values
  if(!is.null(data22$bias)){if(!any(unique(data22$bias)%in%c('low','high','unclear'))) stop("bias values must be either low, high or unclear")}
  if(!is.null(data11$bias)){if(!any(unique(data11$bias)%in%c('low','high','unclear'))) stop("bias values must be either low, high or unclear")}

  if(!is.null(method.bias)){if(is.null(bias)&(method.bias%in%c('adjust1','adjust2'))) stop("bias should be specified if method.bias = 'adjust1' or 'adjust2'")}
  if(!is.null(method.bias)){if(is.null(bias.type)&(method.bias%in%c('adjust1'))) stop("bias.type should be specified if method.bias = 'adjust1'")}
  if(!is.null(data11)){if(any(!unique(data11$design)%in%c('rct','nrs'))) stop("For prt.data, design must be set to either 'rct' ( randomised clinical trials) or 'nrs' ( non-randomised studies)")}
  if(!is.null(data22)){if(any(!unique(data22$design)%in%c('rct','nrs'))) stop("For std.data, design must be set to either 'rct' ( randomised clinical trials) or 'nrs' ( non-randomised studies)")}
  if(!is.null(bias.type)){if(!(bias.type%in%c('add','mult','both')))stop("bias.type need to be set as 'add' (additive), 'mult' (multiplicative) or 'both' (additive and multiplicative)")}

  #if(!(reference %in% c(data11$trt,data22$trt))) stop("Reference treatment is not available in the list of treatments.")
  if(!is.null(data22)){if(!is.numeric(data22$n)|sum(data22$n%%1)!=0|any(data22$n==0)) stop("Sample size must be an integer greater than 0.")}
  if(!is.null(data22)){if(!is.numeric(data22$r)|sum(data22$r%%1)!=0) stop("Outcome must be an integer greater than or equal 0.")}
  if(!is.null(data22)){if(sum(data22$r>data22$n)!=0) stop("Sample size must be greater than number of events.")}
  if(!is.null(data11)){if(!is.numeric(data11$r)|sum(!data11$r%in%c(0,1))!=0) stop('The values of the outcome in prt.data must be either 0 or 1 ')}

  if(!is.null(split.regcoef)){ if(!is.logical(split.regcoef)) stop("split.regcoef is a logical value TRUE/FALSE")}

  # 4. dependency
  if(trt.effect=='random'&!is.null(prior$tau.trt)) warning(" The prior of the heterogeneity between relative treatments parameters is ignored")
  if(reg0.effect=='random'&!is.null(prior$tau.reg0)) warning(" The prior of the heterogeneity between progonostic parameters is ignored")
  if(regw.effect=='random'&!is.null(prior$tau.regw)) warning(" The prior of the heterogeneity between within-study interaction parameters is ignored")
  if(regb.effect=='random'&!is.null(prior$tau.regb)) warning(" The prior of the heterogeneity between between-study interaction parameters is ignored")
  if(bias.effect=='random'&!is.null(prior$tau.gamma)) warning(" The prior of the heterogeneity between bias effect parameters is ignored")
  if(!is.null(data11)&!is.null(unfav)){if(!data11$unfav%in%c(0,1)) stop("The values of 'unfav' should be either 0 or 1")}
  if(!is.null(data22)&!is.null(unfav)){if(!data22$unfav%in%c(0,1)) stop("The values of 'unfav' should be either 0 or 1")}

  # check unfav: unique value 0 per study (repeated for the same treatment)
  if(!is.null(data11$unfav)){
    chk.unfav1 <- data11%>%
      group_by(study)%>%
      group_map(~length(unique(subset(.x,unfav==0, select=c(trt))))!=1)%>%
      unlist()
    if(any(chk.unfav1)) stop("For 'unfav' variable in prt.data, each study could be provided by only a unique 0 for a specific treatment")
  }
  if(!is.null(data22$unfav)){
    chk.unfav2 <- data22%>%
      group_by(study)%>%
      group_map(~length(unique(subset(.x,unfav==0, select=c(trt))))!=1)%>%
      unlist()

    if(any(chk.unfav2)) stop("For 'unfav' variable in std.data, each study could be provided by only a unique 0 for a specific treatment")
  }

  # check unique bias per study
  if(!is.null(data11$bias)){
    chk.bias1 <- data11%>%
      group_by(study)%>%
      group_map(~length(unique(.x$bias))!=1)%>%
      unlist()

    if(any(chk.bias1)) stop("The 'bias' should be a vector of length 2 where the first element is the name of the variable in prt.data and the second for the std.data")
  }
  if(!is.null(data22$bias)){
    chk.bias2 <- data22%>%
      group_by(study)%>%
      group_map(~length(unique(.x$bias))!=1)%>%
      unlist()

    if(any(chk.bias2)) stop("The 'bias' should be a vector of length 2 where the first element is the name of the variable in prt.data and the second for the std.data")
  }


  #====================================
  # jagsdata for IPD

  # pull relevant fields from the data and apply naming convention



  # include/exclude NRS
  if(is.null(method.bias)){
    data1 <- data11
    data2 <- data22
    if(!is.null(data11[data11$design=='nrs',])||!is.null(data22[data22$design=='nrs',])){
      stop('You should specify the method to combine RCT and NRS')
    }
  }else{
    if(method.bias%in%c('naive','adjust1','adjust2')){
      data1 <- data11
      data2 <- data22
    }else if(method.bias=='prior'){
      data1 <- data11[data11$design!='nrs',]
      data2 <- data22[data22$design!='nrs',]
      data1.nrs <- data11[data11$design=='nrs',]
      data2.nrs <- data22[data22$design=='nrs',]
    }else{
      stop("method.bias is either 'naive','prior','adjust1' or 'adjust2'")
    }
  }

  # set a trt key from the two datasets
  trt.df <- data.frame(trt=unique(c(data1$trt,data2$trt)))
  reference <- ifelse(is.null(reference),sort(as.character(trt.df$trt))[1], reference )
  if(!(reference %in% data11$trt)&!(reference %in% data22$trt)) stop("Reference treatment is not present in the list of treatments.")

  trt.key <- trt.df$trt %>% unique %>% sort %>% tibble(trt.ini=.) %>%
    filter(trt.ini!=reference) %>% add_row(trt.ini=reference, .before=1) %>%
    mutate(trt.jags = 1:dim(.)[1])
  # set a study key from the two datasets
  study.df <- data.frame(std.id= unique(c(data1$study,data2$study)))
  study.key <- study.df%>% mutate(study.jags = 1:dim(.)[1])

  if(!is.null(prt.data)){
    #Trt mapping
    #add treatment mapping to data
    data1 %<>% mutate(trt.jags=mapvalues(as.character(trt),
                                         from=trt.key$trt.ini,
                                         to=trt.key$trt.jags,
                                         warn_missing = FALSE)%>%as.integer)
    # Study mapping
    #add study mapping to data
    data1 %<>% mutate(study.jags=mapvalues(study,
                                           from=study.key$std.id,
                                           to=study.key$study.jags,
                                           warn_missing = FALSE)%>% as.integer)

    # create bias_index or x.bias based on RoB and study type RCT or NRS when  method.bias= 'adjust1' or 'adjust2'
    if(!is.null(bias)){
      data1%<>%mutate(bias_index=case_when(
        design=='rct'&bias=='high'~ 1,
        design=='rct'&bias=='low'~ 2,
        design=='nrs'&bias=='high'~ 3,
        design=='nrs'&bias=='low'~ 4,
        bias=='unclear'~ 5
      ))
      bias_index.ipd<- data1%>%group_by(study,bias_index)%>%group_keys()%>%select('bias_index')
      if(!is.null(bias.covariate)){
        # continuous
        if (is.numeric(data1$x.bias) == TRUE) {
          # mean covariate if continuous
          xbias.ipd <- data1%>%
            group_by(study.jags)%>%
            group_map(~mean(.x$x.bias,na.rm = TRUE))%>%unlist()
          # Factor with 2 levels
        } else if (is.factor(data1$x.bias) == TRUE || is.character(data1$x.bias) == TRUE) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data1$x.bias)) > 2)
            stop("crossnma does not currently support bias-regression with categorical variables that have more than two levels.")
          if(length(unique(data1$x.bias)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data1$x.bias) == TRUE)
            data1$x.bias <- as.factor(data1$x.bias)
          data1$x.bias <- as.numeric(data1$x.bias != levels(data1$x.bias)[1])
          data1 <- data1%>%
            group_by(study.jags)%>%
            mutate(x.bias=mean(x.bias,na.rm = TRUE))
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
    if (!is.null(covariate)) {
      # continuous
      if (is.numeric(data1$x1) == TRUE) {
        # mean covariate
        data1 <- data1%>%
          group_by(study.jags)%>%
          mutate(xm1.ipd=mean(x1,na.rm = TRUE))
        # Factor with 2 levels
      } else if (is.factor(data1$x1) == TRUE || is.character(data1$x1) == TRUE ) {
        #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
        if (length(unique(data1$x1)) > 2)
          stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
        if(length(unique(data1$x1)) == 1)
          stop("Covariate should have more than one unique value.")
        if (is.character(data1$x1) == TRUE)
          data1$x1 <- as.factor(data1$x1)
        data1$x1 <- as.numeric(data1$x1 != levels(data1$x1)[1])
        data1 <- data1%>%
          group_by(study.jags)%>%
          mutate(xm1.ipd=mean(x1,na.rm = TRUE))
      } else {stop("Invalid datatype for covariate.")}
      # covariate2
      if(!is.null(data1$x2)){
        # continuous
        if (is.numeric(data1$x2) == TRUE) {
          # mean covariate if continuous
          data1 <- data1%>%
            group_by(study.jags)%>%
            mutate(xm2.ipd=mean(x2,na.rm = TRUE))
          # Factor with 2 levels
        } else if (is.factor(data1$x2) == TRUE || is.character(data1$x2) == TRUE) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data1$x2)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if(length(unique(data1$x2)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data1$x2) == TRUE)
            data1$x2 <- as.factor(data1$x2)
          data1$x2 <- as.numeric(data1$x2 != levels(data1$x2)[1])
          data1 <- data1%>%
            group_by(study.jags)%>%
            mutate(xm2.ipd=mean(x2,na.rm = TRUE))
        } else {stop("Invalid datatype for covariate.")}
      }else {xm2.ipd <- NULL}
      # covariate3
      if(!is.null(data1$x3)){
        # continuous
        if (is.numeric(data1$x3) == TRUE) {
          # mean covariate
          data1 <- data1%>%
            group_by(study.jags)%>%
            mutate(xm3.ipd=mean(x3,na.rm = TRUE))
          # Factor with 2 levels
        } else if (is.factor(data1$x3) == TRUE || is.character(data1$x3) == TRUE) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data1$x3)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if(length(unique(data1$x3)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data1$x3) == TRUE)
            data1$x3 <- as.factor(data1$x3)
          data1$x3 <- as.numeric(data1$x3 != levels(data1$x3)[1])
          data1 <- data1%>%
            group_by(study.jags)%>%
            mutate(xm3.ipd=mean(x3,na.rm = TRUE))
        } else {stop("Invalid datatype for covariate.")}
      }else {xm3.ipd <- NULL}
    } else{
      xm1.ipd <- NULL
      xm2.ipd <- NULL
      xm3.ipd <- NULL
    }

    # create a matrix of treatment per study row
    jagsdata1 <- list()

    # create the matrix of trt index following the values of unfav column (adjust 1&2)
if (method.bias%in%c("adjust1","adjust2")) {
  if(is.null(bias.group)){data1$bias.group <- 1} # Default, make bias adjustment when bias.group is not provided

# From the unfav column create new ref treatment per study
  dd0 <- data1%>%
    group_by(study.jags)%>%
    mutate(ref.trt.std=.data[["trt"]][unfav==0][1])
# For each study, arrange treatments by the new ref
  ns <- length(unique(dd0$study.jags))
  dd1 <-sapply(1:ns,
               function(i){
                 dstd0 <- dd0[dd0$study.jags==unique(dd0$study.jags)[i],]
                 dstd<- dstd0%>%arrange(match(trt,ref.trt.std))
                 }
               ,simplify = FALSE)
  dd2 <- do.call(rbind,dd1)
# create a matrix with the treatment index
  jagsdata1$t.ipd <- dd2 %>%
    select(trt.jags,study.jags)%>% unique()%>%
    group_by(study.jags)%>%
    mutate(arm = row_number())%>%ungroup() %>%
    spread(arm, trt.jags)%>%
    select(-study.jags)%>%
    as.matrix()

    # this returns a message when I run crossnma.model:
    # jagsdata1$t.ipd <-  dd2 %>% group_by(study.jags,trt.jags)%>%
    #   select(trt.jags)%>% unique()%>%
    #   group_by(study.jags)%>%
    #   dplyr::mutate(arm = row_number())%>%ungroup() %>%
    #   spread(arm, trt.jags)%>%
    #   select(-study.jags)%>%
    #   as.matrix()
    # generate JAGS data object
    jagstemp <- data1 %>% select(-c(study,trt,design,bias.group,unfav,bias_index,bias))
    for (v in names(jagstemp)){
      jagsdata1[[v]] <- jagstemp %>% pull(v)
    }
}else{
  jagsdata1$t.ipd <- data1 %>% group_by(study.jags,trt.jags)%>% group_keys()%>%
    group_by(study.jags)%>%
    dplyr::mutate(arm = row_number())%>%ungroup() %>%
    spread(arm, trt.jags)%>%
    select(-study.jags)%>%
    as.matrix()
  # generate JAGS data object
  jagstemp <- data1 %>% select(-c(study,trt,design))
  for (v in names(jagstemp)){
    jagsdata1[[v]] <- jagstemp %>% pull(v)
  }
}

    # create baseline vector
    # data1 %<>% mutate(bl=mapvalues(study,
    #                                from=unique(data1$study),
    #                                to=jagsdata1$t.ipd[,1],
    #                                warn_missing = FALSE) %>% as.integer)





    # add number of treatments, studies, and arms to JAGS data object
    jagsdata1$nt <- trt.key %>% nrow()
    jagsdata1$ns.ipd <- ifelse(!is.null(data1),data1$study%>% unique() %>% length(),0)
    jagsdata1$na.ipd <- data1%>% arrange(study.jags)%>% group_by(study.jags)%>% group_map(~length(unique(.x$trt)))%>%
      unlist()
    jagsdata1$np <- data1 %>% nrow()
    jagsdata1$x.bias <- NULL
    # jagsdata1$bias_index <- NULL
    # jagsdata1$bias <- NULL
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
  if(!is.null(std.data)){
    if(!is.numeric(data2$n)) stop("Sample size must be an integer greater than 0.")
    ifelse(floor(data2$n) != data2$n | data2$n<1, stop("Sample size must be an integer greater than 0."), 1)

    if(!is.numeric(data2$r)) stop("Outcome must be numeric.")


    # add bias_index based on RoB and study design RCT or NRS. when method.bias ='adjust1' or 'adjust2'
    if(!is.null(bias)){
      data2%<>%mutate(bias_index=case_when(
        design=='rct'&bias=='high'~ 1,
        design=='rct'&bias=='low'~ 2,
        design=='nrs'&bias=='high'~ 3,
        design=='nrs'&bias=='low'~ 4,
        bias=='unclear'~ 5
      ))
      bias_index.ad<- data2%>%group_by(study,bias_index)%>%group_keys()%>%select('bias_index')
      if(!is.null(bias.covariate)){
        # continuous
        if (is.numeric(data2$x.bias) == TRUE) {
          # mean covariate if continuous
          xbias.ad <- data2%>%
            group_by(study)%>%
            group_map(~mean(.x$x.bias,na.rm = TRUE))%>%unlist()
          # Factor with 2 levels
        } else if (is.factor(data2$x.bias) == TRUE || is.character(data2$x.bias) == TRUE) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data2$x.bias)) > 2)
            stop("crossnma does not currently support bias-regression with categorical variables that have more than two levels.")
          if(length(unique(data2$x.bias)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data2$x.bias) == TRUE)
            data2$x.bias <- as.factor(data2$x.bias)
          data2$x.bias <- as.numeric(data2$x.bias != levels(data2$x.bias)[1])
          data2 <- data2%>%
            group_by(study)%>%
            mutate(x.bias=mean(x.bias,na.rm = TRUE))
        } else {stop("Invalid datatype for bias covariate.")}

        # bias_index.ad <- NULL
      }else{
      xbias.ad <- NULL
      }

    }else{
      bias_index.ad <- NULL
      xbias.ad <- NULL
    }

    # pre-process the covariate if specified
    if (!is.null(covariate)) {
      # continuous
      if (is.numeric(data2$x1) == TRUE) {
        # mean covariate
        data2 <- data2%>%
          group_by(study)%>%
          mutate(xm1.ad=mean(x1,na.rm = TRUE))
        # Factor with 2 levels
      } else if (is.factor(data2$x1) == TRUE || is.character(data2$x1) == TRUE) {
        #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
        if (length(unique(data2$x1)) > 2)
          stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
        if(length(unique(data2$x1)) == 1)
          stop("Covariate should have more than one unique value.")
        if (is.character(data2$x1) == TRUE)
          data2$x1 <- as.factor(data2$x1)
        data2$x1 <- as.numeric(data2$x1 != levels(data2$x1)[1])
        data2 <- data2%>%
          group_by(study)%>%
          mutate(xm1.ad=mean(x1,na.rm = TRUE))
      } else {stop("Invalid datatype for covariate.")}

      # covariate2
      if(!is.null(data2$x2)){
        # continuous
        if (is.numeric(data2$x2) == TRUE) {
          # mean covariate
          data2 <- data2%>%
            group_by(study)%>%
            mutate(xm2.ad=mean(x2,na.rm = TRUE))
          # Factor with 2 levels
        } else if (is.factor(data2$x2) == TRUE || is.character(data2$x2) == TRUE) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data2$x2)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if(length(unique(data2$x2)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data2$x2) == TRUE)
            data2$x2 <- as.factor(data2$x2)
          data2$x2 <- as.numeric(data2$x2 != levels(data2$x2)[1])
          data2 <- data2%>%
            group_by(study)%>%
            mutate(xm2.ad=mean(x2,na.rm = TRUE))
        } else {stop("Invalid datatype for covariate.")}
      }else {xm2.ad <- NULL}
      # covariate3
      if(!is.null(data2$x3)){
        # continuous
        if (is.numeric(data2$x3) == TRUE) {
          # mean covariate
          data2 <- data2%>%
            group_by(study)%>%
            mutate(xm3.ad=mean(x3,na.rm = TRUE))
          # Factor with 2 levels
        } else if (is.factor(data2$x3) == TRUE || is.character(data2$x3) == TRUE) {
          #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
          if (length(unique(data2$x3)) > 2)
            stop("crossnma does not currently support meta-regression with categorical variables that have more than two levels.")
          if(length(unique(data2$x3)) == 1)
            stop("Covariate should have more than one unique value.")
          if (is.character(data2$x3) == TRUE)
            data2$x3 <- as.factor(data2$x3)
          data2$x3 <- as.numeric(data2$x3 != levels(data2$x3)[1])
          data2 <- data2%>%
            group_by(study)%>%
            mutate(xm3.ad=mean(x3,na.rm = TRUE))
        } else {stop("Invalid datatype for covariate.")}
      }else {xm3.ad <- NULL}
    } else{
      xm1.ad <- NULL
      xm2.ad <- NULL
      xm3.ad <- NULL}


    #generate JAGS data object

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
    # create the matrix of trt index following the values of unfav column (adjust 1&2)
    if (method.bias%in%c("adjust1","adjust2")) {
      if(is.null(bias.group)){data2$bias.group <- 1} # Default, make bias adjustment when bias.group is no provided

      # From the unfav column create new ref treatment per study
      dd0 <- data2%>%
        group_by(study.jags)%>%
        mutate(ref.trt.std=.data[["trt"]][unfav==0])
      # For each study, arrange treatments by the new ref
      ns <- length(unique(dd0$study.jags))
      dd1 <-sapply(1:ns,
                   function(i){
                     dstd0 <- dd0[dd0$study.jags==unique(dd0$study.jags)[i],]
                     dstd<- dstd0%>%arrange(match(trt,ref.trt.std))
                   }
                   ,simplify = FALSE)
      dd2 <- do.call(rbind,dd1)


      # create a matrix with the treatment index
      jagstemp2 <- dd2 %>%arrange(study.jags) %>% group_by(study.jags) %>% dplyr::mutate(arm = row_number()) %>%
        ungroup()%>% dplyr::select(-c(trt,design,bias,ref.trt.std,unfav,bias.group,bias_index))  %>% gather("variable", "value", -study,-study.jags, -arm) %>% spread(arm, value)
      jagsdata2 <- list()
      for (v in unique(jagstemp2$variable)){
        jagsdata2[[v]] <- as.matrix(jagstemp2 %>% filter(variable == v) %>% select(-study,-study.jags, -variable))
      }
    }else{
      jagstemp2 <- data2 %>% arrange(study.jags,trt.jags) %>% group_by(study.jags) %>% dplyr::mutate(arm = row_number()) %>%
        ungroup()%>% select(-c(trt,design,bias))  %>% gather("variable", "value", -study,-study.jags, -arm) %>% spread(arm, value)

      jagsdata2 <- list()
      for (v in unique(jagstemp2$variable)){
        jagsdata2[[v]] <- as.matrix(jagstemp2 %>% filter(variable == v) %>% select(-study,-study.jags, -variable))
      }
    }

    # add number of treatments, studies, and arms to JAGS data object
    jagsdata2$ns.ad <- ifelse(!is.null(data2),data2$study.jags %>% unique()%>%length(),0)
    jagsdata2$na.ad <- data2 %>% group_by(study.jags) %>%dplyr::summarize(n.arms = n()) %>%
      ungroup() %>% select(n.arms) %>% t() %>% as.vector
    # add covariate
    jagsdata2$x1 <- jagsdata2$x2 <- jagsdata2$x3 <- NULL
    jagsdata2$xm1.ad <- unique(data2$xm1.ad)
    jagsdata2$xm2.ad <- unique(data2$xm2.ad)
    jagsdata2$xm3.ad <- unique(data2$xm3.ad)
    jagsdata2$x.bias <- NULL
    # jagsdata2$bias_index <- NULL
    # jagsdata2$bias <- NULL
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
  bmat <- rbind(ifelse(!is.null(data1),list(data1%>% group_by(study.jags)%>%select(bias.group)%>%unique()%>%select(bias.group)), list(NULL))[[1]],
                ifelse(!is.null(data2),list(data2%>% group_by(study.jags)%>%select(bias.group)%>%unique()), list(NULL))[[1]]
                )
  jagsdata$std.in <-bmat$study.jags[bmat$bias.group==1]
  jagsdata$std.act.no <-bmat$study.jags[bmat$bias.group==0]
  jagsdata$std.act.yes <-bmat$study.jags[bmat$bias.group==2]


  #====================================
  # use NRS as prior, jags need to be run for only NRS

  if(method.bias=='prior'){
    # data NRS
    if(!(reference %in% data1.nrs$trt)&!(reference %in% data2.nrs$trt)) stop("Reference treatment should be present in the list of treatments in NRS.")
    trt.df.nrs <- data.frame(trt=unique(c(data1.nrs$trt,as.character(data2.nrs$trt))))

    trt.key.nrs <- trt.df.nrs$trt %>% unique %>% sort %>% tibble(trt.ini=.) %>%
      filter(trt.ini!=reference) %>% add_row(trt.ini=reference, .before=1) %>%
      mutate(trt.jags = 1:dim(.)[1])


    #====================================
    # 1. IPD-NRS

    if(!(reference %in% data1.nrs$trt)) stop("Reference treatment is not present in the list of treatments.")

    #Trt mapping
    data1.nrs %<>% mutate(trt.jags=mapvalues(trt,
                                             from=trt.key.nrs$trt.ini,
                                             to=trt.key.nrs$trt.jags,
                                             warn_missing = FALSE)%>%as.integer)
    data1.nrs %<>% mutate(study.jags=mapvalues(study,
                                               from=unique(study),
                                               to=seq_len(length(unique(study))),
                                               warn_missing = FALSE
    ) %>% as.character %>% as.integer)


    # add a matrix of treatment per study row
    jagsdata.nrs <- list()

    jagsdata.nrs$t.ipd <- data1.nrs %>% group_by(study,trt.jags)%>% group_keys()%>%
      group_by(study)%>%
      dplyr::mutate(arm = row_number())%>%ungroup() %>%
      spread(arm, trt.jags)%>%
      select(-study)%>%
      as.matrix()
    # add baseline vector
    # data1.nrs %<>% mutate(bl=mapvalues(study,
    #                                    from=unique(data1.nrs$study),
    #                                    to=jagsdata.nrs$t.ipd[,1],
    #                                    warn_missing = FALSE) %>% as.integer)

    #generate JAGS data object

    jagstemp.nrs1 <- data1.nrs %>% select(-c(study,trt,design))
    for (v in names(jagstemp.nrs1)){
      jagsdata.nrs[[v]] <- jagstemp.nrs1 %>% pull(v) #%>% as.vector() #%>% select(-trial, -variable)
    }
    #modify BUGS object for the various family/link combinations
    names(jagsdata.nrs)[names(jagsdata.nrs) == "outcome"]  <- "y"
    names(jagsdata.nrs)[names(jagsdata.nrs) == "trt.jags"] <- "trt"
    names(jagsdata.nrs)[names(jagsdata.nrs) == "study.jags"] <- "study"

    #add number of treatments, studies, and arms to BUGS data object
    jagsdata.nrs$nt <- trt.key.nrs %>% nrow()
    jagsdata.nrs$ns.ipd <- ifelse(!is.null(prt.data),data1.nrs$study%>% unique() %>% length(),0)
    jagsdata.nrs$na.ipd <- data1.nrs%>% arrange(study.jags)%>% group_by(study.jags)%>% group_map(~length(unique(.x$trt)))%>%
      unlist()
    jagsdata.nrs$np <- data1.nrs %>% nrow()

    #====================================
    # AD - NRS

    #add treatment mapping to data
    data2.nrs %<>% mutate(trt.jags=mapvalues(trt,
                                             from=trt.key.nrs$trt.ini,
                                             to=trt.key.nrs$trt.jags,
                                             warn_missing = FALSE) %>% as.integer)

    jagstemp.nrs2 <- data2.nrs %>% arrange(study) %>% group_by(study) %>% dplyr::mutate(arm = row_number()) %>%
      ungroup()%>% select(-c(trt,design))  %>% gather("variable", "value", -study, -arm) %>% spread(arm, value)

    for (v in unique(jagstemp.nrs2$variable)){
      jagsdata.nrs[[v]] <- as.matrix(jagstemp.nrs2 %>% filter(variable == v) %>% select(-study, -variable))
    }
    names(jagsdata.nrs)[names(jagsdata.nrs) == "trt.jags"] <- "t.ad"

    #add number of treatments, studies, and arms to JAGS data object
    jagsdata.nrs$ns.ad <- ifelse(!is.null(std.data),data2.nrs$study %>% unique()%>%length(),0)
    jagsdata.nrs$na.ad <- data2.nrs %>% group_by(study) %>%dplyr::summarize(n.arms = n()) %>%
      ungroup() %>% select(n.arms) %>% t() %>% as.vector


    # jags code NRS
    model.nrs <- crossnma.code(ipd = ifelse(nrow(data1.nrs)==0,F,T),
                           ad = ifelse(nrow(data2.nrs)==0,F,T),
                           trt.effect='random',
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

    seeds <- sample(.Machine$integer.max, run.nrs[['n.chains']], replace = FALSE)
    inits <- list()
    for (i in 1:run.nrs[['n.chains']])
      inits[[i]] <- list(.RNG.seed = seeds[i], .RNG.name = "base::Mersenne-Twister")

    jagsmodel.nrs <- jags.model(textConnection(model.nrs),        #Create a connection so JAGS can access the variables
                                jagsdata.nrs,
                                n.chains=run.nrs[['n.chains']],
                                n.adapt=run.nrs[['n.adapt']],
                                inits = inits)
    jagssamples.nrs <- coda.samples(jagsmodel.nrs,
                                    variable.names="d",
                                    n.iter=run.nrs[['n.iter']],
                                    thin=run.nrs[['thin']]
    )

    # # # # # # # # # # # #
    # Output: prior for d's

    # map NRS trt to RCT trt
    trt.key2 <- trt.key %>% mutate(trt.jags.nrs=mapvalues(trt.ini,
                                                          from=trt.key.nrs$trt.ini,
                                                          to=trt.key.nrs$trt.jags,
                                                          warn_missing = FALSE)%>%as.integer)
    d.nrs <- summary(jagssamples.nrs)[[1]][,'Mean']+ifelse(is.null(run.nrs$mean.shift),0,run.nrs$mean.shift)
    prec.nrs <- ifelse(is.null(run.nrs$var.infl),1,run.nrs$var.infl)/(summary(jagssamples.nrs)[[1]][,'SD']^2)

    d.prior.nrs <- ""
    for (i in 2:nrow(trt.key2)) {
      d.prior0 <- paste0("d[",trt.key2$trt.jags[i],
                         "]~dnorm(",
                         ifelse(is.na(d.nrs[trt.key2$trt.jags.nrs[i]])|d.nrs[trt.key2$trt.jags.nrs[i]]==0,0,d.nrs[trt.key2$trt.jags.nrs[i]]),",",
                         ifelse(is.na(prec.nrs[trt.key2$trt.jags.nrs[i]])|prec.nrs[trt.key2$trt.jags.nrs[i]]==Inf,10^-4,prec.nrs[trt.key2$trt.jags.nrs[i]]),
                         ")
                    ")
      d.prior.nrs <- paste0(d.prior.nrs,d.prior0)
    }

  }else {d.prior.nrs <- NULL}

  #=======================
  # jags code
  model <- crossnma.code(ipd = !is.null(prt.data),
                     ad = !is.null(std.data),
                     trt.effect=trt.effect,
                     covariate=covariate,
                     split.regcoef=split.regcoef,
                     reg0.effect=reg0.effect,
                     regb.effect=regb.effect,
                     regw.effect=regw.effect,
                     # reg.effect=reg.effect,
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
                     prior.pi.high.rct=prior[['pi.high.rct']],
                     prior.pi.low.rct=prior[['pi.low.rct']],
                     prior.pi.high.nrs=prior[['pi.high.nrs']],
                     prior.pi.low.nrs=prior[['pi.low.nrs']],
                     method.bias = method.bias,
                     d.prior.nrs=d.prior.nrs)

  crossmodel <- structure(list(jagsmodel=model,
                           data=jagsdata,
                           trt.key=trt.key,
                           study.key=study.key,
                           trt.effect=trt.effect,
                           method.bias=method.bias,
                           covariate=covariate,
                           split.regcoef=split.regcoef,
                           regb.effect=regb.effect,
                           regw.effect=regw.effect,
                           bias.effect=bias.effect,
                           bias.type=bias.type
  ),
  class = "crossnmaModel")
  return(crossmodel)
}

