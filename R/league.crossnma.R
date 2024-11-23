#' League Table
#'
#' @description
#' Produces a league table that contains point estimates of relative
#' effects for all possible pairs of treatments along with 95\%
#' credible intervals obtained with the quantile method.
#'
#' @param x An object created with \code{\link{crossnma}}.
#' @param median A logical indicating whether to use the median
#'   (default) or mean to measure relative treatment effects.
#' @param backtransf A logical indicating whether results should be
#'   back transformed. If \code{backtransf = TRUE}, results for
#'   \code{sm = "OR"} are presented as odds ratios rather than log
#'   odds ratios, for example.
#' @param order A vector of treatment names (character) representing
#'   the order in which to display these treatments.
#' @param cov1.value The participant covariate value of \code{cov1}
#'   for which to report the results. Must be specified for network
#'   meta-regression and when individual participant dataset is used
#'   in the analysis. For dichotomous covariates, a character of the
#'   level (used in the data) should be indicated.
#' @param cov2.value The participant covariate value of \code{cov2}
#'   for which to report the results. Must be specified for network
#'   meta-regression and when individual participant dataset is used
#'   in the analysis. For dichotomous covariates, a character of the
#'   level (used in the data) should be indicated.
#' @param cov3.value The participant covariate value of \code{cov3}
#'   for which to report the results. Must be specified for network
#'   meta-regression and when individual participant dataset is used
#'   in the analysis. For dichotomous covariates, a character of the
#'   level (used in the data) should be indicated.
#' @param digits The number of digits to be used when displaying the
#'   results.
#' @param direction The format to display the league table. Two
#'   options "wide" (default) and "long".
#' @param exp Deprecated argument (replaced by \code{backtransf}).
#' @param \dots Additional arguments (ignored at the moment).
#'
#' @return
#' A league table. Row names indicate comparator treatments.  The
#' table will be displayed in a long or wide formatting.
#'
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch}
#'
#' @seealso \code{\link{crossnma}}
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
#' # Create league tables
#' league(fit)                     #  wide format
#' league(fit, direction = "long") #  long format
#' }
#'
#' @method league crossnma
#' @export


league.crossnma <- function(x,
                            median = TRUE,
                            backtransf = x$model$backtransf,
                            order = NULL,
                            cov1.value=NULL,
                            cov2.value=NULL,
                            cov3.value=NULL,
                            digits = gs("digits"),
                            direction = "wide",
                            exp = backtransf,
                            ...) {

  chkclass(x, "crossnma")
  ##
  chklogical(median)
  if (median)
    central.tdcy <- "median"
  else
    central.tdcy <- "mean"
  ##
  chknumeric(digits, min = 0, length = 1)
  direction <- setchar(direction, c("wide", "long"))
  ##
  missing.backtransf <- missing(backtransf)
  backtransf <-
    deprecated2(backtransf, missing.backtransf, exp, missing(exp))
  chklogical(backtransf)
  if (backtransf & x$model$sm %in% c("MD", "SMD")) {
    if (!missing.backtransf)
      warning("No back transformation of results for (standardised) ",
              "mean differences (argument 'backtransf').")
    backtransf <- FALSE
  }
  ##
  exp <- backtransf
  quant <- x$model$quantiles


  if (!is.null(x$model$covariate) & is.null(cov1.value))
    stop("cov1.value must be specified for network meta-regression")


  ## Bind variables to function
  trt <- Treatment <- Comparator <- cov.ref <- NULL
  #samples <- as.mcmc(x$jagsfit)
  samples <- x$samples

  dmat <-
    do.call(rbind, samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- x$trt.key$trt.ini
  colnames(dmat) <- trt.names
  ##
  if (is.null(order))
    order <- trt.names
  else
    order <- setseq(order, trt.names)


  ##
  ##
  ## In the case of network meta-regression
  ##
  ##
  if (!is.null(x$model$covariate)) {
    ##
    ## (a) If IPD is available
    ##
    if (x$model$data$ns.ipd != 0) {
      bwmat.cov2 <- bwmat.cov3 <- 0
      bbmat.cov2 <- bbmat.cov3 <- 0
      bmat.cov2 <- bmat.cov3 <- 0
      ##
      nc <- length(x$model$covariate)
      nt <- length(trt.names)
      ##
      labs <- x$model$dich.cov.labels
      ##
      if (x$model$split.regcoef) {
        ## betaw
        if (x$model$regw.effect == "independent") {
          bwmat <- do.call(rbind, samples) %>% data.frame() %>%
            select(starts_with("betaw.t_"))
          ## split bwmat by nt_column to generate bwmat.cov for each covariate
          ##
          ## For factor covariate, multiply by 0 or 1 depends on what
          ## value the user indicate in cov1.value
          ##
          bwmat.cov1 <- bwmat[, 1:nt] *
            if (is.numeric(cov1.value))
              cov1.value - x$model$cov.ref[1]
            else
              as.numeric(labs[labs[, 1] == cov1.value, 2])
          ##
          if (nc == 2)
            bwmat.cov2 <-
              bwmat[, (nt + 1):(nt * 2)] *
              if (is.numeric(cov2.value))
                cov2.value - x$model$cov.ref[2]
              else
                as.numeric(labs[labs[, 1] == cov2.value, 2])
          ##
          if (nc == 3)
            bwmat.cov3 <-
              bwmat[, (nt * 2 + 1):(nt * 3)] *
              if (is.numeric(cov3.value))
                cov3.value - x$model$cov.ref[3]
              else
                as.numeric(labs[labs[, 1] == cov3.value, 2])
          ##
          dmat <- dmat + bwmat.cov1 + bwmat.cov2 + bwmat.cov3
        }
        else {
          bwmat <- do.call(rbind, samples) %>% data.frame() %>%
            select(starts_with("bw_"))
          bwmat.cov1 <-
            sweep(cbind(bwmat[, 1]), MARGIN = 2,
                  if (is.numeric(cov1.value))
                    cov1.value - x$model$cov.ref[1]
                  else
                    as.numeric(labs[labs[, 1] == cov1.value, 2]),
                  '*')
          ## repeat the column for each trt to add it to dmat
          bwmat.cov1 <-
            matrix(unlist(rep(bwmat.cov1, each = nt)), ncol = nt, byrow = TRUE)
          ##
          if (nc == 2) {
            bwmat.cov2 <-
              sweep(cbind(bwmat[, 2]), MARGIN = 2,
                    if (is.numeric(cov2.value))
                      cov2.value - x$model$cov.ref[2]
                    else
                      as.numeric(labs[labs[, 1] == cov2.value, 2]),
                    '*')
            bwmat.cov2 <-
              matrix(unlist(rep(bwmat.cov2, each = nt)),
                     ncol = nt, byrow = TRUE)
          }
          ##
          if (nc == 3) {
            bwmat.cov3 <-
              sweep(cbind(bwmat[, 3]), MARGIN = 2,
                    if (is.numeric(cov3.value))
                      cov3.value - x$model$cov.ref[3]
                    else
                      as.numeric(labs[labs[, 1] == cov3.value, 2]),
                    '*')
            bwmat.cov3 <-
              matrix(unlist(rep(bwmat.cov3, each = nt)),
                     ncol = nt, byrow = TRUE)
          }
          ##
          dmat <- dmat + bwmat.cov1 + bwmat.cov2 + bwmat.cov3
        }
        ## betab
        if (x$model$regb.effect == "independent") {
          bbmat <- do.call(rbind, samples) %>% data.frame() %>%
            select(starts_with("betab.t_"))
          ## split bbmat by nt_column to generate bbmat.cov for each
          ## covariate
          stds.mean1 <-
            mean(c(x$model$data$xm1.ad, x$model$data$xm1.ipd), na.rm = TRUE)
          bbmat.cov1 <- bbmat[, 1:nt] *
            if (is.numeric(cov1.value))
              stds.mean1 - x$model$cov.ref[1]
            else
              stds.mean1
          ##
          if (nc == 2) {
            stds.mean2 <-
              mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd), na.rm = TRUE)
            bbmat.cov2 <- bbmat[, (nt + 1):(nt * 2)] *
              if (is.numeric(cov2.value))
                stds.mean2 - x$model$cov.ref[2]
              else
                stds.mean2
          }
          if (nc == 3) {
            stds.mean3 <-
              mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd), na.rm = TRUE)
            bbmat.cov3 <- bbmat[, (nt * 2 + 1):(nt * 3)] *
              if (is.numeric(cov3.value))
                stds.mean3 - x$model$cov.ref[3]
              else
                stds.mean3
          }
          ##
          dmat <- dmat + bbmat.cov1 + bbmat.cov2 + bbmat.cov3
        }
        else {
          stds.mean <-
            c(mean(c(x$model$data$xm1.ad, x$model$data$xm1.ipd),na.rm = TRUE),
              mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd),na.rm = TRUE),
              mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd),na.rm = TRUE))
          ##
          bbmat <- do.call(rbind, samples) %>% data.frame() %>%
            select(starts_with("bb_"))
          ##
          bbmat.cov1 <-
            sweep(cbind(bbmat[, 1]), MARGIN = 2,
                  if (is.numeric(cov1.value))
                    stds.mean[1] - x$model$cov.ref[1]
                  else
                    stds.mean[1],
                  '*')
          bbmat.cov1 <-
            matrix(unlist(rep(bbmat.cov1, each = nt)), ncol = nt, byrow = TRUE)
          ##
          if (nc == 2) {
            bbmat.cov2 <-
              sweep(cbind(bbmat[, 2]), MARGIN = 2,
                    if (is.numeric(cov2.value))
                      stds.mean[2] - x$model$cov.ref[2]
                    else
                      stds.mean[2],
                    '*')
            bbmat.cov2 <-
              matrix(unlist(rep(bbmat.cov2, each = nt)),
                     ncol = nt, byrow = TRUE)
          }
          ##
          if (nc == 3) {
            bbmat.cov3 <-
              sweep(cbind(bbmat[, 3]), MARGIN = 2,
                    if (is.numeric(cov3.value))
                      stds.mean[3] - x$model$cov.ref[3]
                    else
                      stds.mean[3],
                    '*')
            bbmat.cov3 <-
              matrix(unlist(rep(bbmat.cov3, each = nt)),
                     ncol = nt, byrow = TRUE)
          }
          ##
          dmat <- dmat + bbmat.cov1 + bbmat.cov2 + bbmat.cov3
        }
      }
      else {
        if (x$model$regb.effect == "independent" &&
            x$model$regw.effect == "independent") {
          bmat <- do.call(rbind, samples) %>% data.frame() %>%
            select(starts_with("beta.t_"))
          ## For factor covariate, multiply by 0 or 1 depends on what
          ## value the user indicate in cov1.value
          ##
          bmat.cov1 <-
            bmat[, 1:nt] *
            if (is.numeric(cov1.value))
              cov1.value - x$model$cov.ref[1]
            else
              as.numeric(labs[labs[, 1] == cov1.value, 2])
          ##
          if (nc == 2)
            bmat.cov2 <-
              bmat[, (nt + 1):(nt * 2)] *
              if (is.numeric(cov2.value))
                cov2.value - x$model$cov.ref[2]
              else
                as.numeric(labs[labs[, 1] == cov2.value, 2])
          ##
          if (nc == 3)
            bmat.cov3 <-
              bmat[, (nt * 2 + 1):(nt * 3)] *
              if (is.numeric(cov3.value))
                cov3.value - x$model$cov.ref[3]
              else
                as.numeric(labs[labs[, 1] == cov3.value, 2])
          ##
          dmat <- dmat + bmat.cov1 + bmat.cov2 + bmat.cov3
        }
        else {
          bmat <- do.call(rbind, samples) %>% data.frame() %>%
            select(starts_with("b_"))
          bmat.cov1 <-
            sweep(cbind(bmat[, 1]), MARGIN = 2,
                  if (is.numeric(cov1.value))
                    cov1.value - x$model$cov.ref[1]
                  else
                    as.numeric(labs[labs[, 1] == cov1.value, 2]),
                  '*')
          ## repeat the column for each trt to add it to dmat
          bmat.cov1 <-
            matrix(unlist(rep(bmat.cov1, each = nt)), ncol = nt, byrow = TRUE)
          ##
          if (nc == 2) {
            bmat.cov2 <-
              sweep(cbind(bmat[, 2]), MARGIN = 2,
                    if (is.numeric(cov2.value))
                      cov2.value - x$model$cov.ref[2]
                    else
                      as.numeric(labs[labs[, 1] == cov2.value, 2]),
                    '*')
            bmat.cov2 <-
              matrix(unlist(rep(bmat.cov2, each = nt)),
                     ncol = nt, byrow = TRUE)
          }
          ##
          if (nc == 3) {
            bmat.cov3 <-
              sweep(cbind(bmat[, 3]), MARGIN = 2,
                    if (is.numeric(cov3.value))
                      cov3.value - x$model$cov.ref[3]
                    else
                      as.numeric(labs[labs[, 1] == cov3.value, 2]),
                    '*')
            bmat.cov3 <-
              matrix(unlist(rep(bmat.cov3, each = nt)),
                     ncol = nt, byrow = TRUE)
          }
          ##
          dmat + bmat.cov1 + bmat.cov2 + bmat.cov3
        }
      }
    }
    else {
      ##
      ## (b) if only AD is available (i.e., only betab)
      ##
      bbmat.cov2 <- bbmat.cov3 <- 0
      ##
      if (x$model$regb.effect == "independent") {
        bbmat <- do.call(rbind, samples) %>% data.frame() %>%
          select(starts_with("beta.t_"))
        ## split bbmat by nt_column to generate bbmat.cov for each
        ## covariate
        stds.mean1 <- x$model$data$xm1.ad
        bbmat.cov1 <- bbmat[, 1:nt] *
          if (!is.na(x$model$cov.ref[1]))
            stds.mean1 - x$model$cov.ref[1]
          else
            stds.mean1
        ##
        if (nc == 2) {
          stds.mean2 <- x$model$data$xm2.ad
          bbmat.cov2 <-
            bbmat[, (nt + 1):(nt * 2)] *
            if (!is.na(x$model$cov.ref[2]))
              stds.mean2 - x$model$cov.ref[2]
            else
              stds.mean2
        }
        ##
        if (nc == 3) {
          stds.mean3 <- x$model$data$xm3.ad
          bbmat.cov3 <- bbmat[, (nt * 2 + 1):(nt * 3)] *
            if (!is.na(x$model$cov.ref[3]))
              stds.mean3 - x$model$cov.ref[3]
            else
              stds.mean3
        }
        ##
        dmat <- dmat + bbmat.cov1 + bbmat.cov2 + bbmat.cov3
      }
      else {
        stds.mean <-
          c(x$model$data$xm1.ad, x$model$data$xm2.ad, x$model$data$xm3.ad)
        bbmat <- do.call(rbind, samples) %>% data.frame() %>%
          select(starts_with("b_"))
        bbmat.cov1 <-
          sweep(cbind(bbmat[, 1]), MARGIN = 2,
                if (!is.na(x$model$cov.ref[1]))
                  stds.mean[1] - x$model$cov.ref[1]
                else stds.mean[1],
                '*')
        bbmat.cov1 <-
          matrix(unlist(rep(bbmat.cov1, each = nt)), ncol = nt, byrow = TRUE)
        ##
        if (nc == 2) {
          bbmat.cov2 <-
            sweep(cbind(bbmat[, 2]), MARGIN = 2,
                  if (!is.na(x$model$cov.ref[2]))
                    stds.mean[2] - x$model$cov.ref[2]
                  else
                    stds.mean[2],
                  '*')
          bbmat.cov2 <-
            matrix(unlist(rep(bbmat.cov2, each = nt)), ncol = nt, byrow = TRUE)
        }
        ##
        if (nc == 3) {
          bbmat.cov3 <-
            sweep(cbind(bbmat[, 3]), MARGIN = 2,
                  if (!is.na(x$model$cov.ref[3]))
                    stds.mean[3] - x$model$cov.ref[3]
                  else
                    stds.mean[3],
                  '*')
          bbmat.cov3 <-
            matrix(unlist(rep(bbmat.cov3, each = nt)), ncol = nt, byrow = TRUE)
        }
        ##
        dmat <- dmat + bbmat.cov1 + bbmat.cov2 + bbmat.cov3
      }
    }
  }
  ##
  dmat %<>% select(order)
  trt.names <- order
  ##
  dmat %<>% select(order)
  trt.names <- order


  ## # when we adjust for covariate, add a message with the values that we adjusted for
  ## if (!is.null(prt.cov.value)) {
  ##   msg.adjust <-  paste0("The estimates of treatment effect are adjusted for \n
  ##                               participant ",x$model$covariate[1],"= ",cov1.value,
  ##                          "and study mean = ",round(mean(c(x$model$data$xm1.ad, x$model$data$xm1.ipd),na.rm = TRUE),2)
  ##   )
  ##   if (nc == 2) msg.adjust <- paste0(msg.adjust, paste0("\n participant ",x$model$covariate[2],"= ",cov2.value,
  ##                                       "and study mean = ",round(mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd),na.rm = TRUE),2)))
  ##   if (nc == 3) msg.adjust <- paste0(msg.adjust, paste0("\n participant ",x$model$covariate[3],"= ",cov3.value,
  ##                                       "and study mean = ",round(mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd),na.rm = TRUE),2))
  ##                    )
  ## }


  ## Useful functions to compute some statistics
  ##
  calc.report <- function(x, fct = "identity", arg = NULL,
                          trans = "identity") {
    if (is.null(arg))
      eval(call(fct,(call(trans, x))))
    else
      eval(call(fct,(call(trans, x)), arg))
  }
  ##
  exp.mean <- function(x)
    calc.report(x, "mean", trans = "exp")
  exp.median <- function(x)
    calc.report(x, "median", trans = "exp")
  exp.sd <- function(x)
    calc.report(x, "sd", trans = "exp")
  exp.lcl <- function(x)
    calc.report(x, "quantile", arg = quant[1], trans = "exp")
  exp.ucl <- function(x)
    calc.report(x, "quantile", arg = quant[3], trans = "exp")
  ##
  id.mean <- function(x)
    calc.report(x, "mean", trans = "identity")
  id.median <- function(x)
    calc.report(x, "median", trans = "identity")
  id.sd <- function(x)
    calc.report(x, "sd", trans = "identity")
  id.lcl <- function(x)
    calc.report(x, "quantile", arg = quant[1], trans = "identity")
  id.ucl <- function(x)
    calc.report(x, "quantile", arg = quant[3], trans = "identity")
  ##
  colvals <- function(dmat, b.col=1, paste = TRUE) {
    ## Bind variables to function
    key <- value <- trt <- estimate <- lcl <- ucl <- result <- NULL

    base <- colnames(dmat)[b.col]

    dmat2 <- dmat
    new.vars <- paste0(colnames(dmat2), "-", b.col)
    for (i in 1:ncol(dmat)) {
      dmat2[[new.vars[i]]] <- dmat[, i] - dmat[, b.col]
    }
    dmat2 %<>% select(new.vars)
    colnames(dmat2) <- trt.names

    if (central.tdcy == "mean" & exp) {
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = exp.mean)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    }
    else if (central.tdcy == "mean" & exp == FALSE) {
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = id.mean)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    }
    ##
    if (central.tdcy == "median" & exp) {
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = exp.median))%>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    }
    else if (central.tdcy == "median" & exp == FALSE) {
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = id.median)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    }
    ##
    if (exp) {
      tmp.lcl <- dmat2 %>%
        summarise_all(list(lcl = exp.lcl)) %>% gather() %>%
        rename(trt = key, lcl = value) %>%
        mutate(trt = sub("_lcl", "", trt))
      ##
      tmp.ucl <- dmat2 %>%
        summarise_all(list(ucl = exp.ucl)) %>% gather() %>%
        rename(trt = key, ucl = value) %>%
        mutate(trt = sub("_ucl", "", trt))
      ##
      null.value <- 1
    }
    else {
      tmp.lcl <- dmat2 %>%
        summarise_all(list(lcl = id.lcl)) %>% gather() %>%
        rename(trt = key, lcl = value) %>%
        mutate(trt = sub("_lcl", "", trt))
      ##
      tmp.ucl <- dmat2 %>%
        summarise_all(list(ucl = id.ucl)) %>% gather() %>%
        rename(trt = key, ucl = value) %>%
        mutate(trt = sub("_ucl", "", trt))
      ##
      null.value <- 0
    }


    ## Create C-style formatting string from the digits parameter
    ##
    fmt <- paste0("%.", digits, "f")

    if (paste) {
      tmp1 <- left_join(tmp.estimate, tmp.lcl, by = "trt") %>%
        left_join(tmp.ucl, by = "trt") %>%
        mutate(result = paste(sprintf(fmt, estimate),
                              formatCI(round(lcl, digits),
                                       round(ucl, digits)))) %>%
        select(trt, result)
      colnames(tmp1)[2] <-
        as.character(tmp.estimate %>% filter(estimate == null.value) %>%
                     select(trt))
    }
    else {
      tmp1 <- left_join(tmp.estimate, tmp.lcl, by = "trt") %>%
        left_join(tmp.ucl, by = "trt")
      colnames(tmp1)[2] <- central.tdcy
    }
    tmp1
  }


  ##
  ##
  ## Return league table
  ##
  ##
  if (direction == "wide") {
    ##
    ## Wide layout
    ##
    tmp1.list <- list()
    ##
    for (i in 1:ncol(dmat))
      tmp1.list[[i]] <- colvals(dmat, b.col=i)
    ##
    widetable <- suppressMessages(bind_cols(tmp1.list)) %>%
      select(-starts_with("trt")) %>%
      t()
    ##
    colnames(widetable) <- colnames(dmat)
    rownames(widetable) <- colnames(dmat)
    ##
    class(widetable) <- c("league.crossnma", class(widetable))
    attr(widetable, "direction") <- direction
    ##
    for (i in 1:dim(widetable)[1])
      widetable[i, i] <- colnames(widetable)[i]
    ##
    return(widetable)
  }
  else {
    ##
    ## Long layout
    ##
    tmp2.list <- list()
    for (i in 1:ncol(dmat))
      tmp2.list[[i]] <- colvals(dmat, b.col=i, paste=FALSE)
    ##
    longtable <- tmp2.list %>%
      bind_rows() %>%
      mutate(Treatment = trt,
             Comparator =
               rep(trt.names, each = length(trt.names))) %>%
      select(Treatment, Comparator, everything(), -trt)
    ##
    fmt <- paste0("%.", digits, "f")
    ##
    if (median)
      longtable$median <-
        sprintf(fmt, round(longtable$median, digits))
    else
      longtable$mean <-
        sprintf(fmt, round(longtable$mean, digits))
    ##
    longtable$lcl <-
      sprintf(fmt, round(longtable$lcl, digits))
    ##
    longtable$ucl <-
      sprintf(fmt, round(longtable$ucl, digits))
    ##
    longtable <-
      longtable[longtable$Treatment != longtable$Comparator, ,
                drop = FALSE]
    ##
    longtable <- as.matrix(longtable)
    row.names(longtable) <- rep("", nrow(longtable))
    cl <- colnames(longtable)
    cl[cl %in% c("lcl", "ucl")] <- paste0(100 * quant[c(1, 3)], "%")
    colnames(longtable) <- cl
    ##
    class(longtable) <- c("league.crossnma", "matrix")
    attr(longtable, "direction") <- direction
    ##
    return(longtable)
  }


  invisible(NULL)
}





#' @rdname league.crossnma
#' @export league


league <- function(x, ...)
  UseMethod("league")





#' @rdname league.crossnma
#' @method print league.crossnma
#' @export


print.league.crossnma <- function(x, ...) {
  chkclass(x, "league.crossnma")
  ##
  if (attr(x, "direction") == "wide") {
    class(x) <- "matrix"
    rownames(x) <- rep("", nrow(x))
    colnames(x) <- rep("", ncol(x))
    prmatrix(x, quote = FALSE, right = TRUE)
  }
  else if (attr(x, "direction") == "long")
    prmatrix(x, quote = FALSE, right = TRUE, ...)

  invisible(NULL)
}
