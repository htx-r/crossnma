#' League Table
#' @description Produces a league table that contain point estimates of relative effects
#' for all possible pairs of treatments along with 95\% credible intervals obtained with the quantile method.
#' @param x A \code{crossrun} object produced by running \code{crossnma.run()}.
#' @param central.tdcy The statistic to be used to measure the relative treatment effects. The options are "mean" and "median".
#' @param exp If TRUE (default), odds ratios are displayed. If FALSE, log odds ratios will be presented.
#' @param order A vector of treatment names (character) representing the order in which to display these treatments.
#' @param prt.cov1.value  The participant covariate value of \code{cov1} for which to report the results. Must be specified for network meta-regression and when individual participant dataset is used in the analysis.
#' For dichotomous covariates, a character of the level (used in the data) should be indicated.
#' @param prt.cov2.value  The participant covariate value of \code{cov2} for which to report the results. Must be specified for network meta-regression and when individual participant dataset is used in the analysis.
#' For dichotomous covariates, a character of the level (used in the data) should be indicated.
#' @param prt.cov3.value  The participant covariate value of \code{cov3} for which to report the results. Must be specified for network meta-regression and when individual participant dataset is used in the analysis.
#' For dichotomous covariates, a character of the level (used in the data) should be indicated.
#' @param digits The number of digits to be used when displaying the results.
#' @param direction The format to display the league table. Two options "wide" (default) and "long".
#' @param \dots Additional arguments (ignored at the moment).
#'
#' @return \code{table}  A league table. Row names indicate comparator treatments.
#' The table will be displayed in a long or wide formatting.
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
#' #==========================#
#' # Display the league table #
#' #==========================#
#' league(fit,exp = TRUE)                   #  wide formatting
#' league(fit,exp = TRUE, direction="long") #  long formatting
#'
#' @seealso \code{\link{crossnma.run}}
#' 
#' @method league crossnma
#' @export


league.crossnma <- function(x,
                            central.tdcy = "median",
                            exp = FALSE,
                            order = NULL,
                            prt.cov1.value=NULL,
                            prt.cov2.value=NULL,
                            prt.cov3.value=NULL,
                            digits = 2,
                            direction="wide",
                            ...) {

  # Bind variables to function
  trt <- NULL
  Treatment <- NULL
  Comparator <- NULL
  cov.ref <- NULL

  if (class(x) != "crossnma")
    stop("\'x \' must be a valid crossnma object created using the crossnma.run function.")

  if(!is.null(x$model$covariate) & is.null(prt.cov1.value)) stop("prt.cov.value must be specified for network meta-regression")

  dmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- x$trt.key$trt.ini
  colnames(dmat) <- trt.names

if(!is.null(x$model$covariate)){
  if(x$model$data$ns.ipd!=0){ # if we provide IPD with or without AD
    nc <- length(x$model$covariate)
    nt <- length(trt.names)
    if (x$model$split.regcoef){
      #** betaw
      if(x$model$regw.effect=="independent"){
        bwmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("betaw.t_"))
        # split bwmat by nt_column to generate bwmat.cov for each covariate
        bwmat.cov1 <- bwmat[,1:nt]*ifelse(is.numeric(prt.cov1.value),
                                          prt.cov1.value-x$model$cov.ref[1], # for numeric covariate
                                          as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov1.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
        )
        if(nc==2){
          bwmat.cov2 <- bwmat[,(nt+1):(nt*2)]*ifelse(is.numeric(prt.cov2.value),
                                                     prt.cov2.value-x$model$cov.ref[2], # for numeric covariate
                                                     as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov2.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          )
        } else{
          bwmat.cov2 <- 0
        }
        if(nc==3){
          bwmat.cov3 <- bwmat[,(nt*2+1):(nt*3)]*ifelse(is.numeric(prt.cov3.value),
                                                       prt.cov3.value-x$model$cov.ref[3], # for numeric covariate
                                                       as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov3.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          )
        } else{
          bwmat.cov3 <- 0
        }
        dmat <- dmat+bwmat.cov1+bwmat.cov2+bwmat.cov3
      } else{
        bwmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("bw_"))
        bwmat.cov1 <- sweep(cbind(bwmat[,1]),MARGIN = 2,ifelse(is.numeric(prt.cov1.value),
                                                               prt.cov1.value-x$model$cov.ref[1], # for numeric covariate
                                                               as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov1.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
        ),'*') # multiply each column values by prt.cov-cov.ref
        bwmat.cov1 <-matrix(unlist(rep(bwmat.cov1,each=nt)), ncol=nt,byrow = TRUE)       # repeat the column for each trt to add it to dmat
        if(nc==2){
          bwmat.cov2 <- sweep(cbind(bwmat[,2]),MARGIN = 2,ifelse(is.numeric(prt.cov2.value),
                                                                 prt.cov2.value-x$model$cov.ref[2], # for numeric covariate
                                                                 as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov2.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          ),'*')
          bwmat.cov2 <- matrix(unlist(rep(bwmat.cov2,each=nt)), ncol=nt,byrow = TRUE)
        } else{
          bwmat.cov2 <- 0
        }
        if(nc==3){
          bwmat.cov3 <- sweep(cbind(bwmat[,3]),MARGIN = 2,ifelse(is.numeric(prt.cov3.value),
                                                                 prt.cov3.value-x$model$cov.ref[3], # for numeric covariate
                                                                 as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov3.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          ),'*')
          bwmat.cov3 <- matrix(unlist(rep(bwmat.cov3,each=nt)), ncol=nt,byrow = TRUE)
        } else{
          bwmat.cov3 <- 0
        }
        dmat <- dmat+ bwmat.cov1+bwmat.cov2+bwmat.cov3
      }
      # betab
      if(x$model$regb.effect=="independent"){
        bbmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("betab.t_"))
        # split bbmat by nt_column to generate bbmat.cov for each covariate
        stds.mean1 <- mean(c(x$model$data$xm1.ad, x$model$data$xm1.ipd),na.rm = TRUE)
        bbmat.cov1 <- bbmat[,1:nt]*ifelse(is.numeric(prt.cov1.value),stds.mean1-x$model$cov.ref[1],stds.mean1)
        if(nc==2){
          stds.mean2 <- mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd),na.rm = TRUE)
          bbmat.cov2 <- bbmat[,(nt+1):(nt*2)]*ifelse(is.numeric(prt.cov2.value),stds.mean2-x$model$cov.ref[2],stds.mean2)
        } else{
          bbmat.cov2 <- 0
        }
        if(nc==3){
          stds.mean3 <- mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd),na.rm = TRUE)
          bbmat.cov3 <- bbmat[,(nt*2+1):(nt*3)]*ifelse(is.numeric(prt.cov3.value),stds.mean3-x$model$cov.ref[3],stds.mean3)
        } else{
          bbmat.cov3 <- 0
        }
        dmat <- dmat+bbmat.cov1+bbmat.cov2+bbmat.cov3
      } else{
        stds.mean <- c(mean(c(x$model$data$xm1.ad, x$model$data$xm1.ipd),na.rm = TRUE),
                       mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd),na.rm = TRUE),
                       mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd),na.rm = TRUE))

        bbmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("bb_"))
        bbmat.cov1 <- sweep(cbind(bbmat[,1]),MARGIN = 2,ifelse(is.numeric(prt.cov1.value),stds.mean[1]-x$model$cov.ref[1],stds.mean[1]),'*')
        bbmat.cov1 <-matrix(unlist(rep(bbmat.cov1,each=nt)), ncol=nt,byrow = TRUE)
        if(nc==2){
          bbmat.cov2 <- sweep(cbind(bbmat[,2]),MARGIN = 2,ifelse(is.numeric(prt.cov2.value),stds.mean[2]-x$model$cov.ref[2],stds.mean[2]),'*')
          bbmat.cov2 <- matrix(unlist(rep(bbmat.cov2,each=nt)), ncol=nt,byrow = TRUE)
        } else{
          bbmat.cov2 <- 0
        }
        if(nc==3){
          bbmat.cov3 <- sweep(cbind(bbmat[,3]),MARGIN = 2,ifelse(is.numeric(prt.cov3.value),stds.mean[3]-x$model$cov.ref[3],stds.mean[3]),'*')
          bbmat.cov3 <- matrix(unlist(rep(bbmat.cov3,each=nt)), ncol=nt,byrow = TRUE)
        } else{
          bbmat.cov3 <- 0
        }
        dmat <- dmat+ bbmat.cov1+bbmat.cov2+bbmat.cov3
      }

    } else{
      if(x$model$regb.effect=="independent"&&x$model$regw.effect=="independent"){
        bmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("beta.t_"))
        bmat.cov1 <- bmat[,1:nt]*(ifelse(is.numeric(prt.cov1.value),
                                         prt.cov1.value-x$model$cov.ref[1], # for numeric covariate
                                         as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov1.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
        ))
        if(nc==2){
          bmat.cov2 <- bmat[,(nt+1):(nt*2)]*(ifelse(is.numeric(prt.cov2.value),
                                                    prt.cov2.value-x$model$cov.ref[2], # for numeric covariate
                                                    as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov2.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          ))

        } else{
          bmat.cov2 <- 0
        }
        if(nc==3){
          bmat.cov3 <- bmat[,(nt*2+1):(nt*3)]*(ifelse(is.numeric(prt.cov3.value),
                                                      prt.cov3.value-x$model$cov.ref[3], # for numeric covariate
                                                      as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov3.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          ))
        } else{
          bmat.cov3 <- 0
        }
        dmat <- dmat+bmat.cov1+bmat.cov2+bmat.cov3
      } else{
        bmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("b_"))
        bmat.cov1 <- sweep(cbind(bmat[,1]),MARGIN = 2,ifelse(is.numeric(prt.cov1.value),
                                                             prt.cov1.value-x$model$cov.ref[1], # for numeric covariate
                                                             as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov1.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
        ),'*')
        bmat.cov1 <-matrix(unlist(rep(bmat.cov1,each=nt)), ncol=nt,byrow = TRUE)
        if(nc==2){
          bmat.cov2 <- sweep(cbind(bmat[,2]),MARGIN = 2,ifelse(is.numeric(prt.cov2.value),
                                                               prt.cov2.value-x$model$cov.ref[2], # for numeric covariate
                                                               as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov2.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          ),'*')
          bmat.cov2 <- matrix(unlist(rep(bmat.cov2,each=nt)), ncol=nt,byrow = TRUE)
        } else{
          bmat.cov2 <- 0
        }
        if(nc==3){
          bmat.cov3 <- sweep(cbind(bmat[,3]),MARGIN = 2,ifelse(is.numeric(prt.cov3.value),
                                                               prt.cov3.value-x$model$cov.ref[3], # for numeric covariate
                                                               as.numeric(x$model$dich.cov.labels[x$model$dich.cov.labels[,1]==prt.cov3.value,2]) # for factor covariate, multiply by 0 or 1 depends on what value the user indicate in prt.cov.value
          ),'*')
          bmat.cov3 <- matrix(unlist(rep(bmat.cov3,each=nt)), ncol=nt,byrow = TRUE)
        } else{
          bmat.cov3 <- 0
        }
        dmat+ bmat.cov1+bmat.cov2+bmat.cov3
      }
    }
  } else{ # only AD
    # betab
    if(x$model$regb.effect=="independent"){
      bbmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("beta.t_"))
      # split bbmat by nt_column to generate bbmat.cov for each covariate
      stds.mean1 <-x$model$data$xm1.ad
      bbmat.cov1 <- bbmat[,1:nt]*ifelse(!is.na(cov.ref[1]),stds.mean1-x$model$cov.ref[1],stds.mean1)
      if(nc==2){
        stds.mean2 <- mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd),na.rm = TRUE)
        bbmat.cov2 <- bbmat[,(nt+1):(nt*2)]*ifelse(!is.na(cov.ref[2]),stds.mean2-x$model$cov.ref[2],stds.mean2)
      } else{
        bbmat.cov2 <- 0
      }
      if(nc==3){
        stds.mean3 <- mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd),na.rm = TRUE)
        bbmat.cov3 <- bbmat[,(nt*2+1):(nt*3)]*ifelse(!is.na(cov.ref[3]),stds.mean3-x$model$cov.ref[3],stds.mean3)
      } else{
        bbmat.cov3 <- 0
      }
      dmat <- dmat+bbmat.cov1+bbmat.cov2+bbmat.cov3
    } else{
      stds.mean <- c(x$model$data$xm1.ad,
                     x$model$data$xm2.ad,
                     x$model$data$xm3.ad
                     )

      bbmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(starts_with("b_"))
      bbmat.cov1 <- sweep(cbind(bbmat[,1]),MARGIN = 2,ifelse(!is.na(cov.ref[1]),stds.mean[1]-x$model$cov.ref[1],stds.mean[1]),'*')
      bbmat.cov1 <-matrix(unlist(rep(bbmat.cov1,each=nt)), ncol=nt,byrow = TRUE)
      if(nc==2){
        bbmat.cov2 <- sweep(cbind(bbmat[,2]),MARGIN = 2,ifelse(!is.na(cov.ref[2]),stds.mean[2]-x$model$cov.ref[2],stds.mean[2]),'*')
        bbmat.cov2 <- matrix(unlist(rep(bbmat.cov2,each=nt)), ncol=nt,byrow = TRUE)
      } else{
        bbmat.cov2 <- 0
      }
      if(nc==3){
        bbmat.cov3 <- sweep(cbind(bbmat[,3]),MARGIN = 2,ifelse(!is.na(cov.ref[3]),stds.mean[3]-x$model$cov.ref[3],stds.mean[3]),'*')
        bbmat.cov3 <- matrix(unlist(rep(bbmat.cov3,each=nt)), ncol=nt,byrow = TRUE)
      } else{
        bbmat.cov3 <- 0
      }
      dmat <- dmat+ bbmat.cov1+bbmat.cov2+bbmat.cov3
    }
  }
}


  # # when we adjust for covariate, add a message with the values that we adjusted for
  # if(!is.null(prt.cov.value)){
  #   msg.adjust <-  paste0("The estimates of treatment effect are adjusted for \n
  #                               participant ",x$model$covariate[1],"= ",prt.cov1.value,
  #                          "and study mean = ",round(mean(c(x$model$data$xm1.ad, x$model$data$xm1.ipd),na.rm = TRUE),2)
  #   )
  #   if(nc==2) msg.adjust <- paste0(msg.adjust, paste0("\n participant ",x$model$covariate[2],"= ",prt.cov2.value,
  #                                       "and study mean = ",round(mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd),na.rm = TRUE),2)))
  #   if(nc==3) msg.adjust <- paste0(msg.adjust, paste0("\n participant ",x$model$covariate[3],"= ",prt.cov3.value,
  #                                       "and study mean = ",round(mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd),na.rm = TRUE),2))
  #                    )
  # }


  if(is.null(order)){
    order <- trt.names
  }

  dmat %<>% select(order)
  trt.names <- order

  # useful functions to compute some statistics
  calc.report <- function(x, fct="identity", arg=NULL, trans="identity"){
    if(is.null(arg)){
      eval(call(fct,(call(trans, x))))
    } else
      eval(call(fct,(call(trans, x)), arg))
  }
  exp.mean <- function(x) calc.report(x, "mean", trans="exp")
  exp.median <- function(x) calc.report(x, "median", trans="exp")
  exp.sd <- function(x) calc.report(x, "sd", trans="exp")
  exp.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="exp")
  exp.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="exp")

  id.mean <- function(x) calc.report(x, "mean", trans="identity")
  id.median <- function(x) calc.report(x, "median", trans="identity")
  id.sd <- function(x) calc.report(x, "sd", trans="identity")
  id.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="identity")
  id.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="identity")

  colvals <- function(dmat, b.col=1, paste=TRUE) {

    # Bind variables to function
    key <- NULL
    value <- NULL
    trt <- NULL
    estimate <- NULL
    lci <- NULL
    uci <- NULL
    result <- NULL

    base <- colnames(dmat)[b.col]

    dmat2 <- dmat
    new.vars <- paste0(colnames(dmat2), "-", b.col)
    for(i in 1:ncol(dmat)) {
      dmat2[[new.vars[i]]] <- dmat[, i] - dmat[, b.col]
    }
    dmat2 %<>% select(new.vars)
    colnames(dmat2) <- trt.names

    if(central.tdcy=="mean" & exp==TRUE){
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = exp.mean)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    } else if(central.tdcy=="mean"& exp==FALSE){
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = id.mean)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))}
    if(central.tdcy=="median" & exp==TRUE){
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = exp.median))%>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    } else if(central.tdcy=="median"& exp==FALSE){
      tmp.estimate <- dmat2 %>%
        summarise_all(list(estimate = id.median)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    }

    if(exp==TRUE ){
      tmp.lci <- dmat2 %>%
        summarise_all(list(lci = exp.lci)) %>% gather() %>%
        rename(trt = key, lci = value) %>%
        mutate(trt = sub("_lci", "", trt))

      tmp.uci <- dmat2 %>%
        summarise_all(list(uci = exp.uci)) %>% gather() %>%
        rename(trt = key, uci = value) %>%
        mutate(trt = sub("_uci", "", trt))

      null.value <- 1
    } else{
      tmp.lci <- dmat2 %>%
        summarise_all(list(lci = id.lci)) %>% gather() %>%
        rename(trt = key, lci = value) %>%
        mutate(trt = sub("_lci", "", trt))

      tmp.uci <- dmat2 %>%
        summarise_all(list(uci = id.uci)) %>% gather() %>%
        rename(trt = key, uci = value) %>%
        mutate(trt = sub("_uci", "", trt))

      null.value <- 0
    }

    #create C-style formatting string from the digits parameter
    fmt <- paste0("%.", digits, "f")

    if(paste){
      tmp1 <- left_join(tmp.estimate, tmp.lci, by = "trt") %>%
        left_join(tmp.uci, by = "trt") %>%
        mutate(result = paste(sprintf(fmt, estimate),
                              " (",
                              sprintf(fmt, lci),
                              " to ",
                              sprintf(fmt, uci),
                              ")", sep="")) %>%
        select(trt, result)

      colnames(tmp1)[2] <- as.character(tmp.estimate %>% filter(estimate==null.value) %>% select(trt))
    } else{
      tmp1 <- left_join(tmp.estimate, tmp.lci, by = "trt") %>%
        left_join(tmp.uci, by = "trt")

      colnames(tmp1)[2] <- central.tdcy
    }

    return(tmp1)
  }


  #Default layout

  tmp1.list <- list()
  for(i in 1:ncol(dmat)) {
    tmp1.list[[i]] <- colvals(dmat, b.col=i)
  }

  tmp.df <- suppressMessages(bind_cols(tmp1.list)) %>%
    select(-starts_with("trt")) %>%
    t()

  colnames(tmp.df) <- colnames(dmat)
  rownames(tmp.df) <- colnames(dmat)

  for(i in 1:dim(tmp.df)[1]){
    tmp.df[i,i] <- colnames(tmp.df)[i]
  }

  default <- tmp.df

  #long layout

  tmp2.list <- list()
  for(i in 1:ncol(dmat)) {
    tmp2.list[[i]] <- colvals(dmat, b.col=i, paste=FALSE)
  }

  longtable <- tmp2.list %>%
    bind_rows() %>%
    mutate(Treatment = trt,
           Comparator = rep(trt.names, each=length(trt.names))) %>%
    select(Treatment, Comparator, everything(), -trt)

  if(exp==TRUE ){ # & x$link!="identity"
    null.value <- 1
  } else{
    null.value <- 0
  }
if(direction=="wide"){
  tab <- default
} else if(direction=="long"){
tab <- longtable
} else {
  stop("direction takes values 'wide' or 'long' ")
}

  return(tab)

}





#' @rdname league.crossnma
#' @export league


league <- function(x, ...)
  UseMethod("league")
