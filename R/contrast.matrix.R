contrast.matrix <- function(x,
                            median = FALSE,
                            order = NULL,
                            cov1.value = NULL,
                            cov2.value = NULL,
                            cov3.value = NULL,
                            ...) {
  
  chkclass(x, "crossnma")
  ##
  chklogical(median)
  
  if (median)
    central.tdcy <- "median"
  else
    central.tdcy <- "mean"
  
  
  if (!is.null(x$model$covariate) & is.null(cov1.value))
    stop("cov1.value must be specified for network meta-regression")


  ## Bind variables to function
  ##
  trt <- Treatment <- Comparator <- cov.ref <- NULL
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
            c(mean(c(x$model$data$xm1.ad, x$model$data$xm1.ipd), na.rm = TRUE),
              mean(c(x$model$data$xm2.ad, x$model$data$xm2.ipd), na.rm = TRUE),
              mean(c(x$model$data$xm3.ad, x$model$data$xm3.ipd), na.rm = TRUE))
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
          if (!is.na(cov.ref[1]))
            stds.mean1 - x$model$cov.ref[1]
        else
          stds.mean1
        ##
        if (nc == 2) {
          stds.mean2 <- x$model$data$xm2.ad
          bbmat.cov2 <-
            bbmat[, (nt + 1):(nt * 2)] *
            if (!is.na(cov.ref[2]))
              stds.mean2 - x$model$cov.ref[2]
          else
            stds.mean2
        }
        ##
        if (nc == 3) {
          stds.mean3 <- x$model$data$xm3.ad
          bbmat.cov3 <- bbmat[, (nt * 2 + 1):(nt * 3)] *
            if (!is.na(cov.ref[3]))
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
                if (!is.na(cov.ref[1]))
                  stds.mean[1] - x$model$cov.ref[1]
                else stds.mean[1],
                '*')
        bbmat.cov1 <-
          matrix(unlist(rep(bbmat.cov1, each = nt)), ncol = nt, byrow = TRUE)
        ##
        if (nc == 2) {
          bbmat.cov2 <-
            sweep(cbind(bbmat[, 2]), MARGIN = 2,
                  if (!is.na(cov.ref[2]))
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
                  if (!is.na(cov.ref[3]))
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
  
  
  colvals <- function(dmat, b.col = 1, paste = TRUE) {
    ## Bind variables to function
    key <- value <- trt <- estimate <- lcl <- ucl <- result <- NULL
    ##
    dmat2 <- dmat
    new.vars <- paste0(colnames(dmat2), "-", b.col)
    for (i in 1:ncol(dmat)) {
      dmat2[[new.vars[i]]] <- dmat[, i] - dmat[, b.col]
    }
    dmat2 %<>% select(new.vars)
    colnames(dmat2) <- trt.names
    
    if (central.tdcy == "mean") {
      TE <- dmat2 %>%
        summarise_all(list(estimate = mean)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    }
    else {
      TE <- dmat2 %>%
        summarise_all(list(estimate = median)) %>% gather() %>%
        rename(trt = key, estimate = value) %>%
        mutate(trt = sub("_estimate", "", trt))
    }
    TE <- dmat2 %>%
      summarise_all(list(estimate =
                           if (central.tdcy == "median")
                             median
                           else
                             mean)) %>%
      gather() %>%
      rename(trt = key, estimate = value) %>%
      mutate(trt = sub("_estimate", "", trt))
    ##
    sdTE <- dmat2 %>%
      summarise_all(list(sd = sd)) %>% gather() %>%
      rename(trt = key, sd = value) %>%
      mutate(trt = sub("_sd", "", trt))
    ##
    return(list(TE = TE, sdTE = sdTE))
  }
  
  
  ##
  ##
  ## Return two matrices for every relative treatment effect and standard deviation
  ##
  ##

  n <- ncol(dmat)
  ##
  TE.list <- list("list", n)
  ##
  for (i in seq_len(n))
    TE.list[[i]] <- colvals(dmat, b.col = i)[["TE"]]
  ##
  TE <- suppressMessages(bind_cols(TE.list)) %>%
    select(-starts_with("trt")) %>% t()
  ##
  colnames(TE) <- colnames(dmat)
  rownames(TE) <- colnames(dmat)
  ##
  TE[row(TE) == col(TE)] <- 0
  ##
  class(TE) <- "matrix"
  
  sdTE.list <- list("list", n)
  ##
  for (i in seq_len(n))
    sdTE.list[[i]] <- colvals(dmat, b.col = i)[["sdTE"]]
  ##
  sdTE <- suppressMessages(bind_cols(sdTE.list)) %>%
    select(-starts_with("trt")) %>% t()
  ##
  colnames(sdTE) <- colnames(dmat)
  rownames(sdTE) <- colnames(dmat)
  ##
  sdTE[row(sdTE) == col(sdTE)] <- 0
  ##
  class(sdTE) <- "matrix"
  
  
  return(list(TE = t(TE), sdTE = t(sdTE)))
}
