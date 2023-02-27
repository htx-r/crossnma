crossnma.code <- function(ipd = TRUE,
                          ad = TRUE,
                          sm = NULL,
                          max.d = NULL,
                          trt.effect = "random",
                          prior.tau.trt = NULL,
                          ## -------- meta-regression
                          split.regcoef = FALSE,
                          covariate = NULL,
                          reg0.effect = "random",
                          regb.effect = "random",
                          regw.effect = "random",
                          prior.tau.reg0 = NULL,
                          prior.tau.regb = NULL,
                          prior.tau.regw = NULL,
                          ## --------  bias adjustment
                          bias.effect = NULL,
                          bias.type = NULL,
                          bias.covariate = NULL,
                          add.std.in = NULL,
                          add.std.act.no = NULL,
                          add.std.act.yes = NULL,
                          ##
                          prior.tau.gamma = NULL,
                          v = NULL,
                          prior.pi.high.rct = NULL,
                          prior.pi.low.rct = NULL,
                          prior.pi.high.nrs = NULL,
                          prior.pi.low.nrs = NULL,
                          ##
                          method.bias = NULL,
                          d.prior.nrs = NULL # required when method.bias = "prior"
                          ) {


  dnmax <- sprintf(" ~ dnorm(0, (%s * 15)^(-2))", max.d)
  dumax <- sprintf("dunif(0, %s)", max.d)
  ##
  sel.i <- "[j, t.ipd[j, k]]"
  sel.a <- "[j, t.ad[j, k]]"
  
  
  ##
  ##
  ##  Priors
  ##
  ##
  
  ##
  ## Trial baseline effects
  ##
  if (sm == "RR")
    prior.u <-
      paste0("for (j in 1:(ns.ipd + ns.ad)) {",
             "\n  u[j] <- log(p.b[j])",
             "\n  p.b[j] ~ dunif(0, 1)",
             "\n}")
  else
    prior.u <-
      paste0("for (j in 1:(ns.ipd + ns.ad)) {",
             "\n  u[j]", dnmax,
             "\n}")
  
  
  ##
  ## Construct default priors following basic paramaters
  ##
  d.prior <-
    paste0("\nfor (k in 2:nt) {",
           "\n  d[k]", dnmax,
           "\n}")
  ##
  ## Heterogeneity parameters
  ##
  prior.tau.trt <- replaceNULL(prior.tau.trt, dumax)
  prior.tau.reg0 <- replaceNULL(prior.tau.reg0, dumax)
  prior.tau.regb <- replaceNULL(prior.tau.regb, dumax)
  prior.tau.regw <- replaceNULL(prior.tau.regw, dumax)
  prior.tau.gamma <- replaceNULL(prior.tau.gamma, dumax)
  ##
                                        # Bias probabilities (for adjust1 and adjust2)
  ##
  prior.pi.high.rct <- replaceNULL(prior.pi.high.rct, "dbeta(10, 1)")
  prior.pi.low.rct  <- replaceNULL(prior.pi.low.rct, "dbeta(1, 10)")
  prior.pi.high.nrs <- replaceNULL(prior.pi.high.nrs, "dbeta(30, 1)")
  prior.pi.low.nrs  <- replaceNULL(prior.pi.low.nrs, "dbeta(1, 30)")
  
  
  ##
  ##
  ## Meta-regression
  ##
  ##
  
  beta0.prior.ipd <- betab.prior <- betaw.prior.ipd <-
    beta.prior.ipd <- beta.prior.ad <- ""
  mreg.ipd <- mreg.ad <- ""
  betab.consis.ipd <- betaw.consis.ipd <- betab.consis.ad <- ""
  ##
  if (!is.null(covariate)) {
    ##
    ## IPD
    ##
    if (ipd) {
      for (i in 1:length(covariate[[1]])) {
        ## Meta-regression terms
        mreg.ipd0 <-
          paste0(
            " + beta0_", i, "[study[i]] * (x", i, "[i]) + betaw_", i,
            "[study[i], trt[i]] * (x", i,
            "[i] - xm", i, ".ipd[i]) + betab_", i,
            "[study[i], trt[i]] * xm", i, ".ipd[i]")
        ##
        mreg.ipd <- paste0(mreg.ipd, mreg.ipd0)
        ## consistency equations for beta_b and beta_w - up to 3
        betab.consis.ipd0 <-
          paste0("\n    betab_", i, sel.i, " <- betab.t_", i,
                 "[t.ipd[j, k]] - betab.t_", i, "[t.ipd[j, 1]]")
        betaw.consis.ipd0 <-
          paste0("\n    betaw_", i, sel.i, " <- betaw.t_", i,
                 "[t.ipd[j, k]] - betaw.t_", i, "[t.ipd[j, 1]]")
        ##
        betab.consis.ipd <- paste0(betab.consis.ipd, betab.consis.ipd0)
        betaw.consis.ipd <- paste0(betaw.consis.ipd, betaw.consis.ipd0)
      }

      ##
      ## beta0
      ##
      if (reg0.effect == "random") {
        for (i in 1:length(covariate[[1]])) {
          beta0.prior.ipd0 <-
            paste0("\n# Random effect for beta0",
                   "\nfor (j in 1:(ns.ipd)) {",
                   "\n  beta0_", i, "[j] ~ dnorm(b0_", i, ", prec.beta0_", i, ")",
                   "\n}",
                   "\nb0_", i, dnmax,
                   "\nprec.beta0_", i, " <- pow(tau.b0_", i, ", -2)",
                   "\ntau.b0_", i, " ~ ", prior.tau.reg0)
          ##
          beta0.prior.ipd <- paste0(beta0.prior.ipd, beta0.prior.ipd0)
        }
      }
      else if (reg0.effect == "independent") {
        for (i in 1:length(covariate[[1]])) {
          beta0.prior.ipd0 <-
            paste0("\n# Independent effect for beta0",
                   "\nfor (j in 1:(ns.ipd)) {",
                   "\n  beta0_", i, "[j] <- b0_", i, "[j]",
                   "\n  b0_", i, "[j]", dnmax,
                   "\n}")
          ##
          beta0.prior.ipd <- paste0(beta0.prior.ipd, beta0.prior.ipd0)
        }
      }
      else
        stop("The progonostic effect can be assumed either ",
             "'independent' or 'random' across studies")
      
      
      ##
      ## betab and betaw
      ##
      if (!split.regcoef) { # not splitted within and between-study covariate
        if (regb.effect == "independent" || regw.effect == "independent") {
          beta.prior.ipd0 <-
            paste0("\n# Independent effect for beta (within = between)",
                   "\nbeta.t_", i, "[1] <- 0",
                   "\nfor (k in 1:nt) {",
                   "\n  betab.t_", i, "[k] <- beta.t_", i, "[k]",
                   "\n  betaw.t_", i, "[k] <- beta.t_", i, "[k]",
                   "\n}",
                   "\nfor (k in 2:nt) {",
                   "\n  beta.t_", i, "[k]", dnmax,
                   "\n}")
          ##
          beta.prior.ipd <- paste0(beta.prior.ipd, beta.prior.ipd0)
        }
        else if (regb.effect == "random" || regw.effect == "random") {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ipd0 <-
              paste0("\n# Random effects for beta (within = between)",
                     "\nbeta.t_", i, "[1] <- 0",
                     "\nfor (k in 1:nt) {",
                     "\n  betab.t_", i, "[k] <- beta.t_", i, "[k]",
                     "\n  betaw.t_", i, "[k] <- beta.t_", i, "[k]",
                     "\n}",
                     "\nfor (k in 2:nt) {",
                     "\n  beta.t_", i, "[k] ~ dnorm(b_", i, ", prec.beta_", i, ")",
                     "\n}",
                     "\nb_", i, dnmax,
                     "\ntau.b_", i, " ~ ", prior.tau.regw,
                     "\nprec.beta_", i, " <- pow(tau.b_", i, ", -2)")
            ##
            beta.prior.ipd <- paste0(beta.prior.ipd, beta.prior.ipd0)
          }
        }
        else if (regb.effect == "common" & regw.effect == "common") {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ipd0 <-
              paste0("\n# Common effect for beta (within = between)",
                     "\nbetab.t_", i, "[1] <- 0",
                     "\nbetaw.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\n  betab.t_", i, "[k] <- b_", i,
                     "\n  betaw.t_", i, "[k] <- b_", i,
                     "\n}",
                     "\nb_", i, dnmax)
            ##
            beta.prior.ipd <- paste0(beta.prior.ipd, beta.prior.ipd0)
          }
        }
        else
          stop("The regb.effect and regw.effect need to both be assumed ",
               "'random' or 'common' across studies")
      }
      else {
        ## splitted within and between- study covariate between- study
        ## covariate
        if (regb.effect == "independent") {
          betab.prior0 <-
            paste0("\n# Random effect for betab ",
                   "(the between-study covariate effect)",
                   "\nbetab.t_", i, "[1] <- 0",
                   "\nfor (k in 2:nt) {",
                   "\n  betab.t_", i, "[k]", dnmax,
                   "\n}")
          ##
          betab.prior <- paste0(betab.prior, betab.prior0)
        }
        else if (regb.effect == "random") {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
              paste0("\n# Random effect for betab ",
                     "(the between-study covariate effect)",
                     "\nbetab.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\n  betab.t_", i, "[k] ~ dnorm(bb_", i, ", ",
                     "prec.betab_", i, ")",
                     "\n}",
                     "\nbb_", i, dnmax,
                     "\ntau.bb_", i, " ~ ", prior.tau.regb,
                     "\nprec.betab_", i, " <- pow(tau.bb_", i, ", -2)")
            ##
            betab.prior <- paste0(betab.prior, betab.prior0)
          }
        }
        else if (regb.effect == "common") {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
              paste0("\n# Common effect for betab ",
                     "(the between-study covariate effect)",
                     "\nbetab.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\n  betab.t_", i, "[k] <- bb_", i,
                     "\n}",
                     "\nbb_", i, dnmax)
            ##
            betab.prior <- paste0(betab.prior, betab.prior0)
          }
        }
        else
          stop("The between-study covariate effect need to be assumed ",
               "'independent', 'random' or 'common' across studies")


        ##
        ## Within-study covariate
        ##
        if (regw.effect == "independent") {
          betaw.prior.ipd0 <-
            paste0("\n# Random effect for betab ",
                   "(the between-study covariate effect)",
                   "\nbetaw.t_", i, "[1] <- 0",
                   "\nfor (k in 2:nt) {",
                   "\n  betaw.t_", i, "[k]", dnmax,
                   "\n}")
          ##
          betaw.prior.ipd <- paste0(betaw.prior.ipd, betaw.prior.ipd0)
        }
        else if (regw.effect == "random") {
          for (i in 1:length(covariate[[1]])) {
            betaw.prior.ipd0 <-
              paste0("\n# Random effect for betaw ",
                     "(the within-study covariate effect)",
                     "\nbetaw.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\n  betaw.t_", i, "[k] ~ dnorm(bw_", i, ", ",
                     "prec.betaw_", i, ")",
                     "\n}",
                     "\nbw_", i, dnmax,
                     "\nprec.betaw_", i, " <- pow(tau.bw_", i, ", -2)",
                     "\ntau.bw_", i, " ~ ", prior.tau.regw)
            ##
            betaw.prior.ipd <- paste0(betaw.prior.ipd, betaw.prior.ipd0)
          }
        }
        else if (regw.effect == "common") {
          for (i in 1:length(covariate[[1]])) {
            betaw.prior.ipd0 <-
              paste0("\n# Common effect for betaw ",
                     "(the within-study covariate effect)",
                     "betaw.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\n  betaw.t_", i, "[k] <- bw_", i,
                     "\n}",
                     "\nbw_", i, dnmax)
            ##
            betaw.prior.ipd <- paste0(betaw.prior.ipd, betaw.prior.ipd0)
          }
        }
        else
          stop("The within-study covariate effect can be assumed ",
               "'independent', 'random' or 'common' across studies.")
      }
    }
    
    
    ##
    ## AD
    ##
    if (ad) {
      for (i in 1:length(covariate[[1]])) {
        ## meta-regression terms - up to 3
        mreg.ad0 <-
          paste0(" + betab.ad_", i, sel.a, " * xm", i, ".ad[j]")
        betab.consis.ad0 <-
          paste0("\n    betab.ad_", i,
                 sel.a, " <- betab.t_", i,
                 "[t.ad[j, k]] - betab.t_", i, "[t.ad[j, 1]]")
        
        ## consistency equation - up to 3
        mreg.ad <- paste0(mreg.ad, mreg.ad0)
        betab.consis.ad <- paste0(betab.consis.ad, betab.consis.ad0)
      }
      if (!split.regcoef) { # not splitted
        if (regb.effect == "independent" && regw.effect == "independent") {
          beta.prior.ad0 <-
            paste0("\n# Random effect for beta (within = between)",
                   "\nbeta.t_", i, "[1] <- 0",
                   "\nfor (k in 1:nt) {",
                   "\n  betab.t_", i, "[k] <- beta.t_", i, "[k]",
                   "\n}",
                   "\nfor (k in 2:nt) {",
                   "\n  beta.t_", i, "[k]", dnmax,
                   "\n}")
          ##
          beta.prior.ad <- paste0(beta.prior.ad, beta.prior.ad0)
        }
        else if (regb.effect == "random" && regw.effect == "random") {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ad0 <-
              paste0("\n# Random effect for beta (within = between)",
                     "\nbeta.t_", i, "[1] <- 0",
                     "\nfor (k in 1:nt) {",
                     "\n  betab.t_", i, "[k] <- beta.t_", i, "[k]",
                     "\n}",
                     "\nfor (k in 2:nt) {",
                     "\n  beta.t_", i, "[k] ~ dnorm(b_", i, ", prec.beta_",
                     i, ")",
                     "\n}",
                     "\nb_", i, dnmax,
                     "\ntau.b_", i, " ~ ", prior.tau.regb,
                     "\nprec.beta_", i, " <- pow(tau.b_", i, ", -2)")
            ##
            beta.prior.ad <- paste0(beta.prior.ad, beta.prior.ad0)
          }
        }
        else if (regb.effect == "common" & regw.effect == "common") {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ad0 <-
              paste0("\n# Common effect for beta (within = between)",
                     "\nbetab.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\n  betab.t_", i, "[k] <- b_", i,
                     "\n}",
                     "\nb_", i, dnmax)
            ##
            beta.prior.ad <- paste0(beta.prior.ad, beta.prior.ad0)
          }
        }
        else
          stop("The regb.effect and regw.effect need to be assumed both ",
               "'independent', 'random' or 'common' across studies")
      }
      else { # splitted
        betab.prior <- ""
        if (regb.effect == "independent") {
          betab.prior0 <-
            paste0("\n# Random effect for betab ",
                   "(the between-study covariate effect)",
                   "\nbetab.t_", i, "[1] <- 0",
                   "\nfor (k in 2:nt) {",
                   "\n  betab.t_", i, "[k]", dnmax,
                   "\n}")
          ##
          betab.prior <- paste0(betab.prior, betab.prior0)
        }
        else if (regb.effect == "random") {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
              paste0("\n# Random effects for betab ",
                     "(the between-study covariate effect)",
                     "\nbetab.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\nbetab.t_", i, "[k] ~ dnorm(bb_", i, ", prec.betab_", i, ")",
                     "\n}",
                     "\nbb_", i, dnmax,
                     "\ntau.bb_", i, " ~ ", prior.tau.regb,
                     "\nprec.betab_", i, " <- pow(tau.bb_", i, ", -2)")
            ##
            betab.prior <- paste0(betab.prior, betab.prior0)
          }
        }
        else if (regb.effect == "common") {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
              paste0("\n# Random effects for betab ",
                     "(the between-study covariate effect)",
                     "\nbetab.t_", i, "[1] <- 0",
                     "\nfor (k in 2:nt) {",
                     "\nbetab.t_", i, "[k] <- bb_", i,
                     "\n}",
                     "\nbb_", i, dnmax)
            ##
            betab.prior <- paste0(betab.prior, betab.prior0)
          }
        }
        else
          stop("The between-study covariate effect need to be assumed ",
               "'independent', 'random' or 'common' across studies")
      }
    }
  }
  ##
  if (!split.regcoef) {
    if (ipd)
      beta.prior <- beta.prior.ipd
    else
      beta.prior <- beta.prior.ad
  }
  else
    beta.prior <- ""
  
  
  ##
  ## Treatment effect is zero for reference treatment
  ##
  ref.trt.effect.ipd <- ""
  ref.trt.effect.ad <- ""
  ##  
  if (!is.null(covariate)) {
    if (ipd) {
      for (i in 1:length(covariate[[1]])) {
        ## meta-regression terms - Up to 3
        ref.trt.effect.ipd0 <-
          paste0("\n  betaw_", i, "[j, t.ipd[j, 1]] <- 0",
                 "\n  betab_", i, "[j, t.ipd[j, 1]] <- 0")
        ##
        ref.trt.effect.ipd <- paste0(ref.trt.effect.ipd, ref.trt.effect.ipd0)
      }
    }
    if (ad) {
      for (i in 1:length(covariate[[1]])) {
        ## meta-regression terms - up to 3
        ref.trt.effect.ad0 <-
          if (ipd)
            ""
          else
            paste0("\nbetab.ad_", i, "[j, t.ad[j, 1]] <- 0")
        ##
        ref.trt.effect.ad <- paste0(ref.trt.effect.ad, ref.trt.effect.ad0)
      }
    }
  }
  
  
  ##
  ##
  ## Relative treatment effects
  ##
  ##

  if (is.null(v))
    q.prior <- ""
  else
    q.prior <- paste0("q ~ dbeta(", v, ", 1)")
  
  
  if (trt.effect == "random") {
    if (!is.null(v)) {
      theta.effect.ipd <-
        paste0("\n    theta", sel.i, " <- (1 - R[j]) * ",
               "theta.un", sel.i, " + R[j] * theta.adj", sel.i, "",
               "\n    theta.un", sel.i, " ~ ",
               "dnorm(md", sel.i, ", precd", sel.i, ")",
               "\n    theta.adj", sel.i, " ~ ",
               "dnorm(md", sel.i, " + gamma[j], precd", sel.i, " / q)",
               "\n    # Multi-arm correction",
               "\n    md", sel.i, " <- mean[j, k] + sw[j, k]",
               "\n    w[j, k] <- theta", sel.i, " - mean[j, k]",
               "\n    sw[j, k] <- sum(w[j, 1:(k - 1)]) / (k - 1)",
               "\n    precd", sel.i, " <- prec * 2 * (k - 1) / k")
      theta.effect.ad <-
        paste0("\n    theta[j + ns.ipd, t.ad[j, k]] <- ",
               "(1 - R[j + ns.ipd]) * theta.un[j + ns.ipd, t.ad[j, k]] + ",
               "R[j + ns.ipd] * theta.adj[j + ns.ipd, t.ad[j, k]]",
               "\n    theta.un[j + ns.ipd, t.ad[j, k]] ~ ",
               "dnorm(md.ad", sel.a, ", precd.ad", sel.a, ")",
               "\n    theta.adj[j + ns.ipd, t.ad[j, k]] ~ ",
               "\n    dnorm(md.ad", sel.a, " + ",
               "gamma[(j + ns.ipd)], precd.ad", sel.a, " / q)",
               "\n    # Multi-arm correction",
               "\n    md.ad", sel.a, " <- mean.ad[j, k] + sw.ad[j, k]",
               "\n    w.ad[j, k] <- (theta[j + ns.ipd, t.ad[j, k]] - ",
               "mean.ad[j, k])",
               "\n    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)",
               "\n    precd.ad", sel.a, " <- prec * 2 * (k - 1) / k")
      ##
      prior.tau.theta <-
        paste0("\n# Heterogeneity between theta's\n tau ~ ",
               prior.tau.trt, "\n prec <- pow(tau, -2)")
    }
    else {
      theta.effect.ipd <-
        paste0("\n    theta", sel.i, " ~ dnorm(md", sel.i, ", precd", sel.i, ")",
               "\n    # Multi-arm correction",
               "\n    md", sel.i, " <- mean[j, k] + sw[j, k]",
               "\n    w[j, k] <- theta", sel.i, " - mean[j, k]",
               "\n    sw[j, k] <- sum(w[j, 1:(k - 1)]) / (k - 1)",
               "\n    precd", sel.i, " <- prec * 2 * (k - 1) / k")
      ##
      theta.effect.ad <-
        paste0("\n    theta[j + ns.ipd, t.ad[j, k]] ~ ",
               "dnorm(md.ad", sel.a, ", precd.ad", sel.a, ")",
               "\n    # Multi-arm correction",
               "\n    md.ad", sel.a, " <- mean.ad[j, k] + sw.ad[j, k]",
               "\n    w.ad[j, k] <- (theta[j + ns.ipd, t.ad[j, k]] - ",
               "mean.ad[j, k])",
               "\n    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)",
               "\n    precd.ad", sel.a, " <- prec * 2 * (k - 1) / k")
      ##
      prior.tau.theta <-
        paste0("\n# Heterogeneity between theta's\ntau ~ ",
               prior.tau.trt, "\nprec <- pow(tau, -2)")
    }
  }
  else if (trt.effect == "common") {
    theta.effect.ipd <-
      paste0("\n    theta", sel.i, " <- md", sel.i,
             "\n    md", sel.i, " <- mean[j, k]")
    theta.effect.ad <-
      paste0("\n    theta[j + ns.ipd, t.ad[j, k]] <- md.ad[j, t.ad[j, k]]",
             "\n    md.ad[j, t.ad[j, k]] <- mean.ad[j, k]")
    ##
    prior.tau.theta <- ""
  }
  else
    stop("Please indicate the treatment effect model as either ",
         "'random' or 'common' ")
  
  
  ##
  ##
  ## Adjust for NRS
  ##
  ##

  adj.ipd <- ""
  adj.ad <- ""
  adjust.prior <- ""

  if (!is.null(method.bias)) {
    if (method.bias == "adjust1") {
      if (bias.type == "add") {
        if (!is.null(v)) {
          if (trt.effect == "random") {
            adj.ipd <- ""
            adj.ad <- ""
          }
          else {
            adj.ipd <- " + R[study[i]] * gamma[study[i]]"
            adj.ad <- " + R[j + ns.ipd] * gamma[(j + ns.ipd)]"
          }
        }
        else {
          adj.ipd <- " + R[study[i]] * gamma[study[i]]"
          adj.ad <- " + R[j + ns.ipd] * gamma[(j + ns.ipd)]"
        }
      }
      else if (bias.type == "mult") {
        adj.ipd <- " * gamma[study[i]]^R[study[i]]"
        adj.ad <- " * gamma[(j + ns.ipd)]^R[j + ns.ipd]"
      }
      else if (bias.type == "both") {
        adj.ipd <-
          paste0(" * gamma1[study[i]]^R[study[i]] + ",
                 "R[study[i]] * gamma2[study[i]]")
        adj.ad <-
          paste0(" * gamma1[(j + ns.ipd)]^R[j + ns.ipd] + ",
                 "R[j + ns.ipd] * gamma2[(j + ns.ipd)]")
      }
      else
        stop("The bias type should be set as 'add', 'mult' or 'both'.")


      if (bias.type == "both") {
        gamma.effect <- ""
        if (bias.effect == "random") {
          for (i in 1:2) {
            gamma.effect0 <-
              paste0("\n# Random effect for gamma (bias effect)",
                     if (add.std.in)
                       paste0("\nfor (j in std.in) {",
                              "\n  gamma", i, "[j] ~ dnorm(g", i, ", ",
                              "prec.gamma", i, ")",
                              "\n}"),
                     if (add.std.act.no)
                       paste0("\nfor (j in std.act.no) {",
                              "\n  gamma", i, "[j] ~ ",
                              "dnorm(0, prec.gamma", i, ")",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("\nfor (j in std.act.yes) {",
                              "\n  gamma", i, "[j] ~ dnorm(g.act", i, ", ",
                              "prec.gamma", i, ")",
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", i, dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act", i, dnmax),
                     "\nprec.gamma", i, " <- pow(tau.gamma", i, ", -2)",
                     "\ntau.gamma", i, " ~ ", prior.tau.gamma)
            ##
            gamma.effect <- paste0(gamma.effect, gamma.effect0)
          }
        }
        else {
          for (i in 1:2) {
            gamma.effect0 <-
              paste0("\n# Common effect for gamma (bias effect)",
                     if (add.std.in)
                       paste0("\nfor (j in std.in) {",
                              "\n  gamma", i, "[j] <- g", i,
                              "\n}"),
                     if (add.std.act.no)
                       paste0("\nfor (j in std.act.no) {",
                              "\n  gamma", i, "[j] <- 0",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("\nfor (j in std.act.yes) {",
                              "\n  gamma", i, "[j] <- g.act", i,
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", i, dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act", i, dnmax))
            ##
            gamma.effect <- paste0(gamma.effect, gamma.effect0)
          }
          ##
          gamma.effect <- paste0(gamma.effect, "\nprec.gamma <- 0")
          message("Bias effect is assumed common across studies")
        }
      }
      else {
        if (!is.null(v)) { # with v
          gamma.effect <-
            paste0("\n# Common effect for gamma (bias effect)",
                   if (add.std.in)
                     paste0("\nfor (j in std.in) {",
                            "\n  gamma[j] <- g",
                            "\n}"),
                   if (add.std.act.no)
                     paste0("\nfor (j in std.act.no) {",
                            "\n  gamma[j] <- 0",
                            "\n}"),
                   if (add.std.act.yes)
                     paste0("\nfor (j in std.act.yes) {",
                            "\n  gamma[j] <- g.act",
                            "\n}"),
                   if (add.std.in)
                     paste0("\ng" , dnmax),
                   if (add.std.act.yes)
                     paste0("\ng.act", dnmax))
        }
        else { # no v
          if (bias.effect == "random") {
            gamma.effect <-
              paste0("\n# Random effect for gamma (bias effect)",
                     if (add.std.in)
                       paste0("\nfor (j in std.in) {",
                              "\ngamma[j] ~ dnorm(g, prec.gamma)",
                              "\n}"),
                     if (add.std.act.no)
                       paste0("\nfor (j in std.act.no) {",
                              "\ngamma[j] ~ dnorm(0, prec.gamma)",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("\nfor (j in std.act.yes) {",
                              "\ngamma[j] ~ dnorm(g.act, prec.gamma)",
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act" , dnmax),
                     "\ntau.gamma ~ ", prior.tau.gamma,
                     "\nprec.gamma <- pow(tau.gamma, -2)")
          }
          else {
            gamma.effect <-
              paste0("\n# Common effect for gamma (bias effect)",
                     if (add.std.in)
                       paste0("\nfor (j in std.in) {",
                              "\n  gamma[j] <- g",
                              "\n}"),
                     if (add.std.act.no)
                       paste0("\nfor (j in std.act.no) {",
                              "\n  gamma[j] <- 0",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("\nfor (j in std.act.yes) {",
                              "\n  gamma[j] <- g.act",
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act" , dnmax),
                     "\nprec.gamma <- 0"
                     )
            ##
            message("Bias effect is assumed common across studies")
          }
        }
      }
      ##
      ##
      ##
      if (is.null(bias.covariate)) {
        adjust.prior <-
          paste0(gamma.effect,
                 "\n# Bias adjustment",
                 "\nfor (j in 1:(ns.ipd + ns.ad)) {",
                 "\n  R[j] ~ dbern(pi[bias_index[j]])",
                 "\n}",
                 "\npi[1] ~ ", prior.pi.high.rct, " # high RCT",
                 "\npi[2] ~ ", prior.pi.low.rct, " # low RCT",
                 "\npi[3] ~ ", prior.pi.high.nrs, " # high NRS",
                 "\npi[4] ~ ", prior.pi.low.nrs, " # low NRS",
                 "\npi[5] ~ dbeta(1, 1)  # unclear RCT or NRS")
      }
      else {
        adjust.prior <-
          paste0(gamma.effect,
                 "\n# Bias adjustment",
                 "\nfor (j in 1:(ns.ipd + ns.ad)) {",
                 "\n  R[j] ~ dbern(pi[j])",
                 "\n  logit(pi[j]) <- a + b * xbias[j]",
                 "\n}",
                 "\na", dnmax,
                 "\nb", dnmax)
      }
    }
    else if (method.bias == "adjust2") {
      if (is.null(bias.covariate)) { # assign priors for pi[bias_index[j]] based on design and RoB of the study
        if (trt.effect=="random") {
          if (!is.null(v)) {
            theta.effect.ipd <-
              paste0("\n    theta1", sel.i, " ~ dnorm(md", sel.i, ", precd", sel.i, ")",
                     "\n    theta2", sel.i, " ~ dnorm(md", sel.i, " + gamma[j], precd", sel.i, " / q)",
                     "\ntheta", sel.i, " <- (1 - pi[bias_index[j]]) * theta1", sel.i,
                     " + pi[bias_index[j]] * theta2", sel.i,
                     "\n# Multi-arm correction",
                     "\nmd", sel.i, " <- mean[j, k] + sw[j, k]",
                     "\nw[j, k] <- (theta", sel.i, " - mean[j, k])",
                     "\nsw[j, k] <- sum(w[j, 1:(k - 1)]) / (k - 1)",
                     "\nprecd", sel.i, "<- prec * 2 * (k - 1) / k")
            ##
            theta.effect.ad <-
              paste0("\n    theta1[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]], precd.ad[j, t.ad[j, k]])",
                     "\n    theta2[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]] + ",
                     "gamma[j + ns.ipd], precd.ad[j, t.ad[j, k]] / q)",
                     "\n    theta[j + ns.ipd, t.ad[j, k]] <- ",
                     "(1 - pi[bias_index[j]]) * theta1[j + ns.ipd, t.ad[j, k]] + ",
                     "pi[bias_index[j]] * theta2[j + ns.ipd, t.ad[j, k]]",
                     "\n    # Multi-arm correction",
                     "\n    md.ad[j, t.ad[j, k]] <- mean.ad[j, k] + sw.ad[j, k]",
                     "\n    w.ad[j, k] <- (theta[j + ns.ipd, t.ad[j, k]] - mean.ad[j, k])",
                     "\n    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)",
                     "\n    precd.ad[j, t.ad[j, k]] <- prec * 2 * (k - 1) / k")
            ##
            prior.tau.theta <-
              paste0("\n# Heterogeneity between theta's",
                     "\ntau ~", prior.tau.trt,
                     "\nprec <- pow(tau, -2)")
          }
          else {
            theta.effect.ipd <-
              paste0("\n    theta1", sel.i, " ~ dnorm(md", sel.i, ", ",
                     "precd", sel.i, ")",
                     "\n    theta2", sel.i, " ~ dnorm(md", sel.i, " + ",
                     "gamma[j], precd", sel.i, " + (prec.gamma * 2 * (k - 1) / k))",
                     "\n    theta", sel.i, " <- (1 - pi[bias_index[j]]) * ",
                     "theta1", sel.i, " + pi[bias_index[j]] * theta2", sel.i,
                     "\n    # Multi-arm correction",
                     "\n    md", sel.i, "<- mean[j, k] + sw[j, k]",
                     "\n    w[j, k] <- (theta", sel.i, " - mean[j, k])",
                     "\n    sw[j, k] <- sum(w[j, 1:(k - 1)]) / (k - 1)",
                     "\n    precd", sel.i, " <- prec * 2 * (k - 1) / k")
            ##
            theta.effect.ad <-
              paste0("\n    theta1[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]], precd.ad[j, t.ad[j, k]])",
                     "\n    theta2[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]] + ",
                     "gamma[j + ns.ipd], precd.ad[j, t.ad[j, k]] + ",
                     "(prec.gamma * 2 * (k - 1) / k))",
                     "\n    theta[j + ns.ipd, t.ad[j, k]] <- ",
                     "(1 - pi[bias_index[j]]) * ",
                     "theta1[j + ns.ipd, t.ad[j, k]] + pi[bias_index[j]] * ",
                     "theta2[j + ns.ipd, t.ad[j, k]]",
                     "\n    # Multi-arm correction",
                     "\n    md.ad[j, t.ad[j, k]] <- mean.ad[j, k] + sw.ad[j, k]",
                     "\n    w.ad[j, k] <- ",
                     "(theta[j + ns.ipd, t.ad[j, k]] - mean.ad[j, k])",
                     "\n    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)",
                     "\n    precd.ad[j, t.ad[j, k]] <- prec * 2 * (k - 1) / k")
            ##
            prior.tau.theta <-
              paste0("\n# Heterogeneity between theta's",
                     "\ntau ~", prior.tau.trt,
                     "\nprec <- pow(tau, -2)")
          }
        }
        else if (trt.effect == "common") {
          theta.effect.ipd <-
            paste0("\n    theta", sel.i, " <- md", sel.i, " + ",
                   "(pi[bias_index[j]] * gamma[j])",
                   "\n    md", sel.i, " <- mean[j, k]")
          ##
          theta.effect.ad <-
            paste0("\n    theta[j + ns.ipd, t.ad[j, k]] <- ",
                   "md.ad[j, t.ad[j, k]] + (pi[bias_index[j]] * ",
                   "gamma[j + ns.ipd])",
                   "\n    md.ad[j, t.ad[j, k]] <- mean.ad[j, k]")
          ##
          prior.tau.theta <- ""
        }
        else {
          stop("Please indicate the model of treatment effect as either ",
               "'random' or 'common'.")
        }
        
        
        if (!is.null(v)) { # with v
          gamma.effect <-
            paste0("\n# Common effect for gamma (bias effect)\n",
                   if (add.std.in)
                     paste0("for (j in std.in) {",
                            "\n  gamma[j] <- g",
                            "\n}"),
                   if (add.std.act.no)
                     paste0("for (j in std.act.no) {",
                            "\n  gamma[j] <- 0",
                            "\n}"),
                   if (add.std.act.yes)
                     paste0("for (j in std.act.yes) {",
                            "\n  gamma[j] <- g.act",
                            "\n}"),
                   if (add.std.in)
                     paste0("\ng", dnmax),
                   if (add.std.act.yes)
                     paste0("\ng.act", dnmax))
        }
        else { # no v
          if (bias.effect == "random") {
            gamma.effect <-
              paste0("\n# Random effect for gamma (bias effect)",
                     if (add.std.in)
                       paste0("for (j in std.in) {",
                              "\n  gamma[j] ~ dnorm(g, prec.gamma)",
                              "\n}"),
                     if (add.std.act.no)
                       paste0("for (j in std.act.no) {",
                              "\n  gamma[j] ~ dnorm(0, prec.gamma)",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("for (j in std.act.yes) {",
                              "\n  gamma[j] ~ dnorm(g.act, prec.gamma)",
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act", dnmax),
                     if (is.null(v))
                       paste0("\ntau.gamma ~ ", prior.tau.gamma,
                              "\nprec.gamma <- pow(tau.gamma, -2)"))
          }
          else {
            gamma.effect <-
              paste0("\n# Common effect for gamma (bias effect)\n",
                     if (add.std.in)
                       paste0("for (j in std.in) {",
                              "\n  gamma[j] <- g",
                              "\n}"),
                     if (add.std.act.no)
                       paste0("for (j in std.act.no) {",
                              "\n  gamma[j] <- 0",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("for (j in std.act.yes) {",
                              "\n  gamma[j] <- g.act",
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act", dnmax),
                     "\nprec.gamma <- 0")
            ##
            warning("Bias effect is assumed common across studies")
          }
        }
        adjust.prior <-
          paste0(gamma.effect,
                 "\n# Bias adjustment",
                 "\npi[1] ~ ", prior.pi.high.rct, " # high RCT",
                 "\npi[2] ~ ", prior.pi.low.rct, " # low RCT",
                 "\npi[3] ~ ", prior.pi.high.nrs, " # high NRS",
                 "\npi[4] ~ ", prior.pi.low.nrs, " # low NRS",
                 "\npi[5] ~ dbeta(1, 1)  # unclear RCT or NRS")
      }
      else { # estimate pi[j] through a logistic model using the bias covariate
        if (trt.effect == "random") {
          if (!is.null(v)) {
            theta.effect.ipd <-
              paste0("\n    theta1", sel.i, " ~ ",
                     "dnorm(md", sel.i, ", precd", sel.i, ")",
                     "\n    theta2", sel.i, " ~ ",
                     "dnorm(md", sel.i, " + gamma[j], precd", sel.i, " / q)",
                     "\n    theta", sel.i, " <- ",
                     "(1 - pi[j]) * theta1", sel.i, " + pi[j] * theta2", sel.i,
                     "\n    # Multi-arm correction",
                     "\n    md", sel.i, "<- mean[j, k] + sw[j, k]",
                     "\n    w[j, k] <- (theta", sel.i, " - mean[j, k])",
                     "\n    sw[j, k] <- sum(w[j, 1:(k - 1)]) / a(k - 1)",
                     "\n    precd", sel.i, " <- prec * 2 * (k - 1) / k")
            ##
            theta.effect.ad <-
              paste0("\n    theta1[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]], precd.ad[j, t.ad[j, k]])",
                     "\n    theta2[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]] + ",
                     "gamma[j + ns.ipd], precd.ad[j, t.ad[j, k]] / q)",
                     "\n    theta[j + ns.ipd, t.ad[j, k]] <- ",
                     "(1 - pi[j]) * theta1[j + ns.ipd, t.ad[j, k]] + ",
                     "pi[j] * theta2[j + ns.ipd, t.ad[j, k]]",
                     "\n    # Multi-arm correction",
                     "\n    md.ad[j, t.ad[j, k]] <- mean.ad[j, k] + sw.ad[j, k]",
                     "\n    w.ad[j, k] <- (theta[j + ns.ipd, t.ad[j, k]] - mean.ad[j, k])",
                     "\n    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)",
                     "\n    precd.ad[j, t.ad[j, k]] <- prec * 2 * (k - 1) / k")
            ##
            prior.tau.theta <-
              paste0("\n# Heterogeneity between theta's",
                     "\ntau ~", prior.tau.trt,
                     "\nprec <- pow(tau, -2)")
          }
          else {
            theta.effect.ipd <-
              paste0("\n    theta1", sel.i, " ~ ",
                     "dnorm(md", sel.i, ", precd", sel.i, ")",
                     "\n    theta2", sel.i, " ~ ",
                     "dnorm(md", sel.i, " + gamma[j], precd", sel.i, " + ",
                     "(prec.gamma * 2 * (k - 1) / k))",
                     "\n    theta", sel.i, " <- (1 - pi[j]) * theta1",
                     sel.i, " + pi[j] * theta2", sel.i,
                     "\n    # Multi-arm correction",
                     "\n    md", sel.i, "<- mean[j, k] + sw[j, k]",
                     "\n    w[j, k] <- (theta", sel.i, " - mean[j, k])",
                     "\n    sw[j, k] <- sum(w[j, 1:(k - 1)]) / (k - 1)",
                     "\n    precd", sel.i, " <- prec * 2 * (k - 1) / k")
            ##
            theta.effect.ad <-
              paste0("\n    theta1[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]], precd.ad[j, t.ad[j, k]])",
                     "\n    theta2[j + ns.ipd, t.ad[j, k]] ~ ",
                     "dnorm(md.ad[j, t.ad[j, k]] + ",
                     "gamma[j + ns.ipd], precd.ad[j, t.ad[j, k]] + ",
                     "(prec.gamma * 2 * (k - 1) / k))",
                     "\n    theta[j + ns.ipd, t.ad[j, k]] <- ",
                     "(1 - pi[j]) * theta1[j + ns.ipd, t.ad[j, k]] + ",
                     "pi[j] * theta2[j + ns.ipd, t.ad[j, k]]",
                     "\n    # Multi-arm correction",
                     "\n    md.ad[j, t.ad[j, k]] <- mean.ad[j, k] + sw.ad[j, k]",
                     "\n    w.ad[j, k] <- (theta[j + ns.ipd, t.ad[j, k]] - ",
                     "mean.ad[j, k])",
                     "\n    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)",
                     "\n    precd.ad[j, t.ad[j, k]] <- prec * 2 * (k - 1) / k")
            ##
            prior.tau.theta <-
              paste0("\n# Heterogeneity between theta's",
                     "\ntau ~", prior.tau.trt,
                     "\nprec <- pow(tau, -2)")
          }
        }
        else if (trt.effect == "common") {
          theta.effect.ipd <-
            paste0("\n    theta", sel.i, " <- md", sel.i, " + ",
                   "(pi[j] * gamma[j])",
                   "\n    md", sel.i, " <- mean[j, k]")
          ##
          theta.effect.ad <-
            paste0("\n    theta[j + ns.ipd, t.ad[j, k]] <- ",
                   "md.ad[j, t.ad[j, k]] + (pi[j] * gamma[j + ns.ipd])",
                   "\n    md.ad[j, t.ad[j, k]] <- mean.ad[j, k]")
          ##
          prior.tau.theta <- ""
        }
        else
          stop("Please indicate the model of treatment effect as ",
               "either 'random' or 'common'.")
        
        
        if (!is.null(v)) { # with v
          gamma.effect <-
            paste0("\n# Common effect for gamma (bias effect)",
                   if (add.std.in)
                     paste0("\nfor (j in std.in) {",
                            "\n  gamma[j] <- g",
                            "\n}"),
                   if (add.std.act.no)
                     paste0("\nfor (j in std.act.no) {",
                            "\n  gamma[j] <- 0",
                            "\n}"),
                   if (add.std.act.yes)
                     paste0("\nfor (j in std.act.yes) {",
                            "\n  gamma[j] <- g.act",
                            "\n}"),
                   if (add.std.in)
                     paste0("\ng", dnmax),
                   if (add.std.act.yes)
                     paste0("\ng.act", dnmax))
        }
        else { # no v
          if (bias.effect == "random") {
            gamma.effect <-
              paste0("\n# Random effect for gamma (bias effect)",
                     if (add.std.in)
                       paste0("\nfor (j in std.in) {",
                              "\n  gamma[j] ~ dnorm(g, prec.gamma)",
                              "\n}"),
                     if (add.std.act.no)
                       paste0("for (j in std.act.no) {",
                              "\n  gamma[j] ~ dnorm(0, prec.gamma)",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("for (j in std.act.yes) {",
                              "\n  gamma[j] ~ dnorm(g.act, prec.gamma)",
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act", dnmax),
                     if (is.null(v))
                       paste0("\ntau.gamma ~ ", prior.tau.gamma,
                              "\nprec.gamma <- pow(tau.gamma, -2)"))
          }
          else {
            gamma.effect <-
              paste0("\n# Common effect for gamma (bias effect)\n",
                     if (add.std.in)
                       paste0("for (j in std.in) {",
                              "\n  gamma[j] <- g",
                              "\n}"),
                     if (add.std.act.no)
                       paste0("for (j in std.act.no) {",
                              "\n  gamma[j] <- 0",
                              "\n}"),
                     if (add.std.act.yes)
                       paste0("for (j in std.act.yes) {",
                              "\n  gamma[j] <- g.act",
                              "\n}"),
                     if (add.std.in)
                       paste0("\ng", dnmax),
                     if (add.std.act.yes)
                       paste0("\ng.act", dnmax),
                     "\nprec.gamma <- 0")
            ##
            warning("Bias effect is assumed common across studies")
          }
        }
        ##        
        adjust.prior <-
          paste0(gamma.effect,
                 "\n# Bias adjustment",
                 "for (j in 1:(ns.ipd + ns.ad)) {",
                 "\n  logit(pi[j]) <- a + b * xbias[j]",
                 "\n}",
                 "\na", dnmax,
                 "\nb", dnmax)
      }
    }
    else if (method.bias == "prior") {
      d.prior <- d.prior.nrs
    }
    else if (method.bias == "naive" & ipd & ad) {
      message("Both designs are combined naively without acknowledging ",
              "design differences.")
    }
  }
  
  
  ##
  ## Set up the likelihood and the link
  ##
  if (sm == "OR") {
    like.ipd <-
      paste0("\n  # Bernoulli likelihood",
             "\n  y[i] ~ dbern(p[i])")
    link.ipd <-
      "\n  logit(p[i]) <- u[study[i]] + theta[study[i], trt[i]]"
    ##
    like.ad <-
      paste0("\n    # Binomial likelihood of number of events",
             "\n    r[j, k] ~ dbin(pa[j, t.ad[j, k]], n[j, k])")
    link.ad.ref <-
      paste0("\n  # Log odds at reference arm",
             "\n  logit(pa[j, t.ad[j, 1]]) <- u[j + ns.ipd]")
    link.ad <-
      paste("\n    logit(pa[j, t.ad[j, k]]) <-",
            "u[j + ns.ipd] + (theta[j + ns.ipd, t.ad[j, k]])")
  }
  else if (sm == "RR") {
    like.ipd <-
      paste0("\n  # Bernoulli likelihood",
             "\n  y[i] ~ dbern(p[i])")
    link.ipd <-
      "\nlog(p[i]) <- u[study[i]] + theta[study[i], trt[i]]"
    ##
    like.ad <-
      paste0("\n    # binomial likelihood of number of events",
             "\n    r[j, k] ~ dbin(pa[j, t.ad[j, k]], n[j, k])")
    link.ad.ref <-
      paste0("\n  # Log of risk probability at reference arm",
             "\n  log(pa[j, t.ad[j, 1]]) <- u[j + ns.ipd]")
    link.ad <-
      paste0("\n    log(pa[j, t.ad[j, k]]) <- ",
             "u[j + ns.ipd] + (theta[j + ns.ipd, t.ad[j, k]])")
  }
  else if (sm == "MD") {
    like.ipd <-
      paste0("\n  # Normal likelihood",
             "\n  y[i] ~ dnorm(delta[i], ",
             "prec.delta.ipd[study[i], trt.index[i]])")
    link.ipd <-
      "\ndelta[i] <- u[study[i]] + (theta[study[i], trt[i]])"
    ##
    like.ad <-
      paste0("\n    # Normal likelihood",
             "\n    ybar[j, k] ~ dnorm(delta.ad[j, t.ad[j, k]], ",
             "prec.delta.ad[j, k])")
    link.ad.ref <-
      paste0("\n  # Mean outcome at reference arm",
             "\n  delta.ad[j, t.ad[j, 1]] <- u[j + ns.ipd]")
    link.ad <-
      paste0("\n    delta.ad[j, t.ad[j, k]] <- ",
             "u[j + ns.ipd] + (theta[j + ns.ipd, t.ad[j, k]])")
  }
  else if (sm == "SMD") {
    like.ipd <-
      paste0("\n  # Normal likelihood",
             "\n  y[i] ~ dnorm(phi[i], prec.delta.ipd[study[i], trt.index[i]])")
    link.ipd <-
      paste0("\n  phi[i] <- delta[i] * s.pool.ipd[study[i]]",
             "\n  delta[i] <- u[study[i]] + theta[study[i], trt[i]]")
    ##
    like.ad <-
      paste0("\n    # Normal likelihood",
             "\n    ybar[j, k] ~ dnorm(phi.ad[j, t.ad[j, k]], ",
             "prec.delta.ad[j, k])",
             "\n    phi.ad[j, t.ad[j, k]] <- delta.ad[j, t.ad[j, k]] * ",
             "s.pool.ad[j]")
    link.ad.ref <-
      "\n  delta.ad[j, t.ad[j, 1]] <- u[j + ns.ipd] * s.pool.ad[j]"
    link.ad <-
      paste0("\n    delta.ad[j, t.ad[j, k]] <- u[j + ns.ipd] + ",
             "(theta[j + ns.ipd, t.ad[j, k]])")
  }
  
  
  ##
  ##
  ## Combine the code
  ##
  ##
  
  
  ##
  ## IPD part
  ##
  
  ipd.code <- sprintf("

#
# (1) IPD part
#

# Loop through individuals
for (i in 1:np) {%s
  # Link function%s%s%s
}

# Loop through IPD studies
for (j in 1:(ns.ipd)) {
  # At each study reference arm
  w[j, 1] <- 0
  theta[j, t.ipd[j, 1]] <- 0%s
  # Loop through non-reference IPD arms
  for (k in 2:na.ipd[j]) {
    # Synthesize relative treatment effects%s
    # Consistency equation
    mean[j, k] <- d[t.ipd[j, k]] - d[t.ipd[j, 1]]%s%s
  }
}",
like.ipd,           #%s2
link.ipd,           #%s3
adj.ipd,            #%s4
mreg.ipd,           #%s5
ref.trt.effect.ipd, #%s6
theta.effect.ipd,   #%s7
betab.consis.ipd,   #%s8
betaw.consis.ipd)   #%s9
  
  
  ##
  ## AD part
  ##
  ad.code <- sprintf("


#
# (2) AD part
#

# Loop through AD studies
for (j in 1:ns.ad) {
  # Multi-arm correction is zero for reference arm
  w.ad[j, 1] <- 0
  # Treatment effect is zero for reference arm
  theta[j + ns.ipd, t.ad[j, 1]] <- 0%s
  # Loop through AD arms
  for (k in 1:na.ad[j]) {%s
  }
  # Effect in reference arm%s
  # Loop through non-reference AD arms
  for (k in 2:na.ad[j]) {
    # Link function%s%s%s
    # Synthesize relative treatment effects%s
    # Consistency equations
    mean.ad[j, k] <- d[t.ad[j, k]] - d[t.ad[j, 1]]%s
  }
}",
ref.trt.effect.ad,
like.ad,
link.ad.ref,
link.ad,
adj.ad,
mreg.ad,
theta.effect.ad,
betab.consis.ad)
  
  
  ##
  ## Prior part
  ##
  
  prior.code <- sprintf("


#
# (3) Priors
#
%s%s
# The effect is zero at network reference treatment
d[1] <- 0%s%s%s%s%s%s%s",
prior.u,
prior.tau.theta,
d.prior,
beta0.prior.ipd,
betab.prior,
betaw.prior.ipd,
beta.prior,
adjust.prior,
q.prior)
  
  
  ##
  ## Everything
  ##
  
  ad.code <- if (ad) ad.code else ""
  ipd.code <- if (ipd) ipd.code else ""
  code.str <-
    paste0("model {", ipd.code, ad.code, prior.code, "\n\n}\n")
  ##
  code.str
}
