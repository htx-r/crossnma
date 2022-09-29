crossnma.code <- function(ipd = TRUE,
                          ad = TRUE,
                          sm=NULL,
                          trt.effect = 'random',
                          prior.tau.trt = NULL,
                          ## -------- meta-regression
                          split.regcoef = FALSE,
                          covariate = NULL,
                          reg0.effect = 'random',
                          regb.effect = 'random',
                          regw.effect = 'random',
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
                          d.prior.nrs = NULL # required when method.bias='prior'
                          ) {


  ##
  ##
  ## User priors
  ##
  ##

  prior.tau.trt <- replaceNULL(prior.tau.trt, 'dunif(0, 2)')
  prior.tau.reg0 <- replaceNULL(prior.tau.reg0, 'dunif(0, 2)')
  prior.tau.regb <- replaceNULL(prior.tau.regb, 'dunif(0, 2)')
  prior.tau.regw <- replaceNULL(prior.tau.regw, 'dunif(0, 2)')
  prior.tau.gamma<- replaceNULL(prior.tau.gamma, 'dunif(0, 2)')
  prior.pi.high.rct <- replaceNULL(prior.pi.high.rct, 'dbeta(10, 1)')
  prior.pi.low.rct <- replaceNULL(prior.pi.low.rct, 'dbeta(1, 10)')
  prior.pi.high.nrs <- replaceNULL(prior.pi.high.nrs, 'dbeta(30, 1)')
  prior.pi.low.nrs <- replaceNULL(prior.pi.low.nrs, 'dbeta(1, 30)')


  ##
  ##
  ## Meta-regression
  ##
  ##

  metareg.str.ipd <- ""
  betab.consis.ipd <- ""
  betaw.consis.ipd <- ""
  beta0.prior.ipd <- ""
  betab.prior <- ""
  betaw.prior.ipd <- ""
  beta.prior.ipd <- ""
  ##
  beta.prior.ad <- ""
  metareg.str.ad <- ""
  betab.consis.ad <- ""
  ##
  if (!is.null(covariate)) {
    if (ipd) {
      for (i in 1:length(covariate[[1]])) {
        ## meta-regression terms - Up to 3
        metareg.str.ipd0 <-
          paste0(
            " + beta0_", i, "[study[i]] * (x", i, "[i]) + betaw_", i,
            "[study[i], trt[i]] * (x", i,
            "[i] - xm", i, ".ipd[i]) + betab_", i,
            "[study[i], trt[i]] * xm", i, ".ipd[i]")
        ##
        metareg.str.ipd <- paste0(metareg.str.ipd,metareg.str.ipd0)
        ## consistency equations for beta_b and beta_w - Up to 3
        betab.consis.ipd0 <-
          paste0(
            "    betab_", i, "[j, t.ipd[j, k]] <- betab.t_", i,
            "[t.ipd[j, k]] - betab.t_", i, "[t.ipd[j, 1]]")
        betaw.consis.ipd0 <-
          paste0("    betaw_", i, "[j, t.ipd[j, k]] <- betaw.t_", i,
                 "[t.ipd[j, k]] - betaw.t_", i, "[t.ipd[j, 1]]")
        ##
        betab.consis.ipd <- paste0(betab.consis.ipd, "\n", betab.consis.ipd0)
        ##
        betaw.consis.ipd <- paste0(betaw.consis.ipd, "\n", betaw.consis.ipd0)
      }

      ## beta0
      if (reg0.effect == 'random') {
        for (i in 1:length(covariate[[1]])) {
          beta0.prior.ipd0 <-
            paste0("\n\n# Random effect for beta0
for (j in 1:(ns.ipd)) {
  beta0_", i, "[j] ~ dnorm(b0_", i,
", prec.beta0_", i, ")
}
b0_", i, " ~ dnorm(0, .01)\n
prec.beta0_", i, " <- pow(tau.b0_", i, ", -2)
tau.b0_", i, " ~ ", prior.tau.reg0
)
          ##
          beta0.prior.ipd <- paste0(beta0.prior.ipd, beta0.prior.ipd0)
        }
      }
      else if (reg0.effect == 'independent') {
        for (i in 1:length(covariate[[1]])) {
          beta0.prior.ipd0 <-
            paste0("\n\n# Independent effect for beta0
for (j in 1:(ns.ipd)) {
  beta0_", i, "[j] <- b0_", i, "[j]
  b0_", i, "[j] ~ dnorm(0, .01)
}")
          ##
          beta0.prior.ipd <- paste0(beta0.prior.ipd, beta0.prior.ipd0)
        }
      }
      else
        stop("The progonostic effect can be assumed either ",
             "'independent' or 'random' across studies")


      ## betab and betaw

      if (!split.regcoef) { # not splitted within and between- study covariate
        if (regb.effect == 'independent' || regw.effect == 'independent') {
          beta.prior.ipd0 <-
       paste0("\n\n# Random effect for beta (within=between)
beta.t_", i, "[1] <- 0
for (k in 1:nt) {
  betab.t_", i, "[k] <- beta.t_", i, "[k]
  betaw.t_", i, "[k] <- beta.t_", i, "[k]
}
for (k in 2:nt) {
  beta.t_", i, "[k] ~ dnorm(0, 1e-2)
}")
          ##
          beta.prior.ipd <- paste0(beta.prior.ipd, beta.prior.ipd0)
        }
        else if (regb.effect == 'random' || regw.effect == 'random') {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ipd0 <-
       paste0("\n\n# Random effects for beta (within=between)
beta.t_", i, "[1] <- 0
for (k in 1:nt) {
  betab.t_", i, "[k] <- beta.t_", i, "[k]
  betaw.t_", i, "[k] <- beta.t_", i, "[k]
}
for (k in 2:nt) {
  beta.t_", i, "[k] ~ dnorm(b_", i, ", prec.beta_", i, ")
}
b_", i, " ~ dnorm(0, 1e-2)
tau.b_", i, " ~ ", prior.tau.regw,
"\nprec.beta_", i, " <- pow(tau.b_", i, ", -2)")
            ##
            beta.prior.ipd <- paste0(beta.prior.ipd, beta.prior.ipd0)
          }
        }
        else if (regb.effect == 'common' & regw.effect == 'common') {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ipd0 <-
       paste0("\n\n# Common effect for beta (within=between)
betab.t_", i, "[1] <- 0
betaw.t_", i, "[1] <- 0
for (k in 2:nt) {
  betab.t_", i, "[k] <- b_", i, "
  betaw.t_", i, "[k] <- b_", i, "
}
b_", i, " ~ dnorm(0, 1e-2)")
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
        if (regb.effect == 'independent') {
          betab.prior0 <-
            paste0("\n\n# Random effect for betab (the between-study covariate effect)
betab.t_", i, "[1] <- 0
for (k in 2:nt) {
  betab.t_", i, "[k] ~ dnorm(0, 1e-2)
}")
          ##
          betab.prior <- paste0(betab.prior, betab.prior0)
        }
        else if (regb.effect == 'random') {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
              paste0("\n\n# Random effect for betab (the between-study covariate effect)
betab.t_", i, "[1] <- 0
for (k in 2:nt) {
  betab.t_", i, "[k] ~ dnorm(bb_", i, ", prec.betab_", i, ")
}
bb_", i, " ~ dnorm(0, 1e-2)
tau.bb_", i, " ~ ", prior.tau.regb,
    "\nprec.betab_", i, " <- pow(tau.bb_", i, ", -2)")
            ##
            betab.prior <- paste0(betab.prior, betab.prior0)
          }
        }
        else if (regb.effect == 'common') {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
              paste0("\n\n# Common effect for betab (the between-study covariate effect)
betab.t_", i, "[1] <- 0
for (k in 2:nt) {
  betab.t_", i, "[k] <- bb_", i, "
}
bb_", i, " ~ dnorm(0, 1e-2)")
            ##
            betab.prior <- paste0(betab.prior, betab.prior0)
          }
        }
        else
          stop("The between-study covariate effect need to be assumed ",
               "'independent', 'random' or 'common' across studies")


        ## within- study covariate
        if (regw.effect == 'independent') {
          betaw.prior.ipd0 <-
       paste0("\n\n# Random effect for betab (the between-study covariate effect)
betaw.t_", i, "[1] <- 0
for (k in 2:nt) {
  betaw.t_", i, "[k] ~ dnorm(0, 1e-2)
}")
          ##
          betaw.prior.ipd <- paste0(betaw.prior.ipd, betaw.prior.ipd0)
        }
        else if (regw.effect == 'random') {
          for (i in 1:length(covariate[[1]])) {
            betaw.prior.ipd0 <-
       paste0("\n\n# Random effect for betaw (the within-study covariate effect)
betaw.t_", i, "[1] <- 0
for (k in 2:nt) {
  betaw.t_", i, "[k] ~ dnorm(bw_", i, ", prec.betaw_", i, ")
}
bw_", i, " ~ dnorm(0, 1e-2)
prec.betaw_", i, "<- pow(tau.bw_", i, ", -2)",
     "\n  tau.bw_", i, " ~ ", prior.tau.regw)
            ##
          betaw.prior.ipd <- paste0(betaw.prior.ipd, betaw.prior.ipd0)
          }
        }
        else if (regw.effect == 'common') {
          for (i in 1:length(covariate[[1]])) {
            betaw.prior.ipd0 <-
              paste0("\n\n# Common effect for betaw (the within-study covariate effect)
    betaw.t_", i, "[1] <- 0
    for (k in 2:nt) {
      betaw.t_", i, "[k] <- bw_", i, "
    }
    bw_", i, " ~ dnorm(0, 1e-2)
    ")
            ##
            betaw.prior.ipd <- paste0(betaw.prior.ipd, betaw.prior.ipd0)
          }
        }
        else
          stop("The within-study covariate effect can be assumed ",
               "'independent', 'random' or 'common' across studies")
      }
    }


    if (ad) {
      for (i in 1:length(covariate[[1]])) {
        ## meta-regression terms - up to 3
        metareg.str.ad0 <- paste0(" + betab.ad_", i, "[j, t.ad[j, k]] * xm", i, ".ad[j]")
        betab.consis.ad0 <- paste0("    betab.ad_", i, "[j, t.ad[j, k]] <- betab.t_", i, "[t.ad[j, k]] - betab.t_", i, "[t.ad[j, 1]]")

        ## consistency equation - up to 3
        metareg.str.ad <- paste0(metareg.str.ad,metareg.str.ad0)
        betab.consis.ad <- paste0(betab.consis.ad, "\n", betab.consis.ad0)
      }
      if (!split.regcoef) { # not splitted
        if (regb.effect == 'independent' && regw.effect == 'independent') {
          beta.prior.ad0 <-
       paste0("\n\n# Random effect for beta (within=between)

    beta.t_", i, "[1] <- 0
    for (k in 1:nt) {
      betab.t_", i, "[k] <- beta.t_", i, "[k]
    }
    for (k in 2:nt) {
      beta.t_", i, "[k] ~ dnorm(0, 1e-2)
    }")
          ##
          beta.prior.ad <- paste0(beta.prior.ad, beta.prior.ad0)
        }
        else if (regb.effect == 'random' && regw.effect == 'random') {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ad0 <-
       paste0("\n\n# Random effects for beta (within=between)
    beta.t_", i, "[1] <- 0
    for (k in 1:nt) {
      betab.t_", i, "[k] <- beta.t_", i, "[k]
    }
    for (k in 2:nt) {
      beta.t_", i, "[k] ~ dnorm(b_", i, ", prec.beta_", i, ")
    }
    b_", i, " ~ dnorm(0, 1e-2)
    tau.b_", i, " ~ ", prior.tau.regb,
    "\n  prec.beta_", i, " <- pow(tau.b_", i, ", -2)")
            ##
            beta.prior.ad <- paste0(beta.prior.ad, beta.prior.ad0)
          }
        }
        else if (regb.effect == 'common' & regw.effect == 'common') {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ad0 <-
              paste0("\n\n# Common effects for beta (within=between)
     betab.t_", i, "[1] <- 0
     for (k in 2:nt) {
       betab.t_", i, "[k] <- b_", i, "
     }
     b_", i, " ~ dnorm(0, 1e-2)
     ")
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
        if (regb.effect == 'independent') {
          betab.prior0 <-
       paste0("\n\n# Random effect for betab (the between-study covariate effect)
    betab.t_", i, "[1] <- 0
    for (k in 2:nt) {
      betab.t_", i, "[k] ~ dnorm(0, 1e-2)
    }")
          ##
          betab.prior <- paste0(betab.prior, betab.prior0)
        }
        else if (regb.effect == 'random') {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
       paste0("\n\n# Random effects for betab (the between-study covariate effect)
    betab.t_", i, "[1] <- 0
    for (k in 2:nt) {
      betab.t_", i, "[k] ~ dnorm(bb_", i, ", prec.betab_", i, ")
    }
    bb_", i, " ~ dnorm(0, 1e-2)
    tau.bb_", i, " ~ ", prior.tau.regb,
    "\n  prec.betab_", i, " <- pow(tau.bb_", i, ", -2)")
            ##
            betab.prior <- paste0(betab.prior, betab.prior0)
          }
        }
        else if (regb.effect == 'common') {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <-
       paste0("\n\n# Random effects for betab (the between-study covariate effect)
    betab.t_", i, "[1] <- 0
    for (k in 2:nt) {
      betab.t_", i, "[k] <- bb_", i, "
    }
    bb_", i, " ~ dnorm(0, 1e-2)")
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


  ## treatment effect is zero for reference treatment
  ref.trt.effect.ipd <- ""
  ref.trt.effect.ad <- ""

  if (!is.null(covariate)) {
    if (ipd) {
      for (i in 1:length(covariate[[1]])) {
        ## meta-regression terms - Up to 3
        ref.trt.effect.ipd0 <- paste0(
          "\n  betaw_", i, "[j, t.ipd[j, 1]] <- 0\n  betab_", i,
          "[j, t.ipd[j, 1]] <- 0"
        )
        ##
        ref.trt.effect.ipd <- paste0(ref.trt.effect.ipd,ref.trt.effect.ipd0)
      }
    }
    if (ad) {
      for (i in 1:length(covariate[[1]])) {
        ## meta-regression terms - Up to 3
        ref.trt.effect.ad0 <-
          paste0("\n  betab.ad_", i, "[j, t.ad[j, 1]] <- 0"
        )
        ref.trt.effect.ad0 <-
          if (ipd) ""
          else
            ref.trt.effect.ad0
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

  q.prior <- ""

  if (trt.effect == "random") {
    if (!is.null(v)) {
      theta.effect.ipd <- "\n    theta[j, t.ipd[j, k]] <- (1 - R[j]) * theta.un[j, t.ipd[j, k]] + R[j] * theta.adj[j, t.ipd[j, k]]
    theta.un[j, t.ipd[j, k]] ~ dnorm(md[j, t.ipd[j, k]], precd[j, t.ipd[j, k]])
    theta.adj[j, t.ipd[j, k]] ~ dnorm(md[j, t.ipd[j, k]] + gamma[j], precd[j, t.ipd[j, k]] / q)
    # Multi-arm correction
    md[j, t.ipd[j, k]] <- mean[j, k] + sw[j, k]
    w[j, k] <- theta[j, t.ipd[j, k]]  - mean[j, k]
    sw[j, k] <- sum(w[j, 1:(k - 1)]) / (k - 1)
    precd[j, t.ipd[j, k]] <- prec * 2 * (k - 1) / k
  "
      theta.effect.ad <-
        "\n    theta[j + ns.ipd, t.ad[j, k]] <- (1 - R[j + ns.ipd]) * theta.un[j + ns.ipd, t.ad[j, k]] + R[j + ns.ipd] * theta.adj[j + ns.ipd, t.ad[j, k]]
    theta.un[j + ns.ipd, t.ad[j, k]] ~ dnorm(md.ad[j, t.ad[j, k]], precd.ad[j, t.ad[j, k]])
    theta.adj[j + ns.ipd, t.ad[j, k]] ~ dnorm(md.ad[j, t.ad[j, k]] + gamma[(j + ns.ipd)], precd.ad[j, t.ad[j, k]] / q)
    # Multi-arm correction
    md.ad[j, t.ad[j, k]] <- mean.ad[j, k] + sw.ad[j, k]
    w.ad[j, k] <- (theta[j + ns.ipd, t.ad[j, k]]  - mean.ad[j, k])
    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)
    precd.ad[j, t.ad[j, k]] <- prec * 2 * (k - 1) / k
  "
      ##
      q.prior <- paste0("q ~ dbeta(",v, ", 1)")
      ##
      prior.tau.theta <-
        paste0("\n\n# heterogeneity between theta's\ntau ~ ",
               prior.tau.trt, "\nprec <- pow(tau, -2)")
    }
    else {
      theta.effect.ipd <-
        "\n    theta[j, t.ipd[j, k]] ~ dnorm(md[j, t.ipd[j, k]], precd[j, t.ipd[j, k]])
    # Multi-arm correction
    md[j, t.ipd[j, k]] <- mean[j, k] + sw[j, k]
    w[j, k] <- theta[j, t.ipd[j, k]]  - mean[j, k]
    sw[j, k] <- sum(w[j, 1:(k - 1)]) / (k - 1)
    precd[j, t.ipd[j, k]] <- prec * 2 * (k - 1) / k"
      theta.effect.ad <- "\n    theta[j + ns.ipd, t.ad[j, k]] ~ dnorm(md.ad[j, t.ad[j, k]], precd.ad[j, t.ad[j, k]])
    # Multi-arm correction
    md.ad[j, t.ad[j, k]] <- mean.ad[j, k] + sw.ad[j, k]
    w.ad[j, k] <- (theta[j + ns.ipd, t.ad[j, k]]  - mean.ad[j, k])
    sw.ad[j, k] <- sum(w.ad[j, 1:(k - 1)]) / (k - 1)
    precd.ad[j, t.ad[j, k]] <- prec * 2 * (k - 1) / k"
      ##
      prior.tau.theta <-
        paste0("\n\n# heterogeneity between theta's\ntau ~ ",
               prior.tau.trt, "\nprec <- pow(tau, -2)")
    }
  }
  else if (trt.effect == "common") {
    theta.effect.ipd <- "\n    theta[j, t.ipd[j, k]] <- md[j, t.ipd[j, k]]
  md[j, t.ipd[j, k]] <- mean[j, k]"
    theta.effect.ad <- "\n    theta[j + ns.ipd, t.ad[j, k]] <- md.ad[j, t.ad[j, k]]
  md.ad[j, t.ad[j, k]] <- mean.ad[j, k]"
    prior.tau.theta <- ""
  }
  else
    stop("Please indicate the treatment effect model as either ",
         "'random' or 'common' ")
  ##
  d.prior <- "\n  for (k in 2:nt) {
    d[k] ~ dnorm(0, .01)
  }"


  ##
  ##
  ## Adjust for NRS
  ##
  ##

  adjust.str.ipd <- ""
  adjust.str.ad <- ""
  adjust.prior <- ""

  if (!is.null(method.bias)) {
    if (method.bias == "adjust1") {
      if (bias.type == 'add') {
        if (!is.null(v)) {
          if (trt.effect == "random") {
            adjust.str.ipd <- ""
            adjust.str.ad <- ""
          }
          else {
            adjust.str.ipd <- " + R[study[i]] * gamma[study[i]]"
            adjust.str.ad <- " + R[j + ns.ipd] * gamma[(j + ns.ipd)]"
          }
        }
        else {
          adjust.str.ipd <- " + R[study[i]] * gamma[study[i]]"
          adjust.str.ad <- " + R[j + ns.ipd] * gamma[(j + ns.ipd)]"
        }
      }
      else if (bias.type == 'mult') {
        adjust.str.ipd <- "* gamma[study[i]]^R[study[i]]"
        adjust.str.ad <- "* gamma[(j + ns.ipd)]^R[j + ns.ipd]"
      }
      else if (bias.type == 'both') {
        adjust.str.ipd <- "* gamma1[study[i]]^R[study[i]] + R[study[i]] * gamma2[study[i]]"
        adjust.str.ad <- "* gamma1[(j + ns.ipd)]^R[j + ns.ipd] + R[j + ns.ipd] * gamma2[(j + ns.ipd)]"
      }
      else
        stop("The bias type should be set as 'add', 'mult' or 'both'")


      if (bias.type == 'both') {
        gamma.effect <- ""
        if (bias.effect == 'random') {
          for (i in 1:2) {
            gamma.effect0 <-
              paste0("\n\n# Random effect for gamma (bias effect)",
                     ifelse(add.std.in, paste0("\nfor (j in std.in) {\n  gamma", i, "[j] ~ dnorm(g", i, ", prec.gamma", i, ")\n}"), ""),
                     ifelse(add.std.act.no, paste0("\nfor (j in std.act.no) {\n  gamma", i, "[j] ~ dnorm(0, prec.gamma", i, ")\n}"), ""),
                     ifelse(add.std.act.yes, paste0("\nfor (j in std.act.yes) {\n  gamma", i, "[j] ~ dnorm(g.act", i, ", prec.gamma", i, ")\n}"), ""),
                     ifelse(add.std.in, paste0("\ng", i, " ~ dnorm(0, 0.01)"), ""),
                     ifelse(add.std.act.yes, paste0("\ng.act", i, " ~ dnorm(0, 0.01)"), ""),
                     "\nprec.gamma", i, " <- pow(tau.gamma", i, ", -2)",
                     "\ntau.gamma", i, " ~ ", prior.tau.gamma)
            ##
            gamma.effect <- paste0(gamma.effect,gamma.effect0)
          }
        }
        else {
          for (i in 1:2) {
            gamma.effect0 <-
              paste0("\n\n# Common effect for gamma (bias effect)",
                     ifelse(add.std.in, paste0("\nfor (j in std.in) {\n  gamma", i, "[j] <- g", i, "\n  \n}"), ""),
                     ifelse(add.std.act.no, paste0("\nfor (j in std.act.no) {\n  gamma", i, "[j] <- 0\n}"), ""),
                     ifelse(add.std.act.yes, paste0("\nfor (j in std.act.yes) {\n  gamma", i, "[j] <- g.act", i, "\n}"), ""),
                     ifelse(add.std.in, paste0("\ng", i, " ~ dnorm(0, 0.01)"), ""),
                     ifelse(add.std.act.yes, paste0("\ng.act", i, " ~ dnorm(0, 0.01)"), "")
                     )
            ##
            gamma.effect <- paste0(gamma.effect,gamma.effect0)
          }
          ##
          gamma.effect <- paste0(gamma.effect, "\n  prec.gamma <- 0")
          message("Bias effect is assumed common across studies")
        }
      }
      else {
        if (!is.null(v)) { # with v
          gamma.effect <-
            paste0("\n\n# Common effect for gamma (bias effect)",
                   ifelse(add.std.in, "\nfor (j in std.in) {\n  gamma[j] <- g\n}", ""),
                   ifelse(add.std.act.no, "\nfor (j in std.act.no) {\n  gamma[j] <- 0\n}", ""),
                   ifelse(add.std.act.yes, "\nfor (j in std.act.yes) {\n  gamma[j] <- g.act\n}", ""),
                   ifelse(add.std.in, "\ng ~ dnorm(0, 0.01)", ""),
                   ifelse(add.std.act.yes, "\ng.act ~ dnorm(0, 0.01)", "")
                   )
        }
        else { # no v
          if (bias.effect == 'random') {
            gamma.effect <-
              paste0("\n\n# Random effect for gamma (bias effect)",
                     ifelse(add.std.in, "\nfor (j in std.in) {\n  gamma[j] ~ dnorm(g, prec.gamma)\n}", ""),
                     ifelse(add.std.act.no, "\nfor (j in std.act.no) {\n  gamma[j] ~ dnorm(0, prec.gamma)\n}", ""),
                     ifelse(add.std.act.yes, "\nfor (j in std.act.yes) {\n  gamma[j] ~ dnorm(g.act, prec.gamma)\n}", ""),
                     ifelse(add.std.in, "\ng ~ dnorm(0, 0.01)", ""),
                     ifelse(add.std.act.yes, "\ng.act ~ dnorm(0, 0.01)", ""),
                     ifelse(is.null(v), paste0("tau.gamma ~ ", prior.tau.gamma, "\n  prec.gamma <- pow(tau.gamma, -2)"), "")
                     )
          }
          else {
            gamma.effect <-
              paste0("\n\n# Common effect for gamma (bias effect)",
                     ifelse(add.std.in, "\nfor (j in std.in) {\n  gamma[j] <- g\n}", ""),
                     ifelse(add.std.act.no, "\nfor (j in std.act.no) {\n  gamma[j] <- 0\n}", ""),
                     ifelse(add.std.act.yes, "\nfor (j in std.act.yes) {\n  gamma[j] <- g.act\n}", ""),
                     ifelse(add.std.in, "\ng ~ dnorm(0, 0.01)", ""),
                     ifelse(add.std.act.yes, "\ng.act ~ dnorm(0, 0.01)", ""),
                     "\nprec.gamma <- 0"
                     )

            message("Bias effect is assumed common across studies")
          }
        }
      }
      if (is.null(bias.covariate)) {
        adjust.prior <-
          paste0(gamma.effect, "\n\n# Bias adjustment
for (j in 1:(ns.ipd + ns.ad)) {
  R[j] ~ dbern(pi[bias_index[j]])
}
pi[1] ~ ", prior.pi.high.rct,  # high RCT
"
pi[2] ~ ", prior.pi.low.rct, # low RCT
"
pi[3] ~ ", prior.pi.high.nrs, # high NRS
"
pi[4] ~ ", prior.pi.low.nrs,  # low NRS
"
pi[5] ~ dbeta(1, 1)")  # unclear RCT or NRS")
      }
      else {
        adjust.prior <-
          paste0(gamma.effect, "\n\n# Bias adjustment
for (j in 1:(ns.ipd + ns.ad)) {
  R[j] ~ dbern(pi[j])
  logit(pi[j]) <- a + b * xbias[j]
}
a ~ dnorm(0, 1e-2)
b ~ dnorm(0, 1e-2)")
      }
    }else if(method.bias=="adjust2"){
      if(is.null(bias.covariate)){ # assign priors for pi[bias_index[j]] based on design and RoB
        if(trt.effect=="random"){
          if(!is.null(v)){
            theta.effect.ipd <- "
        theta1[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]],precd[j,t.ipd[j,k]])
        theta2[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]]+gamma[j],precd[j,t.ipd[j,k]]/q)
        theta[j,t.ipd[j,k]] <- (1-pi[bias_index[j]])*theta1[j,t.ipd[j,k]]+pi[bias_index[j]]*theta2[j,t.ipd[j,k]]
        # multi-arm correction
      md[j,t.ipd[j,k]]<- mean[j,k] + sw[j,k]
      w[j,k]<- (theta[j,t.ipd[j,k]]  - mean[j,k])
      sw[j,k]<- sum(w[j,1:(k-1)])/(k-1)
      precd[j,t.ipd[j,k]]<- prec *2*(k-1)/k
        "

            theta.effect.ad <- "
        theta1[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]])
      theta2[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]]+gamma[j+ns.ipd],precd.ad[j,t.ad[j,k]]/q)
      theta[j+ns.ipd,t.ad[j,k]] <- (1-pi[bias_index[j]])*theta1[j+ns.ipd,t.ad[j,k]]+pi[bias_index[j]]*theta2[j+ns.ipd,t.ad[j,k]]
        # multi-arm correction
      md.ad[j,t.ad[j,k]]<- mean.ad[j,k] + sw.ad[j,k]
      w.ad[j,k]<- (theta[j+ns.ipd,t.ad[j,k]]  - mean.ad[j,k])
      sw.ad[j,k]<- sum(w.ad[j,1:(k-1)])/(k-1)
      precd.ad[j,t.ad[j,k]]<- prec *2*(k-1)/k
        "
            prior.tau.theta <- paste0("# heterogeneity between theta's
                                  tau ~",prior.tau.trt,
                                      "\n prec<- pow(tau,-2)")
            q.prior <- paste0("q~dbeta(",v,",1)")
          } else{


            # theta[j,t.ipd[j,k]] ~ dnormmix(c(md[j,t.ipd[j,k]],md[j,t.ipd[j,k]]+gamma[j]),c(precd[j,t.ipd[j,k]],precd[j,t.ipd[j,k]]+prec.gamma), c(1-pi[bias_index[j]],pi[bias_index[j]]))
            # theta[j+ns.ipd,t.ad[j,k]] ~ dnormmix(c(md.ad[j,t.ad[j,k]],md.ad[j,t.ad[j,k]]+gamma[j+ns.ipd]),c(precd.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]]+prec.gamma), c(1-pi[bias_index[j]],pi[bias_index[j]]))

            theta.effect.ipd <- "
        theta1[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]],precd[j,t.ipd[j,k]])
        theta2[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]]+gamma[j],precd[j,t.ipd[j,k]]+(prec.gamma*2*(k-1)/k))
        theta[j,t.ipd[j,k]] <- (1-pi[bias_index[j]])*theta1[j,t.ipd[j,k]]+pi[bias_index[j]]*theta2[j,t.ipd[j,k]]
        # multi-arm correction
      md[j,t.ipd[j,k]]<- mean[j,k] + sw[j,k]
      w[j,k]<- (theta[j,t.ipd[j,k]]  - mean[j,k])
      sw[j,k]<- sum(w[j,1:(k-1)])/(k-1)
      precd[j,t.ipd[j,k]]<- prec *2*(k-1)/k
        "

            theta.effect.ad <- "
        theta1[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]])
      theta2[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]]+gamma[j+ns.ipd],precd.ad[j,t.ad[j,k]]+(prec.gamma*2*(k-1)/k))
      theta[j+ns.ipd,t.ad[j,k]] <- (1-pi[bias_index[j]])*theta1[j+ns.ipd,t.ad[j,k]]+pi[bias_index[j]]*theta2[j+ns.ipd,t.ad[j,k]]
        # multi-arm correction
      md.ad[j,t.ad[j,k]]<- mean.ad[j,k] + sw.ad[j,k]
      w.ad[j,k]<- (theta[j+ns.ipd,t.ad[j,k]]  - mean.ad[j,k])
      sw.ad[j,k]<- sum(w.ad[j,1:(k-1)])/(k-1)
      precd.ad[j,t.ad[j,k]]<- prec *2*(k-1)/k
        "
            prior.tau.theta <- paste0("# heterogeneity between theta's
                                  tau ~",prior.tau.trt,
                                      "\n prec<- pow(tau,-2)")
          }
        }else if(trt.effect=="common"){
          theta.effect.ipd <- "theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]+(pi[bias_index[j]]*gamma[j])
        md[j,t.ipd[j,k]]<- mean[j,k]"
          theta.effect.ad <- "theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]+(pi[bias_index[j]]*gamma[j+ns.ipd])
        md.ad[j,t.ad[j,k]]<- mean.ad[j,k]"
          prior.tau.theta <- " "
        }else{
          stop("Please indicate the model of treatment effect as either 'random' or 'common' ")

        }

        if(!is.null(v)){ # with v
          gamma.effect <- paste0("# Common effect for gamma (bias effect)\n",
                                 ifelse(add.std.in,"for (j in std.in) {gamma[j]<-g}\n",""),
                                 ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]<-0}\n",""),
                                 ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]<-g.act}\n",""),
                                 ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                                 ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n","")
          )
        } else{ # no v
          if(bias.effect=='random'){
            gamma.effect <- paste0("# Random effect for gamma (bias effect)\n",
                                   ifelse(add.std.in,"for (j in std.in) {gamma[j]~dnorm(g,prec.gamma)}\n",""),
                                   ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]~dnorm(0,prec.gamma)}\n",""),
                                   ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]~dnorm(g.act,prec.gamma)}\n",""),
                                   ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                                   ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n",""),
                                   ifelse(is.null(v),paste0("tau.gamma~",prior.tau.gamma,"\n prec.gamma <- pow(tau.gamma,-2)"),"")

            )
          }else{
            gamma.effect <- paste0("# Common effect for gamma (bias effect)\n",
                                   ifelse(add.std.in,"for (j in std.in) {gamma[j]<-g}\n",""),
                                   ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]<-0}\n",""),
                                   ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]<-g.act}\n",""),
                                   ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                                   ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n",""),
                                   "prec.gamma <- 0"
            )
            warning("Bias effect is assumed common across studies")
          }
        }
        adjust.prior <- paste0(gamma.effect,"
                      # bias adjustment
                      pi[1]~",prior.pi.high.rct,  # high RCT
                               "
                      pi[2]~",prior.pi.low.rct, # low RCT
                               "
                      pi[3]~",prior.pi.high.nrs, # high NRS
                               "
                      pi[4]~",prior.pi.low.nrs,  # low NRS
                               "
                      pi[5]~ dbeta(1,1)"
        )# unclear RCT or NRS")
      }else{ # estimate pi[j] through a logistic model using the bias covariates
        if(trt.effect=="random"){
          if(!is.null(v)){
            theta.effect.ipd <- "
        theta1[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]],precd[j,t.ipd[j,k]])
        theta2[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]]+gamma[j],precd[j,t.ipd[j,k]]/q)
        theta[j,t.ipd[j,k]] <- (1-pi[j])*theta1[j,t.ipd[j,k]]+pi[j]*theta2[j,t.ipd[j,k]]
        # multi-arm correction
      md[j,t.ipd[j,k]]<- mean[j,k] + sw[j,k]
      w[j,k]<- (theta[j,t.ipd[j,k]]  - mean[j,k])
      sw[j,k]<- sum(w[j,1:(k-1)])/(k-1)
      precd[j,t.ipd[j,k]]<- prec *2*(k-1)/k
        "

            theta.effect.ad <- "
        theta1[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]])
      theta2[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]]+gamma[j+ns.ipd],precd.ad[j,t.ad[j,k]]/q)
      theta[j+ns.ipd,t.ad[j,k]] <- (1-pi[j])*theta1[j+ns.ipd,t.ad[j,k]]+pi[j]*theta2[j+ns.ipd,t.ad[j,k]]
        # multi-arm correction
      md.ad[j,t.ad[j,k]]<- mean.ad[j,k] + sw.ad[j,k]
      w.ad[j,k]<- (theta[j+ns.ipd,t.ad[j,k]]  - mean.ad[j,k])
      sw.ad[j,k]<- sum(w.ad[j,1:(k-1)])/(k-1)
      precd.ad[j,t.ad[j,k]]<- prec *2*(k-1)/k
        "
            prior.tau.theta <- paste0("# heterogeneity between theta's
                                  tau ~",prior.tau.trt,
                                      "\n prec<- pow(tau,-2)")
            q.prior <- paste0("q~dbeta(",v,",1)")
          } else{


            # theta[j,t.ipd[j,k]] ~ dnormmix(c(md[j,t.ipd[j,k]],md[j,t.ipd[j,k]]+gamma[j]),c(precd[j,t.ipd[j,k]],precd[j,t.ipd[j,k]]+prec.gamma), c(1-pi[j],pi[j]))
            # theta[j+ns.ipd,t.ad[j,k]] ~ dnormmix(c(md.ad[j,t.ad[j,k]],md.ad[j,t.ad[j,k]]+gamma[j+ns.ipd]),c(precd.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]]+prec.gamma), c(1-pi[j],pi[j]))

            theta.effect.ipd <- "
        theta1[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]],precd[j,t.ipd[j,k]])
        theta2[j,t.ipd[j,k]] ~dnorm(md[j,t.ipd[j,k]]+gamma[j],precd[j,t.ipd[j,k]]+(prec.gamma*2*(k-1)/k))
        theta[j,t.ipd[j,k]] <- (1-pi[j])*theta1[j,t.ipd[j,k]]+pi[j]*theta2[j,t.ipd[j,k]]
        # multi-arm correction
      md[j,t.ipd[j,k]]<- mean[j,k] + sw[j,k]
      w[j,k]<- (theta[j,t.ipd[j,k]]  - mean[j,k])
      sw[j,k]<- sum(w[j,1:(k-1)])/(k-1)
      precd[j,t.ipd[j,k]]<- prec *2*(k-1)/k
        "

            theta.effect.ad <- "
        theta1[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]])
      theta2[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]]+gamma[j+ns.ipd],precd.ad[j,t.ad[j,k]]+(prec.gamma*2*(k-1)/k))
      theta[j+ns.ipd,t.ad[j,k]] <- (1-pi[j])*theta1[j+ns.ipd,t.ad[j,k]]+pi[j]*theta2[j+ns.ipd,t.ad[j,k]]
        # multi-arm correction
      md.ad[j,t.ad[j,k]]<- mean.ad[j,k] + sw.ad[j,k]
      w.ad[j,k]<- (theta[j+ns.ipd,t.ad[j,k]]  - mean.ad[j,k])
      sw.ad[j,k]<- sum(w.ad[j,1:(k-1)])/(k-1)
      precd.ad[j,t.ad[j,k]]<- prec *2*(k-1)/k
        "
            prior.tau.theta <- paste0("# heterogeneity between theta's
                                  tau ~",prior.tau.trt,
                                      "\n prec<- pow(tau,-2)")
          }
        }else if(trt.effect=="common"){
          theta.effect.ipd <- "theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]+(pi[j]*gamma[j])
        md[j,t.ipd[j,k]]<- mean[j,k]"
          theta.effect.ad <- "theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]+(pi[j]*gamma[j+ns.ipd])
        md.ad[j,t.ad[j,k]]<- mean.ad[j,k]"
          prior.tau.theta <- " "
        }else{
          stop("Please indicate the model of treatment effect as either 'random' or 'common' ")

        }

        if(!is.null(v)){ # with v
          gamma.effect <- paste0("# Common effect for gamma (bias effect)\n",
                                 ifelse(add.std.in,"for (j in std.in) {gamma[j]<-g}\n",""),
                                 ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]<-0}\n",""),
                                 ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]<-g.act}\n",""),
                                 ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                                 ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n","")
          )
        } else{ # no v
          if(bias.effect=='random'){
            gamma.effect <- paste0("# Random effect for gamma (bias effect)\n",
                                   ifelse(add.std.in,"for (j in std.in) {gamma[j]~dnorm(g,prec.gamma)}\n",""),
                                   ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]~dnorm(0,prec.gamma)}\n",""),
                                   ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]~dnorm(g.act,prec.gamma)}\n",""),
                                   ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                                   ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n",""),
                                   ifelse(is.null(v),paste0("tau.gamma~",prior.tau.gamma,"\n prec.gamma <- pow(tau.gamma,-2)"),"")

            )
          }else{
            gamma.effect <- paste0("# Common effect for gamma (bias effect)\n",
                                   ifelse(add.std.in,"for (j in std.in) {gamma[j]<-g}\n",""),
                                   ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]<-0}\n",""),
                                   ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]<-g.act}\n",""),
                                   ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                                   ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n",""),
                                   "prec.gamma <- 0"
            )
            warning("Bias effect is assumed common across studies")
          }
        }

        adjust.prior <- paste0(gamma.effect,"
                      # bias adjustment
                      for (j in 1:(ns.ipd+ns.ad)) {logit(pi[j]) <- a+b*xbias[j]}
                               a~dnorm(0,1e-2)
                               b~dnorm(0,1e-2)"
        )
      }
    }
    else if (method.bias == "prior") {
      d.prior <- d.prior.nrs
    }
    else if (method.bias == 'naive') {
      message("Both designs are combined naively without acknowledging ",
              "design differences")
    }
  }
  else
    message("The data is analyzed assuming the studies have the same design")

  # prior for trial baseline
if(sm=="RR"){
  prior.u <- "for (j in 1:(ns.ipd + ns.ad)) {
  u[j] <- log(p.b[j])
  p.b[j] ~ dunif(0, 1)}"
} else {
  prior.u <-"for (j in 1:(ns.ipd + ns.ad)) {u[j] ~ dnorm(0, .01)}"
}


  #! #-------------------------------------------------#
  #------  Set up the likelihood and the link  ------#
  #-------------------------------------------------#
  if(sm=="OR"){
  like.str.ipd <- " y[i]~dbern(p[i]) # bernoulli likelihood"
  link.str.ipd <- "logit(p[i]) <- u[study[i]]+theta[study[i],trt[i]]"
  like.str.ad <- " r[j,k] ~ dbin(pa[j,t.ad[j,k]],n[j,k]) # binomial likelihood of number of events"
  link.str.ad.ref <- "logit(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm"
  link.str.ad <-"logit(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]])"

  }  else if(sm=="RR"){
    like.str.ipd <- " y[i]~dbern(p[i]) "
    link.str.ipd <- "log(p[i]) <- u[study[i]]+theta[study[i],trt[i]]"
    like.str.ad <- " r[j,k] ~ dbin(pa[j,t.ad[j,k]],n[j,k]) # binomial likelihood of number of events"
    link.str.ad.ref <- "log(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm"
    link.str.ad <-"log(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]])"

  } else if(sm=="MD") {
  like.str.ipd <- " y[i]~dnorm(delta[i],prec.delta.ipd[study[i],trt[i]]) "
  link.str.ipd <- "delta[i] <- u[study[i]]+(theta[study[i],trt[i]])
  prec.delta.ipd[study[i],trt[i]] <- pow(sd[study[i],trt[i]],-2)"

  like.str.ad <- " ybar[j,k] ~ dnorm(delta.ad[j,t.ad[j,k]],prec.delta.ad[j,k])
  prec.delta.ad[j,k] <- pow(se[j,k],-2)"
  link.str.ad.ref <- "delta.ad[j,t.ad[j,1]] <- u[j+ns.ipd]"
  link.str.ad <-"delta.ad[j,t.ad[j,k]] <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]])"

  }else if(sm=="SMD") {
  like.str.ipd <- " y[i]~dnorm(phi[i],prec.delta.ad[study[i],trt[i]])"
  link.str.ipd <- "phi[i] <- delta[i]*s.pool.ipd[study[i]]
  delta[i] <- u[study[i]]+theta[study[i],trt[i]]
  prec.delta.ad[study[i],trt[i]] <- pow(sd[study[i],trt[i]],-2)"

  like.str.ad <- " ybar[j,k] ~ dnorm(phi.ad[j,t.ad[j,k]],prec.delta.ad[j,k])
  phi.ad[j,t.ad[j,k]] <- delta.ad[j,t.ad[j,k]]*s.pool.ad[j]
  prec.delta.ad[j,k] <- pow(se[j,k],-2)"
  link.str.ad.ref <- "delta.ad[j,t.ad[j,1]] <- u[j+ns.ipd]*s.pool.ad[j] "
  link.str.ad <-"delta.ad[j,t.ad[j,k]] <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]])"
  }

  ##
  ##
  ## Combine the code
  ##
  ##


  #-------------------------------#
  #------ IPD part


  ipd.code <- sprintf("
#
# (1) IPD part
#

# Loop through individuals
for (i in 1:np) {
  #  likelihood
  %s
  # link function
  %s %s %s
}

# Loop through IPD studies
for (j in 1:(ns.ipd)) {
  # Multi-arm correction is zero for reference arm
  w[j, 1] <- 0
  # Treatment effect is zero for reference arm
  theta[j, t.ipd[j, 1]] <- 0%s
  # Loop through non-referent IPD arms
  for (k in 2:na.ipd[j]) {
    # Synthesize relative treatment effects
    %s
    # Consistency equation
    mean[j, k] <- d[t.ipd[j, k]] - d[t.ipd[j, 1]]%s%s
  }
}",like.str.ipd, link.str.ipd,adjust.str.ipd, metareg.str.ipd, ref.trt.effect.ipd, theta.effect.ipd, betab.consis.ipd, betaw.consis.ipd)


  ad.code <- sprintf("

#
# (2) AD part
#

# Loop through AD studies
for (j in 1:ns.ad) {
  # Multi-arm correction is zero for referent arm
  w.ad[j, 1] <- 0
  # Treatment effect is zero for referent arm
  theta[j + ns.ipd, t.ad[j, 1]] <- 0%s
  # Loop through AD arms
  for (k in 1:na.ad[j]) {
    # likelihood
    %s
  }
  # effect in referent arm
  %s
  # Loop through non-referent AD arms
  for (k in 2:na.ad[j]) {
    # link function
    %s %s %s
    # Synthesize relative treatment effects
    %s
    # Consistency equations
    mean.ad[j, k] <- d[t.ad[j, k]] - d[t.ad[j, 1]]%s
  }
}",ref.trt.effect.ad,like.str.ad, link.str.ad.ref,link.str.ad, adjust.str.ad,  metareg.str.ad, theta.effect.ad, betab.consis.ad)


  prior.code <- sprintf("

#
# (3) Priors
#

%s
%s

# effect is zero at network reference treatment
d[1] <- 0%s

%s%s%s%s%s%s",
prior.u,prior.tau.theta, d.prior, beta0.prior.ipd, betab.prior, betaw.prior.ipd, beta.prior, adjust.prior, q.prior)

  ad.code <- ifelse(ad, ad.code, "")
  ipd.code <- ifelse(ipd, ipd.code, "")
  code.str <- paste0('model {', ipd.code, ad.code, prior.code, '\n}\n')

  return(code.str)
}
