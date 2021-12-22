crossnma.code <- function(ipd = T,
                         ad = T,
                         trt.effect='random',
                         prior.tau.trt=NULL,
                         # -------- meta-regression
                         split.regcoef =F,
                         covariate=NULL,

                         reg0.effect='random',
                         regb.effect='random',
                         regw.effect='random',

                         prior.tau.reg0=NULL,
                         prior.tau.regb=NULL,
                         prior.tau.regw=NULL,
                         # --------  bias adjustment
                         bias.effect=NULL,
                         bias.type=NULL,
                         bias.covariate=NULL,
                         add.std.in=NULL,
                         add.std.act.no=NULL,
                         add.std.act.yes=NULL,

                         prior.tau.gamma=NULL,

                         prior.pi.high.rct=NULL,
                         prior.pi.low.rct=NULL,
                         prior.pi.high.nrs=NULL,
                         prior.pi.low.nrs=NULL,

                         method.bias = NULL,
                         d.prior.nrs=NULL # required when method.bias='prior'
) {

  #-----------------------------------------#
  #------ User priors ------#
  #-----------------------------------------#
  prior.tau.trt <- ifelse(is.null(prior.tau.trt),'dunif(0,2)',prior.tau.trt)
  prior.tau.reg0 <- ifelse(is.null(prior.tau.reg0),'dunif(0,2)',prior.tau.reg0)
  prior.tau.regb <- ifelse(is.null(prior.tau.regb),'dunif(0,2)',prior.tau.regb)
  prior.tau.regw <- ifelse(is.null(prior.tau.regw),'dunif(0,2)',prior.tau.regw)
  prior.tau.gamma<- ifelse(is.null(prior.tau.gamma),'dunif(0,2)',prior.tau.gamma )
  prior.pi.high.rct <- ifelse(is.null(prior.pi.high.rct),'dbeta(10,1)',prior.pi.high.rct )
  prior.pi.low.rct <- ifelse(is.null(prior.pi.low.rct),'dbeta(1,10)',prior.pi.low.rct )
  prior.pi.high.nrs <- ifelse(is.null(prior.pi.high.nrs),'dbeta(30,1)',prior.pi.high.nrs)
  prior.pi.low.nrs <- ifelse(is.null(prior.pi.low.nrs),'dbeta(1,30)',prior.pi.low.nrs)

  #-----------------------------------------#
  #------ meta regression ------#
  #-----------------------------------------#
  metareg.str.ipd <- ""
  betab.consis.ipd <- ""
  betaw.consis.ipd <- ""
  beta0.prior.ipd <- ""
  betab.prior <- ""
  betaw.prior.ipd <- ""
  beta.prior.ipd <- ""

  beta.prior.ad <- ""
  metareg.str.ad <- ""
  betab.consis.ad <- ""

  if (!is.null(covariate)) {
    if(ipd) {
      for (i in 1:length(covariate[[1]])) {
        # meta-regression terms - Up to 3
        metareg.str.ipd0 <- paste0(
          "+beta0_",i,"[study[i]]*(x",i,"[i])+betaw_",i,"[study[i],trt[i]]*(x",i,"[i]-xm",i,".ipd[i])+betab_",i,"[study[i],trt[i]]*xm",i,".ipd[i]"
        )
        metareg.str.ipd <- paste0(metareg.str.ipd,metareg.str.ipd0)
        # consistency equations for beta_b and beta_w - Up to 3
        betab.consis.ipd0 <- paste0("betab_",i,"[j,t.ipd[j,k]] <- betab.t_",i,"[t.ipd[j,k]] - betab.t_",i,"[t.ipd[j,1]]")
        betaw.consis.ipd0 <- paste0("betaw_",i,"[j,t.ipd[j,k]] <- betaw.t_",i,"[t.ipd[j,k]] - betaw.t_",i,"[t.ipd[j,1]]")
        betab.consis.ipd <- paste0(betab.consis.ipd," \n ",betab.consis.ipd0)
        betaw.consis.ipd <- paste0(betaw.consis.ipd," \n ",betaw.consis.ipd0)
      }

      # beta0
      if(reg0.effect=='random'){
        for (i in 1:length(covariate[[1]])) {
          beta0.prior.ipd0 <- paste0("
                              # Random effect for beta0
                              for (j in 1:(ns.ipd)) {
                              beta0_",i,"[j] ~ dnorm(b0_",i,",prec.beta0_",i,")
                              }\n
                              b0_",i," ~ dnorm(0,.01)\n
                              prec.beta0_",i," <- pow(tau.b0_",i,", -2)\n
                              tau.b0_",i,"~",prior.tau.reg0
          )
          beta0.prior.ipd <- paste0(beta0.prior.ipd,beta0.prior.ipd0)
        }
      }else if(reg0.effect=='independent') {
        for (i in 1:length(covariate[[1]])) {
          beta0.prior.ipd0 <- paste0("
                              # Independent effect for beta0
                              for (j in 1:(ns.ipd)) {
                              beta0_",i,"[j] <- b0_",i,"[j] \n
                              b0_",i,"[j]  ~ dnorm(0,.01)
                              }")
          beta0.prior.ipd <- paste0(beta0.prior.ipd,beta0.prior.ipd0)
        }
      }else{
        stop("The progonostic effect can be assumed either 'independent' or 'random' across studies")
      }

      # betab and betaw

      if(!split.regcoef) { # not splitted within and between- study covariate
        if(regb.effect=='random'||regw.effect=='random'){
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ipd0 <- paste0("
           # Random effects for beta (within=between)
           beta.t_",i,"[1] <- 0 \n
           for(k in 1:nt){
           betab.t_",i,"[k] <- beta.t_",i,"[k] \n
           betaw.t_",i,"[k] <- beta.t_",i,"[k]
           } \n
           for (k in 2:nt){
           beta.t_",i,"[k]~dnorm(b_",i,",prec.beta_",i,")
           } \n
           b_",i,"~dnorm(0,1e-2) \n
           tau.b_",i,"~",prior.tau.regw,
           "\n prec.beta_",i," <- pow(tau.b_",i,",-2)
                            ")
            beta.prior.ipd <- paste0(beta.prior.ipd,beta.prior.ipd0)
          }
        }else if(regb.effect=='common'&regw.effect=='common') {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ipd0 <- paste0("
            # Common effect for beta (within=between)
            betab.t_",i,"[1] <- 0 \n
                              betaw.t_",i,"[1] <- 0 \n
                              for (k in 2:nt) {
                              betab.t_",i,"[k]<-b_",i," \n
                              betaw.t_",i,"[k]<-b_",i,"
          } \n
          b_",i,"~dnorm(0,1e-2)
                                     ")
            beta.prior.ipd <- paste0(beta.prior.ipd,beta.prior.ipd0)
          }
        }else{
          stop("The covariate effect need to be assumed 'random' or 'common' across studies")

        }
      } else{ # splitted within and between- study covariate
        # between- study covariate
        if(regb.effect=='random'){
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <- paste0("
            # Random effect for betab (the between-study covariate effect)
            betab.t_",i,"[1] <- 0 \n
                     for (k in 2:nt){
                     betab.t_",i,"[k]~dnorm(bb_",i,",prec.betab_",i,")
                     } \n
                    bb_",i,"~dnorm(0,1e-2) \n
                    tau.bb_",i,"~",prior.tau.regb,
                                   "\n prec.betab_",i," <- pow(tau.bb_",i,",-2)
                                 ")
            betab.prior <- paste0(betab.prior,betab.prior0)
          }
        }else if(regb.effect=='common') {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <- paste0("
            # Common effect for betab (the between-study covariate effect)
            betab.t_",i,"[1] <- 0 \n
            for (k in 2:nt) {
            betab.t_",i,"[k]<-bb_",i,"
            } \n
            bb_",i,"~dnorm(0,1e-2)
                                   ")
            betab.prior <- paste0(betab.prior,betab.prior0)
          }
        }else{
          stop("The between-study covariate effect need to be assumed 'random' or 'common' across studies")
        }

        # within- study covariate
        if(regw.effect=='random'){
          for (i in 1:length(covariate[[1]])) {
            betaw.prior.ipd0 <- paste0("
            # Random effect for betaw (the within-study covariate effect)
            betaw.t_",i,"[1] <- 0 \n
            for (k in 2:nt){
            betaw.t_",i,"[k]~dnorm(bw_",i,",prec.betaw_",i,")
            } \n
            bw_",i,"~dnorm(0,1e-2) \n
            prec.betaw_",i,"<- pow(tau.bw_",i,",-2)",                          "
            \n tau.bw_",i,"~",prior.tau.regw)
            betaw.prior.ipd <- paste0(betaw.prior.ipd,betaw.prior.ipd0)
          }
        }else if(regw.effect=='common') {
          for (i in 1:length(covariate[[1]])) {
            betaw.prior.ipd0 <- paste0("
          # Common effect for betaw (the within-study covariate effect)
          betaw.t_",i,"[1] <- 0 \n
          for (k in 2:nt){
          betaw.t_",i,"[k]<-bw_",i,"
          } \n
          bw_",i,"~dnorm(0,1e-2)
                                       ")
          betaw.prior.ipd <- paste0(betaw.prior.ipd,betaw.prior.ipd0)
          }
        }else{
          stop("The within-study covariate effect can be assumed 'random' or 'common' across studies")

        }
      }
    }

    if(ad){
      for (i in 1:length(covariate[[1]])) {
        # meta-regression terms - up to 3
        metareg.str.ad0 <- paste0("+betab.ad_",i,"[j,t.ad[j,k]]*xm",i,".ad[j]")
        betab.consis.ad0 <- paste0("betab.ad_",i,"[j,t.ad[j,k]] <- betab.t_",i,"[t.ad[j,k]] - betab.t_",i,"[t.ad[j,1]]")

        # consistency equation - up to 3
        metareg.str.ad <- paste0(metareg.str.ad,metareg.str.ad0)
        betab.consis.ad <- paste0(betab.consis.ad," \n ",betab.consis.ad0)
      }
      if(!split.regcoef) { # not splitted
        if(regb.effect=='random'||regw.effect=='random'){
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ad0 <- paste0("
            # Random effects for beta (within=between)
            beta.t_",i,"[1] <- 0 \n
           for(k in 1:nt){
           betab.t_",i,"[k] <- beta.t_",i,"[k]
           } \n
           for (k in 2:nt){
           beta.t_",i,"[k]~dnorm(b_",i,",prec.beta_",i,")
           } \n
           b_",i,"~dnorm(0,1e-2) \n
           tau.b_",i,"~",prior.tau.regb,
           "\n prec.beta_",i," <- pow(tau.b_",i,",-2)
           ")
            beta.prior.ad <- paste0(beta.prior.ad,beta.prior.ad0)
          }
        }else if(regb.effect=='common'&regw.effect=='common') {
          for (i in 1:length(covariate[[1]])) {
            beta.prior.ad0 <- paste0("
            # Common effects for beta (within=between)
            betab.t_",i,"[1] <- 0 \n
            for (k in 2:nt) {
            betab.t_",i,"[k]<-b_",i,"
            } \n
            b_",i,"~dnorm(0,1e-2)
                                    ")
            beta.prior.ad <- paste0(beta.prior.ad,beta.prior.ad0)
          }
        }else{
          stop("The between-study covariate effect need to be assumed 'random' or 'common' across studies")

        }
      } else { # splitted
        betab.prior <- ""
        if(regb.effect=='random'){
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <- paste0("
            # Random effects for betab (the between-study covariate effect)
            betab.t_",i,"[1] <- 0 \n
            for (k in 2:nt){
            betab.t_",i,"[k]~dnorm(bb_",i,",prec.betab_",i,")
            } \n
            bb_",i,"~dnorm(0,1e-2) \n
            tau.bb_",i,"~",prior.tau.regb,
            "\n prec.betab_",i," <- pow(tau.bb_",i,",-2)
            ")
            betab.prior <- paste0(betab.prior,betab.prior0)
          }
        }else if(regb.effect=='common') {
          for (i in 1:length(covariate[[1]])) {
            betab.prior0 <- paste0("
            # Random effects for betab (the between-study covariate effect)
            betab.t_",i,"[1] <- 0 \n
            for (k in 2:nt){
            betab.t_",i,"[k] <- bb_",i,"
            } \n
            bb_",i,"~dnorm(0,1e-2)
                                ")
            betab.prior <- paste0(betab.prior,betab.prior0)
          }

        }else{
          stop("The between-study covariate effect need to be assumed 'random' or 'common' across studies")

        }
      }
    }
  }


  if(!split.regcoef){
    if(ipd){
    beta.prior <- beta.prior.ipd
    }else{
      beta.prior <-beta.prior.ad
    }
  }else {
    beta.prior <- ""
  }

  # treatment effect is zero for reference treatment
  ref.trt.effect.ipd <- ""
  ref.trt.effect.ad <- ""

  if(!is.null(covariate)){
  if(ipd){
    for (i in 1:length(covariate[[1]])) {
      # meta-regression terms - Up to 3
      ref.trt.effect.ipd0 <- paste0(
        "betaw_",i,"[j,t.ipd[j,1]] <- 0 \n","betab_",i,"[j,t.ipd[j,1]] <- 0 \n"
        )
      ref.trt.effect.ipd <- paste0(ref.trt.effect.ipd,ref.trt.effect.ipd0)
    }
  }
  if(ad){
    for (i in 1:length(covariate[[1]])) {
      # meta-regression terms - Up to 3
      ref.trt.effect.ad0 <- paste0(
        "betab_",i,"[j,t.ipd[j,1]] <- 0 \n"
      )
      ref.trt.effect.ad0 <- ifelse(ipd,"", ref.trt.effect.ad0)
      ref.trt.effect.ad <- paste0(ref.trt.effect.ad,ref.trt.effect.ad0)
    }
  }
    }

    # Relative treatment effect
  if(trt.effect=="random"){
    theta.effect.ipd <- "
    theta[j,t.ipd[j,k]] ~ dnorm(md[j,t.ipd[j,k]],precd[j,t.ipd[j,k]])
  # multi-arm correction
      md[j,t.ipd[j,k]]<- mean[j,k] + sw[j,k]
      w[j,k]<- (theta[j,t.ipd[j,k]]  - mean[j,k])
      sw[j,k]<- sum(w[j,1:(k-1)])/(k-1)
      precd[j,t.ipd[j,k]]<- prec *2*(k-1)/k"
    theta.effect.ad <- "theta[j+ns.ipd,t.ad[j,k]] ~ dnorm(md.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]])
  # multi-arm correction
      md.ad[j,t.ad[j,k]]<- mean.ad[j,k] + sw.ad[j,k]
      w.ad[j,k]<- (theta[j+ns.ipd,t.ad[j,k]]  - mean.ad[j,k])
      sw.ad[j,k]<- sum(w.ad[j,1:(k-1)])/(k-1)
      precd.ad[j,t.ad[j,k]]<- prec *2*(k-1)/k"
    prior.tau.theta <- paste0("
    # heterogeneity between theta's
                              tau ~",prior.tau.trt,
                            "\n prec<- pow(tau,-2)")
  }else if(trt.effect=="common"){
    theta.effect.ipd <- "
    theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]
    md[j,t.ipd[j,k]]<- mean[j,k]"
    theta.effect.ad <- "
    theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]
    md.ad[j,t.ad[j,k]]<- mean.ad[j,k]"
    prior.tau.theta <- ""

  }else{
    stop("Please indicate the treatment effect model as either 'random' or 'common' ")
  }

  d.prior <- "for(k in 2:nt) {d[k] ~ dnorm(0,.01)}"

  #-------------------------------#
  #------ adjust for NRS ------#
  #-------------------------------#

  adjust.str.ipd <- ""
  adjust.str.ad <- ""
  adjust.prior <- ""
  if(!is.null(method.bias)){
    if(method.bias=="adjust1"){
      if(bias.type=='add'){
        adjust.str.ipd <- "+R[study[i]]*gamma[study[i]]"
        adjust.str.ad <- "+R[j+ns.ipd]*gamma[(j+ns.ipd)]"
      } else if(bias.type=='mult'){
        adjust.str.ipd <- "*gamma[study[i]]^R[study[i]]"
        adjust.str.ad <- "*gamma[(j+ns.ipd)]^R[j+ns.ipd]"
      } else if(bias.type=='both') {
        adjust.str.ipd <- "*gamma1[study[i]]^R[study[i]]+R[study[i]]*gamma2[study[i]]"
        adjust.str.ad <- "*gamma1[(j+ns.ipd)]^R[j+ns.ipd]+R[j+ns.ipd]*gamma2[(j+ns.ipd)]"
      } else{
        stop("The bias type should be set as 'add', 'mult' or 'both'")
      }
      if(bias.type=='both'){
        gamma.effect <- ""
        if(bias.effect=='random'){
          for (i in 1:2) {
            gamma.effect0 <- paste0("# Random effect for gamma (bias effect)\n",
                                    ifelse(add.std.in,paste0("for (j in std.in) {gamma",i,"[j]~dnorm(g",i,",prec.gamma",i,")}\n"),""),
                                    ifelse(add.std.act.no,paste0("for (j in std.act.no) {gamma",i,"[j]~dnorm(0,prec.gamma",i,")}\n"),""),
                                    ifelse(add.std.act.yes,paste0("for (j in std.act.yes) {gamma",i,"[j]~dnorm(g.act",i,",prec.gamma",i,")}\n"),""),
                                    ifelse(add.std.in,paste0("g",i,"~dnorm(0, 0.01)\n"),""),
                                    ifelse(add.std.act.yes,paste0("g.act",i,"~dnorm(0, 0.01)\n"),""),
                                    "prec.gamma",i," <- pow(tau.gamma",i,",-2)",
                                    "\n tau.gamma",i,"~",prior.tau.gamma)
            gamma.effect <- paste0(gamma.effect,gamma.effect0)
          }
        }else{
          for (i in 1:2) {
            gamma.effect0 <- paste0("# Common effect for gamma (bias effect)\n",
                                    ifelse(add.std.in,paste0("for (j in std.in) {gamma",i,"[j]<-g",i,"}\n"),""),
                                    ifelse(add.std.act.no,paste0("for (j in std.act.no) {gamma",i,"[j]<-0}\n"),""),
                                    ifelse(add.std.act.yes,paste0("for (j in std.act.yes) {gamma",i,"[j]<-g.act",i,"}\n"),""),
                                    ifelse(add.std.in,paste0("g",i,"~dnorm(0, 0.01)\n"),""),
                                    ifelse(add.std.act.yes,paste0("g.act",i,"~dnorm(0, 0.01)\n"),"")
            )
            gamma.effect <- paste0(gamma.effect,gamma.effect0)
          }
          gamma.effect <- paste0(gamma.effect,"\n prec.gamma <- 0")
          warning("Bias effect is assumed common across studies")
        }
      }else {
        if(bias.effect=='random'){
          gamma.effect <- paste0("# Random effect for gamma (bias effect)\n",
                                  ifelse(add.std.in,"for (j in std.in) {gamma[j]~dnorm(g,prec.gamma)}\n",""),
                                  ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]~dnorm(0,prec.gamma)}\n",""),
                                  ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]~dnorm(g.act,prec.gamma)}\n",""),
                                  ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                                  ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n",""),
                                  "tau.gamma~",prior.tau.gamma,
                                  "\n prec.gamma <- pow(tau.gamma,-2)"
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
      if(is.null(bias.covariate)){
        adjust.prior <- paste0(gamma.effect,"
                      # bias adjustment
                      for (j in 1:(ns.ipd+ns.ad)) {R[j]~dbern(pi[bias_index[j]])}
                      pi[1]~",prior.pi.high.rct,  # high RCT
                               "
                      pi[2]~",prior.pi.low.rct, # low RCT
                               "
                      pi[3]~",prior.pi.high.nrs, # high NRS
                               "
                      pi[4]~",prior.pi.low.nrs,  # low NRS
                               "
                      pi[5]~dbeta(1,1)")  # unclear RCT or NRS")
      }else{
        adjust.prior <- paste0(gamma.effect,"
                      # bias adjustment
                      for (j in 1:(ns.ipd+ns.ad)) {R[j]~dbern(pi[j])
                      logit(pi[j]) <- a+b*xbias[j]}
                               a~dnorm(0,1e-2)
                               b~dnorm(0,1e-2)"
        )
      }

    }else if(method.bias=="adjust2"){
      if(trt.effect=="random"){
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
      }else if(trt.effect=="common"){
        theta.effect.ipd <- "theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]+(pi[bias_index[j]]*gamma[j])
        md[j,t.ipd[j,k]]<- mean[j,k]"
        theta.effect.ad <- "theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]+(pi[bias_index[j]]*gamma[j+ns.ipd])
        md.ad[j,t.ad[j,k]]<- mean.ad[j,k]"
        prior.tau.theta <- " "
      }else{
        stop("Please indicate the model of treatment effect as either 'random' or 'common' ")

      }

      if(bias.effect=='random'){
        gamma.effect <- paste0("# Random effect for gamma (bias effect)\n",
                               ifelse(add.std.in,"for (j in std.in) {gamma[j]~dnorm(g,prec.gamma)}\n",""),
                               ifelse(add.std.act.no,"for (j in std.act.no) {gamma[j]~dnorm(0,prec.gamma)}\n",""),
                               ifelse(add.std.act.yes,"for (j in std.act.yes) {gamma[j]~dnorm(g.act,prec.gamma)}\n",""),
                               ifelse(add.std.in,"g~dnorm(0, 0.01)\n",""),
                               ifelse(add.std.act.yes,"g.act~dnorm(0, 0.01)\n",""),
                               "tau.gamma~",prior.tau.gamma,
                               "\n prec.gamma <- pow(tau.gamma,-2)"
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
      if(is.null(bias.covariate)){
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
      }else{
        adjust.prior <- paste0(gamma.effect,"
                      # bias adjustment
                      for (j in 1:(ns.ipd+ns.ad)) {logit(pi[j]) <- a+b*xbias[j]}
                               a~dnorm(0,1e-2)
                               b~dnorm(0,1e-2)"
        )
      }
    }else if(method.bias=="prior"){
      d.prior <- d.prior.nrs
    }else if(method.bias=='naive'){
      warning("Both designs are combined naively without acknowledging the design differences")
    }
  }else{
    message("The data is analyzed assuming the studies has the same design")
  }


  #### Combine the code

  #-------------------------------#
  #------ IPD part


  ipd.code <- sprintf("
  #*** IPD part
  for (i in 1:np) {  # loop through individuals
    y[i]~dbern(p[i]) # binomial likelihood
    logit(p[i]) <- u[study[i]]+theta[study[i],trt[i]]%s%s # logistic transformation - to estimate Odds Ratio (OR)
  }

  for(j in 1:(ns.ipd)){      # loop through IPD studies
     w[j,1]<- 0              # multi-arm correction is zero for reference arm
     theta[j,t.ipd[j,1]]<- 0 # treatment effect is zero for reference arm
     %s
     for (k in 2:na.ipd[j]){ # loop through non-referent IPD arms

      # Synthesize relative treatment effects
      %s

      # consistency equation
      mean[j,k] <-d[t.ipd[j,k]] - d[t.ipd[j,1]]
      %s
      %s
     }
    }
", adjust.str.ipd,metareg.str.ipd,ref.trt.effect.ipd,theta.effect.ipd, betab.consis.ipd, betaw.consis.ipd)

  #-------------------------------#
  #------ AD part

  ad.code <- sprintf("
  #*** AD part
  for (j in 1:ns.ad) {            # loop through AD studies
    w.ad[j,1]<- 0                 # multi-arm correction is zero for referent arm
    theta[j+ns.ipd,t.ad[j,1]]<- 0 # treatment effect is zero for referent arm
    %s
    for (k in 1:na.ad[j]) {                 # loop through AD arms
      r[j,k] ~ dbin(pa[j,t.ad[j,k]],n[j,k]) # binomial likelihood of number of events
          }
    logit(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm

    for (k in 2:na.ad[j]){                  # loop through non-referent AD arms
      logit(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]])%s%s # logistic transformation Log Odds Ratio

      # Synthesize relative treatment effects
      %s
      # consistency equations
      mean.ad[j,k] <-d[t.ad[j,k]] - d[t.ad[j,1]]
      %s
    }
}",ref.trt.effect.ad,adjust.str.ad, metareg.str.ad,theta.effect.ad, betab.consis.ad)

  #-------------------------------#
  #------ priors
  prior.code <- sprintf("
  #*** PRIORS
  for (j in 1:(ns.ipd+ns.ad)) {u[j] ~ dnorm(0,.01)} # log-odds in referent arm
  %s
  d[1] <- 0 # treatment effect is zero for reference treatment
  %s
  # Compute OR for each comparison
  for(i in 1:(nt-1)) {
    for (j in (i+1):nt) {
      OR[j,i]<- exp(d[j] - d[i])
      LOR[j,i]<- d[j] - d[i]}
      }
  %s
  %s
  %s
  %s
  %s
  ",prior.tau.theta,d.prior,beta0.prior.ipd, betab.prior, betaw.prior.ipd,beta.prior,adjust.prior)
  ad.code <- ifelse(ad, ad.code, "")
  ipd.code <- ifelse(ipd, ipd.code, "")
  code.str <- paste0('model {',ipd.code, ad.code, prior.code,'}')
  return(code.str)
}

