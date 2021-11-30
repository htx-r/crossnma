rhat <- function(x,
                           confidence = 0.95,
                           transform = FALSE,
                           autoburnin = FALSE,
                           multivariate = TRUE
){
  x <- as.mcmc.list(x)
  if (nchain(x) < 2)
    stop("You need at least two chains")
  if (autoburnin && start(x) < end(x)/2)
    x <- window(x, start = end(x)/2 + 1)
  Niter <- niter(x)
  Nchain <- nchain(x)
  Nvar <- nvar(x)
  xnames <- varnames(x)
  if (transform)
    x <- gelman.transform(x)
  x <- lapply(x, as.matrix)
  S2 <- array(sapply(x, var, simplify = TRUE),
              dim = c(Nvar, Nvar, Nchain)
  )
  W <- apply(S2, c(1, 2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE),
                 nrow = Nvar, ncol = Nchain)
  B <- Niter * var(t(xbar))
  if (Nvar > 1 && multivariate) {  #ph-edits
    # CW <- chol(W)
    #    #This is W^-1*B.
    # emax <- eigen(
    #  backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), transpose = TRUE),
    # symmetric = TRUE, only.values = TRUE)$values[1]
    emax <- 1
    mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
  }  else {
    mpsrf <- NULL
  }

  w <- diag(W)
  b <- diag(B)
  s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, var)/Nchain
  var.b <- (2 * b^2)/(Nchain - 1)
  cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 *
                                    muhat * var(t(s2), t(xbar)))
  V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b +
              2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) *
    R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  out <- list(psrf = psrf, mpsrf = mpsrf, B = B, W = W) #added ph
  class(out) <- "gelman.diag"
  return( out )
}
