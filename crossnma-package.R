pdf(file = "RpLoTs-crossnma.pdf", width = 12, height = 9)


sortstr <- function(x){
  if (is.list(x))
    res <- str(x[sort(names(x))], list.len = 999)
  else
    res <- x
  res
}


sortunclass <- function(x){
  if (is.list(x))
    res <- unclass(x[sort(names(x))])
  else
    res <- x
  res
}

library(crossnma)
packageDescription("meta")
packageDescription("netmeta")
packageDescription("crossnma")


set.seed(19091909)


## crossnma.Rd

head(ipddata) # participant-level data
stddata # study-level data

mod <- crossnma.model(treat, id, relapse, n, design,level.ma=0.95,
  prt.data = ipddata, std.data = stddata,
  reference = "A", trt.effect = "random", method.bias = "naive")

set.seed(1909)
fit <- crossnma(mod, n.burnin = 10, n.iter = 50,
  n.thin = 1, n.chains = 3)

summary(fit)
plot(fit)


## crossnma.model.Rd

mod

summary(mod)





## heatplot.crossnma.Rd

heatplot(fit)





## league.crossnma.Rd

# Create league tables
league1 <- league(fit)                     #  wide format
league2 <- league(fit, direction = "long") #  long format
league1
league2





## netconnection.crossnma.Rd

nc1 <- netconnection(fit)
nc1





## netconnection.crossnma.Rd

nc2 <- netconnection(mod)
nc2





## netgraph.crossnma.Rd

ng1 <- netgraph(fit, plastic = FALSE, cex.points = 7, adj = 0.5)
ng1





## netgraph.crossnma.Rd

ng2 <- netgraph(mod, plastic = FALSE, cex.points = 7, adj = 0.5)
ng2





## plot.crossnma.Rd

plot(fit)





## print.crossnma.Rd
## print.crossnma.model.Rd
## print.summary.crossnma.Rd
## print.summary.crossnma.model.Rd





## summary.crossnma.Rd

sfit <- summary(fit)
sfit





##
##
## Vignette
##
##

set.seed(1910)
settings.meta(digits = 3)
cilayout("(", " to ")

mod1 <- crossnma.model(treat, id, relapse, n, design,
  prt.data = ipddata, std.data = stddata,
  trt.effect = "random",
  #---------- bias adjustment ----------
  method.bias = "naive",
  #---------- assign a prior ----------
  prior.tau.trt = "dunif(0, 3)"
  )
summary(mod1)
##
jagsfit1 <- crossnma(mod1, n.iter = 5000, n.burnin = 2000, n.thin = 1)
jagsfit1
summary(jagsfit1, backtransf = FALSE)
par(mar = rep(2, 4), mfrow = c(2, 3))
plot(jagsfit1)

mod2 <- crossnma.model(treat, id, relapse, n, design,
  prt.data = ipddata, std.data = stddata,
  trt.effect = "random",
  #---------- bias adjustment ----------
  method.bias = "naive",
  #----------  meta-regression ----------
  cov1 = age,
  split.regcoef = FALSE
  )
summary(mod2)
##
jagsfit2 <- crossnma(mod2, n.iter = 5000, n.burnin = 2000, n.thin = 1)
summary(jagsfit2, backtransf = FALSE)
league(jagsfit2, cov1.value = 38, digits = 2)
league(jagsfit2, cov1.value = 38, digits = 2, direction = "long")


mod3 <- crossnma.model(treat, id, relapse, n, design,
  prt.data = ipddata, std.data = stddata,
  reference = "D", trt.effect = "random",
  #----------  meta-regression ----------
  cov1 = age,
  split.regcoef = FALSE,
  #---------- bias adjustment ----------
  method.bias = "prior",
  run.nrs.trt.effect= "common",
  run.nrs.var.infl = 0.6, run.nrs.mean.shift = 0,
  run.nrs.n.iter = 10000, run.nrs.n.burnin = 4000,
  run.nrs.n.thin = 1, run.nrs.n.chains = 2
  )
summary(mod3)
##
jagsfit3 <- crossnma(mod3, n.iter = 5000, n.burnin = 2000, n.thin = 1)
heatplot(jagsfit3, cov1.value = 38,
  size = 6, size.trt = 20, size.axis = 12)


mod4 <- crossnma.model(treat, id, relapse, n, design,
  prt.data = ipddata, std.data = stddata,
  trt.effect = "random",
  #---------- bias adjustment ----------
  method.bias = "adjust1",
  bias.type = "add",
  bias.effect = "common",
  bias = rob,
  unfav = unfavored,
  bias.group = bias.group,
  bias.covariate = year
)
summary(mod4)
##
jagsfit4 <- crossnma(mod4, n.iter = 5000, n.burnin = 2000, n.thin = 1)
summary(jagsfit4, backtransf = FALSE)


mod5 <- crossnma.model(treat, id, relapse, n, design,
  prt.data = ipddata, std.data = stddata,
  trt.effect = "random",
  #---------- bias adjustment ----------
  method.bias = "adjust2",
  bias.type = "add",
  bias = rob,
  unfav = unfavored,
  bias.group = bias.group
)
summary(mod5)
##
jagsfit5 <- crossnma(mod5, n.iter = 5000, n.burnin = 2000, n.thin = 1)
summary(jagsfit5, backtransf = FALSE)


mod6 <- crossnma.model(treat, id, relapse, n, design,
  prt.data = ipddata, std.data = stddata,
  trt.effect = "random",level.ma=0.95,
  #---------- bias adjustment ----------
  method.bias = "naive",
  #----------  meta-regression ----------
  cov1 = age, cov2 = year,
  split.regcoef = FALSE
  )
summary(mod6)
##
jagsfit6 <- crossnma(mod6, n.iter = 5000, n.burnin = 2000, n.thin = 1)
summary(jagsfit6, backtransf = FALSE)





##
##
## Save R objects
##
##

save.image("crossnma-package.rda")





##
##
## Structure of R objects
##
##

sortstr(mod)
sortstr(fit)
sortstr(sfit)
sortstr(league1)
sortstr(league2)
sortstr(ng1)
sortstr(ng2)
##
sortstr(mod1)
sortstr(jagsfit1)
sortstr(mod2)
sortstr(jagsfit2)
sortstr(mod3)
sortstr(jagsfit3)
sortstr(mod4)
sortstr(jagsfit4)
sortstr(mod5)
sortstr(jagsfit5)
sortstr(mod6)
sortstr(jagsfit6)





##
##
## Unclass R objects
##
##

sortunclass(mod)
sortunclass(fit)
sortunclass(league1)
sortunclass(league2)
sortunclass(nc1)
sortunclass(nc2)
sortunclass(ng1)
sortunclass(ng2)
sortunclass(sfit)
##
sortunclass(mod1)
sortunclass(jagsfit1)
sortunclass(mod2)
sortunclass(jagsfit2)
sortunclass(mod3)
sortunclass(jagsfit3)
sortunclass(mod4)
sortunclass(jagsfit4)
sortunclass(mod5)
sortunclass(jagsfit5)
sortunclass(mod6)
sortunclass(jagsfit6)
