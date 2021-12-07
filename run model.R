devtools::install_github("htx-r/crossnma",force = TRUE)
library(crossnma)

#-------- MCMC settings --------#
n.adapt = 20
n.iter=100
n.burnin = 40
thin=1
n.chains=2

mod3 <- crossnma.model(prt.data=NULL,
                       std.data=std.data,
                       trt='trt',
                       study='study',
                       outcome='outcome',
                       n='n',
                       design='design',
                       reference='A',
                       trt.effect='common',
                       #---------- bias adjustment ----------
                       method.bias='adjust2',
                       bias='bias',
                       bias.type='add',
                       bias.effect='common',
                       unfav="unfav",
                       bias.group="bias.group"
)


# run jags
jagsfit3 <- crossnma.run(model=mod3,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor=c('LOR'))

# output
summary(jagsfit3,expo=F)
coda::traceplot(jagsfit3$samples)



mod5 <- crossnma.model(prt.data=prt.data,
                   std.data=std.data,
                   trt='trt',
                   study='study',
                   outcome='outcome',
                   n='n',
                   design='design',
                   reference='A',
                   trt.effect='random',
                   #---------- bias adjustment ----------
                   method.bias='adjust2',
                   bias='bias',
                   bias.type='add',
                   unfav='unfav',
                   bias.group="bias.group",
                   bias.effect='common'
                   )
# run jags
jagsfit5 <- crossnma.run(model=mod5,
                       n.adapt = 50,
                        n.iter=500,
                        n.burnin = 200,
                        thin=1,
                        n.chains=2)
