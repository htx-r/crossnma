#'@export
netplot <- function(prt.data=prt.data,
                    std.data=std.data){
prt.data.ad0 <- sapply(1:length(unique(prt.data$study)), function(i){
  with(prt.data,
       data.frame(
         study=length(unique(std.data$study))+unique(prt.data$study)[i],
         trt=unique(prt.data[prt.data$study==unique(prt.data$study)[i],]$trt),
         outcome=sum(prt.data[prt.data$study==unique(prt.data$study)[i],]$outcome),
         n =nrow(prt.data[prt.data$study==unique(prt.data$study)[i],]),
         design=unique(prt.data[prt.data$study==unique(prt.data$study)[i],]$design)
       ))
}, simplify = F)
prt.data.ad <- do.call(rbind,prt.data.ad0)
all.data <- rbind(std.data[c('study','trt','outcome','n','design')],prt.data.ad)
data.pair = suppressWarnings(pairwise(treat = trt, event = outcome, n = n,
                          data = all.data, studlab = study, sm = "OR",warn = FALSE))
net1 = suppressWarnings(netmeta(data.pair,warn = FALSE))

## ------------------------------------------------------------------------
return(netgraph(net1,offset = 0.1)) # ,col=1:net1$m col=as.numeric(all.data$design)

}



