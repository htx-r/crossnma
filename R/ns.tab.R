ns.tab <- function(prt.data=prt.data,
                   std.data=std.data){
  m <- matrix(c(ns.ipd.rct=length(unique(prt.data[prt.data$design=='rct',]$study)),
                ns.ipd.nrs=length(unique(prt.data[prt.data$design=='nrs',]$study)),
                ns.ad.rct=length(unique(std.data[std.data$design=='rct',]$study)),
                ns.ad.nrs=length(unique(std.data[std.data$design=='nrs',]$study)))
              , 2,2)
  row.names(m) <- c('RCT','NRS')
  colnames(m) <- c('IPD','AD')
  return(m)
}
