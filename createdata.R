library('dplyr')
library(magrittr)
library(plyr)
# manipulate IPD data
prt.data <- read.csv("IPD for analysis")
#** add unfav and bias.group columns
# IPD
# myprt.data$unfav <- myprt.data$bias.group <- NA
# study_plac <- unique(myprt.data$STUDYID[myprt.data$TRT01A=="Placebo"]) # index of studies with placebo
# for(i in 1:length(study_plac)){
#   myprt.data[myprt.data$STUDYID==study_plac[i]&myprt.data$TRT01A=="Placebo",]["unfav"] <- 0
#   myprt.data[myprt.data$STUDYID==study_plac[i]&myprt.data$TRT01A!="Placebo",]["unfav"]  <- 1
#   myprt.data[myprt.data$STUDYID==study_plac[i],]["bias.group"]<- 1
# }
#
# # NRS
# myprt.data[myprt.data$STUDYID=="NRS"&myprt.data$TRT01A=="Glatiramer acetate",]["unfav"] <- 0
# myprt.data[myprt.data$STUDYID=="NRS"&myprt.data$TRT01A!="Glatiramer acetate",]["unfav"]  <- 1
# myprt.data[myprt.data$STUDYID=="NRS",]["bias.group"]<- 1
# prt.data <- myprt.data
newdata <- list()
i=1
for (i in 1:4) {
prt.data.t <- prt.data[prt.data$STUDYID==levels(prt.data$STUDYID)[i],]
prt.data.t$trt <- as.factor(as.character(prt.data.t$TRT01A))
prt.data.t$SEX <- sample(0:1, nrow(prt.data.t), replace=T,prob=c(0.60,0.40))
#---- 1. Run a simple logistic regression + 2. Estimate the parameters
m <- glm(RELAPSE2year~trt+AGE+SEX,
         data=prt.data.t,
         family = binomial(link = "logit"))

#---- 3. change the age
agenew <- round(prt.data.t$AGE+rnorm(1,mean=0,sd=1.5))

# #---- 4. change the EDSS, round and bw 0 and 10
# edssnew <- prt.data.t$EDSSBL+rnorm(1,mean=0,sd=0.5)

#---- 4. sex
sex <- prt.data.t$SEX
#---- 5. find the probability to relapse
p <- predict.glm(m, data.frame(trt=prt.data.t$trt,
                          AGE=agenew,
                          SEX=sex),type = "response")

relapse <- rbinom(nrow(prt.data.t),1,p)
newdata0 <- data.frame(study=i,
                           outcome=relapse,
                           trt=prt.data.t$trt,
                           design=ifelse(prt.data.t$design=='nrs','nrs','rct'),
                           age=agenew,
                           sex=sex,
                           bias=prt.data.t$bias
                           # unfav=prt.data.t$unfav,
                           # bias.group=prt.data.t$bias.group
                           )
if(prt.data.t$design=='nrs'){
  newdata[[i]] <- newdata0[1:150,]
}else{
  newdata[[i]] <- newdata0[sort(sample(1:nrow(newdata0),round(nrow(newdata0)/2))),]
}

}
newdata2[sort(sample(1:nrow(newdata2),round(nrow(newdata2)/2))),]
newdata2 <- do.call(rbind,newdata)

newdata2 %<>% mutate(trt.new=mapvalues(as.character(trt),
                                     from=c('Placebo', 'Dimethyl fumarate', 'Glatiramer acetate','Natalizumab' ),
                                     to=c('A', 'B','C', 'D'),
                                     warn_missing = FALSE))

newdata2$trt <- newdata2$trt.new
newdata2 <- newdata2[,-which(names(newdata2) %in% c('trt.new'))]


newdata2 %<>% mutate(year=mapvalues(as.character(study),
                                       from=1:4,
                                       to=c(2002, 2007,20012, 2015),
                                       warn_missing = FALSE)%>%as.numeric())

names(newdata2)
prt.data <- newdata2
table(prt.data$bias.group)
prt.data$unfav <- prt.data$bias.group <- NA
study_plac <- unique(prt.data$study[prt.data$trt=="A"]) # index of studies with placebo

for(i in 1:length(study_plac)){
  prt.data[prt.data$study==study_plac[i]&prt.data$trt=="A",]["unfav"] <- 0
  prt.data[prt.data$study==study_plac[i]&prt.data$trt!="A",]["unfav"]  <- 1
  prt.data[prt.data$study==study_plac[i],]["bias.group"]<- 1
}

# NRS
prt.data[prt.data$design=="nrs"&prt.data$trt=="B",]["unfav"] <- 0
prt.data[prt.data$design=="nrs"&prt.data$trt!="B",]["unfav"]  <- 1
prt.data[prt.data$design=="nrs",]["bias.group"]<- 1

# prt.data <- newdata2[,c("STUDYID","RELAPSE2year","TRT01A","design",  "AGE","SEX","bias","unfav","bias.group")]
# names(prt.data) <- c("study","r","treat","design","age","bias","unfav","bias.group")

# manipulate aggregate data
std.data0 <- read.csv("AD for analysis")

# AD
std.data0$unfav <- std.data0$bias.group <- NA

for(i in 1:length(unique(std.data0$study))){
  std.data0[std.data0$study==unique(std.data0$study)[i]&std.data0$treat=="Placebo",]["unfav"] <- 0
  std.data0[std.data0$study==unique(std.data0$study)[i]&std.data0$treat!="Placebo",]["unfav"]  <- 1
  std.data0[std.data0$study==unique(std.data0$study)[i],]["bias.group"]<- 1
}
std.data <- std.data0[,c('study','r','n','treat','design','age','sex','bias','unfav','bias.group','year')]
names(std.data) <- c('study','outcome','n','trt','design','age','sex','bias','unfav','bias.group','year')
std.data %<>% mutate(trt.new=mapvalues(as.character(trt),
                                       from=c('Placebo', 'Dimethyl fumarate', 'Glatiramer acetate','Natalizumab' ),
                                       to=c('A', 'B','C', 'D'),
                                       warn_missing = FALSE))

std.data$trt <- std.data$trt.new
std.data <- std.data[,-which(names(std.data) %in% c('trt.new'))]


std.data$year <-  mapvalues(as.character(std.data$study),
                                    from=unique(std.data$study),
                                    to=c(2010,2015),
                                    warn_missing = FALSE)%>%as.numeric()
std.data$study <- rep(1:length(unique(std.data$study)),each=2)


usethis::use_data(prt.data,std.data)
std.data
# create all Vignette related files
#usethis::use_vignette("gnma")
#devtools::document()


