##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s)
##
document("../crossnma")


##
## (3) Build R package and PDF file with help pages
##
build("../crossnma")
build_manual("../crossnma")


##
## (4) Install R package
##
install("../crossnma")
