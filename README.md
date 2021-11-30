# crossnma
crossnma (CROSS-design and CROSS-Format synthesis using Network Meta-Analysis) is an R Package that allows for synthesize the data that can be in various formats; individual-participants (IPD) or aggregated (AD) and can come from non-randomized or randomized studies (NRS, RCT). The package implements Bayesian models in network meta-analysis and meta-regression through JAGS software.
# Installation
The package can be installed in Rstudio directly from GitHub using the R package `devtools` 
```
# Install first devtool
install.packages("devtools")

# Then install crosnma
devtools::install_github("TasnimHamza/crosnma")
```
For the bias-adjustment 2 method, we need further to install the R package `rjags` to load the mixture of normals module `load.module("mix")`.
```
# Install the rjags
install.packages("rjags")
library("rjags")
# To be able to load the mix module 
load.module("mix")
```
# How to use crosnma?
There are two steps to run the NMA/NMR model. The first step is to create a JAGS model using `crosnma.model()` which produces the JAGS code and the data. In the second step, the output of that function will be used in `crosnma.run()` to run the MCMC (Markov chain Monte Carlo) through JAGS.

In the vignettes [here](), we illustrate how to use `crosnma` through several examples.
# How to cite crosnma?

