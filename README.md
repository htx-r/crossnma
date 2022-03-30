# crossnma
crossnma (CROSS-design and CROSS-format synthesis using Network Meta-Analysis) is an R Package that allows for synthesizing the data that can be in various formats; individual-participants (IPD) or aggregate (AD) and can come from non-randomized or randomized studies (NRS, RCT). The package implements Bayesian models in network meta-analysis and network meta-regression through JAGS software.
# Installation
The package can be installed in Rstudio directly from GitHub using the R package `remotes` 
```
# Install first remotes
install.packages("remotes")

# Then install crossnma
remotes::install_github("htx-r/crossnma")
```

# How to use crossnma?
There are two steps to run the NMA/NMR model. The first step is to create a JAGS model using `crossnma.model()` which produces the JAGS code and the data. In the second step, the output of that function will be used in `crossnma.run()` to run the MCMC (Markov chain Monte Carlo) through JAGS.

We illustrate how to use `crossnma` through several examples in the vignette:
```
vignette("crossnma", package = "crossnma")
```


# How to cite crossnma?

```
citation(package = "crossnma")
```
