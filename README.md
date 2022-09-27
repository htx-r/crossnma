# crossnma: Cross-Design and Cross-Format Synthesis using Network Meta-Analysis and Network Meta-Regression
Official Git repository of R package **crossnma**

[![CRAN Version](http://www.r-pkg.org/badges/version/crossnma)](https://cran.r-project.org/package=crossnma)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/crossnma)](http://cranlogs.r-pkg.org/badges/crossnma)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/crossnma)](http://cranlogs.r-pkg.org/badges/grand-total/crossnma)


## Description

**crossnma** is an R package that allows for synthesizing the data
from randomized or non-randomized studies coming from
individual-participants or aggregate data. The package implements
Bayesian models in network meta-analysis and network meta-regression
through JAGS software.


### References

[Hamza T, Chalkou K, Pellegrini F, et al. (2022), *Synthesizing cross-design evidence and cross-format data using network meta-regression*, arXiv:2203.06350](https://www.doi.org/10.48550/arXiv.2203.06350)


## Installation

### Current stable [![CRAN Version](http://www.r-pkg.org/badges/version/crossnma)](https://cran.r-project.org/package=crossnma) release:
```r
install.packages("crossnma")
```

### Current beta / GitHub release:

Installation using R package
[**remotes**](https://cran.r-project.org/package=remotes):
```r
install.packages("remotes")
remotes::install_github("htx-r/crossnma")
```


## How to use crossnma?

There are two steps to conduct a network meta-analysis or
meta-regression. The first step is to create a JAGS model using
`crossnma.model()` which produces the JAGS code and transforms the
data to the JAGS format. In the second step, the output of that
function will be used in `crossnma()` to run the MCMC (Markov chain
Monte Carlo) through JAGS.

We illustrate how to use **crossnma** through several examples in the vignette:

```
vignette("crossnma", package = "crossnma")
```


## How to cite crossnma?

```
citation(package = "crossnma")
```


### Bug Reports:

```r
bug.report(package = "crossnma")
```

The bug.report function is not supported in RStudio. Please send an
email to Tasnim Hamza <tasnim.hamza@ispm.unibe.ch> if you use RStudio.

You can also report bugs on GitHub under [Issues](https://github.com/htx-r/crossnma/issues/).
