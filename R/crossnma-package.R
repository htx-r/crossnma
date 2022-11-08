#' crossnma: An R package for synthesizing cross-design evidence and
#' cross-format data using Bayesian methods in network meta-analysis
#' and network meta-regression
#'
#' @description
#' An R package \bold{crossnma} for performing (network) meta-analysis
#' and (network) meta-regression (allows including up to 3 covariates)
#' of individual participant data and aggregate data or combination of
#' both. Each format can come from randomized controlled trials or
#' non-randomized studies. Estimates are generated in a Bayesian
#' framework using JAGS. The implemented models are described by Hamza
#' et al. 2022 \doi{10.48550/arXiv.2203.06350}.
#'
#' @details
#' The evidence in network meta-analysis (NMA) typically comes from
#' randomized controlled trials (RCT) where aggregate data (AD) are
#' extracted from published reports.  Retrieving individual
#' participant data (IPD) allows considering participant covariates to
#' explain some of the heterogeneity/inconsistency in the network and
#' identify effect modifiers.  Additionally, evidence from
#' non-randomized studies (NRS) reflects the reality in clinical
#' practice and bridges the efficacy-effectiveness gap. The
#' cross-NMA/NMR model is a Bayesian suite for evidence synthesis
#' which extends and integrates four different approaches that combine
#' RCT and NRS evidence into a three-level hierarchical model for the
#' synthesis of IPD and AD.  The four approaches account for
#' differences in the design and risk of bias in the RCT and NRS
#' evidence.  These four approaches variously ignoring differences in
#' risk of bias, using NRS to construct penalized treatment effect
#' priors and bias-adjustment models that control the contribution of
#' information from high risk of bias studies in two different ways.
#'
#' Further details:
#' \itemize{
#' \item To have a list of all R functions available in
#'   \bold{crossnma} type \code{help(package = "crossnma")}
#' \item The R command \code{citation("crossnma")} shows how to cite
#'   \bold{crossnma} in publications.
#' \item To report problems and bugs send an email to
#'   \email{tasnim.hamza@@ispm.unibe.ch}
#' \item The development version of \bold{crossnma} is available on
#'   GitHub \url{https://github.com/htx-r/crossnma}.
#'}
#' @name crossnma-package
#'
#' @docType package
#'
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch},
#' Guido Schwarzer \email{sc@@imbi.uni-freiburg.de},
#' Georgia Salanti \email{georgia.salanti@@ispm.unibe.ch}
#'
#' @references
#' Saramago P, Sutton AJ, Cooper NJ, Manca A (2012):
#' Mixed treatment comparisons using aggregate and individual
#' participant level data.
#' \emph{Statistics in Medicine},
#' \bold{10;31(28)}, 3516-36
#'
#' Dias, Sofia, N. J. Welton, V. C. C. Marinho, G. Salanti, J.P.T
#' Higgins, and A. E. Ades (2010):
#' Estimation and Adjustment of Bias in Randomized Evidence by Using
#' Mixed Treatment Comparison Meta-Analysis.
#' \emph{Journal of the Royal Statistical Society},
#' \bold{173}, 613-29
#'
#' Plummer, Martyn. (2003):
#' JAGS: A Program for Analysis of Bayesian Graphical Models Using
#' Gibbs Sampling.
#'
#' Tramacere, Del Giovane, I, and G Filippini (2015):
#' Immunomodulators and Immunosuppressants for Relapsing‐remitting
#' Multiple Sclerosis: A Network Meta‐analysis.
#' \emph{Cochrane Database of Systematic Reviews, no. 9.} John Wiley &
#' Sons, Ltd. \doi{10.1002/14651858.CD011381.pub2}.
#'
#' Verde, Pablo Emilio. (2020):
#' A Bias-Corrected Meta-Analysis Model for Combining, Studies of
#' Different Types and Quality.
#' \emph{Biometrical Journal},
#' \doi{10.1002/bimj.201900376}
#'
#' @keywords package
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr "%<>%"
#' @importFrom plyr mapvalues
#' @importFrom rlang quo
#' @importFrom graphics plot
#' @importFrom stats end qf start update var window
#' @importFrom utils packageDescription
#' @import dplyr
#' @import ggplot2
#' @import rjags
#' @import tidyr
#' @import coda
#' @import netmeta
#' @import purrr


NULL
