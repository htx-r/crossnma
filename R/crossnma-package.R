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
#' et al. (2023).
#'
#' @details
#' The evidence in network meta-analysis (NMA) typically comes from
#' randomized controlled trials (RCT) where aggregate data (AD) are
#' extracted from published reports. Retrieving individual participant
#' data (IPD) allows considering participant covariates to explain
#' some of the heterogeneity/inconsistency in the network and identify
#' effect modifiers. Additionally, evidence from non-randomized
#' studies (NRS) reflects the reality in clinical practice and bridges
#' the efficacy-effectiveness gap. The cross-NMA/NMR model is a
#' Bayesian suite for evidence synthesis which extends and integrates
#' four different approaches that combine RCT and NRS evidence into a
#' three-level hierarchical model for the synthesis of IPD and AD. The
#' four approaches account for differences in the design and risk of
#' bias in the RCT and NRS evidence. These four approaches variously
#' ignoring differences in risk of bias, using NRS to construct
#' penalized treatment effect priors and bias-adjustment models that
#' control the contribution of information from high risk of bias
#' studies in two different ways.
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
#' @author Tasnim Hamza \email{tasnim.hamza@@ispm.unibe.ch},
#' Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de},
#' Georgia Salanti \email{georgia.salanti@@ispm.unibe.ch}
#'
#' @references
#' Dias S, Welton NJ, Marinho VCC et al. (2010):
#' Estimation and adjustment of bias in randomized evidence by using
#' mixed treatment comparison meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{173}, 613-29
#'
#' Hamza T, Chalkou K, Pellegrini F et al. (2023):
#' Synthesizing cross-design evidence and cross-format data using
#' network meta-regression.
#' \emph{Research Synthesis Methods},
#' \bold{14}, 283-300
#' 
#' Plummer M (2003):
#' JAGS: A program for analysis of Bayesian graphical models using
#' Gibbs sampling
#'
#' Saramago P, Sutton AJ, Cooper NJ, Manca A (2012):
#' Mixed treatment comparisons using aggregate and individual
#' participant level data.
#' \emph{Statistics in Medicine},
#' \bold{31}, 3516-36
#'
#' Tramacere I, Del Giovane C, Salanti G et al. (2015):
#' Immunomodulators and immunosuppressants for relapsing-remitting
#' multiple sclerosis: a network meta-analysis.
#' \emph{Cochrane Database of Systematic Reviews}, \bold{9},
#' John Wiley & Sons, Ltd. \doi{10.1002/14651858.CD011381.pub2}
#'
#' Verde PE (2021):
#' A bias-corrected meta-analysis model for combining studies of
#' different types and quality.
#' \emph{Biometrical Journal},
#' \bold{63}, 406-22
#'
#' @keywords package
#'
#' @importFrom graphics plot
#' @importFrom ggplot2 aes element_blank element_text geom_text
#'   geom_tile ggplot scale_fill_gradient2 scale_x_discrete
#'   scale_y_discrete theme theme_dark
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr .data arrange bind_cols bind_rows case_when do
#'   everything filter group_by group_keys group_map group_size
#'   left_join mutate n pull rename row_number select starts_with
#'   summarize summarise_all ungroup
#' @importFrom plyr mapvalues
#' @importFrom tidyr as_tibble gather spread unnest
#' @importFrom rlang quo
#' @importFrom stats start end qf sd var update window
#' @importFrom utils packageDescription combn
#' @importFrom meta gs
#' @importFrom netmeta pairwise netmeta netgraph heatplot
#'   netconnection
#' @importFrom rjags jags.model coda.samples
#' @importFrom coda as.mcmc.list effectiveSize nchain niter nvar
#'   traceplot varnames

"_PACKAGE"

NULL
