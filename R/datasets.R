#' Simulated individual participant dataset.
#'
#' @description
#' A dataset containing 1944 participants who are treated in four
#' different treatments: A, B, C and D. The dataset includes four
#' studies. The outcome is binary. There are 10 attributes on
#' individual level.
#'
#' @format A data frame with 1944 rows and 10 variables:
#' \describe{
#'   \item{id}{numeric, study identifier}
#'   \item{relapse}{binary data, respond indicator, 0=no relapse and
#'     1=relapse}
#'   \item{treat}{character, indicating the assigned treatment to each
#'     participant}
#'   \item{design}{character, design of the study, either 'rct' or
#'     'nrs'}
#'   \item{age}{numeric, age of the participant}
#'   \item{sex}{binary data, sex of the participant, 0=female and
#'     1=male}
#'   \item{rob}{character, the risk of bias of the study, 'low',
#'     'high','unclear'}
#'   \item{unfavored}{numeric, the indicator of the unfavored
#'     treatment in each study, values are 0 or 1}
#'   \item{bias.group}{numeric, the bias effect of the study, 1 = if
#'     the study has inactive treatment and adjust for bias effect, 2=
#'     if the study has active treatments and it is assumed another
#'     bias effect, 0=no bias adjustment}
#'   \item{year}{numeric, the year study published}
#' }


"ipddata"





#' Simulated aggregate dataset.
#'
#' @description
#' The dataset includes two randomized-controlled trials (RCTs),
#' comparing treatments A and C. The outcome is binary represented as
#' the number of participants with at least one relapse.
#'
#' @format A data frame with 4 rows and 11 variables:
#' \describe{
#'   \item{id}{numeric, study identifier}
#'   \item{n}{numeric, the sample size}
#'   \item{relapse}{numeric, the number of relapses}
#'   \item{treat}{character, indicating the assigned treatment to
#'     participants in each study arm}
#'   \item{design}{character, design of the study, either 'rct' or
#'     'nrs'}
#'   \item{age}{numeric, the mean age of participants in each study}
#'   \item{sex}{numeric, the proportion of females on each study}
#'   \item{rob}{character, the risk of bias of the study, 'low',
#'     'high','unclear'}
#'   \item{unfavored}{numeric, the indicator of the unfavored
#'     treatment in each study, values are 0 or 1}
#'   \item{bias.group}{numeric, the bias effect of the study, 1 =
#'     study has inactive treatment and adjust for bias effect, 2=
#'     study has active treatments and another adjustment for bias
#'     effect, 0=no bias adjustment}
#'   \item{year}{numeric, the year published of the study}
#' }


"stddata"
