#' A simulated individual participant data.
#'
#' A dataset containing 2950 participants who are treated in 4 different treatments: A, B, C and D.
#' The dataset includes 4 different studies. The outcome is binary. There are 10 attributes on individual level.
#'
#' @format A data frame with 2950 rows and 10 variables:
#' \describe{
#'   \item{study}{numeric, study identifier}
#'   \item{outcome}{binary data, respond indicator, 0=no relapse and 1=relapse}
#'   \item{trt}{character, indicating the assigned treatment to each participant}
#'   \item{design}{character, design of the study, either 'rct' or 'nrs'}
#'   \item{age}{numeric, age of the participant}
#'   \item{sex}{binary data, sex of the participant, 0=Female and 1=Male}
#'   \item{bias}{character, the risk of bias of the study, 'low', 'high','unclear'}
#'   \item{unfav}{numeric, the indicator of the unfavored treatment in each study, values are 0 or 1}
#'   \item{bias.group}{numeric, the bias effect of the study, 1 = if the study has inactive treatment and adjust for bias effect, 2= if the study has active treatments and it is assumed another bias effect, 0=no bias adjustment}
#'   \item{year}{numeric, the year study published}
#'
#' }
#'
"prt.data"


#' A simulated aggregate data.
#'
#' @description The dataset includes 2 randomized-controlled trials (RCTs), comparing 2 treatments: A and C.
#' The outcome is binary represented as the number of participants with at least one relapse.
#'
#' @format A data frame with 4 rows and 11 variables:
#' \describe{
#'   \item{study}{numeric, study identifier}
#'   \item{outcome}{numeric, the number of responders}
#'   \item{n}{numeric, the sample size}
#'   \item{trt}{character, indicating the assigned treatment to participants in each study arm}
#'   \item{design}{character, design of the study, either 'rct' or 'nrs'}
#'   \item{age}{numeric, the mean age of participants in each study}
#'   \item{sex}{numeric, the proportion of females on each study}
#'   \item{bias}{character, the risk of bias of the study, 'low', 'high','unclear'}
#'   \item{unfav}{numeric, the indicator of the unfavored treatment in each study, values are 0 or 1}
#'   \item{bias.group}{numeric, the bias effect of the study, 1 = study has inactive treatment and adjust for bias effect, 2= study has active treatments and another adjustment for bias effect, 0=no bias adjustment}
#'   \item{year}{numeric, the year published of the study}
#' }
"std.data"
