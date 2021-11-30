#' A simulated individual participant data.
#'
#' A dataset containing 2950 participants who are treated in 4 different treatments: A, B, C and D.
#' The dataset includes 4 different studies. The outcome is binary. There are 10 attributes on individual level.
#'
#' @format A data frame with 2950 rows and 8 variables:
#' \describe{
#'   \item{study}{numeric, study identifier}
#'   \item{outcome}{binary data, respond indicator, 0=respond and 1=no-respond}
#'   \item{trt}{character, indicating the assigned treatment to each participant}
#'   \item{design}{character, design of the study, either 'rct' or 'nrs'}
#'   \item{age}{numeric, age of the participant}
#'   \item{sex}{binary data, sex of the participant, 0=Female and 1=Male}
#'   \item{bias}{character, the risk of bias of the study, 'low', 'high','unclear'}
#'   \item{unfav}{numeric, the indicator of the unfavoured treatment of the study, 0, 1}
#'   \item{bias.group}{numeric, the bias effect of the study, 1 = study has inactive treatment and adjust for bias effect, 2= study has active treatments and another adjustment for bias effect, 0=no bias adjustment}
#'   \item{year}{numeric, the year published of the study}
#'
#' }
#'
"prt.data"


#' A simulated aggregate data.
#' The dataset includes 2 randomised-controlled trials (RCTs), comparing 2 treatments: A, C.
#' The outcome is binary that is represented as the number of participants with a respond.
#'
#'#' @format A data frame with 4 rows and 11 variables:
#' \describe{
#'   \item{study}{numeric, study identifier}
#'   \item{outcome}{numeric, the number of responders}
#'   \item{n}{numeric, the sample size}
#'   \item{trt}{character, indicating the assigned treatment to participants in each study arm}
#'   \item{design}{character, design of the study, either 'rct' or 'nrs'}
#'   \item{age}{numeric, the mean age of participants in each study}
#'   \item{sex}{numeric, the proportion of males on each study}
#'   \item{bias}{character, the risk of bias of the study, 'low', 'high','unclear'}
#'   \item{unfav}{numeric, the indicator of the unfavoured treatment of the study, 0, 1}
#'   \item{bias.group}{numeric, the bias effect of the study, 1 = study has inactive treatment and adjust for bias effect, 2= study has active treatments and another adjustment for bias effect, 0=no bias adjustment}
#'   \item{year}{numeric, the year published of the study}
#'
#' }
"std.data"
