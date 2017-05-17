### wally-package.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 28 2017 (11:13) 
## Version: 
## Last-Updated: May 16 2017 (08:11) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' 
#' divat data
#'
#' Extracted data from a french population based cohort (DIVAT cohort). The dataset includes 
#' followup information on kidney graft failure  outcome and predicted 5-year risks based on 
#' based on the subject specific information which includes age, gender,
#' cardiovascular and diabetes histories, monitoring of the evolution of the kidney function
#' measured via serum creatinine and relevant characteristics of his or her kidney donor.
#' Graft failure is defined as either death with functioning kidney graft or return to dialysis.
#' The prediction model from which the predictions have been computed has been previously fitted
#' using an independent training sample from the DIVAT data. Details about data and modeling can
#' be found in Fournier et al. (2016).
#' 
#' @name divat
#' @docType data
#' @format A subsample consisting of 1300 observations on the following 3 variables.
#' \describe{ \item{pi}{5-year risk prediction of kidney graft failure.}
#' \item{status}{0=censored, 1=kidney graft failure}
#' \item{time}{time to event (i.e., time to kidney graft failure or loss of follow-up)}}
#' @references
#' Fournier, M. C., Foucher, Y., Blanche, P., Buron, F., Giral, M., & Dantan, E. (2016).
#' A joint model for longitudinal and time-to-event data to better assess the specific
#' role of donor and recipient factors on long-term kidney transplantation outcomes.
#' European journal of epidemiology, 31(5), 469-479.
#' 
#' @keywords datasets
#' @examples
#' 
#' data(divat)
NULL

#' threecity data
#'
#' Extracted data from a french population based cohort (Three-City cohort). The dataset includes 
#' followup information on dementia outcome and predicted 5-year risks based on 
#' based on the subject specific information which includes age, gender,
#' education level and cognitive decline measured by a psychometric test
#' (Mini Mental State Examination). The prediction model from which the
#' predictions have been computed has been fitted on independent training
#' data from the Paquid cohort, another french population based cohort with similar design (see Reference Blanche et al. 2015 for details) .
#' 
#' @name threecity
#' @docType data
#' @format A subsample consisting of 2000 observations on the following 3 variables.
#' \describe{ \item{pi}{5-year absolute risk predictions of dementia.}
#' \item{status}{0=censored, 1=dementia, 2=death dementia free}
#' \item{time}{time to event (i.e., time to
#' either dementia, death dementia free or loss of follow-up)}}
#' @references
#' Blanche, P., Proust-Lima, C., Loubere, L., Berr, C., Dartigues, J. F., Jacqmin-Gadda, H. (2015).
#' Quantifying and comparing dynamic
#' predictive accuracy of joint models for longitudinal marker and
#' time-to-event in presence of censoring and competing risks. 
#' Biometrics, 71(1), 102-113.
#' 
#' @source
#'
#' Web-appendix of Blanche et al. (2015).
#' 
#' @keywords datasets
#' @examples
#' 
#' data(threecity)
#'
#' 
#' @importFrom prodlim Hist jackknife prodlim sindex
#' @importFrom stats rmultinom
#' @importFrom grDevices col2rgb gray
#' @importFrom graphics bxp  abline axis box legend lines mtext par plot points segments text title polygon par boxplot
#' @importFrom utils capture.output find head select.list tail
#' @importFrom stats as.formula coef delete.response drop.terms family formula get_all_vars glm median model.frame model.matrix model.response na.fail na.omit optim pnorm predict qnorm quantile rbinom reformulate rexp runif sd setNames smooth terms terms.formula time uniroot update update.formula var wilcox.test
NULL




######################################################################
### wally-package.R ends here
