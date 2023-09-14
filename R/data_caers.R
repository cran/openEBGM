#' Dietary supplement reports and products
#'
#' A dataset for dietary supplement adverse event reports from 2012 containing
#' CAERS product and adverse event reports as reported to the FDA. This
#' particular dataset contains only products which were reported to be dietary
#' supplements (industry code 54) reported in the year 2012. and includes 2874
#' unique product names and 1328 unique adverse events. There are a total of
#' 3356 unique reports. In addition, there is also one stratification variable,
#' indicating whether the patient is male or female
#'
#' @format A data frame with 20156 rows and 4 variables:
#' \describe{
#'   \item{\code{id}}{Identification number}
#'   \item{\code{var1}}{Name of the product}
#'   \item{\code{var2}}{Name of the symptom/event category}
#'   \item{\code{strat1}}{Gender of the patient associated with report}
#' }
#'
#'
#' @details Further details about the data can be found using the links below.
#'
#' @source CFSAN Adverse Event Reporting System (FDA Center for Food Safety and
#'   Nutrition)
#' @source \url{https://www.fda.gov/food/compliance-enforcement-food}
#' @source
#'   \url{https://www.fda.gov/media/97035/download}
"caers"
