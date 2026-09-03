#' Dietary supplement reports
#'
#' A small subset of data from the Human Foods Complaint System (HFCS) database
#' (formerly CAERS). This dataset contains only reports received in 2025 for
#' dietary supplement products (industry/product code 54). 
#'
#' @format A data frame with 15,204 rows and 4 columns:
#' \describe{
#'   \item{\code{id}}{The unique number that identifies each case.}
#'   \item{\code{strat1}}{The sex of the individual reported to have experienced
#'     the adverse event.}
#'   \item{\code{var1}}{The verbatim brands and/or product names indicated to
#'     have been used by the consumer reported to have experienced the adverse
#'     event. An adverse event report may specify consumption of a single
#'     product or multiple products.}
#'   \item{\code{var2}}{The symptom(s) experienced by the injured consumer as
#'     specified by the reporter and coded by FDA according to the Medical Data
#'     Dictionary for Regulatory Activities (MedDRA).
#'     Each adverse event report may indicate one or more symptoms.}
#' }
#'
#' @details Further details about the data can be found using the links below.
#'
#' @source Human Foods Complaint System (HFCS)
#' @source \url{https://www.fda.gov/food/compliance-enforcement-food/human-foods-complaint-system-hfcs}
#' @source
#'   \url{https://www.fda.gov/media/97035/download}
"caers"
