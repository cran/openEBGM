#' Construct an openEBGM object
#'
#' \code{ebScores} calculates EBGM scores as well as the quantiles from the
#' posterior distribution and returns an object of class openEBGM.
#'
#' @param processed A data frame resulting from running \code{\link{processRaw}}.
#' @param hyper_estimate A list resulting from running \code{\link{autoHyper}}.
#' @param quantiles Either a numeric vector of desired quantiles to be
#'   calculated from the posterior distribution or NULL for no calculation of
#'   quantiles.
#' @param digits A whole number scalar specifying how many decimal places to
#'   use for rounding \emph{EBGM} and the quantiles scores.
#' @return An openEBGM object (class S3) containing:
#'   \itemize{
#'     \item{\emph{data}: }{A data frame containing the results (scores, etc.).}
#'     \item{\emph{hyper_parameters}: }{A list containing the hyperparameter
#'                                      estimation results.}
#'     \item{\emph{quantiles}: }{The chosen percentiles.}
#'   }
#'
#' @details This function takes the processed data as well as the hyperparameter
#'   estimates and instantiates an object of class openEBGM. This object then
#'   contains additional calculations, such as the EBGM score, and the quantiles
#'   that are supplied by the quantiles parameter at the time of calling the
#'   function.
#'
#' @details The function allows for the choice of an arbitrary amount of
#'   quantiles or no quantiles at all to be calculated. This may be helpful for
#'   large datasets, or when the EBGM score is the only metric of interest.
#'
#' @examples
#' theta_init <- data.frame(alpha1 = c(0.2, 0.1),
#'                          beta1  = c(0.1, 0.1),
#'                          alpha2 = c(2,   10),
#'                          beta2  = c(4,   10),
#'                          p      = c(1/3, 0.2)
#'                          )
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, bin_size = 100, keep_pts = 100)
#' squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 20)
#' suppressWarnings(
#'   hypers <- autoHyper(data = squashed, theta_init = theta_init)
#' )
#' obj <- ebScores(processed = proc, hyper_estimate = hypers,
#'                 quantiles = c(10, 90))
#'
#' @keywords openEBGM
#' @export

ebScores <- function(processed, hyper_estimate, quantiles = c(5, 95),
                     digits = 2) {

  .checkInputs_ebScores(processed, hyper_estimate, quantiles, digits)

  theta <- hyper_estimate$estimates
  N <- processed$N
  E <- processed$E
  q_n <- Qn(theta_hat = theta, N = N, E = E)
  EBGM <- ebgm(theta_hat = theta, N = N, E = E, qn = q_n, digits = digits)
  #Check to see if EBGM is already in the data, if so, don't duplicate column
  EBGM_sort <- EBGM[order(EBGM, decreasing = TRUE)]
  compare_func <- function(num_vec, EBGM_sort) {
    if(class(num_vec) == "numeric") {
      num_vec <- num_vec[order(num_vec, decreasing = TRUE)]
      if(identical(EBGM_sort, num_vec)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  }
  compare_vars <- unlist(lapply(processed, compare_func, EBGM_sort = EBGM))
  if(any(compare_vars)) {
    names(processed)[compare_vars] <- "EBGM"
  }
  ##End check
  if(!is.null(quantiles)) {
    quant_calcs <- lapply(quantiles, quantBisect, theta_hat = theta,
                          N = N, E = E, qn = q_n, digits = digits)
    quant_df <- do.call(cbind, quant_calcs)
    quant_df <- as.data.frame(quant_df)
    quantiles <- ifelse(nchar(quantiles) == 1, paste("0", quantiles, sep = ""), quantiles)
    names(quant_df) <- paste("QUANT_", quantiles, sep = "")
  }
  ebout <- list()
  class(ebout) <- "openEBGM"
  if(any(grepl("EBGM", names(processed)))) {
    if(!is.null(quantiles)) {
      data_df <- cbind(processed, quant_df)
    } else {
      data_df <- processed
    }
  } else {
    if(!is.null(quantiles)) {
      data_df <- cbind(processed, EBGM, quant_df)
    } else {
      data_df <- cbind(processed, EBGM)
    }
  }
  ebout$data <- data_df
  ebout$hyper_parameters <- hyper_estimate
  ebout$quantiles <- quantiles
  ebout
}
