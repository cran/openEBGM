#' Squash data for hyperparameter estimation
#'
#' \code{squashData} squashes data by binning expected counts, \emph{E}, for a
#'   given actual count, \emph{N}, using bin means as the expected counts for
#'   the reduced data set. The squashed points are weighted by bin size. Data
#'   can be squashed to reduce computational burden (see DuMouchel et al.,
#'   2001) when estimating the hyperparameters.
#'
#' @param data A data frame (typically from \code{\link{processRaw}} or a
#'   previous call to \code{\link{squashData}}) containing columns named
#'   \emph{N}, \emph{E}, and (possibly) \emph{weight}. Can contain additional
#'   columns, which will be ignored.
#' @param count A non-negative scalar whole number for the count size, \emph{N},
#'   used for binning
#' @param bin_size A scalar whole number (>= 2)
#' @param keep_pts A nonnegative scalar whole number for number of points with
#'   the largest expected counts to leave unsquashed. Used to help prevent
#'   \dQuote{oversquashing}.
#' @param min_bin A positive scalar whole number for the minimum number of bins
#'   needed. Used to help prevent \dQuote{oversquashing}.
#' @param min_pts A positive scalar whole number for the minimum number of
#'   original (unsquashed) points needed for squashing. Used to help prevent
#'   \dQuote{oversquashing}.
#' @return A data frame with column names \emph{N}, \emph{E}, and
#'   \emph{weight} containing the reduced data set.
#'
#' @details Can be used iteratively (count = 1, then 2, etc.).
#' @details The \emph{N} column in \code{data} will be coerced using
#'   \code{\link{as.integer}}, and \emph{E} will be coerced using
#'   \code{\link{as.numeric}}. Missing data are not allowed.
#' @details Since the distribution of expected counts, \emph{E}, tends to be
#'   skewed to the right, the largest \emph{E}s are not squashed by default.
#'   This behavior can be changed by setting the \code{keep_pts} argument to
#'   zero (0); however, this is not recommended. Squashing the largest \emph{E}s
#'   could result in a large loss of information, so it is recommended to use a
#'   value of one (100) or more for \code{keep_pts}.
#' @details Values for \code{keep_pts}, \code{min_bin}, and \code{min_pts}
#'   should typically be at least as large as the default values.
#' @examples
#' set.seed(483726)
#' dat <- data.frame(var1 = letters[1:26], var2 = LETTERS[1:26],
#'                   N = c(rep(0, 11), rep(1, 10), rep(2, 4), rep(3, 1)),
#'                   E = round(abs(c(rnorm(11, 0), rnorm(10, 1), rnorm(4, 2),
#'                             rnorm(1, 3))), 3)
#'                   )
#' (zeroes <- squashData(dat, count = 0, bin_size = 3, keep_pts = 1,
#'                       min_bin = 2, min_pts = 2))
#' (ones <- squashData(zeroes, bin_size = 2, keep_pts = 1,
#'                     min_bin = 2, min_pts = 2))
#' (twos <- squashData(ones, count = 2, bin_size = 2, keep_pts = 1,
#'                     min_bin = 2, min_pts = 2))
#'
#' squashData(zeroes, bin_size = 2, keep_pts = 0,
#'            min_bin = 2, min_pts = 2)
#' squashData(zeroes, bin_size = 2, keep_pts = 1,
#'            min_bin = 2, min_pts = 2)
#' squashData(zeroes, bin_size = 2, keep_pts = 2,
#'            min_bin = 2, min_pts = 2)
#' squashData(zeroes, bin_size = 2, keep_pts = 3,
#'            min_bin = 2, min_pts = 2)
#'
#' @references DuMouchel W, Pregibon D (2001). "Empirical Bayes Screening for
#'   Multi-item Associations." In \emph{Proceedings of the Seventh ACM SIGKDD
#'   International Conference on Knowledge Discovery and Data Mining}, KDD '01,
#'   pp. 67-76. ACM, New York, NY, USA. ISBN 1-58113-391-X.
#' @keywords openEBGM
#' @seealso \code{\link{processRaw}} for data preparation
#' @import data.table
#' @export
squashData <- function(data, count = 1, bin_size = 50, keep_pts = 100,
                       min_bin = 50, min_pts = 500) {

  #Check inputs & coerce to data table
  data <- .checkInputs_squashData(data, count, bin_size, keep_pts, min_bin,
                                  min_pts)

  hold         <- data[N != count, .(N, E, weight)]  #not squashed
  maybe_squash <- data[N == count, .(N, E, weight)]
  data.table::setorder(maybe_squash, E)  #need similar Es in same bin
  num_maybe_squash <- nrow(maybe_squash)

  #Largest values are most variable, so keep (i.e., do not squash)
  num_squash <- num_maybe_squash - keep_pts
  if (num_squash < bin_size) {
    stop("reduce 'bin_size' or 'keep_pts'")
  }
  squash     <- maybe_squash[1:num_squash, ]
  keep       <- maybe_squash[(num_squash + 1):num_maybe_squash, ]
  num_remain <- num_squash %% bin_size  #partial bin count

  #Create bin indices and squash points
  if (num_remain == 0) {
    num_bins <- num_squash / bin_size
    squash[, bin_index := rep(1:num_bins, each = bin_size)]
    squash <- squash[, j = list(E = mean(E, na.rm = TRUE)),
                     by = .(N, bin_index)]
    squash[, bin_index := NULL]
    squash[, weight := rep.int(bin_size, num_bins)]
  } else {
    num_bins_full    <- floor(num_squash / bin_size)
    bin_full_index   <- rep(1:num_bins_full, each = bin_size)
    bin_remain_index <- rep.int(num_bins_full + 1, num_remain)
    squash[, bin_index := c(bin_full_index, bin_remain_index)]
    squash <- squash[, j = list(E = mean(E, na.rm = TRUE)),
                     by = .(N, bin_index)]
    squash[, bin_index := NULL]
    weights_full <- rep.int(bin_size, num_bins_full)
    squash[, weight := c(weights_full, num_remain)]
  }

  if (keep_pts == 0) {
    results <- data.table::rbindlist(list(hold, squash))
  } else {
    results <- data.table::rbindlist(list(hold, squash, keep))
  }

  data.table::setorder(results, N)
  row.names(results) <- NULL
  data.table::setDF(results)
  results
}

#Hack to trick 'R CMD check'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "weight", "bin_index"))
}
