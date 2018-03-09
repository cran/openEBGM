#' Squash data for hyperparameter estimation
#'
#' \code{squashData} squashes data by binning expected counts, \emph{E}, for a
#'   given actual count, \emph{N}, using bin means as the expected counts for
#'   the reduced data set. The squashed points are weighted by bin size. Data
#'   can be squashed to reduce computational burden (see DuMouchel et al.,
#'   2001) when estimating the hyperparameters.
#'
#' @param data A data frame containing columns named \emph{N}, \emph{E},
#'   and (possibly) \emph{weight}. Can contain additional columns, which will
#'   be ignored.
#' @param count A non-negative scalar whole number for the count size, \emph{N},
#'   used for binning
#' @param bin_size A scalar whole number (>= 2)
#' @param keep_bins A nonnegative scalar whole number for number of bins of the
#'   largest expected counts to leave unsquashed. Used to help prevent
#'   \dQuote{oversquashing}.
#' @param min_bin A positive scalar whole number for the minimum number of bins
#'   needed. Used to help prevent \dQuote{oversquashing}.
#' @param min_pts A positive scalar whole number for the minimum number of
#'   points needed for squashing. Used to help prevent \dQuote{oversquashing}.
#' @return A data frame with column names \emph{N}, \emph{E}, and
#'   \emph{weight} containing the reduced data set.
#'
#' @details Can be used iteratively (count = 1, then 2, etc.).
#' @details Typically, \code{\link{processRaw}} is used to create the data frame
#'   supplied to the \code{data} argument.
#' @details The \emph{N} column in \code{data} will be coerced using
#'   \code{\link{as.integer}}, and \emph{E} will be coerced using
#'   \code{\link{as.numeric}}. Missing data are not allowed.
#' @details Since the distribution of expected counts, \emph{E}, tends to be
#'   skewed to the right, the largest \emph{E}s are not squashed by default.
#'   This behavior can be changed by setting the \code{keep_bins} argument to
#'   zero (0); however, this is not recommended. Squashing the largest \emph{E}s
#'   could result in a large loss of information, so it is recommended to use a
#'   value of one (1) or more for \code{keep_bins}.
#' @details Values for \code{keep_bins}, \code{min_bin}, and \code{min_pts}
#'   should typically be at least as large as the default values.

#' @examples
#' set.seed(483726)
#' dat <- data.frame(var1 = letters[1:26], var2 = LETTERS[1:26],
#'                   N = c(rep(0, 11), rep(1, 10), rep(2, 4), rep(3, 1)),
#'                   E = round(abs(c(rnorm(11, 0), rnorm(10, 1), rnorm(4, 2),
#'                             rnorm(1, 3))), 3)
#'                   )
#' (zeroes <- squashData(dat, count = 0, bin_size = 3, keep_bins = 1,
#'                       min_bin = 2, min_pts = 2))
#' (ones <- squashData(zeroes, bin_size = 2, keep_bins = 1,
#'                     min_bin = 2, min_pts = 2))
#' (twos <- squashData(ones, count = 2, bin_size = 2, keep_bins = 1,
#'                     min_bin = 2, min_pts = 2))
#'
#' squashData(dat, count = 0, bin_size = 3, keep_bins = 0,
#'            min_bin = 2, min_pts = 2)
#' squashData(dat, count = 0, bin_size = 3, keep_bins = 1,
#'            min_bin = 2, min_pts = 2)
#' squashData(dat, count = 0, bin_size = 3, keep_bins = 2,
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
squashData <- function(data, count = 1, bin_size = 50, keep_bins = 2,
                       min_bin = 50, min_pts = 500) {

  #Check inputs & coerce to data table
  data <- .checkInputs_squashData(data, count, bin_size, keep_bins, min_bin,
                                  min_pts)

  hold <- data[N != count, .(N, E, weight)]  #not squashed
  maybe_squash <- data[N == count, .(N, E, weight)]  #squash most points
  maybe_squash <- maybe_squash[order(E), ]  #need similar Es in same bin

  #Largest values are most variable, so keep (i.e., do not squash)
  num_squash <- nrow(maybe_squash) - (keep_bins * bin_size)
  if (num_squash < bin_size) {
    stop("reduce 'bin_size' or 'keep_bins'")
  }
  squash    <- maybe_squash[1:num_squash, ]
  keep      <- maybe_squash[(num_squash + 1):nrow(maybe_squash), ]
  remainder <- nrow(squash) %% bin_size

  #Create bins and squash points
  if (remainder == 0) {
    num_bins <- nrow(squash) / bin_size
    bins   <- sort(rep.int(1:num_bins, bin_size))
    squash <- squash[, bins := bins]
    squash <- squash[, mean(E, na.rm = TRUE), by = .(N, bins)]
    weight <- rep(bin_size, num_bins)
    squash <- squash[, weight := weight]
  } else {
    num_bins <- (nrow(squash) / bin_size) + 1
    bins   <- c(sort(rep.int(1:(num_bins - 1), bin_size)),
                rep.int(num_bins, remainder))
    squash <- squash[, bins := bins]
    squash <- squash[, mean(E, na.rm = TRUE), by = .(N, bins)]
    weight <- c(rep(bin_size, num_bins - 1), remainder)
    squash <- squash[, weight := weight]
  }
  squash <- squash[, bins := NULL]

  if (keep_bins == 0) {
    results <- data.table::rbindlist(list(hold, squash))
  } else {
    results <- data.table::rbindlist(list(hold, squash, keep))
  }

  results <- results[order(N), ]
  row.names(results) <- NULL
  as.data.frame(results)
}

#Hack to trick 'R CMD check'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "weight", "bins"))
}
