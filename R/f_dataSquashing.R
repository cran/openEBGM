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
#'   value of 100 or more for \code{keep_pts}.
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
#' @seealso \code{\link{processRaw}} for data preparation and
#'   \code{\link{autoSquash}} for automatically squashing an entire data set in
#'   one function call
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


#' Automated data squashing
#'
#' \code{autoSquash} squashes data by calling \code{\link{squashData}} once for
#'   each count (\emph{N}), removing the need to repeatedly squash the same data
#'   set.
#'
#' @param data A data frame (typically from \code{\link{processRaw}}) containing
#'   columns named \emph{N}, \emph{E}, and (possibly) \emph{weight}. Can contain
#'   additional columns, which will be ignored.
#' @param keep_pts A vector of whole numbers for the number of points to leave
#'   unsquashed for each count (\emph{N}). See the 'Details' section.
#' @param cut_offs A vector of whole numbers for the cutoff values of unsquashed
#'   data used to determine how many "super points" to end up with after
#'   squashing each count (\emph{N}). See the 'Details' section.
#' @param num_super_pts A vector of whole numbers for the number of
#'   "super points" to end up with after squashing each count (\emph{N}). Length
#'   must be 1 more than length of \code{cut_offs}. See the 'Details' section.
#' @return A data frame with column names \emph{N}, \emph{E}, and
#'   \emph{weight} containing the reduced data set.
#'
#' @details See \code{\link{squashData}} for details on squashing a given
#'   count (\emph{N}).
#' @details The elements in \code{keep_pts} determine how many points are left
#'   unsquashed for each count (\emph{N}). The first element in \code{keep_pts}
#'   is used for the smallest \emph{N} (usually 1). Each successive element is
#'   used for each successive \emph{N}. Once the last element is reached, it is
#'   used for all other \emph{N}.
#' @details For counts that are squashed, \code{cut_offs} and
#'   \code{num_super_pts} determine how the points are squashed. For instance,
#'   by default, if a given \emph{N} contains less than 500 points to be
#'   squashed, then those points are squashed to 50 "super points".
#'
#' @examples
#' data(caers)
#' proc <- processRaw(caers)
#' table(proc$N)
#'
#' squash1 <- autoSquash(proc)
#' ftable(squash1[, c("N", "weight")])
#'
#' squash2 <- autoSquash(proc, keep_pts = c(50, 5))
#' ftable(squash2[, c("N", "weight")])
#'
#' squash3 <- autoSquash(proc, keep_pts = 100,
#'                       cut_offs = c(250, 500),
#'                       num_super_pts = c(20, 60, 125)
#' )
#' ftable(squash3[, c("N", "weight")])
#'
#' @references DuMouchel W, Pregibon D (2001). "Empirical Bayes Screening for
#'   Multi-item Associations." In \emph{Proceedings of the Seventh ACM SIGKDD
#'   International Conference on Knowledge Discovery and Data Mining}, KDD '01,
#'   pp. 67-76. ACM, New York, NY, USA. ISBN 1-58113-391-X.
#' @keywords openEBGM
#' @seealso \code{\link{processRaw}} for data preparation and
#'   \code{\link{squashData}} for squashing individual counts
#' @export
autoSquash <-
  function(data, keep_pts = c(100, 75, 50, 25),
           cut_offs = c(500, 1e+03, 1e+04, 1e+05, 5e+05, 1e+06, 5e+06),
           num_super_pts = c(50, 75, 150, 500, 750, 1e+03, 2e+03, 5e+03)) {
  #Automatically squashes data for all counts

  .checkInputs_autoSquash(data, keep_pts, cut_offs, num_super_pts)

  len_keep_pts  <- length(keep_pts)
  last_keep_pts <- keep_pts[len_keep_pts]
  len_cut_offs  <- length(cut_offs)
  last_num_super_pts <- num_super_pts[length(num_super_pts)]

  numKeepPts <- function(count_position) {
    #Finds number of points to keep unsquashed for a given count
    if (count_position <= len_keep_pts) {
      return(keep_pts[count_position])
    } else {
      return(last_keep_pts)
    }
  }

  numSuperPts <- function(num_original, count_position) {
    #Finds number of "super points" for a given count
    keep_pts_val  <- numKeepPts(count_position)
    pts_to_squash <- num_original - keep_pts_val
    cut_off_index <- which(pts_to_squash < c(cut_offs, Inf))[1]
    num_super_pts[cut_off_index]
  }

  binSize <- function(num_original, count_position) {
    #Finds bin size for a given count
    keep_pts_val      <- numKeepPts(count_position)
    num_super_pts_val <- numSuperPts(num_original, count_position)
    bin_size <- (num_original - keep_pts_val) / num_super_pts_val
    max(floor(bin_size), 2L)
  }

  data$weight <- 1L  #in case no squashing occurs
  count_table <- table(data$N)
  count_table <- count_table[count_table > max(keep_pts) + 1]
  counts      <- names(count_table)
  for (i in counts) {
    num_pts    <- count_table[[i]]
    i_position <- match(i, counts)
    num_supers <- numSuperPts(num_pts, i_position)
    bin_size   <- binSize(num_pts, i_position)
    keep_pts_value <- numKeepPts(i_position)
    data <- squashData(data, count = i, bin_size = bin_size,
                       keep_pts = keep_pts_value, min_bin = 1, min_pts = 1)
  }
  data
}
