#' Process raw data
#'
#' \code{processRaw} finds the actual and expected counts using the methodology
#' described by DuMouchel (1999); however, an adjustment is applied to expected
#' counts to prevent double-counting (i.e., using unique marginal ID counts
#' instead of contingency table marginal totals). Also calculates the relative
#' reporting ratio (\emph{RR}) and the proportional reporting ratio
#' (\emph{PRR}).
#'
#' @param data A data frame containing columns named: \emph{id}, \emph{var1},
#'   and \emph{var2}. Possibly includes columns for stratification variables
#'   identified by the substring \emph{'strat'} (e.g. \emph{strat_age}). Other
#'   columns will be ignored.
#' @param stratify A logical scalar (TRUE or FALSE) specifying if stratification
#'   is to be used.
#' @param zeroes A logical scalar specifying if zero counts should be included.
#'   Using zero counts is only recommended for small data sets because it will
#'   dramatically increase run time.
#' @param digits A whole number scalar specifying how many decimal places to
#'   use for rounding \emph{RR} and \emph{PRR}.
#' @param max_cats A whole number scalar specifying the maximum number of
#'   categories to allow in any given stratification variable. Used to help
#'   prevent a situation where the user forgets to categorize a continuous
#'   variable, such as age.
#' @return A data frame with actual counts (\emph{N}), expected counts
#'   (\emph{E}), relative reporting ratio (\emph{RR}), and proportional
#'   reporting ratio (\emph{PRR}) for \emph{var1-var2} pairs.
#'
#' @details An \emph{id} column must be included in \code{data}. If your data
#'   set does not include IDs, make a column of unique IDs using \code{df$id <-
#'   1:nrow(df)}. However, unique IDs should only be constructed if the cells in
#'   the contingency table are mutually exclusive. For instance, unique IDs for
#'   each row in \code{data} are not appropriate with CAERS data since a given
#'   report can include multiple products and/or adverse events.
#' @details Stratification variables are identified by searching for the
#'   substring \emph{'strat'}. Only variables containing \emph{'strat'} (case
#'   sensitive) will be used as stratification variables. \emph{PRR}
#'   calculations ignore stratification, but \emph{E} and \emph{RR} calculations
#'   are affected. A warning will be displayed if any stratum contains less than
#'   50 unique IDs.
#' @details If a \emph{PRR} calculation results in division by zero, \code{Inf}
#'   is returned.
#' @section Warnings:
#'   Use of the \code{zeroes = TRUE} option will result in a considerable
#'   increase in runtime. Using zero counts is not recommended if the contingency
#'   table is moderate or large in size (~500K cells or larger). However, using
#'   zeroes could be useful if the optimization algorithm fails to converge
#'   when estimating hyperparameters.
#' @section Warnings:
#' Any columns in \code{data} containing the substring \emph{'strat'} in the
#'   column name will be considered stratification variables, so verify that you
#'   do not have any extraneous columns with that substring.
#'
#' @examples
#' dat <- data.frame(
#'   var1 = c("product_A", rep("product_B", 3), "product_C",
#'            rep("product_A", 2), rep("product_B", 2), "product_C"),
#'   var2 = c("event_1", rep("event_2", 2), rep("event_3", 2),
#'            "event_2", rep("event_3", 3), "event_1"),
#'   strat1 = c(rep("Male", 5), rep("Female", 3), rep("Male", 2)),
#'   strat2 = c(rep("age_cat1", 5), rep("age_cat1", 3), rep("age_cat2", 2))
#' )
#' dat$id <- 1:nrow(dat)
#' processRaw(dat, digits = 3)
#' suppressWarnings(
#'   processRaw(dat, stratify = TRUE, digits = 3)
#' )
#' processRaw(dat, zeroes = TRUE, digits = 3)
#' suppressWarnings(
#'   processRaw(dat, stratify = TRUE, zeroes = TRUE, digits = 3)
#' )
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @keywords openEBGM
#' @import data.table
#' @export
processRaw <- function(data, stratify = FALSE, zeroes = FALSE, digits = 2,
                       max_cats = 10) {

  .checkInputs_processRaw(data, stratify, zeroes)

  data <- data.table::as.data.table(data)
  data <- unique(data)

  #Actual var1/var2 combination counts
  actual   <- data[, .countUnique(id), by = .(var1, var2)]
  actual$N <- actual$V1; actual$V1 <- NULL
  if (zeroes == TRUE) {
    data.table::setkeyv(actual, c("var1", "var2"))
    actual <- actual[data.table::CJ(unique(var1), unique(var2))]
    actual <- actual[is.na(N), N := 0L]
  }

  #Unstratified marginal counts
  v1_marg      <- data[, .countUnique(id), by = .(var1)]
  v1_marg$N_v1 <- v1_marg$V1; v1_marg$V1 <- NULL
  v2_marg      <- data[, .countUnique(id), by = .(var2)]
  v2_marg$N_v2 <- v2_marg$V1; v2_marg$V1 <- NULL

  #Expected counts
  if (stratify == FALSE) {
    counts       <- merge(actual, v1_marg, by = "var1",
                          all.x = zeroes, sort = FALSE)
    counts       <- merge(counts, v2_marg, by = "var2",
                          all.x = zeroes, sort = FALSE)
    counts$N_tot <- .countUnique(data$id)
    counts$E     <- (counts$N_v1 * counts$N_v2) / counts$N_tot
  } else {
    data <- .checkStrata_processRaw(data, max_cats)  #adds 'stratum' column
    v1_marg_str          <- data[, .countUnique(id), by = .(var1, stratum)]
    v1_marg_str$N_v1_str <- v1_marg_str$V1; v1_marg_str$V1 <- NULL
    v2_marg_str          <- data[, .countUnique(id), by = .(var2, stratum)]
    v2_marg_str$N_v2_str <- v2_marg_str$V1; v2_marg_str$V1 <- NULL
    strat_tot            <- data[, .countUnique(id), by = .(stratum)]
    strat_tot$N_tot_str  <- strat_tot$V1; strat_tot$V1 <- NULL
    counts <- merge(actual, v1_marg_str, by = "var1",
                    all.x = zeroes, sort = FALSE, allow.cartesian = TRUE)
    counts <- counts[is.na(N_v1_str), N_v1_str := 0L]
    counts <- merge(counts, v2_marg_str, by = c("var2", "stratum"),
                    all.x = zeroes, sort = FALSE)
    counts <- counts[is.na(N_v2_str), N_v2_str := 0L]
    counts <- merge(counts, strat_tot, by = "stratum",
                    all.x = zeroes, sort = FALSE)
    counts$E <- (counts$N_v1_str * counts$N_v2_str) / counts$N_tot_str
    counts <- counts[, sum(E, na.rm = TRUE), by = .(var1, var2)]
    counts$E <- counts$V1; counts$V1 <- NULL
    counts <- merge(actual, counts, by = c("var1", "var2"),
                    all.x = zeroes, sort = FALSE)
    counts <- merge(counts, v1_marg, by = "var1", all.x = zeroes, sort = FALSE)
    counts <- merge(counts, v2_marg, by = "var2", all.x = zeroes, sort = FALSE)
    counts$N_tot <- .countUnique(data$id)
  }

  #Add disproportionality measures
  counts$RR  <- round(counts$N / counts$E, digits)
  counts$RR  <- ifelse(is.nan(counts$RR), 0, counts$RR)
  PRR_num    <- counts$N / counts$N_v1
  PRR_den    <- (counts$N_v2 - counts$N) / (counts$N_tot - counts$N_v1)
  counts$PRR <- round(PRR_num / PRR_den, digits)
  counts$PRR <- ifelse(is.nan(counts$PRR), Inf, counts$PRR)
  vars        <- c("var1", "var2", "N", "E", "RR", "PRR")
  counts$var1 <- as.factor(counts$var1)
  counts$var2 <- as.factor(counts$var2)
  counts      <- counts[order(var1, var2), vars, with = FALSE]
  as.data.frame(counts)
}

#Hack to trick 'R CMD check'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "id", "var1", "var2", "N", "stratum",
                           "N_v1_str", "N_v2_str", "E"))
}
