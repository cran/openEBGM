#Various helper functions

.isMissing_str <- function(x) {
  #Determines if any elements in a string vector are missing.
  any(is.na(x) | x == "")
}

.isMissing_num <- function(x) {
  #Determines if any elements in a numeric vector are missing.
  any(is.na(x) | is.nan(x) | is.infinite(x))
}

.countUnique <- function(x) {
  #Counts the number of unique elements in a vector.
  length(unique(x))
}
