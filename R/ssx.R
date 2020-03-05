#' @name ssx
#' @title Sum of squares
#' @param X numeric matrix
#' @details Not exported; exists for testing C code.
#' @return vector with the row sum of squares
#' @useDynLib bcTSNE

ssx <- function(X) {
  N <- NROW(X)
  res <- .C("ssx", as.numeric(X), N, NCOL(X), rep(0.0, N))
  res[[4]]
}