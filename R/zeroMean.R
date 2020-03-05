#' @name zeroMean
#' @title Subtract the column means from X
#' @param X numeric matrix
#' @details Not exported; exists for testing C code.
#' @return numeric matrix with column means subtracted
#' @useDynLib bcTSNE

zeroMean <- function(X) {
  N <- NROW(X)
  D <- NCOL(X)
  res <- .C("zeroMean", as.numeric(X), N, D)
  matrix(res[[1]], N, D)
}