#' @name sqdist
#' @title Calculate squared Euclidean distance
#' @param X numeric matrix
#' @details Not exported; exists for testing C code.
#' @return numeric squared distance matrix
#' @useDynLib bcTSNE

sqdist <- function(X) {
  N <- NROW(X)
  res <- .C("sqdist", as.numeric(X), N, NCOL(X), rep(0.0, N*N))
  matrix(res[[4]], N, N)
}