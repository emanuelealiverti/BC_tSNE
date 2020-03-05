#' @name ols
#' @title Ordinary least squares, solves B = AX for X.
#' @param A numeric matrix
#' @param B numeric matrix
#' @details Not exported; exists for testing C code.
#' @return numeric matrix (X)
#' @useDynLib bcTSNE

ols <- function(A, B) {
  M <- NROW(A)
  N <- NCOL(A)
  NRHS <- NCOL(B)
  res <- .C("ols", as.numeric(A), M, N, as.numeric(B), NRHS, matrix(0, N, NRHS))
  list(W = res[[6]])
}