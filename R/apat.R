#' @name apat
#' @title A + t(A)
#' @param A numeric matrix
#' @details Not exported; exists for testing C code.
#' @return numeric matrix (A + t(A))
#' @useDynLib bcTSNE

apat <- function(A) {
  N <- NROW(A)
  res <- .C("apat", as.numeric(A), N)
  matrix(res[[1]], N, N)
}