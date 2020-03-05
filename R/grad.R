#' @name grad
#' @title Calculate t-SNE gradient
#' @param Y numeric matrix, lower dimension embedding
#' @param pval numeric matrix, input data p-values
#' @param Z numeric covariate matrix
#' @details Not exported; exists for testing C code.
#' @return numeric matrix, t-SNE gradient
#' @useDynLib bcTSNE

grad <- function(Y, pval, Z) {
  stopifnot(as.numeric(Y) && as.numeric(pval) && as.numeric(Z))
  N <- NROW(Y)
  J <- NCOL(Y)
  K <- NCOL(Z)
  res <- .C("grad", Y, pval, N, J, Z, K, matrix(0.0, N, J))
  matrix(res[[7]], N, J)
}