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


#' @name calcPvals
#' @title Calculate t-SNE p-values based on a distance matrix
#' @param D numeric matrix, distance matrix
#' @param perplexity numeric of length 1, t-SNE perplexity 
#' @details Not exported; exists for testing C code.
#' @return numeric matrix of p-values based on the given perplexity
#' @useDynLib bcTSNE

calcPvals <- function(D, perplexity = 30) {
  N <- NROW(D)
  stopifnot(is.numeric(perplexity) && is.numeric(D)) 
  res <- .C("calcPvals", as.numeric(D), N, log(perplexity), rep(0.0, N*N))
  matrix(res[[4]], N, N)
}


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