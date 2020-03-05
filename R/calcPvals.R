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


