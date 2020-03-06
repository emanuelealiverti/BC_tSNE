#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <Rmath.h> 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h> 
#define LAPACK_COL_MAJOR  102

void ssx(double *X, int *N, int *D, double *SSX) {
  
  /* X is an array (1d column-major representation of 2d array); N is number of 
   * rows in the 2d X array; D is number of columns. SSX is the row-wise 
   * sum of squares. 
   */
  
  for (int n = 0; n < N[0]; n++) {
    for (int d = 0; d < D[0]; d++ ) {
      SSX[n] += X[d*N[0] + n] * X[d*N[0] + n];
    }
  }
  
}

void sqdist(double *X, int *N, int *D, double *dist) {
  
  /* X is an array (1d column-major representation of 2d array); N is number of 
   * rows in the 2d X array; D is number of columns. dist is the 1d column-major
   * representation of the squared euclidian distance matrix of X (N x N). 
   */
  
  double *SSX = (double*) calloc(N[0], sizeof(double)); // sum of squares X
  if (SSX == NULL) error("Memory allocation failed.");  
  ssx(X, N, D, SSX);
  for (int d = 0; d < N[0]; d++) {
    for (int n = 0; n < N[0]; n++ ) {
      dist[d*N[0] + n] += SSX[n] + SSX[d];
    }
  }
  free(SSX);
  double a = -2.0;
  double b =  1.0;
  // a * op(A) %*% op(B) + b * C; op(A) (M x K), op(B) (K x N), C (M x N);
  // For mat M, LDM is the first dimension (LDM, *) of non-op matrix
  // TA gives op for A; TB gives op for B; "N" means none, "T" is transpose
  //        dgemm(TA,  TB,  M,  N,  K,  a,  A, LDA, B, LDB, b,  C,    LDC)
  F77_CALL(dgemm)("N", "T", N,  N,  D,  &a, X, N,   X, N,   &b, dist, N);
  
}

void apat(double *X, int *N) {
  // X + t(X); X: N x N matrix
  for (int d = 0; d < N[0]; d++) {
    for (int n = d; n < N[0]; n++) {
      X[d*N[0] + n] += X[n*N[0] + d];
      X[n*N[0] + d] = X[d*N[0] + n];
    }
  }
}

void calcPvals(double *dist, int *N, double *lp, double *pval) {
  
  /* dist is the distance matrix in column-major form; N is the row number; 
   * lp is the log(perplexity); pval is the p-value matrix in column-major
   * form
   */
  
  for(int d = 0; d < N[0]; d++) {
    bool found = false;
    double beta = 1.0;
    double minBeta = -DBL_MAX;
    double maxBeta =  DBL_MAX;
    double tol = 1e-5;
    double sumPval;
    int iter = 0;
    while (!found && iter < 200) {
      // Compute Gaussian kernel row
      for (int n = 0; n < N[0]; n++) {
        pval[d*N[0] + n] = exp(-beta*dist[d*N[0] + n]);
      }
      pval[d*N[0] + d] = 0.0;
      // Compute entropy of current row
      sumPval = 0.0;
      for (int n = 0; n < N[0]; n++) sumPval += pval[d*N[0] + n];
      double H = 0.0;
      for (int n = 0; n < N[0]; n++) {
        H += beta*dist[d*N[0] + n]*pval[d*N[0] + n];
      }
      H = (H/sumPval) + log(sumPval);
      // Evaluate whether the entropy is within the tolerance level
      double Hdiff = H - lp[0];
      if (fabs(Hdiff) < tol) {
        found = true;
      } else {
        if (Hdiff > 0) {
          minBeta = beta;
          if (maxBeta == DBL_MAX || maxBeta == -DBL_MAX) beta *= 2.0;
          else beta = (beta + maxBeta) / 2.0;
        } else {
          maxBeta = beta;
          if (minBeta == -DBL_MAX || minBeta == DBL_MAX) beta /= 2.0;
          else beta = (beta + minBeta) / 2.0;
        }
      }
      iter++;
    }
    // Column normalize P
    for (int n = 0; n < N[0]; n++) pval[d*N[0] + n] /= sumPval;
  }
  
  // Make pvals symmetric; normalize so all valuess sum to 1
  apat(pval, N);
  double sumPval = 0;
  for (int i = 0; i < N[0]*N[0]; i++) sumPval += pval[i];
  for (int i = 0; i < N[0]*N[0]; i++) pval[i] /= sumPval;
  
}

void ols(double *A, int *M, int *N, double *B, int *NRHS, double *W) {
  
  // Solves B = AW for W. A (M x N), B (M x NRHS), W (N x NRHS)
  
  int info;
  int lwork = -1; // put in "query mode"
  double worksize;
  /* copes for dY and Z, which will be destroyed during fitting; dYcopy will
   * contain the beta values -- it's confusing.
   */
  double *Ac = (double*) malloc(M[0]*N[0]*sizeof(double));
  double *Bc = (double*) malloc(M[0]*NRHS[0]*sizeof(double));
  if (Ac == NULL || Bc == NULL) error("Memory allocation failed.");
  for (int i = 0; i < M[0]*N[0];    i++) Ac[i] = A[i];
  for (int i = 0; i < M[0]*NRHS[0]; i++) Bc[i] = B[i];
  // Solves B = AX for X; A (M x N), B (M x NRHS)
  //        dgels(TA,  M, N, NRHS, A,  LDA, B,  LDB, WORK,      LWORK,  INFO)
  F77_CALL(dgels)("N", M, N, NRHS, Ac, M,   Bc, M,   &worksize, &lwork, &info);
  int nwork = (int) worksize;
  double *work = (double*) malloc(nwork*sizeof(double)); 
  if (work == NULL) error("Memory allocation failed.");
  // Solve the system of equations.
  F77_CALL(dgels)("N", M, N, NRHS, Ac, M,   Bc, M,   work,      &nwork, &info);
  if (info != 0) error("OLS step to correct for batch failed.");
  // Get the beta values out of dYc -- we need the first K rows and J columns
  for (int i = 0; i < NRHS[0]; i++) {
    for (int j = 0; j < N[0]; j++) {
      W[i*N[0] + j] = Bc[i*M[0] + j];
    }
  }
  free(Ac);
  free(Bc);
}

void grad(double *Y, double *pval, int *N, int *J, 
          double *Z, int *K, double *dY) {
  
  /* Y is an array (1d column-major representation of 2d array) with the 
   * reduced dimension embeddings of X. N is number of rows in the 2d Y array; 
   * J is number of columns. Z is the N x K model matrix (1d column major) 
   * with covariate information. dY is the gradient. pval is the p-values 
   * calculated on X.   
   */
  
  // Re-zero gradient
  for (int i = 0; i < N[0]*J[0]; i++) dY[i] = 0.0;
  
  // Compute the squared euclidian distances (dist)
  double *dist = (double*) calloc(N[0]*N[0], sizeof(double)); // Y distances
  if (dist == NULL) error("Memory allocation failed.");
  sqdist(Y, N, J, dist);
  
  // Compute q-values
  double *qval = (double*) calloc(N[0]*N[0], sizeof(double)); // q-values 
  if (qval == NULL) error("Memory allocation failed.");
  for (int d = 0; d < N[0]; d++) { // iterate columns
    for (int n = 0; n < N[0]; n++) { // iterate rows
      if (d == n) {
        continue;
      } else {
        qval[d * N[0] + n] = 1.0/(1.0 + dist[d * N[0] + n]);
      }
    }
  }
  double sumQ = 0.0;
  for (int i = 0; i < N[0]*N[0]; i++) sumQ += qval[i];
  
  for (int d = 0; d < N[0]; d++) {
    for (int n = 0; n < N[0]; n++) {
      if (n == d) continue;
      double mult = (pval[d*N[0] + n] - qval[d*N[0] + n]/sumQ)*qval[d*N[0] + n];
      for (int j = 0; j < J[0]; j++) {
        dY[j*N[0] + d] += ((Y[j*N[0] + d] - Y[j*N[0] + n])*mult);
      }
    }
  }
  free(dist);
  free(qval);
  
  /* Perform OLS dY ~ ZW; weights in (K x J) W
   * http://www.stat.cmu.edu/~hseltman/711/CfromR/bootReg2.cpp 
   */ 
  double *W = (double*) malloc(K[0]*J[0]*sizeof(double)); // beta/weight values
  if (W == NULL) error("Memory allocation failed.");
  // ols(double *A, int *M, int *N, double *B, int *NRHS, double *W)
  ols(Z, N, K, dY, J, W);
  
  // calculate dYhat & subtract from the gradient
  double *dYhat = (double*) calloc(N[0]*J[0], sizeof(double)); 
  if (dYhat == NULL) error("Memory allocation failed.");
  double a = 1.0;
  double b = 0.0;
  // a * op(A) %*% op(B) + b * C; op(A) (M x K), op(B) (K x N), C (M x N);
  // For mat M, LDM is the first dimension (LDM, *) of non-op matrix
  // TA gives op for A; TB gives op for B; "N" means none, "T" is transpose
  //        dgemm(TA,  TB,  M,  N,  K,  a,   A, LDA, B, LDB, b,  C,     LDC)
  F77_CALL(dgemm)("N", "N", N,  J,  K,  &a,  Z, N,   W, K,   &b, dYhat, N);
  for (int i = 0; i < N[0]*J[0]; i++) dY[i] -= dYhat[i];
  free(dYhat);
  free(W);
}

void zeroMean(double *X, int *N, int *D) {
  
  /* X (N x D) is a column-major rep of 2d array; subtracts the column mean
   * from each column.
   */
  
  double *mean = (double*) calloc(D[0], sizeof(double));
  if (mean == NULL) error("Memory allocation failed.");
  for(int d = 0; d < D[0]; d++) {
    for(int n = 0; n < N[0]; n++) mean[d] += X[d*N[0] + n];
    mean[d] /= (double) N[0];
    for(int n = 0; n < N[0]; n++) X[d*N[0] + n] -= mean[d];
  }
  free(mean); 
}

void bctsne(double *X, int *N, int *D, double *Z, int *K,
            int *J, double *lp, double *Y, int *maxIter) {
  
  /* X is an array (1d column-major representation of 2d array) with the values 
   * to perform t-SNE; all transofrmations (e.g. PCA) are performed prior to the 
   * algorithm. N is number of rows in the 2d X array; D is number of columns. 
   * J is output dimension. lp is log perplexity. Z is the column-major
   * representation of the N x K covariate model matrix.  
   */
  
  // Compute the squared euclidian distances (dist) of X
  double *dist = (double*) calloc(N[0]*N[0], sizeof(double)); // distances
  if (dist == NULL) error("Memory allocation failed.");  
  sqdist(X, N, D, dist);
  
  // Compute p-values
  double *pval = (double*) calloc(N[0]*N[0], sizeof(double)); // p-values 
  double pvalExag = 4.0;
  if (pval == NULL) error("Memory allocation failed.");  
  calcPvals(dist, N, lp, pval);
  free(dist);
  // Exaggerate p-values and adjust very small values.
  for (int i = 0; i < N[0]*N[0]; i++) {
    pval[i] *= pvalExag; 
    if (pval[i] < 1e-12) pval[i] = 1e-12; 
  }
  
  // Initialize arrays for computing output, Y
  double *dY    = (double*) calloc(N[0]*J[0], sizeof(double));
  double *iY    = (double*) calloc(N[0]*J[0], sizeof(double));
  double *gains = (double*) malloc(N[0]*J[0]*sizeof(double));
  if (dY == NULL || iY == NULL || gains == NULL) {
    error("Memory allocation failed.");
  }  
  GetRNGstate(); // get .Random.seed
  for (int i = 0; i < N[0]*J[0]; i++) {
    Y[i] = rnorm(0.0, 1.0); // initialize Y with rnorm
    gains[i] = 1.0;
  }
  PutRNGstate(); // set .Random.seed
  
  double mom = 0.5; // initial momentum
  double eta = 200.0; // 
  double minGain = 0.01; // Coerce any smaller gains to 0.01
  
  // Iterate and update Y
  for (int iter = 0; iter < maxIter[0]; iter++) {
    
    if (iter == 250) mom = 0.8;
    if (iter == 100) {
      for (int i = 0; i < N[0]*N[0]; i++) {
        pval[i] /= pvalExag; 
        if (pval[i] < 1e-12) pval[i] = 1e-12; 
      }
    }
    
    // Calculate dY
    grad(Y, pval, N, J, Z, K, dY);
    for (int i = 0; i < N[0]*J[0]; i++) {
      gains[i] = dY[i] > 0 && iY[i] > 0 ? gains[i]*0.8 : gains[i] + 0.2;
      if (gains[i] < minGain) gains[i] = minGain;
      iY[i] = mom*iY[i] - eta*(gains[i]*dY[i]);
      Y[i] += iY[i];
    }
    
    // Re-center Y
    zeroMean(Y, N, J);
    
  }
  
  free(dY);
  free(iY);
  free(gains);
  free(pval);
  
}



