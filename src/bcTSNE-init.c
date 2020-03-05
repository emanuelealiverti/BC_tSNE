#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "Rdefines.h"

void ssx(double *X, int *N, int *D, double *SSX);
void sqdist(double *X, int *N, int *D, double *dist);
void calcPvals(double *dist, int *N, double *lp, double *pval);
void apat(double *X, int *N);
void grad(double *Y, double *pval, int *N, int *J, 
          double *Z, int *K, double *dY);
void ols(double *A, int *M, int *N, double *B, int *NRHS, double *W);
void zeroMean(double *X, int *N, int *D);
void bctsne(double *X, int *N, int *D, double *Z, int *K,
            int *J, double *lp, double *Y, int *maxIter);

static const R_CMethodDef cMethods[] = {
  {"ssx",       (DL_FUNC) &ssx,         4},
  {"sqdist",    (DL_FUNC) &sqdist,      4},
  {"calcPvals", (DL_FUNC) &calcPvals,   4},
  {"apat",      (DL_FUNC) &apat,        2},
  {"grad",      (DL_FUNC) &grad,        7},
  {"ols",       (DL_FUNC) &ols,         6},
  {"zeroMean",  (DL_FUNC) &zeroMean,    3},
  {"bctsne",    (DL_FUNC) &bctsne,      9},
  {NULL, NULL, 0}
};

void R_init_bcTSNE(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE); 
  R_forceSymbols(dll, TRUE);
}
