#include <iomanip>
#include <string>
#include <limits>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

///////////////////////////////////////////////////////////////////
//u index 
///////////////////////////////////////////////////////////////////
extern "C" {
  SEXP mkUIndx(SEXP n_r, SEXP m_r, SEXP nnIndx_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP nnIndxLU_r, SEXP searchType_r);
}

///////////////////////////////////////////////////////////////////
//Brute force 
///////////////////////////////////////////////////////////////////
extern "C" {
  SEXP mkNNIndx(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP nnIndx_r, SEXP nnDist_r, SEXP nnIndxLU_r, SEXP nThreads_r);
}

///////////////////////////////////////////////////////////////////
//code book
///////////////////////////////////////////////////////////////////
double dmi(double *x, double *c, int inc);

double dei(double *x, double *c, int inc);

void fastNN(int m, int n, double *coords, int ui, double *u, int *sIndx, int *rSIndx, double *rSNNDist);

extern "C" {
  SEXP mkNNIndxCB(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP nnIndx_r, SEXP nnDist_r, SEXP nnIndxLU_r, SEXP nThreads_r);
}
