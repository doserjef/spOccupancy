#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP PGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pocc_r, SEXP pdet_r, 
             SEXP J_r, SEXP K_r, SEXP betaStarting_r, SEXP alphaStarting_r, 
	     SEXP zStarting_r, SEXP zLongIndx_r, SEXP muBeta_r, 
	     SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

}
