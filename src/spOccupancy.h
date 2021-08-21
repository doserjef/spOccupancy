#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP PGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pocc_r, SEXP pdet_r, 
             SEXP J_r, SEXP K_r, SEXP betaStarting_r, SEXP alphaStarting_r, 
	     SEXP zStarting_r, SEXP zLongIndx_r, SEXP muBeta_r, 
	     SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

  SEXP spPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coordsD_r, SEXP pocc_r, SEXP pdet_r, 
	       SEXP J_r, SEXP K_r, SEXP betaStarting_r, SEXP alphaStarting_r, 
	       SEXP zStarting_r, SEXP wStarting_r, SEXP phiStarting_r, 
	       SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
	       SEXP zLongIndx_r, SEXP muBeta_r, SEXP muAlpha_r, 
	       SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP phiA_r, SEXP phiB_r, 
	       SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
	       SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	       SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	       SEXP nReport_r);

  SEXP spPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, SEXP pocc_r, SEXP pdet_r, 
	           SEXP J_r, SEXP K_r, SEXP m_r, SEXP nnIndx_r, 
		   SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		   SEXP betaStarting_r, SEXP alphaStarting_r, 
	           SEXP zStarting_r, SEXP phiStarting_r,
	           SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
	           SEXP zLongIndx_r, SEXP muBeta_r, SEXP muAlpha_r, 
	           SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP phiA_r, SEXP phiB_r, 
	           SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		   SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	           SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	           SEXP nReport_r);

  SEXP msPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pocc_r, SEXP pdet_r, 
	       SEXP J_r, SEXP K_r, SEXP N_r, SEXP betaStarting_r, 
	       SEXP alphaStarting_r, SEXP zStarting_r, SEXP betaCommStarting_r, 
	       SEXP alphaCommStarting_r, SEXP tauBetaStarting_r, SEXP tauAlphaStarting_r, 
	       SEXP zLongIndx_r, SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	       SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauBetaA_r, 
	       SEXP tauBetaB_r, SEXP tauAlphaA_r, SEXP tauAlphaB_r, 
	       SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);


}
