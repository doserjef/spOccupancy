#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP PGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pocc_r, SEXP pdet_r, 
             SEXP J_r, SEXP K_r, SEXP betaStarting_r, SEXP alphaStarting_r, 
	     SEXP zStarting_r, SEXP zLongIndx_r, SEXP muBeta_r, 
	     SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

  SEXP PGOccREDet(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XpRE_r, SEXP lambdaP_r, 
		  SEXP pocc_r, SEXP pdet_r, SEXP pDetRE_r,  
                  SEXP J_r, SEXP K_r, SEXP nDetRE_r, SEXP nDetRELong_r, 
		  SEXP betaStarting_r, SEXP alphaStarting_r, 
	          SEXP sigmaSqPStarting_r, SEXP alphaStarStarting_r, 
	          SEXP zStarting_r, SEXP zLongIndx_r, 
		  SEXP alphaStarIndx_r, SEXP muBeta_r, SEXP muAlpha_r, 
		  SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP sigmaSqPA_r, 
	          SEXP sigmaSqPB_r, SEXP nSamples_r, SEXP nThreads_r, 
		  SEXP verbose_r, SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, 
		  SEXP nPost_r);

  SEXP PGOccREOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, 
		  SEXP lambdaPsi_r,
		  SEXP pocc_r, SEXP pdet_r, SEXP pOccRE_r, 
                  SEXP J_r, SEXP K_r, SEXP nOccRE_r, 
		  SEXP nOccRELong_r,
		  SEXP betaStarting_r, SEXP alphaStarting_r, 
	          SEXP sigmaSqPsiStarting_r, 
		  SEXP betaStarStarting_r, 
	          SEXP zStarting_r, SEXP zLongIndx_r, 
		  SEXP betaStarIndx_r, SEXP muBeta_r, 
	          SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
	          SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r,
	          SEXP nSamples_r, SEXP nThreads_r, 
		  SEXP verbose_r, SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, 
		  SEXP nPost_r);

  SEXP PGOccREBoth(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
		   SEXP lambdaPsi_r, SEXP lambdaP_r, 
		   SEXP pocc_r, SEXP pdet_r, SEXP pOccRE_r, SEXP pDetRE_r,  
                   SEXP J_r, SEXP K_r, SEXP nOccRE_r, SEXP nDetRE_r, 
		   SEXP nOccRELong_r, SEXP nDetRELong_r, 
		   SEXP betaStarting_r, SEXP alphaStarting_r, 
	           SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
		   SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
	           SEXP zStarting_r, SEXP zLongIndx_r, 
		   SEXP betaStarIndx_r, SEXP alphaStarIndx_r, SEXP muBeta_r, 
	           SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
	           SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, SEXP sigmaSqPA_r, 
	           SEXP sigmaSqPB_r, SEXP nSamples_r, SEXP nThreads_r, 
		   SEXP verbose_r, SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, 
		   SEXP nPost_r);

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

  SEXP spMsPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coordsD_r, SEXP pocc_r, SEXP pdet_r, 
	         SEXP J_r, SEXP K_r, SEXP N_r, SEXP betaStarting_r, 
	         SEXP alphaStarting_r, SEXP zStarting_r, SEXP betaCommStarting_r, 
	         SEXP alphaCommStarting_r, SEXP tauBetaStarting_r, SEXP tauAlphaStarting_r, 
		 SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		 SEXP nuStarting_r, SEXP zLongIndx_r, SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	         SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauBetaA_r, 
	         SEXP tauBetaB_r, SEXP tauAlphaA_r, SEXP tauAlphaB_r, SEXP phiA_r, 
		 SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, 
		 SEXP nuA_r, SEXP nuB_r, SEXP tuning_r, 
		 SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, 
	         SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);
 
  SEXP spMsPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, SEXP pocc_r, SEXP pdet_r, 
	             SEXP J_r, SEXP K_r, SEXP N_r, SEXP m_r, SEXP nnIndx_r, 
		     SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
		     SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		     SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		     SEXP tauBetaStarting_r, SEXP tauAlphaStarting_r, 
		     SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		     SEXP nuStarting_r, SEXP zLongIndx_r, SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	             SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauBetaA_r, 
	             SEXP tauBetaB_r, SEXP tauAlphaA_r, SEXP tauAlphaB_r, SEXP phiA_r, 
		     SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		     SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		     SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

  SEXP intPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pOcc_r, SEXP pDet_r, SEXP pDetLong_r, 
	        SEXP J_r, SEXP JLong_r, SEXP K_r, SEXP nObs_r, SEXP nObsLong_r, SEXP nData_r, 
		SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
		SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP sigmaAlpha_r, 
		SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

  SEXP spIntPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coordsD_r, 
		  SEXP pOcc_r, SEXP pDet_r, SEXP pDetLong_r, 
	          SEXP J_r, SEXP JLong_r, SEXP K_r, SEXP nObs_r, SEXP nObsLong_r, SEXP nData_r, 
		  SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		  SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		  SEXP nuStarting_r, SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
		  SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP sigmaAlpha_r, 
		  SEXP phiA_r, SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, 
		  SEXP nuA_r, SEXP nuB_r, SEXP tuning_r, 
		  SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		  SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

  SEXP spIntPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, 
		      SEXP pOcc_r, SEXP pDet_r, SEXP pDetLong_r, 
	              SEXP J_r, SEXP JLong_r, SEXP K_r, SEXP nObs_r, SEXP nObsLong_r, SEXP nData_r, 
	              SEXP m_r, SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP uIndx_r, 
		      SEXP uIndxLU_r, SEXP uiIndx_r,
		      SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		      SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		      SEXP nuStarting_r, SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
		      SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP sigmaAlpha_r, 
		      SEXP phiA_r, SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, 
		      SEXP nuA_r, SEXP nuB_r, SEXP tuning_r, SEXP covModel_r, 
		      SEXP nBatch_r, SEXP batchLength_r, 
		      SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

}
