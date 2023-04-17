#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP PGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
             SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
             SEXP betaStarting_r, SEXP alphaStarting_r, 
             SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
             SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
             SEXP zStarting_r, SEXP zLongIndx_r, 
             SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	     SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
	     SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
             SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, SEXP sigmaSqPA_r, 
             SEXP sigmaSqPB_r, SEXP nSamples_r, SEXP nThreads_r, 
             SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
	     SEXP chainInfo_r);

  SEXP spPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coordsD_r, SEXP XRE_r, SEXP XpRE_r,
	       SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
	       SEXP betaStarting_r, SEXP alphaStarting_r, SEXP sigmaSqPsiStarting_r,
	       SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
	       SEXP zStarting_r, SEXP wStarting_r, SEXP phiStarting_r, 
	       SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
	       SEXP zLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	       SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, SEXP muBeta_r, SEXP muAlpha_r, 
	       SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP phiA_r, SEXP phiB_r, 
	       SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
	       SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
	       SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	       SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	       SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	       SEXP nReport_r, SEXP samplesInfo_r, SEXP chainInfo_r, SEXP fixedSigmaSq_r, 
	       SEXP sigmaSqIG_r);

  SEXP spPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r, SEXP XpRE_r,
	           SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
		   SEXP m_r, SEXP nnIndx_r, 
		   SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		   SEXP betaStarting_r, SEXP alphaStarting_r, SEXP sigmaSqPsiStarting_r,
		   SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
	           SEXP zStarting_r, SEXP wStarting_r, SEXP phiStarting_r, 
	           SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
	           SEXP zLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		   SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, SEXP muBeta_r, SEXP muAlpha_r, 
	           SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP phiA_r, SEXP phiB_r, 
	           SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		   SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		   SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	           SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	           SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	           SEXP nReport_r, SEXP samplesInfo_r, SEXP chainInfo_r, SEXP fixedParams_r, 
		   SEXP sigmaSqIG_r);

  SEXP spPGOccPredict(SEXP J_r, SEXP pOcc_r, SEXP X0_r, SEXP q_r, 
		      SEXP obsD_r, SEXP obsPredD_r, SEXP betaSamples_r, 
		      SEXP thetaSamples_r, SEXP wSamples_r, 
		      SEXP betaStarSiteSamples_r,
		      SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, 
		      SEXP verbose_r, SEXP nReport_r);

  SEXP spPGOccNNGPPredict(SEXP coords_r, SEXP J_r, 
		          SEXP pOcc_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, 
			  SEXP q_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			  SEXP thetaSamples_r, SEXP wSamples_r, 
			  SEXP betaStarSiteSamples_r, SEXP nSamples_r, 
			  SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			  SEXP nReport_r);

  SEXP msPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
	       SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
	       SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
	       SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
	       SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
	       SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
	       SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
	       SEXP zLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	       SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
	       SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	       SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	       SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, 
	       SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
	       SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	       SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
	       SEXP samplesInfo_r, SEXP chainInfo_r);

  SEXP spMsPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coordsD_r, 
		 SEXP XRE_r, SEXP XpRE_r, SEXP consts_r, SEXP K_r, 
		 SEXP nOccRELong_r, SEXP nDetRELong_r, 
		 SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		 SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		 SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		 SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		 SEXP nuStarting_r, SEXP sigmaSqPsiStarting_r, 
		 SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, 
		 SEXP alphaStarStarting_r, SEXP zLongIndx_r, 
		 SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		 SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
		 SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	         SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	         SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP phiA_r, 
		 SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		 SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r,
		 SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		 SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		 SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		 SEXP samplesInfo_r, SEXP chainInfo_r, SEXP fixedSigmaSq_r);
 
  SEXP spMsPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, 
		     SEXP XRE_r, SEXP XpRE_r, SEXP consts_r, SEXP K_r, 
		     SEXP nOccRELong_r, SEXP nDetRELong_r, SEXP m_r, SEXP nnIndx_r, 
		     SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
		     SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		     SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		     SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		     SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		     SEXP nuStarting_r, SEXP sigmaSqPsiStarting_r, 
		     SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, 
		     SEXP alphaStarStarting_r, SEXP zLongIndx_r, 
		     SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		     SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
		     SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	             SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	             SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP phiA_r, 
		     SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		     SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r,
		     SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		     SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		     SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		     SEXP samplesInfo_r, SEXP chainInfo_r, SEXP fixedSigmaSq_r);

  SEXP spMsPGOccPredict(SEXP J_r, SEXP N_r, SEXP pOcc_r, SEXP X0_r, SEXP q_r, 
		      SEXP obsD_r, SEXP obsPredD_r, SEXP betaSamples_r, 
		      SEXP thetaSamples_r, SEXP wSamples_r, 
		      SEXP betaStarSiteSamples_r,
		      SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, 
		      SEXP verbose_r, SEXP nReport_r);

  SEXP spMsPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP N_r, 
		            SEXP pOcc_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, 
			    SEXP q_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			    SEXP thetaSamples_r, SEXP wSamples_r, 
			    SEXP betaStarSiteSamples_r, SEXP nSamples_r, 
			    SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			    SEXP nReport_r);

  SEXP intPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pOcc_r, SEXP pDet_r, SEXP pDetLong_r, 
	        SEXP J_r, SEXP JLong_r, SEXP K_r, SEXP nObs_r, SEXP nObsLong_r, SEXP nData_r, 
		SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
		SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP sigmaAlpha_r, 
		SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, SEXP currChain_r, 
		SEXP nChain_r);

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
		  SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		  SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, SEXP currChain_r, SEXP nChain_r, 
		  SEXP fixedSigmaSq_r, SEXP sigmaSqIG_r);

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
		      SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
      	              SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, SEXP currChain_r, SEXP nChain_r, 
		      SEXP fixedSigmaSq_r, SEXP sigmaSqIG_r);

  SEXP lfMsPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
		 SEXP consts_r, SEXP K_r,
		 SEXP nOccRELong_r, SEXP nDetRELong_r, 
		 SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		 SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		 SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		 SEXP lambdaStarting_r, SEXP sigmaSqPsiStarting_r, 
		 SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, 
		 SEXP alphaStarStarting_r, SEXP zLongIndx_r, 
		 SEXP betaStarIndx_r, SEXP betaLevelIndx_r, SEXP alphaStarIndx_r, 
		 SEXP alphaLevelIndx_r, SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	         SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	         SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, 
		 SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r,
		 SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		 SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		 SEXP samplesInfo_r, SEXP chainInfo_r);

  SEXP sfMsPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r, SEXP XpRE_r, 
		     SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
		     SEXP m_r, SEXP nnIndx_r, 
		     SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
		     SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		     SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		     SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		     SEXP phiStarting_r, SEXP lambdaStarting_r, SEXP nuStarting_r, 
		     SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
		     SEXP betaStarStarting_r, SEXP alphaStarStarting_r, SEXP zLongIndx_r,
		     SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		     SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
		     SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	             SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	             SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP phiA_r, 
		     SEXP phiB_r, SEXP nuA_r, SEXP nuB_r, 
		     SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		     SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		     SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		     SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		     SEXP samplesInfo_r, SEXP chainInfo_r, SEXP tauSqIG_r);

  SEXP sfMsPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP N_r, SEXP q_r,
		            SEXP pOcc_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, 
			    SEXP JStr_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			    SEXP thetaSamples_r, SEXP lambdaSamples_r, 
			    SEXP wSamples_r, SEXP betaStarSiteSamples_r, 
			    SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			    SEXP nReport_r);

  SEXP lfJSDM(SEXP y_r, SEXP X_r, SEXP XRE_r, SEXP consts_r, SEXP nOccRELong_r, 
		 SEXP betaStarting_r, SEXP betaCommStarting_r, 
		 SEXP tauSqBetaStarting_r, SEXP lambdaStarting_r, 
		 SEXP sigmaSqPsiStarting_r, SEXP betaStarStarting_r, 
		 SEXP betaStarIndx_r, SEXP betaLevelIndx_r,  
		 SEXP muBetaComm_r, SEXP SigmaBetaComm_r, 
		 SEXP tauSqBetaA_r, SEXP tauSqBetaB_r, 
		 SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r,
		 SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		 SEXP samplesInfo_r, SEXP chainInfo_r);

  SEXP sfJSDMNNGP(SEXP y_r, SEXP X_r, SEXP coords_r, SEXP XRE_r, 
		  SEXP consts_r, SEXP nOccRELong_r, SEXP m_r, SEXP nnIndx_r, 
		  SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
		  SEXP betaStarting_r, SEXP betaCommStarting_r, SEXP tauSqBetaStarting_r, 
		  SEXP phiStarting_r, SEXP lambdaStarting_r, SEXP nuStarting_r, 
		  SEXP sigmaSqPsiStarting_r, SEXP betaStarStarting_r, SEXP wStarting_r,
		  SEXP betaStarIndx_r, SEXP betaLevelIndx_r, SEXP muBetaComm_r, 
	          SEXP SigmaBetaComm_r, SEXP tauSqBetaA_r, SEXP tauSqBetaB_r, SEXP phiA_r, 
		  SEXP phiB_r, SEXP nuA_r, SEXP nuB_r, 
		  SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		  SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		  SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		  SEXP samplesInfo_r, SEXP chainInfo_r, SEXP monitors_r);

  SEXP tPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, SEXP consts_r, 
	      SEXP nOccRELong_r, SEXP nDetRELong_r, 
	      SEXP betaStarting_r, SEXP alphaStarting_r, 
	      SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r,
	      SEXP betaStarStarting_r, SEXP alphaStarStarting_r, SEXP zStarting_r, 
	      SEXP zLongIndx_r, SEXP zYearIndx_r, SEXP zDatIndx_r, 
	      SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	      SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
	      SEXP muBeta_r, SEXP SigmaBeta_r, 
	      SEXP muAlpha_r, SEXP SigmaAlpha_r, 
	      SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
	      SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	      SEXP ar1_r, SEXP ar1Vals_r, SEXP tuning_r,
	      SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r,
	      SEXP nThreads_r, SEXP verbose_r, 
	      SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, 
	      SEXP currChain_r, SEXP nChain_r);

  SEXP stPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r, 
		   SEXP XpRE_r, SEXP consts_r, 
	           SEXP nOccRELong_r, SEXP nDetRELong_r, 
		   SEXP m_r, SEXP nnIndx_r, SEXP nnIndxLU_r, 
		   SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		   SEXP betaStarting_r, SEXP alphaStarting_r, 
		   SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r,
		   SEXP betaStarStarting_r, SEXP alphaStarStarting_r,
	           SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
		   SEXP wStarting_r, SEXP zStarting_r, 
	           SEXP zLongIndx_r, SEXP zYearIndx_r, SEXP zDatIndx_r, 
		   SEXP zLongSiteIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		   SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
		   SEXP muBeta_r, SEXP SigmaBeta_r, 
		   SEXP muAlpha_r, SEXP SigmaAlpha_r, 
		   SEXP phiA_r, SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r,
		   SEXP nuA_r, SEXP nuB_r, SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		   SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		   SEXP ar1_r, SEXP ar1Vals_r,
		   SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	           SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	           SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, 
		   SEXP currChain_r, SEXP nChain_r, SEXP sigmaSqIG_r);

  SEXP stPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP nYearsMax_r,
		          SEXP pOcc_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, 
			  SEXP q_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			  SEXP thetaSamples_r, SEXP wSamples_r, 
			  SEXP betaStarSiteSamples_r, SEXP etaSamples_r, SEXP nSamples_r, 
			  SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			  SEXP nReport_r);

  SEXP svcPGBinomNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP coords_r, SEXP XRE_r, 
	            SEXP consts_r, SEXP weights_r, SEXP nRELong_r, SEXP m_r, SEXP nnIndx_r, 
		    SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		    SEXP betaStarting_r, SEXP sigmaSqPsiStarting_r,
		    SEXP betaStarStarting_r,  
	            SEXP wStarting_r, SEXP phiStarting_r, 
	            SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
	            SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		    SEXP muBeta_r, SEXP SigmaBeta_r, SEXP phiA_r, SEXP phiB_r, 
	            SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		    SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
	            SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	            SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	            SEXP nReport_r, SEXP samplesInfo_r, SEXP chainInfo_r, SEXP fixedParams_r, 
		    SEXP sigmaSqIG_r);

  SEXP svcPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP pOcc_r, SEXP pTilde_r, 
		           SEXP m_r, SEXP X0_r, SEXP Xw0_r, SEXP coords0_r, 
			   SEXP weights0_r,
			   SEXP JStr_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			   SEXP thetaSamples_r, SEXP wSamples_r, 
			   SEXP betaStarSiteSamples_r, SEXP nSamples_r, 
			   SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			   SEXP nReport_r);

  SEXP svcPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r, SEXP XpRE_r,
	            SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
		    SEXP m_r, SEXP nnIndx_r, 
		    SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		    SEXP betaStarting_r, SEXP alphaStarting_r, SEXP sigmaSqPsiStarting_r,
		    SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
	            SEXP zStarting_r, SEXP wStarting_r, SEXP phiStarting_r, 
	            SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
	            SEXP zLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		    SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, SEXP muBeta_r, SEXP muAlpha_r, 
	            SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP phiA_r, SEXP phiB_r, 
	            SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		    SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		    SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	            SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	            SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	            SEXP nReport_r, SEXP samplesInfo_r, SEXP chainInfo_r, 
		    SEXP fixedParams_r, SEXP sigmaSqIG_r);

  SEXP svcTPGBinomNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP coords_r, SEXP XRE_r, 
		       SEXP consts_r, SEXP weights_r, SEXP nRELong_r, 
		       SEXP m_r, SEXP nnIndx_r, SEXP nnIndxLU_r, 
		       SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		       SEXP betaStarting_r, SEXP sigmaSqPsiStarting_r, SEXP betaStarStarting_r, 
	               SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
		       SEXP wStarting_r, SEXP zYearIndx_r, SEXP zDatIndx_r, 
		       SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		       SEXP muBeta_r, SEXP SigmaBeta_r, 
		       SEXP phiA_r, SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r,
		       SEXP nuA_r, SEXP nuB_r, SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		       SEXP ar1_r, SEXP ar1Vals_r,
		       SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	               SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	               SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, 
		       SEXP currChain_r, SEXP nChain_r, SEXP sigmaSqIG_r);

  SEXP svcTPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP nYearsMax_r,
		            SEXP pOcc_r, SEXP pTilde_r, SEXP m_r, SEXP X0_r, SEXP Xw0_r, 
			    SEXP coords0_r, SEXP weights0_r, 
			    SEXP q_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			    SEXP thetaSamples_r, SEXP wSamples_r, 
			    SEXP betaStarSiteSamples_r, SEXP etaSamples_r, SEXP nSamples_r, 
			    SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			    SEXP nReport_r);

  SEXP svcTPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r, 
		     SEXP XpRE_r, SEXP consts_r, 
	             SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
		     SEXP m_r, SEXP nnIndx_r, SEXP nnIndxLU_r, 
		     SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		     SEXP betaStarting_r, SEXP alphaStarting_r, 
		     SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r,
		     SEXP betaStarStarting_r, SEXP alphaStarStarting_r,
	             SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
		     SEXP wStarting_r, SEXP zStarting_r, 
	             SEXP zLongIndx_r, SEXP zYearIndx_r, SEXP zDatIndx_r, 
		     SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		     SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
		     SEXP muBeta_r, SEXP SigmaBeta_r, 
		     SEXP muAlpha_r, SEXP SigmaAlpha_r, 
		     SEXP phiA_r, SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r,
		     SEXP nuA_r, SEXP nuB_r, SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		     SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		     SEXP ar1_r, SEXP ar1Vals_r,
		     SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	             SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	             SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, 
		     SEXP currChain_r, SEXP nChain_r, SEXP sigmaSqIG_r);

  SEXP intMsPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, 
	       SEXP consts_r, SEXP nOccRELong_r, SEXP pDetLong_r, 
	       SEXP nObsLong_r, SEXP NLong_r,
	       SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
	       SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
	       SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
	       SEXP sigmaSqPsiStarting_r, SEXP betaStarStarting_r, 
	       SEXP zLongIndx_r, SEXP alphaCommIndx_r, SEXP dataIndx_r, 
	       SEXP alphaSpIndx_r, SEXP alphaDatIndx_r, SEXP alphaPDetIndx_r,
	       SEXP spLongIndx_r, SEXP obsLongIndx_r, SEXP spSiteIndx_r,
	       SEXP spDatLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	       SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	       SEXP SigmaBetaComm_r, SEXP sigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	       SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, 
	       SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
	       SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
	       SEXP samplesInfo_r, SEXP chainInfo_r);

  SEXP postHocLM(SEXP y_r, SEXP X_r, SEXP XRE_r, 
                 SEXP consts_r, SEXP nRELong_r, 
                 SEXP betaStarting_r, SEXP tauSqStarting_r,  
                 SEXP sigmaSqStarting_r, SEXP betaStarStarting_r, 
                 SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	         SEXP muBeta_r, SEXP SigmaBeta_r, 
		 SEXP tauSqA_r, SEXP tauSqB_r, 
                 SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nSamples_r, 
                 SEXP verbose_r, SEXP nReport_r, 
	         SEXP chainInfo_r);

  SEXP svcMsPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP N_r, SEXP q_r,
		             SEXP pOcc_r, SEXP pTilde_r, SEXP m_r, 
			     SEXP X0_r, SEXP Xw0_r, SEXP coords0_r, 
			     SEXP JStr_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			     SEXP thetaSamples_r, SEXP lambdaSamples_r, 
			     SEXP wSamples_r, SEXP betaStarSiteSamples_r, 
			     SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, 
			     SEXP verbose_r, SEXP nReport_r);

  SEXP svcMsPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP Xp_r, SEXP coords_r, 
                      SEXP XRE_r, SEXP XpRE_r, SEXP rangeInd_r, SEXP consts_r, 
		      SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
		      SEXP m_r, SEXP nnIndx_r, 
		      SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
		      SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		      SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		      SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		      SEXP phiStarting_r, SEXP lambdaStarting_r, SEXP nuStarting_r, 
		      SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
		      SEXP betaStarStarting_r, SEXP alphaStarStarting_r, SEXP zLongIndx_r,
		      SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		      SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
		      SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	              SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	              SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP phiA_r, 
		      SEXP phiB_r, SEXP nuA_r, SEXP nuB_r, 
		      SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		      SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		      SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		      SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		      SEXP samplesInfo_r, SEXP chainInfo_r);

}
