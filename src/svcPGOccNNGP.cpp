#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

void updateBFSVC(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';

  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuUnifb));
  int threadID = 0;
  double e;
  int mm = m*m;

#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
    for(i = 0; i < n; i++){
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
      if(i > 0){
	for(k = 0; k < nnIndxLU[n+i]; k++){
	  e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
	  c[m*threadID+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
	    C[mm*threadID+l*nnIndxLU[n+i]+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  }
	}
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){Rf_error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){Rf_error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
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
		    SEXP fixedParams_r, SEXP sigmaSqIG_r, SEXP gridIndx_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, l, ll, ii, k, s, r, q, info, nProtect=0;
    int status = 0; // For AMCMC. 
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";

    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *Xp = REAL(Xp_r);
    // Order: covariate, site
    double *Xw = REAL(Xw_r);
    int *XpRE = INTEGER(XpRE_r); 
    int *XRE = INTEGER(XRE_r); 
    int m = INTEGER(m_r)[0]; 
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1]; 
    int pOcc = INTEGER(consts_r)[2];
    int pOccRE = INTEGER(consts_r)[3];
    int nOccRE = INTEGER(consts_r)[4];
    int pDet = INTEGER(consts_r)[5];
    int pDetRE = INTEGER(consts_r)[6];
    int nDetRE = INTEGER(consts_r)[7];
    int pTilde = INTEGER(consts_r)[8];
    int Jw = INTEGER(consts_r)[9];
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    int ppTilde = pTilde * pTilde;
    int JwpTilde = Jw * pTilde;
    // Priors
    double *muBeta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(muBeta_r), &inc, muBeta, &inc);
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *SigmaBetaInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *SigmaAlphaInv = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlpha_r), &inc, SigmaAlphaInv, &inc);
    double *phiA = REAL(phiA_r);
    double *phiB = REAL(phiB_r); 
    double *nuA = REAL(nuA_r); 
    double *nuB = REAL(nuB_r); 
    double *sigmaSqA = REAL(sigmaSqA_r); 
    double *sigmaSqB = REAL(sigmaSqB_r); 
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double *tuning = REAL(tuning_r); 
    double *coords = REAL(coords_r);
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    double *K = REAL(K_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    double acceptRate = REAL(acceptRate_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int *fixedParams = INTEGER(fixedParams_r);
    int sigmaSqIG = INTEGER(sigmaSqIG_r)[0];
    int thinIndx = 0; 
    int sPost = 0; 
    int *gridIndx = INTEGER(gridIndx_r);

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    /**********************************************************************
     * Print Information 
     * *******************************************************************/
    if(verbose){
      if (currChain == 1) {
        Rprintf("----------------------------------------\n");
        Rprintf("\tModel description\n");
        Rprintf("----------------------------------------\n");
        Rprintf("NNGP Occupancy model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
	Rprintf("Number of spatially-varying coefficients: %i \n", pTilde);
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
        Rprintf("Using %i nearest neighbors.\n\n", m);
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
        Rprintf("Adaptive Metropolis with target acceptance rate: %.1f\n", 100*acceptRate);
      }
      Rprintf("----------------------------------------\n");
      Rprintf("\tChain %i\n", currChain);
      Rprintf("----------------------------------------\n");
      Rprintf("Sampling ... \n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&nOccRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Spatial processes
    double *w = (double *) R_alloc(JwpTilde, sizeof(double));   
    F77_NAME(dcopy)(&JwpTilde, REAL(wStarting_r), &inc, w, &inc);
    // Spatial variance
    double *sigmaSq = (double *) R_alloc(pTilde, sizeof(double)); 
    F77_NAME(dcopy)(&pTilde, REAL(sigmaSqStarting_r), &inc, sigmaSq, &inc); 
    // Spatial range parameter
    double *phi = (double *) R_alloc(pTilde, sizeof(double)); 
    F77_NAME(dcopy)(&pTilde, REAL(phiStarting_r), &inc, phi, &inc); 
    // Spatial smoothing parameter for Matern
    double *nu = (double *) R_alloc(pTilde, sizeof(double)); 
    F77_NAME(dcopy)(&pTilde, REAL(nuStarting_r), &inc, nu, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(zStarting_r), &inc, z, &inc);
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(J, sizeof(double)); zeros(omegaOcc, J);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); zeros(kappaOcc, J);
    double *zStar = (double *) R_alloc(J, sizeof(double)); zeros(zStar, J);
    
    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOcc * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDet * nPost);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), J * nPost);
    SEXP wSamples_r; 
    PROTECT(wSamples_r = Rf_allocMatrix(REALSXP, JwpTilde, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), JwpTilde * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), J * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = Rf_allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = Rf_allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetRE * nPost);
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = Rf_allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPsiSamples_r), pOccRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nOccRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nOccRE * nPost);
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), J * nPost);
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int JpOcc = J * pOcc; 
    int JpOccRE = J * pOccRE; 
    int nObspDet = nObs * pDet;
    int nObspDetRE = nObs * pDetRE;
    int jj, kk;
    double tmp_0, tmp_02; 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_one = (double *) R_alloc(1, sizeof(double)); 
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = zero; 
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_pTilde = (double *) R_alloc(pTilde, sizeof(double));
    double *tmp_ppTilde = (double *) R_alloc(ppTilde, sizeof(double));
   
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *yWAIC = (double *) R_alloc(J, sizeof(double)); zeros(yWAIC, J);
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *piProdWAIC = (double *) R_alloc(J, sizeof(double)); 
    ones(piProdWAIC, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J);

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, 
                   SigmaAlphaInvMuAlpha, &inc FCONE);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(J, sizeof(double)); 
    zeros(betaStarSites, J); 
    int *betaStarLongIndx = (int *) R_alloc(JpOccRE, sizeof(int));
    // Initial sums
    for (j = 0; j < J; j++) {
      for (l = 0; l < pOccRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nOccRE);
        betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]];
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(alphaStarObs, nObs); 
    int *alphaStarLongIndx = (int *) R_alloc(nObspDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarLongIndx[l * nObs + i] = which(XpRE[l * nObs + i], alphaLevelIndx, nDetRE);
        alphaStarObs[i] += alphaStar[alphaStarLongIndx[l * nObs + i]];
      }
    }
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pOccRE, sizeof(int)); 
    for (l = 0; l < pOccRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nOccRE); 
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int)); 
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE); 
    }

    /**********************************************************************
     * Set up spatial stuff and MH stuff
     * *******************************************************************/
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }  
    int nThetapTilde = nTheta * pTilde;
    double *accept = (double *) R_alloc(nThetapTilde, sizeof(double)); zeros(accept, nThetapTilde); 
    double *theta = (double *) R_alloc(nThetapTilde, sizeof(double));
    double logPostCurr = 0.0, logPostCand = 0.0;
    double logDet;  
    double phiCand = 0.0, nuCand = 0.0, sigmaSqCand = 0.0;  
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = Rf_allocMatrix(REALSXP, nThetapTilde, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nThetapTilde * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = Rf_allocMatrix(REALSXP, nThetapTilde, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nThetapTilde * nBatch);
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = Rf_allocMatrix(REALSXP, nThetapTilde, nPost)); nProtect++; 
    zeros(REAL(thetaSamples_r), nThetapTilde * nPost);
    double b, e, aij, aa; 
    double *a = (double *) R_alloc(pTilde, sizeof(double));
    double *v = (double *) R_alloc(pTilde, sizeof(double));
    double *mu = (double *) R_alloc(pTilde, sizeof(double));
    double *var = (double *) R_alloc(ppTilde, sizeof(double)); zeros(var, ppTilde);
    double *ff = (double *) R_alloc(pTilde, sizeof(double));
    double *gg = (double *) R_alloc(pTilde, sizeof(double));
    // Initiate spatial values
    for (i = 0; i < pTilde; i++) {
      theta[sigmaSqIndx * pTilde + i] = sigmaSq[i]; 
      theta[phiIndx * pTilde + i] = phi[i]; 
      if (corName == "matern") {
        theta[nuIndx * pTilde + i] = nu[i]; 
      } 
    } // i
    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(Jw-m-1)*m);

    // For NNGP
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * pTilde, sizeof(double));
    double *F = (double *) R_alloc(Jw * pTilde, sizeof(double));
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(Jw, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads * pTilde, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads * pTilde, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(pTilde*sizeBK, sizeof(double));

    // Initiate B and F for each SVC
    for (i = 0; i < pTilde; i++) {
    updateBFSVC(&B[i * nIndx], &F[i*Jw], &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, Jw, m, theta[sigmaSqIndx * pTilde + i], theta[phiIndx * pTilde + i], nu[i], covModel, &bk[i * sizeBK], nuB[i]);
    }
    // Spatial process sums for each site
    double *wSites = (double *) R_alloc(J, sizeof(double));
    // For each location, multiply w x Xw
    for (j = 0; j < J; j++) {
      wSites[j] = 0.0;
      for (ll = 0; ll < pTilde; ll++) {
        wSites[j] += w[gridIndx[j] * pTilde + ll] * Xw[ll * J + j];
      }
    }

    GetRNGstate();
  
    /**********************************************************************
     * Begin Sampler 
     * *******************************************************************/
    for (s = 0, q = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, q++) {
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
        for (j = 0; j < J; j++) {
          omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + wSites[j] + betaStarSites[j]);
        } // j
        /********************************************************************
         *Update Detection Auxiliary Variables 
         *******************************************************************/
        // Note that all of the variables are sampled, but only those at 
        // locations with z[j] == 1 actually effect the results. 
        if (nObs == J) {
          for (i = 0; i < nObs; i++) {
            if (z[zLongIndx[i]] == 1.0) {
              omegaDet[i] = rpg(K[i], F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i]);
	    }
          } // i
        } else {
          for (i = 0; i < nObs; i++) {
            if (z[zLongIndx[i]] == 1.0) {
              omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i]);
	    }
          } // i
        }
             
        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          kappaOcc[j] = z[j] - 1.0 / 2.0; 
          tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * (wSites[j] + betaStarSites[j]); 
	  zStar[j] = kappaOcc[j] / omegaOcc[j];
        } // j
        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
        for (j = 0; j < pOcc; j++) {
          tmp_pOcc[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
        for(j = 0; j < J; j++){
          for(i = 0; i < pOcc; i++){
            tmp_JpOcc[i*J+j] = X[i*J+j]*omegaOcc[j];
          }
        }

        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
        for (j = 0; j < ppOcc; j++) {
          tmp_ppOcc[j] += SigmaBetaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri here failed\n");}
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
	if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
        
        /********************************************************************
         *Update Detection Regression Coefficients
         *******************************************************************/
        // /********************************
        //  * Compute b.alpha
        //  *******************************/
        // First multiply kappDet * the current occupied values, such that values go 
        // to 0 if they z == 0 and values go to kappaDet if z == 1
        if (nObs == J) {
          for (i = 0; i < nObs; i++) {
            kappaDet[i] = (y[i] - K[i]/2.0) * z[zLongIndx[i]];
            tmp_nObs[i] = kappaDet[i] - omegaDet[i] * alphaStarObs[i]; 
            tmp_nObs[i] *= z[zLongIndx[i]]; 
          } // i
        } else {
          for (i = 0; i < nObs; i++) {
            kappaDet[i] = (y[i] - 1.0/2.0) * z[zLongIndx[i]];
            tmp_nObs[i] = kappaDet[i] - omegaDet[i] * alphaStarObs[i]; 
            tmp_nObs[i] *= z[zLongIndx[i]]; 
          } // i
        }
        
        F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
        for (j = 0; j < pDet; j++) {
          tmp_pDet[j] += SigmaAlphaInvMuAlpha[j]; 
        } // j

        /********************************
         * Compute A.alpha
         * *****************************/
        for (j = 0; j < nObs; j++) {
          for (i = 0; i < pDet; i++) {
            tmp_nObspDet[i*nObs + j] = Xp[i * nObs + j] * omegaDet[j] * z[zLongIndx[j]];
          } // i
        } // j

        F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

        for (j = 0; j < ppDet; j++) {
          tmp_ppDet[j] += SigmaAlphaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf A.alpha failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri A.alpha failed\n");}
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        mvrnorm(alpha, tmp_pDet2, tmp_ppDet, pDet);

        /********************************************************************
         *Update Occupancy random effects variance
         *******************************************************************/
        for (l = 0; l < pOccRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nOccRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqPsi[l] = rigamma(sigmaSqPsiA[l] + nOccRELong[l] / 2.0, sigmaSqPsiB[l] + tmp_0); 
        }

        /********************************************************************
         *Update Detection random effects variance
         *******************************************************************/
        for (l = 0; l < pDetRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nDetRELong[l], &alphaStar[alphaStarStart[l]], &inc, &alphaStar[alphaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] / 2.0, sigmaSqPB[l] + tmp_0); 
        }

        /********************************************************************
         *Update Occupancy random effects
         *******************************************************************/
        if (pOccRE > 0) {
          // Update each individual random effect one by one. 
          for (l = 0; l < nOccRE; l++) {
            /********************************
             * Compute b.beta.star
             *******************************/
            zeros(tmp_one, inc);
            tmp_0 = 0.0;	      
            // Only allow information to come from when XRE == betaLevelIndx[l]. 
            // aka information only comes from the sites with any given level 
            // of a random effect. 
            for (j = 0; j < J; j++) {
              if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l]) {
                tmp_02 = 0.0;
                for (ll = 0; ll < pOccRE; ll++) {
                  tmp_02 += betaStar[betaStarLongIndx[ll * J + j]];
	        } 
                tmp_one[0] += kappaOcc[j] - (F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + 
          		    tmp_02 - betaStar[l] + wSites[j]) * omegaOcc[j];
                tmp_0 += omegaOcc[j];
              }
            }
            /********************************
             * Compute A.beta.star
             *******************************/
            tmp_0 += 1.0 / sigmaSqPsi[betaStarIndx[l]]; 
            tmp_0 = 1.0 / tmp_0; 
            betaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
          }
        
          // Update the RE sums for the current species
          zeros(betaStarSites, J);
          for (j = 0; j < J; j++) {
            for (l = 0; l < pOccRE; l++) {
              betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]];
            }
          }
        }

        /********************************************************************
         *Update Detection random effects
         *******************************************************************/
        if (pDetRE > 0) {
          // Update each individual random effect one by one. 
          for (l = 0; l < nDetRE; l++) {
            /********************************
             * Compute b.alpha.star
             *******************************/
            // Only allow information to come from when z[r] == 1 and XpRE == alphaLevelIndx[l]
            zeros(tmp_one, inc);
            tmp_0 = 0.0;
            for (i = 0; i < nObs; i++) {
              if ((z[zLongIndx[i]] == 1.0) && (XpRE[alphaStarIndx[l] * nObs + i] == alphaLevelIndx[l])) {
                tmp_02 = 0.0;
                for (ll = 0; ll < pDetRE; ll++) {
                  tmp_02 += alphaStar[alphaStarLongIndx[ll * nObs + i]];
	        } 
                tmp_one[0] += kappaDet[i] - (F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + tmp_02 - alphaStar[l]) * omegaDet[i];
        	      tmp_0 += omegaDet[i];
              }
            }
            /********************************
             * Compute A.alpha.star
             *******************************/
            tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
            tmp_0 = 1.0 / tmp_0; 
            alphaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
          }
          zeros(alphaStarObs, nObs); 
          // Update the RE sums for the current species
          for (i = 0; i < nObs; i++) {
            for (l = 0; l < pDetRE; l++) {
              alphaStarObs[i] += alphaStar[alphaStarLongIndx[l * nObs + i]]; 
            }
          }
        }

        /********************************************************************
         *Update w (spatial random effects)
         *******************************************************************/
	// Update B and F for all svcs
        // for (ll = 0; ll < pTilde; ll++) {
        //   updateBFSVC(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], 
	//       	&C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, 
	//       	J, m, theta[sigmaSqIndx * pTilde + ll], 
	//       	theta[phiIndx * pTilde + ll], nu[ll], covModel, 
	//       	&bk[ll * sizeBK], nuB[ll]);
        // }

	for (ii = 0; ii < Jw; ii++) {

          for (ll = 0; ll < pTilde; ll++) { // row

            a[ll] = 0; 
	    v[ll] = 0; 

	    if (uIndxLU[Jw + ii] > 0){ // is ii a neighbor for anybody
	      for (j = 0; j < uIndxLU[Jw+ii]; j++){ // how many locations have ii as a neighbor
	        b = 0;
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[ii]+j]; // jj is the index of the jth location who has ii as a neighbor
	        for(k = 0; k < nnIndxLU[Jw+jj]; k++){ // these are the neighbors of the jjth location
	          kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	          if(kk != ii){ //if the neighbor of jj is not ii
	    	    b += B[ll*nIndx + nnIndxLU[jj]+k]*w[kk * pTilde + ll]; //covariance between jj and kk and the random effect of kk
	          }
	        } // k
	        aij = w[jj * pTilde + ll] - b;
	        a[ll] += B[ll*nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]]*aij/F[ll*Jw + jj];
	        v[ll] += pow(B[ll * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]],2)/F[ll * Jw + jj];
	      } // j
	    }
	    
	    e = 0;
	    for(j = 0; j < nnIndxLU[Jw+ii]; j++){
	      e += B[ll * nIndx + nnIndxLU[ii]+j]*w[nnIndx[nnIndxLU[ii]+j] * pTilde + ll];
	    }

	    ff[ll] = 1.0 / F[ll * Jw + ii];
	    gg[ll] = e / F[ll * Jw + ii];
	  } // ll
            
	  zeros(tmp_pTilde, pTilde);
	  zeros(tmp_ppTilde, ppTilde);
	  for (j = 0; j < J; j++) {
            if (gridIndx[j] == ii) {
              for (ll = 0; ll < pTilde; ll++) {
                // tmp_pTilde = X_tilde' %*% omega_beta  
                tmp_pTilde[ll] = Xw[ll * J + j] * omegaOcc[j];
	        // Compute tmp_pTilde %*% t(Xw)
	        for (k = 0; k < pTilde; k++) { // column
                  tmp_ppTilde[ll * pTilde + k] += tmp_pTilde[ll] * 
	                                         Xw[k * J + j];
	        } // k
	      }
	    }
	  }
	  // var
	  F77_NAME(dcopy)(&ppTilde, tmp_ppTilde, &inc, var, &inc);
	  for (k = 0; k < pTilde; k++) {
            var[k * pTilde + k] += ff[k] + v[k]; 
          } // k
	  F77_NAME(dpotrf)(lower, &pTilde, var, &pTilde, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotrf var failed\n");}
	  F77_NAME(dpotri)(lower, &pTilde, var, &pTilde, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri var failed\n");}

	  // mu
	  zeros(mu, pTilde);
	  for (j = 0; j < J; j++) {
            if (gridIndx[j] == ii) {
	      for (k = 0; k < pTilde; k++) {
                mu[k] += (zStar[j] - F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) - betaStarSites[j]) * omegaOcc[j] * Xw[k * J + j];
              } // k
	    }
	  }
	  for (k = 0; k < pTilde; k++) {
            mu[k] += gg[k] + a[k];
	  }
	  F77_NAME(dsymv)(lower, &pTilde, &one, var, &pTilde, mu, &inc, &zero, tmp_pTilde, &inc FCONE);

	  F77_NAME(dpotrf)(lower, &pTilde, var, &pTilde, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf var 2 failed\n");}

	  mvrnorm(&w[ii * pTilde], tmp_pTilde, var, pTilde);
        } // ii
	

	// Compute Xw %*% w = wSites. 
	// This calculation is correct (confirmed April 27)
	for (j = 0; j < J; j++) {
          wSites[j] = 0.0;
	  for (ll = 0; ll < pTilde; ll++) {
            wSites[j] += w[gridIndx[j] * pTilde + ll] * Xw[ll * J + j];
	  }
	}

        /********************************************************************
         *Update spatial covariance parameters
         *******************************************************************/
	for (ll = 0; ll < pTilde; ll++) {
          /******************************************************************
           *Update sigmaSq
           *****************************************************************/
          if (sigmaSqIG) {
            aa = 0;
            if (!fixedParams[3]) {
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:aa, logDet)
#endif
              for (j = 0; j < Jw; j++){
                if(nnIndxLU[Jw+j] > 0){
                  e = 0;
                  for(i = 0; i < nnIndxLU[Jw+j]; i++){
                    e += B[ll * nIndx + nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i] * pTilde + ll];
                  }
                  b = w[j * pTilde + ll] - e;
                }else{
                  b = w[j * pTilde + ll];
                }	
                aa += b*b/F[ll * Jw + j];
              }

	      theta[sigmaSqIndx * pTilde + ll] = rigamma(sigmaSqA[ll] + Jw / 2.0, sigmaSqB[ll] + 0.5 * aa * theta[sigmaSqIndx * pTilde + ll]); 
	    }
	  }
      
          /******************************************************************
           *Update phi (and nu if matern)
           *****************************************************************/
          // Current
	  if (!fixedParams[2] || !fixedParams[3]) {
            if (corName == "matern"){ 
	      nu[ll] = theta[nuIndx * pTilde + ll];
       	    }
            updateBFSVC(&B[ll * nIndx], &F[ll*Jw], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, Jw, m, theta[sigmaSqIndx * pTilde + ll], theta[phiIndx * pTilde + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[ll]);
	  }
          aa = 0;
          logDet = 0;

	  if (!fixedParams[2]) {
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
            for (j = 0; j < Jw; j++){
              if (nnIndxLU[Jw+j] > 0){
                e = 0;
                for (ii = 0; ii < nnIndxLU[Jw+j]; ii++){
                  e += B[ll * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * pTilde + ll];
                }
                b = w[j * pTilde + ll] - e;
              } else{
                b = w[j * pTilde + ll];
              }	
              aa += b*b/F[ll * Jw + j];
              logDet += log(F[ll * Jw + j]);
            }
      
            logPostCurr = -0.5 * logDet - 0.5 * aa;
            logPostCurr += log(theta[phiIndx * pTilde + ll] - phiA[ll]) + log(phiB[ll] - theta[phiIndx * pTilde + ll]); 
            if(corName == "matern"){
       	      logPostCurr += log(theta[nuIndx * pTilde + ll] - nuA[ll]) + log(nuB[ll] - theta[nuIndx * pTilde + ll]); 
            }
	    if (sigmaSqIG == 0) {
              logPostCurr += log(theta[sigmaSqIndx * pTilde + ll] - sigmaSqA[ll]) +
	                     log(sigmaSqB[ll] - theta[sigmaSqIndx * pTilde + ll]);
	    }
            
            // Candidate
            phiCand = logitInv(rnorm(logit(theta[phiIndx * pTilde + ll], phiA[ll], phiB[ll]), exp(tuning[phiIndx * pTilde + ll])), phiA[ll], phiB[ll]);
            if (corName == "matern"){
      	      nuCand = logitInv(rnorm(logit(theta[nuIndx * pTilde + ll], nuA[ll], nuB[ll]), exp(tuning[nuIndx * pTilde + ll])), nuA[ll], nuB[ll]);
            }
	  if (sigmaSqIG == 0) {
            sigmaSqCand = logitInv(rnorm(logit(theta[sigmaSqIndx * pTilde + ll], sigmaSqA[ll], sigmaSqB[ll]),
				         exp(tuning[sigmaSqIndx * pTilde + ll])), sigmaSqA[ll], sigmaSqB[ll]);
	  }
     
	  if (sigmaSqIG) { 
            updateBFSVC(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, Jw, m, theta[sigmaSqIndx * pTilde + ll], phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
	  } else {
            updateBFSVC(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, Jw, m, sigmaSqCand, phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
	  }
      
            aa = 0;
            logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
            for (j = 0; j < Jw; j++){
              if (nnIndxLU[Jw+j] > 0){
                e = 0;
                for (ii = 0; ii < nnIndxLU[Jw+j]; ii++){
                  e += BCand[nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * pTilde + ll];
                }
                b = w[j * pTilde + ll] - e;
              } else{
                b = w[j * pTilde + ll];
                }	
                aa += b*b/FCand[j];
                logDet += log(FCand[j]);
            }
            
            logPostCand = -0.5*logDet - 0.5*aa;      
            logPostCand += log(phiCand - phiA[ll]) + log(phiB[ll] - phiCand); 
            if (corName == "matern"){
              logPostCand += log(nuCand - nuA[ll]) + log(nuB[ll] - nuCand); 
            }
	    if (sigmaSqIG == 0) {
              logPostCand += log(sigmaSqCand - sigmaSqA[ll]) + log(sigmaSqB[ll] - sigmaSqCand);
	    }

            if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {

              F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[ll * nIndx], &inc);
              F77_NAME(dcopy)(&Jw, FCand, &inc, &F[ll * Jw], &inc);
              
	      theta[phiIndx * pTilde + ll] = phiCand;
              accept[phiIndx * pTilde + ll]++;
              if (corName == "matern") {
                nu[ll] = nuCand; 
	        theta[nuIndx * pTilde + ll] = nu[ll]; 
                accept[nuIndx * pTilde + ll]++; 
              }
	      if (sigmaSqIG == 0) {
                theta[sigmaSqIndx * pTilde + ll] = sigmaSqCand;
	        accept[sigmaSqIndx * pTilde + ll]++;
	      }
            }
	  }
	} // ll

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
        if (nObs == J) {
          for (i = 0; i < nObs; i++) {
            detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
            psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + wSites[zLongIndx[i]] + betaStarSites[zLongIndx[i]], zero, one); 
            piProd[zLongIndx[i]] = pow(1.0 - detProb[i], K[i]);
	    piProdWAIC[zLongIndx[i]] *= pow(detProb[i], y[i]);
	    piProdWAIC[zLongIndx[i]] *= pow(1.0 - detProb[i], K[i] - y[i]);
            ySum[zLongIndx[i]] = y[i]; 	
          } // i
        } else {
          for (i = 0; i < nObs; i++) {
            detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
            if (tmp_J[zLongIndx[i]] == 0) {
              psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + wSites[zLongIndx[i]] + betaStarSites[zLongIndx[i]], zero, one); 
            }
            piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
	    piProdWAIC[zLongIndx[i]] *= pow(detProb[i], y[i]);
	    piProdWAIC[zLongIndx[i]] *= pow(1.0 - detProb[i], 1 - y[i]);
            ySum[zLongIndx[i]] += y[i]; 	
            tmp_J[zLongIndx[i]]++;
          } // i
        }
        // Compute occupancy probability and the integrated likelihood for WAIC
        for (j = 0; j < J; j++) {
          psiNum = psi[j] * piProd[j]; 
          if (ySum[j] == zero) {
            z[j] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j])));          
            yWAIC[j] = (1.0 - psi[j]) + psi[j] * piProdWAIC[j]; 
          } else {
            z[j] = one; 
            yWAIC[j] = psi[j] * piProdWAIC[j];
          }
          piProd[j] = one;
          piProdWAIC[j] = one;
          ySum[j] = zero; 
          tmp_J[j] = 0; 
        } // j

        /********************************************************************
         *Save samples
         *******************************************************************/
	if (q >= nBurn) {
          thinIndx++; 
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
            F77_NAME(dcopy)(&JwpTilde, w, &inc, &REAL(wSamples_r)[sPost*JwpTilde], &inc); 
	    F77_NAME(dcopy)(&nThetapTilde, theta, &inc, 
			    &REAL(thetaSamples_r)[sPost*nThetapTilde], &inc); 
	    F77_NAME(dcopy)(&J, z, &inc, &REAL(zSamples_r)[sPost*J], &inc); 
            if (pOccRE > 0) {
              F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, 
                  	    &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
              F77_NAME(dcopy)(&nOccRE, betaStar, &inc, 
                  	    &REAL(betaStarSamples_r)[sPost*nOccRE], &inc);
            }
            if (pDetRE > 0) {
              F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, 
                  	    &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
              F77_NAME(dcopy)(&nDetRE, alphaStar, &inc, 
                  	    &REAL(alphaStarSamples_r)[sPost*nDetRE], &inc);
            }
            F77_NAME(dcopy)(&J, yWAIC, &inc, 
        		    &REAL(likeSamples_r)[sPost*J], &inc);
	    sPost++; 
	    thinIndx = 0; 
	  }
	}

        R_CheckUserInterrupt();
      } // r (end batch)

      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (ll = 0; ll < pTilde; ll++) {
        for (j = 0; j < nTheta; j++) {
          REAL(acceptSamples_r)[s * nThetapTilde + j * pTilde + ll] = accept[j * pTilde + ll]/batchLength; 
          REAL(tuningSamples_r)[s * nThetapTilde + j * pTilde + ll] = tuning[j * pTilde + ll]; 
          if (accept[j * pTilde + ll] / batchLength > acceptRate) {
            tuning[j * pTilde + ll] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          } else{
              tuning[j * pTilde + ll] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
            }
          accept[j * pTilde + ll] = 0;
        } // j
      } // ll
      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tCoefficient\tParameter\tAcceptance\tTuning\n");	  
          for (ll = 0; ll < pTilde; ll++) {
            Rprintf("\t%i\t\tphi\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nThetapTilde + phiIndx * pTilde + ll], exp(tuning[phiIndx * pTilde + ll]));
	    if (corName == "matern") {
              Rprintf("\t%i\t\tnu\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nThetapTilde + nuIndx * pTilde + ll], exp(tuning[nuIndx * pTilde + ll]));
	    }
	    if (sigmaSqIG == 0) {
	      Rprintf("\t%i\t\tsigmaSq\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nThetapTilde + sigmaSqIndx * pTilde + ll], exp(tuning[sigmaSqIndx * pTilde + ll]));
	    }
          } // ll
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;        
    } // s (sample loop)
    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }


    // This is necessary when generating random numbers in C.     
    PutRNGstate();

    //make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 9;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 5, wSamples_r); 
    SET_VECTOR_ELT(result_r, 6, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 7, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 8, likeSamples_r); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 9, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 10, alphaStarSamples_r);
    }
    if (pOccRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 11;
      } else {
        tmp_0 = 9;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }

    // Rf_mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 10, Rf_mkChar("alpha.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, Rf_mkChar("beta.star.samples")); 
    }
   
    // Set the names of the output list.  
    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}



