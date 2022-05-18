#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

void updateBF1MsRE(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){
    
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
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
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
		     SEXP samplesInfo_r, SEXP chainInfo_r, SEXP sigmaSqInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, l, s, g, t, q, r, info, nProtect=0, ii;    

    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    // Sorted by visit, then by site, then by species. 
    // (e.g., visit 1, site 1, sp 1, v1, s1, sp2, 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *coords = REAL(coords_r); 
    int *XpRE = INTEGER(XpRE_r); 
    int *XRE = INTEGER(XRE_r);
    int m = INTEGER(m_r)[0]; 
    double *Xp = REAL(Xp_r);
    // Load constants
    int N = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int nObs = INTEGER(consts_r)[2]; 
    int pOcc = INTEGER(consts_r)[3];
    int pOccRE = INTEGER(consts_r)[4];
    int nOccRE = INTEGER(consts_r)[5];
    int pDet = INTEGER(consts_r)[6];
    int pDetRE = INTEGER(consts_r)[7];
    int nDetRE = INTEGER(consts_r)[8];
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *muAlphaComm = REAL(muAlphaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *SigmaAlphaCommInv = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlphaComm_r), &inc, SigmaAlphaCommInv, &inc);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *tauSqAlphaA = REAL(tauSqAlphaA_r); 
    double *tauSqAlphaB = REAL(tauSqAlphaB_r); 
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
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nOccRELong = INTEGER(nOccRELong_r);
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
    int status = 0; 
    int thinIndx = 0; 
    int fixedSigmaSq = INTEGER(sigmaSqInfo_r)[0];
    int sigmaSqIG = INTEGER(sigmaSqInfo_r)[1];
    int sPost = 0; 

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
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
        Rprintf("NNGP Multispecies Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
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
    }

    /**********************************************************************
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pOccN = pOcc * N; 
    int pDetN = pDet * N; 
    int nObsN = nObs * N; 
    int nDetREN = nDetRE * N; 
    int nOccREN = nOccRE * N; 
    int JN = J * N;
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    int jj, kk;
    double tmp_0; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_beta = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_alpha = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = 0; 
    }
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Community level
    double *betaComm = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    double *alphaComm = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaCommStarting_r), &inc, alphaComm, &inc);
    double *tauSqAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(tauSqAlphaStarting_r), &inc, tauSqAlpha, &inc);
    // Species level
    double *beta = (double *) R_alloc(pOccN, sizeof(double));   
    F77_NAME(dcopy)(&pOccN, REAL(betaStarting_r), &inc, beta, &inc);
    // Occurrence random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Detection covariates
    double *alpha = (double *) R_alloc(pDetN, sizeof(double));   
    F77_NAME(dcopy)(&pDetN, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREN, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREN, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Latent occurrence random effects
    double *betaStar = (double *) R_alloc(nOccREN, sizeof(double)); 
    F77_NAME(dcopy)(&nOccREN, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Spatial random effects
    double *w = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(wStarting_r), &inc, w, &inc);
    // Spatial variance
    double *sigmaSq = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, REAL(sigmaSqStarting_r), &inc, sigmaSq, &inc); 
    // Spatial range parameter
    double *phi = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, REAL(phiStarting_r), &inc, phi, &inc); 
    // Spatial smoothing parameter for Matern
    double *nu = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, REAL(nuStarting_r), &inc, nu, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
    // Only need to set aside J locations since you don't save these 
    // for each species
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(J, sizeof(double)); zeros(omegaOcc, J);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); zeros(kappaOcc, J);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++; 
    SEXP tauSqAlphaSamples_r; 
    PROTECT(tauSqAlphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++; 
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOccN, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDetN, nPost)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    // Spatial random effects
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    SEXP sigmaSqPSamples_r;
    SEXP alphaStarSamples_r;
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetREN, nPost)); nProtect++;
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nOccREN, nPost)); nProtect++;
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++;

    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
   
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObsN, sizeof(double)); 
    double *yWAIC = (double *) R_alloc(JN, sizeof(double)); zeros(yWAIC, JN);
    double *psi = (double *) R_alloc(JN, sizeof(double)); 
    zeros(psi, JN); 
    double *piProd = (double *) R_alloc(JN, sizeof(double)); 
    ones(piProd, JN); 
    double *piProdWAIC = (double *) R_alloc(JN, sizeof(double)); 
    ones(piProdWAIC, JN); 
    double *ySum = (double *) R_alloc(JN, sizeof(double)); zeros(ySum, JN);

    // For normal priors
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, 
                   SigmaAlphaCommInvMuAlpha, &inc FCONE);
    // Put community level variances in a pOcc x POcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level variances in a pDet x pDet matrix. 
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauSqAlpha[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

    /**********************************************************************
     * Prep for random effects (if they exist)
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JN, sizeof(double)); 
    zeros(betaStarSites, JN); 
    // Initial sums (initiate with the first species)
    for (i = 0; i < N; i++) {
      for (j = 0; j < J; j++) {
        for (l = 0; l < pOccRE; l++) {
          betaStarSites[i * J + j] += betaStar[i * nOccRE + which(XRE[l * J + j], betaLevelIndx, nOccRE)];
        }
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObsN, sizeof(double)); 
    zeros(alphaStarObs, nObsN); 
    // Get sums of the current REs for each site/visit combo for all species
    for (i = 0; i < N; i++) {
      for (r = 0; r < nObs; r++) {
        for (l = 0; l < pDetRE; l++) {
          alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE)];
        }
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
     Set up spatial stuff
     * *******************************************************************/
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }
    int nThetaN = nTheta * N; 
    double *theta = (double *) R_alloc(nThetaN, sizeof(double));
    for (i = 0; i < N; i++) {
      theta[sigmaSqIndx * N + i] = sigmaSq[i]; 
      theta[phiIndx * N + i] = phi[i]; 
      if (corName == "matern") {
        theta[nuIndx * N + i] = nu[i]; 
      } 
    } // i
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nThetaN, nPost)); nProtect++; 
    // For NNGP
    double a, v, b, e, mu, var, aij; 

    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP. Create a copy of these for each species. Increases storage 
    // space that is needed, but reduces amount of computations. 
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * N, sizeof(double));
    double *F = (double *) R_alloc(J * N, sizeof(double));
    // Only need one of these. 
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads*N, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads*N, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(N*sizeBK, sizeof(double));

    // Initiate B and F for each species
    for (i = 0; i < N; i++) {
    updateBF1MsRE(&B[i * nIndx], &F[i*J], &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * N + i], theta[phiIndx * N + i], nu[i], covModel, &bk[i * sizeBK], nuB[0]);
    }

    /**********************************************************************
     Set up stuff for Adaptive MH and other misc
     * *******************************************************************/
    double logPostCurr = 0.0, logPostCand = 0.0; 
    logPostCurr = R_NegInf; 
    double *accept = (double *) R_alloc(nThetaN, sizeof(double)); zeros(accept, nThetaN); 
    double phiCand = 0.0, nuCand = 0.0, sigmaSqCand = 0.0; 
    double logDet; 
    // MCMC info if desired
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 

    GetRNGstate();

    /**********************************************************************
     Start sampling
     *********************************************************************/
    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {

        /********************************************************************
         Update Community level Occupancy Coefficients
         *******************************************************************/
        /********************************
         Compute b.beta.comm
         *******************************/
        zeros(tmp_pOcc, pOcc); 
        for (i = 0; i < N; i++) {
          F77_NAME(dgemv)(ytran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, &beta[i], &N, &one, tmp_pOcc, &inc FCONE); 
        } // i
        for (q = 0; q < pOcc; q++) {
          tmp_pOcc[q] += SigmaBetaCommInvMuBeta[q];  
        } // j

        /********************************
         Compute A.beta.comm
         *******************************/
        for (q = 0; q < ppOcc; q++) {
          tmp_ppOcc[q] = SigmaBetaCommInv[q] + N * TauBetaInv[q]; 
        }
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotri ABetaComm failed\n");}
        // A.beta.inv %*% b.beta
        // 1 * tmp_ppOcc * tmp_pOcc + 0 * tmp_pOcc2  = tmp_pOcc2
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        // Computes cholesky of tmp_pp again stored back in tmp_ppOcc. This chol(A.beta.inv)
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
        // Args: destination, mu, cholesky of the inverse covariance matrix, dimension
        mvrnorm(betaComm, tmp_pOcc2, tmp_ppOcc, pOcc);
        /********************************************************************
         Update Community level Detection Coefficients
         *******************************************************************/
        /********************************
         * Compute b.alpha.comm
         *******************************/
         zeros(tmp_pDet, pDet); 
         for (i = 0; i < N; i++) {
           F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, TauAlphaInv, &pDet, &alpha[i], &N, &one, tmp_pDet, &inc FCONE); 
         } // i
         for (q = 0; q < pDet; q++) {
           tmp_pDet[q] += SigmaAlphaCommInvMuAlpha[q];  
         } // j
        /********************************
         * Compute A.alpha.comm
         *******************************/
        for (q = 0; q < ppDet; q++) {
          tmp_ppDet[q] = SigmaAlphaCommInv[q] + N * TauAlphaInv[q]; 
        }
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotri AAlphaComm failed\n");}
        // A.alpha.inv %*% b.alpha
        // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2  = tmp_pDet2
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
        // Computes cholesky of tmp_pp again stored back in tmp_ppDet. This chol(A.alpha.inv)
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        // Args: destination, mu, cholesky of the inverse covariance matrix, dimension
        mvrnorm(alphaComm, tmp_pDet2, tmp_ppDet, pDet);

        /********************************************************************
         Update Community Occupancy Variance Parameter
        ********************************************************************/
        for (q = 0; q < pOcc; q++) {
          tmp_0 = 0.0;  
          for (i = 0; i < N; i++) {
            tmp_0 += (beta[q * N + i] - betaComm[q]) * (beta[q * N + i] - betaComm[q]);
          } // i
          tmp_0 *= 0.5;
          tauSqBeta[q] = rigamma(tauSqBetaA[q] + N / 2.0, tauSqBetaB[q] + tmp_0); 
        } // q
        // This is correct, nothing wrong here. 
        for (q = 0; q < pOcc; q++) {
          TauBetaInv[q * pOcc + q] = tauSqBeta[q]; 
        } // i
        F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
        /********************************************************************
         Update Community Detection Variance Parameter
        ********************************************************************/
        for (q = 0; q < pDet; q++) {
          tmp_0 = 0.0;  
          for (i = 0; i < N; i++) {
            tmp_0 += (alpha[q * N + i] - alphaComm[q]) * (alpha[q * N + i] - alphaComm[q]);
          } // i
          tmp_0 *= 0.5;
          tauSqAlpha[q] = rigamma(tauSqAlphaA[q] + N / 2.0, tauSqAlphaB[q] + tmp_0); 
        } // q
        for (q = 0; q < pDet; q++) {
          TauAlphaInv[q * pDet + q] = tauSqAlpha[q]; 
        } // i
        F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
        F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

        /********************************************************************
         *Update Occupancy random effects variance
         *******************************************************************/
        for (l = 0; l < pOccRE; l++) {
          tmp_0 = 0.0; 
          for (i = 0; i < N; i++) {
            tmp_0 += F77_NAME(ddot)(&nOccRELong[l], &betaStar[i*nOccRE + betaStarStart[l]], &inc, &betaStar[i*nOccRE + betaStarStart[l]], &inc); 
          }
          tmp_0 *= 0.5; 
          sigmaSqPsi[l] = rigamma(sigmaSqPsiA[l] + nOccRELong[l] * N / 2.0, sigmaSqPsiB[l] + tmp_0);
        }

        /********************************************************************
         *Update Detection random effects variance
         *******************************************************************/
        for (l = 0; l < pDetRE; l++) {
          tmp_0 = 0.0; 
          for (i = 0; i < N; i++) {
            tmp_0 += F77_NAME(ddot)(&nDetRELong[l], &alphaStar[i*nDetRE + alphaStarStart[l]], &inc, &alphaStar[i*nDetRE + alphaStarStart[l]], &inc); 
          }
          tmp_0 *= 0.5; 
          sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] * N / 2.0, sigmaSqPB[l] + tmp_0);
        }

        /********************************************************************
         *Update Species-Specific Regression Parameters
         *******************************************************************/
        for (i = 0; i < N; i++) {  
          /********************************************************************
           *Update Occupancy Auxiliary Variables 
           *******************************************************************/
          for (j = 0; j < J; j++) {
            omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + w[j * N + i] + betaStarSites[i * J + j]);
          } // j
          /********************************************************************
           *Update Detection Auxiliary Variables 
           *******************************************************************/
          // Note that all of the variables are sampled, but only those at 
          // locations with z[j] == 1 actually effect the results. 
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              if (z[zLongIndx[r] * N + i] == 1.0) {
                omegaDet[r] = rpg(K[r], F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
	      }
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              if (z[zLongIndx[r] * N + i] == 1.0) {
                omegaDet[r] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
	      }
            } // r
          }
             
          /********************************************************************
           *Update Occupancy Regression Coefficients
           *******************************************************************/
          for (j = 0; j < J; j++) {
            kappaOcc[j] = z[j * N + i] - 1.0 / 2.0; 
	    tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * (w[j * N + i] + betaStarSites[i * J + j]); 
          } // j
          /********************************
           * Compute b.beta
           *******************************/
          // t(X) * tmp_J1 + 0 * tmp_pOcc = tmp_pOcc. 
          // dgemv computes linear combinations of different variables. 
          F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
          // TauBetaInv %*% betaComm + tmp_pOcc = tmp_pOcc
          F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, betaComm, &inc, &one, tmp_pOcc, &inc FCONE); 

          /********************************
           * Compute A.beta
           * *****************************/
          // t(X) %*% diag(omegaOcc)
          for(j = 0; j < J; j++){
            for(q = 0; q < pOcc; q++){
              tmp_JpOcc[q*J+j] = X[q*J+j]*omegaOcc[j];
            }
          }
          // This finishes off A.beta
          // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
          F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
          for (q = 0; q < ppOcc; q++) {
            tmp_ppOcc[q] += TauBetaInv[q]; 
          } // j
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
	  if(info != 0){error("c++ error: dpotrf here failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_beta, tmp_pOcc2, tmp_ppOcc, pOcc);
          // Can eventually get rid of this and change order of beta. 
          for (q = 0; q < pOcc; q++) {
            beta[q * N + i] = tmp_beta[q]; 
          }
        
          /********************************************************************
           *Update Detection Regression Coefficients
           *******************************************************************/
          /********************************
           * Compute b.alpha
           *******************************/
          // First multiply kappDet * the current occupied values, such that values go 
          // to 0 if they z == 0 and values go to kappaDet if z == 1
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              kappaDet[r] = (y[r * N + i] - K[r]/2.0) * z[zLongIndx[r] * N + i];
              tmp_nObs[r] = kappaDet[r] - omegaDet[r] * alphaStarObs[i * nObs + r]; 
              tmp_nObs[r] *= z[zLongIndx[r] * N + i]; 
            } // r
          } else { 
            for (r = 0; r < nObs; r++) {
              kappaDet[r] = (y[r * N + i] - 1.0/2.0) * z[zLongIndx[r] * N + i];
              tmp_nObs[r] = kappaDet[r] - omegaDet[r] * alphaStarObs[i * nObs + r]; 
              tmp_nObs[r] *= z[zLongIndx[r] * N + i]; 
            } // r
          }
          
          F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE);  
          F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, TauAlphaInv, &pDet, alphaComm, &inc, &one, tmp_pDet, &inc FCONE); 

          /********************************
           * Compute A.alpha
           * *****************************/
          for (r = 0; r < nObs; r++) {
            for (q = 0; q < pDet; q++) {
              tmp_nObspDet[q*nObs + r] = Xp[q * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i];
            } // i
          } // j

          // This finishes off A.alpha
          // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
          F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

          for (q = 0; q < ppDet; q++) {
            tmp_ppDet[q] += TauAlphaInv[q]; 
          } // q
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
          F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
          // A.alpha.inv %*% b.alpha
          // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2 
          // (which is currently nothing) = tmp_pDet2
          F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
          // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf here failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDet);
          for (q = 0; q < pDet; q++) {
            alpha[q * N + i] = tmp_alpha[q];
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
                  tmp_one[0] += kappaOcc[j] - (F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + betaStarSites[i * J + j] - betaStar[i * nOccRE + l] + w[j * N + i]) * omegaOcc[j];
	          tmp_0 += omegaOcc[j];
	        }
              }
              /********************************
               * Compute A.beta.star
               *******************************/
              tmp_0 += 1.0 / sigmaSqPsi[betaStarIndx[l]]; 
              tmp_0 = 1.0 / tmp_0; 
              betaStar[i * nOccRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
            }

            // Update the RE sums for the current species
            zeros(&betaStarSites[i * J], J);
            for (j = 0; j < J; j++) {
              for (l = 0; l < pOccRE; l++) {
                betaStarSites[i * J + j] += betaStar[i * nOccRE + which(XRE[l * J + j], betaLevelIndx, nOccRE)];
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
              for (r = 0; r < nObs; r++) {
                if ((z[zLongIndx[r] * N + i] == 1.0) && (XpRE[alphaStarIndx[l] * nObs + r] == alphaLevelIndx[l])) {
                  tmp_one[0] += kappaDet[r] - (F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r] - alphaStar[i * nDetRE + l]) * omegaDet[r];
	  	  tmp_0 += omegaDet[r];
	        }
	      }
              /********************************
               * Compute A.alpha.star
               *******************************/
              tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
              tmp_0 = 1.0 / tmp_0; 
              alphaStar[i * nDetRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
            }
            zeros(&alphaStarObs[i * nObs], nObs); 
            // Update the RE sums for the current species
            for (r = 0; r < nObs; r++) {
              for (l = 0; l < pDetRE; l++) {
                alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE)]; 
              }
            }
          }

          /********************************************************************
           *Update w (spatial random effects)
           *******************************************************************/
	  for (ii = 0; ii < J; ii++ ) {
            a = 0;
	    v = 0;
	    if (uIndxLU[J + ii] > 0){ // is ii a neighbor for anybody
	      for (j = 0; j < uIndxLU[J+ii]; j++){ // how many locations have ii as a neighbor
	        b = 0;
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[ii]+j]; // jj is the index of the jth location who has ii as a neighbor
	        for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	          kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	          if(kk != ii){ //if the neighbor of jj is not ii
	    	    b += B[i*nIndx + nnIndxLU[jj]+k]*w[kk * N + i]; //covariance between jj and kk and the random effect of kk
	          }
	        }
	        aij = w[jj * N + i] - b;
	        a += B[i*nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]]*aij/F[i*J + jj];
	        v += pow(B[i * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]],2)/F[i * J + jj];
	      }
	    }
	    
	    e = 0;
	    for(j = 0; j < nnIndxLU[J+ii]; j++){
	      e += B[i * nIndx + nnIndxLU[ii]+j]*w[nnIndx[nnIndxLU[ii]+j] * N + i];
	    }
	    
	    mu = (kappaOcc[ii] / omegaOcc[ii] - F77_NAME(ddot)(&pOcc, &X[ii], &J, &beta[i], &N)- betaStarSites[i * J + ii])*omegaOcc[ii] + e/F[i*J + ii] + a;
	    
	    var = 1.0/(omegaOcc[ii] + 1.0/F[i * J + ii] + v);
	    
	    w[ii * N + i] = rnorm(mu*var, sqrt(var));

          } // ii

          /********************************************************************
           *Update sigmaSq
           *******************************************************************/
          if (!fixedSigmaSq) {
            if (sigmaSqIG) {
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:a, logDet)
#endif

              for (j = 0; j < J; j++){
                if(nnIndxLU[J+j] > 0){
                  e = 0;
                  for(ii = 0; ii < nnIndxLU[J+j]; ii++){
                    e += B[i * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * N + i];
                  }
                  b = w[j * N + i] - e;
                }else{
                  b = w[j * N + i];
                }	
                a += b*b/F[i * J + j];
              }

	      theta[sigmaSqIndx * N + i] = rigamma(sigmaSqA[i] + J / 2.0, sigmaSqB[i] + 0.5 * a * theta[sigmaSqIndx * N + i]); 
	    }
	  }

          /********************************************************************
           *Update phi (and nu if matern)
           *******************************************************************/
          // Current
          if (corName == "matern"){ 
	    nu[i] = theta[nuIndx * N + i];
       	  }
          updateBF1MsRE(&B[i * nIndx], &F[i*J], &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * N + i], theta[phiIndx * N + i], nu[i], covModel, &bk[i * sizeBK], nuB[i]);
          a = 0;
          logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:a, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += B[i * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * N + i];
              }
              b = w[j * N + i] - e;
            } else{
              b = w[j * N + i];
            }	
            a += b*b/F[i * J + j];
            logDet += log(F[i * J + j]);
          }
      
          logPostCurr = -0.5 * logDet - 0.5 * a;
          logPostCurr += log(theta[phiIndx * N + i] - phiA[i]) + log(phiB[i] - theta[phiIndx * N + i]); 
          if(corName == "matern"){
       	    logPostCurr += log(theta[nuIndx * N + i] - nuA[i]) + log(nuB[i] - theta[nuIndx * N + i]); 
          }
	  if (sigmaSqIG == 0) {
            logPostCurr += log(theta[sigmaSqIndx * N + i] - sigmaSqA[i]) + 
		           log(sigmaSqB[i] - theta[sigmaSqIndx * N + i]);
	  }
          
          // Candidate
          phiCand = logitInv(rnorm(logit(theta[phiIndx * N + i], phiA[i], phiB[i]), exp(tuning[phiIndx * N + i])), phiA[i], phiB[i]);
          if (corName == "matern"){
      	    nuCand = logitInv(rnorm(logit(theta[nuIndx * N + i], nuA[i], nuB[i]), exp(tuning[nuIndx * N + i])), nuA[i], nuB[i]);
          }
	  if (sigmaSqIG == 0) {
	    sigmaSqCand = logitInv(rnorm(logit(theta[sigmaSqIndx * N + i], sigmaSqA[i], sigmaSqB[i]), 
	  			 exp(tuning[sigmaSqIndx * N + i])), sigmaSqA[i], sigmaSqB[i]); 
	  }
     
	  if (sigmaSqIG) { 
            updateBF1MsRE(BCand, FCand, &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * N + i], phiCand, nuCand, covModel, &bk[i * sizeBK], nuB[i]);
	  } else {
            updateBF1MsRE(BCand, FCand, &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, sigmaSqCand, phiCand, nuCand, covModel, &bk[i * sizeBK], nuB[i]);
	  }
      
          a = 0;
          logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:a, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += BCand[nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * N + i];
              }
              b = w[j * N + i] - e;
            } else{
              b = w[j * N + i];
              }	
              a += b*b/FCand[j];
              logDet += log(FCand[j]);
          }
          
          logPostCand = -0.5*logDet - 0.5*a;      
          logPostCand += log(phiCand - phiA[i]) + log(phiB[i] - phiCand); 
          if (corName == "matern"){
            logPostCand += log(nuCand - nuA[i]) + log(nuB[i] - nuCand); 
          }
	  if (sigmaSqIG == 0) {
            logPostCand += log(sigmaSqCand - sigmaSqA[i]) + log(sigmaSqB[i] - sigmaSqCand);
	  }

          if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {

            F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[i * nIndx], &inc);
            F77_NAME(dcopy)(&J, FCand, &inc, &F[i * J], &inc);
            
	    theta[phiIndx * N + i] = phiCand;
            accept[phiIndx * N + i]++;
            if (corName == "matern") {
              nu[i] = nuCand; 
	      theta[nuIndx * N + i] = nu[i]; 
              accept[nuIndx * N + i]++; 
            }
	    if (sigmaSqIG == 0) {
              theta[sigmaSqIndx * N + i] = sigmaSqCand;
	      accept[sigmaSqIndx * N + i]++;
	    }
          }

          /********************************************************************
           *Update Latent Occupancy
           *******************************************************************/
          // Compute detection probability 
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
              psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + w[zLongIndx[r] * N + i] + betaStarSites[i * J + zLongIndx[r]], zero, one); 
              piProd[zLongIndx[r] * N + i] = pow(1.0 - detProb[i * nObs + r], K[r]);
	      piProdWAIC[zLongIndx[r] * N + i] *= pow(detProb[i * nObs + r], y[r * N + i]);
	      piProdWAIC[zLongIndx[r] * N + i] *= pow(1.0 - detProb[i * nObs + r], K[r] - y[r * N + i]);
              ySum[zLongIndx[r] * N + i] = y[r * N + i]; 
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
              if (tmp_J[zLongIndx[r]] == 0) {
                psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + w[zLongIndx[r] * N + i] + betaStarSites[i * J + zLongIndx[r]], zero, one); 
              }
              piProd[zLongIndx[r] * N + i] *= (1.0 - detProb[i * nObs + r]);
	      piProdWAIC[zLongIndx[r] * N + i] *= pow(detProb[i * nObs + r], y[r * N + i]);
	      piProdWAIC[zLongIndx[r] * N + i] *= pow(1.0 - detProb[i * nObs + r], 
	          	                            1.0 - y[r * N + i]);
              ySum[zLongIndx[r] * N + i] += y[r * N + i]; 
              tmp_J[zLongIndx[r]]++;
            } // r
          }
          // Compute occupancy probability and the integrated likelihood for WAIC
          for (j = 0; j < J; j++) {
            psiNum = psi[j * N + i] * piProd[j * N + i]; 
            if (ySum[j * N + i] == zero) {
              z[j * N + i] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j * N + i])));
              yWAIC[j * N + i] = (1.0 - psi[j * N + i]) + psi[j * N + i] * piProdWAIC[j * N + i];
            } else {
              z[j * N + i] = one; 
              yWAIC[j * N + i] = psi[j * N + i] * piProdWAIC[j * N + i]; 
            }
            // Reset variables
            piProd[j * N + i] = one;
	    piProdWAIC[j * N + i] = one;
            ySum[j * N + i] = zero; 
            tmp_J[j] = 0; 
          } // j
        } // i

        /********************************************************************
         *Save samples
         *******************************************************************/
	if (g >= nBurn) {
	  thinIndx++;
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, alphaComm, &inc, &REAL(alphaCommSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&pOcc, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, tauSqAlpha, &inc, &REAL(tauSqAlphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&pOccN, beta, &inc, &REAL(betaSamples_r)[sPost*pOccN], &inc); 
            F77_NAME(dcopy)(&pDetN, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDetN], &inc); 
            F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&JN, w, &inc, &REAL(wSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&nThetaN, theta, &inc, &REAL(thetaSamples_r)[sPost*nThetaN], &inc);
	    if (pDetRE > 0) {
              F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, 
			      &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
              F77_NAME(dcopy)(&nDetREN, alphaStar, &inc, 
			      &REAL(alphaStarSamples_r)[sPost*nDetREN], &inc);
	    }
	    if (pOccRE > 0) {
              F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, 
			      &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
              F77_NAME(dcopy)(&nOccREN, betaStar, &inc, 
			      &REAL(betaStarSamples_r)[sPost*nOccREN], &inc);
	    }
            F77_NAME(dcopy)(&JN, yWAIC, &inc, &REAL(likeSamples_r)[sPost*JN], &inc); 
	    sPost++; 
	    thinIndx = 0; 
	  }
	}
        R_CheckUserInterrupt();
      } // t (end batch)

      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (i = 0; i < N; i++) {
        for (k = 0; k < nTheta; k++) {
          REAL(acceptSamples_r)[s * nThetaN + k * N + i] = accept[k * N + i]/batchLength; 
          REAL(tuningSamples_r)[s * nThetaN + k * N + i] = tuning[k * N + i]; 
	  // Rprintf("accept[k * N + i]: %f\n", accept[k * N + i]); 
          if (accept[k * N + i] / batchLength > acceptRate) {
            tuning[k * N + i] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          } else{
              tuning[k * N + i] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
            }
          accept[k * N + i] = 0.0;
        } // k
      } // i

      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
	  Rprintf("\tSpecies\t\tParameter\tAcceptance\tTuning\n");	  
	  for (i = 0; i < N; i++) {
	    Rprintf("\t%i\t\tphi\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nThetaN + phiIndx * N + i], exp(tuning[phiIndx * N + i]));
	    if (corName == "matern") {
	      Rprintf("\t%i\t\tnu\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nThetaN + nuIndx * N + i], exp(tuning[nuIndx * N + i]));
	    }
	    if (sigmaSqIG == 0) {
	      Rprintf("\t%i\t\tsigmaSq\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nThetaN + sigmaSqIndx * N + i], exp(tuning[sigmaSqIndx * N + i]));
	    }
	  } // i
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      } 
      status++;        

    } // all batches
    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }
  
    // This is necessary when generating random numbers in C.     
    PutRNGstate();

    // make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 13;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauSqAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, zSamples_r);
    SET_VECTOR_ELT(result_r, 7, psiSamples_r);
    SET_VECTOR_ELT(result_r, 8, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 9, wSamples_r); 
    SET_VECTOR_ELT(result_r, 10, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 11, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 12, likeSamples_r); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 13, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 14, alphaStarSamples_r);
    }
    if (pOccRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 15;
      } else {
        tmp_0 = 13;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }

    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 10, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 11, mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 12, mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 13, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 14, mkChar("alpha.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, mkChar("beta.star.samples")); 
    }
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


