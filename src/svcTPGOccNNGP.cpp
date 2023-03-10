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

void updateBFSVCT(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

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
		     SEXP currChain_r, SEXP nChain_r, SEXP sigmaSqIG_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, r, ll, rr, ii, t, info, nProtect=0, l;
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
    // Order: covariate, site
    double *Xw = REAL(Xw_r);
    double *Xp = REAL(Xp_r);
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
    int nYearsMax = INTEGER(consts_r)[8];
    int pTilde = INTEGER(consts_r)[9];
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc;
    int ppTilde = pTilde * pTilde;
    int JpTilde = J * pTilde;
    /**********************************
     * Priors
     * *******************************/
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
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *zYearIndx = INTEGER(zYearIndx_r); 
    int *zDatIndx = INTEGER(zDatIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    /**********************************
     * AR1 Parameters 
     * *******************************/
    int ar1 = INTEGER(ar1_r)[0];
    double rhoA = REAL(ar1Vals_r)[0];
    double rhoB = REAL(ar1Vals_r)[1];
    double sigmaSqTA = REAL(ar1Vals_r)[2];
    double sigmaSqTB = REAL(ar1Vals_r)[3];
    int sigmaSqIG = INTEGER(sigmaSqIG_r)[0];
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nThin = INTEGER(nThin_r)[0]; 
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    int currChain = INTEGER(currChain_r)[0];
    int nChain = INTEGER(nChain_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int thinIndx = 0; 
    int sPost = 0; 

    // Some constants
    int JnYears = J * nYearsMax;
    int nnYears = nYearsMax * nYearsMax;

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
        Rprintf("Spatial NNGP Trend Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i years.\n\n", J, nYearsMax);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
	Rprintf("Number of spatially-varying coefficients: %i \n", pTilde);
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
        Rprintf("Using %i nearest neighbors.\n\n", m);
	if (ar1) {
          Rprintf("Using an AR(1) temporal autocorrelation matrix.\n\n");
	}
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
    // Occurrence fixed effects
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&nOccRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Detection fixed effects
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Spatial processes
    double *w = (double *) R_alloc(JpTilde, sizeof(double));
    F77_NAME(dcopy)(&JpTilde, REAL(wStarting_r), &inc, w, &inc); 
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
    double *z = (double *) R_alloc(JnYears, sizeof(double));   
    F77_NAME(dcopy)(&JnYears, REAL(zStarting_r), &inc, z, &inc);
    // PG auxiliary variables 
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(JnYears, sizeof(double)); zeros(omegaOcc, JnYears);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(JnYears, sizeof(double)); zeros(kappaOcc, JnYears);
    double *zStar = (double *) R_alloc(JnYears, sizeof(double)); zeros(zStar, JnYears);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, JnYears, nPost)); nProtect++; 
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, JpTilde, nPost)); nProtect++; 
    SEXP etaSamples_r; 
    if (ar1) {
      PROTECT(etaSamples_r = allocMatrix(REALSXP, nYearsMax, nPost)); nProtect++; 
    }
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, JnYears, nPost)); nProtect++; 
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nOccRE, nPost)); nProtect++;
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = allocMatrix(REALSXP, JnYears, nPost)); nProtect++;
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int JnYearspOccRE = J * nYearsMax * pOccRE; 
    int nObspDet = nObs * pDet;
    int nObspDetRE = nObs * pDetRE;
    int jj, kk;
    double tmp_0 = 0.0;
    double tmp_02;
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_ppOcc2 = (double *) R_alloc(ppOcc, sizeof(double));
    zeros(tmp_ppOcc2, ppOcc);
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_one = (double *) R_alloc(1, sizeof(double)); 
    int *tmp_JnYearsInt = (int *) R_alloc(JnYears, sizeof(int));
    for (j = 0; j < JnYears; j++) {
      tmp_JnYearsInt[j] = zero; 
    }
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J1, J);
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double));
    zeros(tmp_nObs, nObs);
    double *tmp_JnYears = (double *) R_alloc(JnYears, sizeof(double));
    zeros(tmp_JnYears, JnYears);
    double *tmp_JnYearspOcc = (double *) R_alloc(JnYears * pOcc, sizeof(double));
    zeros(tmp_JnYearspOcc, JnYears * pOcc);
    double *tmp_pTilde = (double *) R_alloc(pTilde, sizeof(double));
    zeros(tmp_pTilde, pTilde);
    double * tmp_ppTilde = (double *) R_alloc(ppTilde, sizeof(double));
    zeros(tmp_ppTilde, ppTilde);

    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
    double *psi = (double *) R_alloc(JnYears, sizeof(double)); 
    zeros(psi, JnYears); 
    double *piProd = (double *) R_alloc(JnYears, sizeof(double)); 
    ones(piProd, JnYears); 
    double *yWAIC = (double *) R_alloc(JnYears, sizeof(double)); ones(yWAIC, JnYears);
    double *ySum = (double *) R_alloc(JnYears, sizeof(double)); zeros(ySum, JnYears);
    double *piProdWAIC = (double *) R_alloc(JnYears, sizeof(double)); 
    ones(piProdWAIC, JnYears); 

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, 
                   SigmaAlphaInvMuAlpha, &inc FCONE);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site/year-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JnYears, sizeof(double)); 
    zeros(betaStarSites, JnYears); 
    int *betaStarLongIndx = (int *) R_alloc(JnYearspOccRE, sizeof(int));
    // Initial sums
    for (t = 0; t < nYearsMax; t++) {
      for (j = 0; j < J; j++) {
        for (l = 0; l < pOccRE; l++) {
          betaStarLongIndx[l * JnYears + t * J + j] = which(XRE[l * JnYears + t * J + j], 
                                                            betaLevelIndx, nOccRE);
          betaStarSites[t * J + j] += betaStar[betaStarLongIndx[l * JnYears + t * J + j]];
        } // l
      } // j
    } // t
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
    int nTheta, sigmaSqIndx, phiIndx, nuIndx, sigmaSqTIndx, rhoIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }  
    int nThetapTilde = nTheta * pTilde;
    int nThetaSave = nThetapTilde;
    if (ar1 == 1) {
      nThetaSave += 2;
      sigmaSqTIndx = nThetapTilde;
      rhoIndx = sigmaSqTIndx + 1;
    }
    double *accept = (double *) R_alloc(nThetaSave, sizeof(double)); zeros(accept, nThetaSave); 
    double *theta = (double *) R_alloc(nThetaSave, sizeof(double));
    double logPostCurr = 0.0, logPostCand = 0.0;
    double logDet;  
    double phiCand = 0.0, nuCand = 0.0, rhoCand = 0.0, sigmaSqCand = 0.0;  
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nThetaSave, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nThetaSave, nBatch)); nProtect++; 
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nThetaSave, nPost)); nProtect++; 
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
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * pTilde, sizeof(double));
    double *F = (double *) R_alloc(J * pTilde, sizeof(double));
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads * pTilde, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads * pTilde, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(pTilde*sizeBK, sizeof(double));

    // Initiate B and F for each SVC
    for (i = 0; i < pTilde; i++) {
    updateBFSVCT(&B[i * nIndx], &F[i*J], &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * pTilde + i], theta[phiIndx * pTilde + i], nu[i], covModel, &bk[i * sizeBK], nuB[i]);
    }
    // Spatial process sums for each site/year. 
    double *wSites = (double *) R_alloc(JnYears, sizeof(double)); zeros(wSites, JnYears);
    // For each location, multiply w x Xw
    for (t = 0; t < nYearsMax; t++) {
      wSites[t * J + j] = 0.0;
      for (j = 0; j < J; j++) {
        for (ll = 0; ll < pTilde; ll++) {
          wSites[t * J + j] += w[j * pTilde + ll] * Xw[ll * JnYears + t * J + j];
        }
      }
    }
    /**********************************************************************
     * Set up AR1 stuff
     *********************************************************************/
    double rho; 
    if (ar1) {
      theta[sigmaSqTIndx] = REAL(ar1Vals_r)[5];
      theta[rhoIndx] = REAL(ar1Vals_r)[4];
      rho = theta[rhoIndx];
    }
    double *SigmaEta = (double *) R_alloc(nnYears, sizeof(double));
    double *SigmaEtaCand = (double *) R_alloc(nnYears, sizeof(double));
    double *tmp_nYearsMax = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nYearsMax2 = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nnYears = (double *) R_alloc(nnYears, sizeof(double));
    if (ar1) {
      AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
      clearUT(SigmaEta, nYearsMax);
      F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){error("c++ error: Cholesky failed in initial time covariance matrix\n");}
      F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){error("c++ error: Cholesky inverse failed in initial time covariance matrix\n");}
    }
    double *eta = (double *) R_alloc(nYearsMax, sizeof(double)); zeros(eta, nYearsMax);
    // For sigmaSqT sampler
    double aSigmaSqTPost = 0.5 * nYearsMax + sigmaSqTA;
    double bSigmaSqTPost = 0.0;
    double *etaTRInv = (double *) R_alloc(nYearsMax, sizeof(double));

    GetRNGstate();

    /**********************************************************************
     * Begin Sampler 
     * *******************************************************************/
    for (s = 0, rr = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, rr++) {
        
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
	zeros(tmp_JnYears, JnYears);
	for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            // Only calculate omegaOcc when there are observations at that 
	    // site/year combo. Otherwise, you're just wasting time. 
	    if (zDatIndx[t * J + j] == 1) { 
              tmp_JnYears[t * J + j] = F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, 
	         	                      beta, &inc);
	      omegaOcc[t * J + j] = rpg(1.0, tmp_JnYears[t * J + j] + wSites[t * J + j] + betaStarSites[t * J + j] + eta[t]);
	      // Update kappa values along the way. 
              kappaOcc[t * J + j] = z[t * J + j] - 1.0 / 2.0; 
	      zStar[t * J + j] = kappaOcc[t * J + j] / omegaOcc[t * J + j];
	    }
	  } // j 
	} // t
	zeros(kappaDet, nObs);
        /********************************************************************
         *Update Detection Auxiliary Variables 
         *******************************************************************/
	// Only need to sample locations where z[i] == 1.0; 
	for (t = 0; t < nYearsMax; t++) {
          for (i = 0; i < nObs; i++) {
	    // If current observation is from the current year. 
            if (zYearIndx[zLongIndx[i]] == t) {
              // Reset
	      omegaDet[i] = 0.0;
	      kappaDet[i] = 0.0;
              // If the site is occupied and data were collected at that site. 
              if ((z[zLongIndx[i]] == 1.0) && (zDatIndx[zLongIndx[i]] == 1)) {
                omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i]);
	        // Update kappa values along the way. 
                kappaDet[i] = (y[i] - 1.0 / 2.0);
	      }
	    }
          } // i
	} // t

        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        zeros(tmp_JnYears, JnYears);
	for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            if (zDatIndx[t * J + j] == 1) {
              tmp_JnYears[t * J + j] = kappaOcc[t * J + j] - omegaOcc[t * J + j] * (wSites[t * J + j] + betaStarSites[t * J + j] + eta[t]); 
	    }
	  } // j
	} // t
        /********************************
         * Compute b.beta
         *******************************/
	// This is fine, because the elements in tmp_JnYears corresponding
	// to unobserve site/time locations is set to 0 and not changed. 
        F77_NAME(dgemv)(ytran, &JnYears, &pOcc, &one, X, &JnYears, 
			tmp_JnYears, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
        for (j = 0; j < pOcc; j++) {
          tmp_pOcc[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
	// This is fine, because omegaOcc == 0 for the site/year combos
	// that don't have any observations at them, which will cause this
	// whole product to go to 0.  
	for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            if (zDatIndx[t * J + j] == 1) {
              for(i = 0; i < pOcc; i++){
                tmp_JnYearspOcc[i * JnYears + t * J + j] = X[i * JnYears + t * J + j] * omegaOcc[t * J + j];
              } // i
	    }
	  } // j
	} // t

        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &JnYears, &one, X, &JnYears, tmp_JnYearspOcc, &JnYears, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
        for (j = 0; j < ppOcc; j++) {
          tmp_ppOcc[j] += SigmaBetaInv[j]; 
        } // j


        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf A.beta failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotri A.beta failed\n");}
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
	if(info != 0){error("c++ error: dpotrf A.beta2 failed\n");}
        mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);

	/********************************************************************
         *Update Detection covariates
         *******************************************************************/
        /********************************
         * Compute b.alpha
         *******************************/
	// Only calculate when z == 1 and the site/year combo has an observation. 
	zeros(tmp_nObs, nObs);
        for (i = 0; i < nObs; i++) {
          if ((z[zLongIndx[i]] == 1.0) && (zDatIndx[zLongIndx[i]] == 1)) {
            tmp_nObs[i] = kappaDet[i] - omegaDet[i] * alphaStarObs[i]; 
            tmp_nObs[i] *= z[zLongIndx[i]]; 
	  }
        } // i
        
        // This is fine, since tmp_nObs is set to 0 for situations without any data.	
        F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
        for (j = 0; j < pDet; j++) {
          tmp_pDet[j] += SigmaAlphaInvMuAlpha[j]; 
        } // j

        /********************************
         * Compute A.alpha
         * *****************************/
	// This is fine, since omegaDet is zero when there are no observed data
	// at that site/year combo. 
        for (j = 0; j < nObs; j++) {
          for (i = 0; i < pDet; i++) {
            tmp_nObspDet[i * nObs + j] = Xp[i * nObs + j] * omegaDet[j] * z[zLongIndx[j]];
          } // i
        } // j

        F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

        for (j = 0; j < ppDet; j++) {
          tmp_ppDet[j] += SigmaAlphaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
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
	    for (t = 0; t < nYearsMax; t++) {
              for (j = 0; j < J; j++) {
                if (XRE[betaStarIndx[l] * JnYears + t * J + j] == betaLevelIndx[l]) {
                   if (zDatIndx[t * J + j] == 1) {
                    tmp_02 = 0.0;
                    for (ll = 0; ll < pOccRE; ll++) {
                      tmp_02 += betaStar[betaStarLongIndx[ll * JnYears + t * J + j]];
	            } 
                    tmp_one[0] += kappaOcc[t * J + j] - (F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, beta, &inc) + 
                  	          tmp_02 - betaStar[l] + wSites[t * J + j] + eta[t]) * omegaOcc[t * J + j];
                    tmp_0 += omegaOcc[t * J + j];
		  }
                }
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
          zeros(betaStarSites, JnYears);
	  for (t = 0; t < nYearsMax; t++) {
            for (j = 0; j < J; j++) {
              for (l = 0; l < pOccRE; l++) {
                betaStarSites[t * J + j] += betaStar[betaStarLongIndx[l * JnYears + t * J + j]];
              }
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
	for (ii = 0; ii < J; ii++ ) {
          zeros(tmp_ppTilde, ppTilde);
          for (ll = 0; ll < pTilde; ll++) {
            for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + ii] == 1) {
                tmp_0 = Xw[ll * JnYears + t * J + ii] * omegaOcc[t * J + ii];
	        for (k = 0; k < pTilde; k++) {
                  tmp_ppTilde[ll * pTilde + k] += tmp_0 * Xw[k * JnYears + t * J + ii];
	        }
	      }
	    } // t
            a[ll] = 0.0;
	    v[ll] = 0.0;
	    if (uIndxLU[J + ii] > 0){ // is i a neighbor for anybody
	      for (j = 0; j < uIndxLU[J+ii]; j++){ // how many locations have ii as a neighbor
	        b = 0.0;
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[ii]+j]; // jj is the index of the jth location who has ii as a neighbor
	        for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	          kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	          if(kk != ii){ //if the neighbor of jj is not ii
	    	    b += B[ll * nIndx + nnIndxLU[jj]+k]*w[kk * pTilde + ll]; //covariance between jj and kk and the random effect of kk
	          }
	        }
	        aij = w[jj * pTilde + ll] - b;
	        a[ll] += B[ll * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]]*aij/F[ll * J + jj];
	        v[ll] += pow(B[ll * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]],2)/F[ll * J + jj];
	      }
	    }
	    e = 0.0;
	    for(j = 0; j < nnIndxLU[J+ii]; j++){
	      e += B[ll * nIndx + nnIndxLU[ii]+j]*w[nnIndx[nnIndxLU[ii]+j] * pTilde + ll];
	    }
	    ff[ll] = 1.0 / F[ll * J + ii];
	    gg[ll] = e / F[ll * J + ii];
	  } // ll (svc)

	  // var
	  F77_NAME(dcopy)(&ppTilde, tmp_ppTilde, &inc, var, &inc);
	  for (k = 0; k < pTilde; k++) {
            var[k * pTilde + k] += ff[k] + v[k]; 
          } // k
	  F77_NAME(dpotrf)(lower, &pTilde, var, &pTilde, &info FCONE);
          if(info != 0){error("c++ error: dpotrf var failed\n");}
	  F77_NAME(dpotri)(lower, &pTilde, var, &pTilde, &info FCONE);
          if(info != 0){error("c++ error: dpotri var failed\n");}

	  // mu
	  for (k = 0; k < pTilde; k++) {
            mu[k] = 0.0;
	    for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + ii] == 1) {
                mu[k] += (zStar[t * J + ii] - F77_NAME(ddot)(&pOcc, &X[t * J + ii], &JnYears, beta, &inc) - betaStarSites[t * J + ii] - eta[t]) * omegaOcc[t * J + ii] * Xw[k * JnYears + t * J + ii];
	      }
	    } // t
	    mu[k] += gg[k] + a[k];
	  } // k
	  F77_NAME(dsymv)(lower, &pTilde, &one, var, &pTilde, mu, &inc, &zero, tmp_pTilde, &inc FCONE);

	  F77_NAME(dpotrf)(lower, &pTilde, var, &pTilde, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf var 2 failed\n");}

	  mvrnorm(&w[ii * pTilde], tmp_pTilde, var, pTilde);
        } // ii (site)

	// Compute wSites. 
        for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            wSites[t * J + j] = 0.0;
            for (ll = 0; ll < pTilde; ll++) {
              wSites[t * J + j] += w[j * pTilde + ll] * Xw[ll * JnYears + t * J + j];
	    }
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
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:aa, logDet)
#endif
            for (j = 0; j < J; j++){
              if(nnIndxLU[J+j] > 0){
                e = 0;
                for(i = 0; i < nnIndxLU[J+j]; i++){
                  e += B[ll * nIndx + nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i] * pTilde + ll];
                }
                b = w[j * pTilde + ll] - e;
              }else{
                b = w[j * pTilde + ll];
              }	
              aa += b*b/F[ll * J + j];
            }

	    theta[sigmaSqIndx * pTilde + ll] = rigamma(sigmaSqA[ll] + J / 2.0, sigmaSqB[ll] + 0.5 * aa * theta[sigmaSqIndx * pTilde + ll]); 
	  }

          /******************************************************************
           *Update phi (and nu if matern)
           *****************************************************************/
          // Current
          if (corName == "matern"){ 
	    nu[ll] = theta[nuIndx * pTilde + ll];
       	  }
          updateBFSVCT(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * pTilde + ll], theta[phiIndx * pTilde + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[ll]);
          aa = 0;
          logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += B[ll * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * pTilde + ll];
              }
              b = w[j * pTilde + ll] - e;
            } else{
              b = w[j * pTilde + ll];
            }	
            aa += b*b/F[ll * J + j];
            logDet += log(F[ll * J + j]);
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
	  // Rprintf("SVC: %i, phiCand: %f\n", ll, phiCand);
	  // Rprintf("SVC: %i, tuning: %f\n", ll, exp(tuning[phiIndx * pTilde + ll]));
          if (corName == "matern"){
      	    nuCand = logitInv(rnorm(logit(theta[nuIndx * pTilde + ll], nuA[ll], nuB[ll]), exp(tuning[nuIndx * pTilde + ll])), nuA[ll], nuB[ll]);
          }
	  if (sigmaSqIG == 0) {
            sigmaSqCand = logitInv(rnorm(logit(theta[sigmaSqIndx * pTilde + ll], sigmaSqA[ll], sigmaSqB[ll]),
				         exp(tuning[sigmaSqIndx * pTilde + ll])), sigmaSqA[ll], sigmaSqB[ll]);
	  }
     
	  if (sigmaSqIG) { 
            updateBFSVCT(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * pTilde + ll], phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
	  } else {
            updateBFSVCT(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, sigmaSqCand, phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
	  }
      
          aa = 0;
          logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
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
            F77_NAME(dcopy)(&J, FCand, &inc, &F[ll * J], &inc);
            
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
	} // ll
	if (ar1) {
          /********************************************************************
           *Update sigmaSqT
           *******************************************************************/
	  // Form correlation matrix. 
          AR1(nYearsMax, theta[rhoIndx], 1.0, SigmaEta);
	  clearUT(SigmaEta, nYearsMax);
	  F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	  F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){error("c++ error: Cholesky inverse failed in covariance matrix\n");}
	  fillUTri(SigmaEta, nYearsMax);
	  // Compute t(eta) %*% SigmaEta^-1 %*% eta
          for (t = 0; t < nYearsMax; t++) {
            etaTRInv[t] = F77_NAME(ddot)(&nYearsMax, &SigmaEta[t], &nYearsMax, 
	  		               eta, &inc); 
	  }
          bSigmaSqTPost = F77_NAME(ddot)(&nYearsMax, etaTRInv, &inc, eta, &inc);	
	  bSigmaSqTPost /= 2.0;
	  bSigmaSqTPost += sigmaSqTB;
	  theta[sigmaSqTIndx] = rigamma(aSigmaSqTPost, bSigmaSqTPost);
          
          /********************************************************************
           *Update rho
           *******************************************************************/
	  rho = theta[rhoIndx];
	  rhoCand = logitInv(rnorm(logit(rho, rhoA, rhoB), exp(tuning[rhoIndx])), rhoA, rhoB); 
	  theta[rhoIndx] = rhoCand; 

	  // Construct proposal covariance matrix. 
          AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEtaCand);
	  clearUT(SigmaEtaCand, nYearsMax);

          /********************************
           * Proposal
           *******************************/
	  // Invert SigmaEtaCand and log det cov. 
          logPostCand = 0.0;
	  F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
	  if(info != 0){error("c++ error: Cholesky failed in proposal covariance matrix\n");}
	  // Get log of the determinant of the covariance matrix. 
	  for (k = 0; k < nYearsMax; k++) {
	    logPostCand += 2.0 * log(SigmaEtaCand[k*nYearsMax+k]);
	  } // k
	  F77_NAME(dpotri)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
	  if(info != 0){error("c++ error: Cholesky inverse failed in proposal covariance matrix\n");}
          logPostCand = 0.0; 
	  // Jacobian and Uniform prior. 
	  logPostCand += log(rhoCand - rhoA) + log(rhoB - rhoCand); 
	  F77_NAME(dsymv)(lower, &nYearsMax, &one,  SigmaEtaCand, &nYearsMax, eta, &inc, &zero, tmp_nYearsMax, &inc FCONE);
	  logPostCand += -0.5*logPostCand-0.5*F77_NAME(ddot)(&nYearsMax, eta, &inc, tmp_nYearsMax, &inc);
          /********************************
           * Current
           *******************************/
	  theta[rhoIndx] = rho; 
          AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
	  clearUT(SigmaEta, nYearsMax);
          logPostCurr = 0.0;
	  F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	  for (k = 0; k < nYearsMax; k++) {
	    logPostCurr += 2.0 * log(SigmaEta[k*nYearsMax+k]);
	  } // k
	  F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){error("c++ error: Cholesky inverse failed in covariance matrix\n");}
          logPostCurr = 0.0; 
	  logPostCurr += log(rho - rhoA) + log(rhoB - rho); 
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &nYearsMax, &one, SigmaEta, &nYearsMax, eta, &inc, &zero, 
	  		tmp_nYearsMax, &inc FCONE);
	  logPostCurr += -0.5*logPostCurr-0.5*F77_NAME(ddot)(&nYearsMax, eta, &inc, tmp_nYearsMax, &inc);

	  // MH Accept/Reject
	  if (runif(0.0, 1.0) <= exp(logPostCand - logPostCurr)) {
            theta[rhoIndx] = rhoCand;
            accept[rhoIndx]++;
	    F77_NAME(dcopy)(&nnYears, SigmaEtaCand, &inc, SigmaEta, &inc); 
          }
          /********************************************************************
           *Update eta 
           *******************************************************************/
          /********************************
           * Compute b.w
           *******************************/
	  zeros(tmp_nYearsMax, nYearsMax);
          for (j = 0; j < J; j++) {
            for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + j] == 1) {
                tmp_nYearsMax[t] += kappaOcc[t * J + j] - omegaOcc[t * J + j] * (F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, beta, &inc) + betaStarSites[t * J + j] + wSites[t * J + j]);
	      }
	    }
          }
          /********************************
           * Compute A.w
           *******************************/
	  // Copy inverse covariance matrix into tmp_JJ
	  F77_NAME(dcopy)(&nnYears, SigmaEta, &inc, tmp_nnYears, &inc); 
          for (j = 0; j < J; j++) {
	    for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + j] == 1) {
                tmp_nnYears[t * nYearsMax + t] += omegaOcc[t * J + j];
	      }
	    } // t
	  } // j

          // Cholesky of A.eta
          F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf on A.eta failed\n");}
	  // Inverse of A.eta
          F77_NAME(dpotri)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
          if(info != 0){error("c++ error: dpotri on A.eta failed\n");}
          // A.eta.inv %*% b.eta. Stored in tmp_
          F77_NAME(dsymv)(lower, &nYearsMax, &one, tmp_nnYears, &nYearsMax, 
	  		tmp_nYearsMax, &inc, &zero, tmp_nYearsMax2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
	  if(info != 0){error("c++ error: dpotrf on A.eta failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(eta, tmp_nYearsMax2, tmp_nnYears, nYearsMax);
	}

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
        for (i = 0; i < nObs; i++) {
          detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
          if (tmp_JnYearsInt[zLongIndx[i]] == 0) {
            psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &JnYears, 
	  			                      beta, &inc) + 
	    			         wSites[zLongIndx[i]] + betaStarSites[zLongIndx[i]] + 
					 eta[zYearIndx[zLongIndx[i]]], 
				         zero, one); 
          }
          piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
	  piProdWAIC[zLongIndx[i]] *= pow(detProb[i], y[i]);
	  piProdWAIC[zLongIndx[i]] *= pow(1.0 - detProb[i], 1 - y[i]);
          ySum[zLongIndx[i]] += y[i]; 	
          tmp_JnYearsInt[zLongIndx[i]]++;
        } // i
	for (t = 0; t < nYearsMax; t++) {
          // Compute occupancy probability 
          for (j = 0; j < J; j++) {
            psiNum = psi[t * J + j] * piProd[t * J + j]; 
	    // If the site j was sampled in year t
	    if (zDatIndx[t * J + j] == 1) {
              if (ySum[t * J + j] == zero) {
                z[t * J + j] = rbinom(one, psiNum / (psiNum + (1.0 - psi[t * J + j])));
	        yWAIC[t * J + j] = (1.0 - psi[t * J + j]) + psi[t * J + j] * piProdWAIC[t * J + j];	
              } else {
                z[t * J + j] = one; 
                yWAIC[t * J + j] = psi[t * J + j] * piProdWAIC[t * J + j];
              }
	    } else {
              psi[t * J + j] = logitInv(F77_NAME(ddot)(&pOcc, &X[t * J + j], 
					&JnYears, beta, &inc) + wSites[t * J + j] + 
			                betaStarSites[t * J + j] + eta[t], zero, one); 
              z[t * J + j] = rbinom(one, psi[t * J + j]);
	    }
            piProd[t * J + j] = one;
	    piProdWAIC[t * J + j] = one;
            ySum[t * J + j] = zero; 
            tmp_JnYearsInt[t * J + j] = 0; 
          } // j
	} // t


        /********************************************************************
         *Save samples
         *******************************************************************/
	if (rr >= nBurn) {
          thinIndx++; 
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&JnYears, psi, &inc, &REAL(psiSamples_r)[sPost*JnYears], &inc); 
            F77_NAME(dcopy)(&JpTilde, w, &inc, &REAL(wSamples_r)[sPost*JpTilde], &inc); 
	    F77_NAME(dcopy)(&nThetaSave, theta, &inc, 
			    &REAL(thetaSamples_r)[sPost*nThetaSave], &inc); 
	    F77_NAME(dcopy)(&JnYears, z, &inc, &REAL(zSamples_r)[sPost*JnYears], &inc); 
	    if (ar1) {
	      F77_NAME(dcopy)(&nYearsMax, eta, &inc, &REAL(etaSamples_r)[sPost*nYearsMax], &inc);
	    }
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
            F77_NAME(dcopy)(&JnYears, yWAIC, &inc, 
        		    &REAL(likeSamples_r)[sPost*JnYears], &inc);
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
      // Adjust tuning for AR1 parameters
      if (ar1 == 1) {
        // Only need to do rho since sigmaSqT is not an MH parameter 
        REAL(acceptSamples_r)[s * rhoIndx] = accept[rhoIndx]/batchLength; 
        REAL(tuningSamples_r)[s * rhoIndx] = tuning[j * rhoIndx]; 
        if (accept[rhoIndx] / batchLength > acceptRate) {
          tuning[rhoIndx] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        } else{
            tuning[rhoIndx] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          }
        accept[rhoIndx] = 0;
      }
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
	  if (ar1) {
	    Rprintf("\t1\t\trho\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * rhoIndx], exp(tuning[rhoIndx]));
	  }
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
 
    //make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 9;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }
    if (ar1) {
      nResultListObjs += 1;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

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
    int ar1Ind = 0;
    if (ar1) {
      if (pOccRE > 0) {
        ar1Ind = tmp_0 + 2;
      } else if (pDetRE > 0) {
        ar1Ind = 11; 
      } else {
        ar1Ind = 9;
      }
      SET_VECTOR_ELT(result_r, ar1Ind, etaSamples_r);
    }
    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 9, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 10, mkChar("alpha.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, mkChar("beta.star.samples")); 
    }
    if (ar1) {
      SET_VECTOR_ELT(resultName_r, ar1Ind, mkChar("eta.samples")); 
    }
    
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

