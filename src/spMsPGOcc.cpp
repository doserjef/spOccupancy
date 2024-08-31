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

extern "C" {
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
                 SEXP samplesInfo_r, SEXP chainInfo_r, SEXP sigmaSqInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, l, s, a, b, q, r, ll, info, nProtect=0;
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
    double *coordsD = REAL(coordsD_r); 
    int *XpRE = INTEGER(XpRE_r); 
    int *XRE = INTEGER(XRE_r);
    // Xp is sorted by parameter then site, then visit (parameter 1 site 1, p1 site 2, etc) 
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
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nOccRELong = INTEGER(nOccRELong_r);
    double *K = REAL(K_r); 
    double *phiAccept = (double *) R_alloc(N, sizeof(double)); 
    zeros(phiAccept, N); 
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
        Rprintf("Spatial Multi-species Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
#ifdef _OPENMP
        Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
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
    int JpOccRE = J * pOccRE; 
    int nObspDetRE = nObs * pDetRE;
    int JJ = J * J; 
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_beta = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_alpha = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    double * tmp_JJ = (double *) R_alloc(JJ, sizeof(double)); 
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
    PROTECT(betaCommSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), pOcc * nPost);
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaCommSamples_r), pDet * nPost);
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++; 
    zeros(REAL(tauSqBetaSamples_r), pOcc * nPost);
    SEXP tauSqAlphaSamples_r; 
    PROTECT(tauSqAlphaSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++; 
    zeros(REAL(tauSqAlphaSamples_r), pDet * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOccN, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOccN * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, pDetN, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDetN * nPost);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), JN * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), JN * nPost);
    // Spatial parameters
    SEXP wSamples_r; 
    PROTECT(wSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), JN * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r;
    SEXP alphaStarSamples_r;
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = Rf_allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = Rf_allocMatrix(REALSXP, nDetREN, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetREN * nPost);
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = Rf_allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPsiSamples_r), pOccRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nOccREN, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nOccREN * nPost);
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), JN * nPost);
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
   
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObsN, sizeof(double)); 
    double *psi = (double *) R_alloc(JN, sizeof(double)); 
    zeros(psi, JN); 
    double *yWAIC = (double *) R_alloc(JN, sizeof(double)); zeros(yWAIC, JN);
    double *piProd = (double *) R_alloc(JN, sizeof(double)); 
    ones(piProd, JN); 
    double *piProdWAIC = (double *) R_alloc(JN, sizeof(double)); 
    ones(piProdWAIC, JN); 
    double *ySum = (double *) R_alloc(JN, sizeof(double)); zeros(ySum, JN);

    // For normal priors
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    // dgemv computes linear combinations of different variables. 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, 
                   SigmaAlphaCommInvMuAlpha, &inc FCONE);
    // Put community level variances in a pOcc x POcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = 1.0 / tauSqBeta[i]; 
    } // i
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = 1.0 / tauSqAlpha[i]; 
    } // i

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JN, sizeof(double)); 
    zeros(betaStarSites, JN); 
    int *betaStarLongIndx = (int *) R_alloc(JpOccRE, sizeof(int));
    // Initial sums (initiate with the first species)
    for (j = 0; j < J; j++) {
      for (l = 0; l < pOccRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nOccRE);
        for (i = 0; i < N; i++) {
          betaStarSites[i * J + j] += betaStar[i * nOccRE + betaStarLongIndx[l * J + j]];
        }
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObsN, sizeof(double)); 
    zeros(alphaStarObs, nObsN); 
    int *alphaStarLongIndx = (int *) R_alloc(nObspDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo for all species
    for (r = 0; r < nObs; r++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarLongIndx[l * nObs + r] = which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE);
        for (i = 0; i < N; i++) {
          alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + alphaStarLongIndx[l * nObs + r]];
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
    double *currTheta = (double *) R_alloc(nTheta, sizeof(double)); 
    for (i = 0; i < N; i++) {
      theta[sigmaSqIndx * N + i] = sigmaSq[i]; 
      theta[phiIndx * N + i] = phi[i]; 
      if (corName == "matern") {
        theta[nuIndx * N + i] = nu[i]; 
      } 
    } // i
    // Initiate currTheta with first species
    for (q = 0; q < nTheta; q++) {
      currTheta[q] = theta[q * N];   
    }
    // Currently copying over the covariance matrix for each species. This
    // will require an additional inverse of the matrix, but only requires 
    // storage of JJ instead of JJN. 
    double *C = (double *) R_alloc(JJ, sizeof(double)); 
    double *CCand = (double *) R_alloc(JJ, sizeof(double));
    double *R = (double *) R_alloc(JJ, sizeof(double)); 
    double *tmp_JD = (double *) R_alloc(J, sizeof(double));
    double *tmp_JD2 = (double *) R_alloc(J, sizeof(double));
    // Get spatial correlation matrix for first species
    if (sigmaSqIG) {
      spCorLT(coordsD, J, currTheta, corName, R); 
    }
    // Get spatial covariance matrix 
    spCovLT(coordsD, J, currTheta, corName, C); 
    F77_NAME(dpotrf)(lower, &J, C, &J, &info FCONE); 
    if(info != 0){Rf_error("c++ error: Cholesky failed in initial covariance matrix\n");}
    F77_NAME(dpotri)(lower, &J, C, &J, &info FCONE); 
    if(info != 0){Rf_error("c++ error: Cholesky inverse failed in initial covariance matrix\n");}
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = Rf_allocMatrix(REALSXP, nThetaN, nPost)); nProtect++; 
    zeros(REAL(thetaSamples_r), nThetaN * nPost);

    /**********************************************************************
     Set up stuff for Adaptive MH and other misc
     * *******************************************************************/
    double logMHRatio, logPostCurr = 0.0, logPostCand = 0.0, detCand = 0.0, detCurr = 0.0; 
    logPostCurr = R_NegInf; 
    double *accept = (double *) R_alloc(nThetaN, sizeof(double)); zeros(accept, nThetaN); 
    double phiCand = 0.0, nuCand = 0.0, sigmaSqCand = 0.0; 
    // For sigmaSq Sampler
    double aSigmaSqPost = 0.0; 
    double bSigmaSqPost = 0.0; 
    // MCMC info if desired
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = Rf_allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nThetaN * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = Rf_allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nThetaN * nBatch);
    double *wTRInv = (double *) R_alloc(J, sizeof(double)); 

    GetRNGstate();

    for (s = 0, a = 0; s < nBatch; s++) {
      for (b = 0; b < batchLength; b++, a++) {

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
        if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri ABetaComm failed\n");}
        // A.beta.inv %*% b.beta
        // 1 * tmp_ppOcc * tmp_pOcc + 0 * tmp_pOcc2  = tmp_pOcc2
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        // Computes cholesky of tmp_pp again stored back in tmp_ppOcc. This chol(A.beta.inv)
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
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
        if(info != 0){Rf_error("c++ error: dpotrf AAlphaComm failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri AAlphaComm failed\n");}
        // A.alpha.inv %*% b.alpha
        // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2  = tmp_pDet2
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
        // Computes cholesky of tmp_pp again stored back in tmp_ppDet. This chol(A.alpha.inv)
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf AAlphaComm failed\n");}
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
        for (q = 0; q < pOcc; q++) {
          TauBetaInv[q * pOcc + q] = tauSqBeta[q]; 
        } // i
        F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}
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
        if(info != 0){Rf_error("c++ error: dpotrf TauAlphaInv failed\n");}
        F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri TauAlphaInv failed\n");}

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
          if(info != 0){Rf_error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_beta, tmp_pOcc2, tmp_ppOcc, pOcc);
          // Can eventually get rid of this and change order of beta. 
          for (q = 0; q < pOcc; q++) {
            beta[q * N + i] = tmp_beta[q]; 
            // Rprintf("beta[%i]: %f\n", q * N + i, tmp_beta[q]); 
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
            // Rprintf("omegaDet[%i]: %f\n", r, omegaDet[r]); 
            for (q = 0; q < pDet; q++) {
              tmp_nObspDet[q*nObs + r] = Xp[q * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i];
            } // i
          } // j

          // This finishes off A.alpha
          // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
          F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

          for (q = 0; q < ppDet; q++) {
            tmp_ppDet[q] += TauAlphaInv[q]; 
            // Rprintf("TauAlphaInv: %f\n", TauAlphaInv[q]); 
          } // q
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf A.alpha failed\n");}
          F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri A.alpha failed\n");}
          // A.alpha.inv %*% b.alpha
          // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2 
          // (which is currently nothing) = tmp_pDet2
          F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
          // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDet);
          for (q = 0; q < pDet; q++) {
            alpha[q * N + i] = tmp_alpha[q];
          }

	  /********************************************************************
           *Update sigmaSq
           *******************************************************************/
	  // Update the current theta parameters
	  for (q = 0; q < nTheta; q++) {
            currTheta[q] = theta[q * N + i];
          }
	  // Still need to do this, even when sigmaSq is fixed based on how you 
	  // coded this up. 
	  // Get inverse correlation matrix
          spCorLT(coordsD, J, currTheta, corName, R); 
	  fillUTri(R, J); 
	  // Save R in tmp_JJ
          F77_NAME(dcopy)(&JJ, R, &inc, tmp_JJ, &inc); 
          F77_NAME(dpotrf)(lower, &J, R, &J, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky failed in correlation matrix\n");}
          F77_NAME(dpotri)(lower, &J, R, &J, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky inverse failed in correlation matrix\n");}
	  fillUTri(R, J); 

	  if (!fixedSigmaSq) {
            if (sigmaSqIG) {
	      // t(w) %*% R^-1
	      // Definitely a better way to do this. 
	      for (j = 0; j < J; j++) {
                wTRInv[j] = F77_NAME(ddot)(&J, &R[j], &J, &w[i], &N);  
              } // j
	      // wTRInv %*% w
	      bSigmaSqPost = F77_NAME(ddot)(&J, wTRInv, &inc, &w[i], &N); 
	      bSigmaSqPost /= 2.0; 
	      bSigmaSqPost += sigmaSqB[i]; 
	      aSigmaSqPost = 0.5 * J + sigmaSqA[i]; 
	      theta[sigmaSqIndx * N + i] = rigamma(aSigmaSqPost, bSigmaSqPost); 
	      currTheta[sigmaSqIndx] = theta[sigmaSqIndx * N + i]; 
	    }
	  } 
	  // Get inverse covariance matrix from correlation matrix
	  for (j = 0; j < JJ; j++) {
            C[j] = 1.0 / currTheta[sigmaSqIndx] * R[j]; 
	    tmp_JJ[j] = currTheta[sigmaSqIndx] * tmp_JJ[j]; 
          }
	  

          /********************************************************************
           *Update phi (and nu if matern)
           *******************************************************************/
	  if (corName == "matern") {
            nu[i] = currTheta[nuIndx]; 
	    nuCand = logitInv(rnorm(logit(currTheta[nuIndx], nuA[i], nuB[i]), exp(tuning[nuIndx * N + i])), nuA[i], nuB[i]); 
          }
	  phi[i] = currTheta[phiIndx]; 
	  phiCand = logitInv(rnorm(logit(currTheta[phiIndx], phiA[i], phiB[i]), exp(tuning[phiIndx * N + i])), phiA[i], phiB[i]);
	  currTheta[phiIndx] = phiCand; 
	  if (corName == "matern") {
	    currTheta[nuIndx] = nuCand; 
          }
	  if (sigmaSqIG == 0) {
	    sigmaSq[i] = currTheta[sigmaSqIndx]; 
	    sigmaSqCand = logitInv(rnorm(logit(sigmaSq[i], sigmaSqA[i], sigmaSqB[i]), 
	  			 exp(tuning[sigmaSqIndx * N + i])), sigmaSqA[i], sigmaSqB[i]); 
	    currTheta[sigmaSqIndx] = sigmaSqCand; 
	  }

	  // Construct proposal covariance matrix (stored in CCand). 
	  spCovLT(coordsD, J, currTheta, corName, CCand); 

          /********************************
           * Proposal
           *******************************/
	  // Invert CCand and log det cov. 
          detCand = 0.0;
	  F77_NAME(dpotrf)(lower, &J, CCand, &J, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
	  // Get log of the determinant of the covariance matrix. 
	  for (k = 0; k < J; k++) {
            // Rprintf("detCand Value: %f\n", 2.0 * log(CCand[k * J + k])); 
	    detCand += 2.0 * log(CCand[k*J+k]);
	  } // k
	  F77_NAME(dpotri)(lower, &J, CCand, &J, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky inverse failed in covariance matrix\n");}
          logPostCand = 0.0; 
	  // Jacobian and Uniform prior. 
	  logPostCand += log(phiCand - phiA[i]) + log(phiB[i] - phiCand); 
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &J, &one,  CCand, &J, &w[i], &N, &zero, tmp_JD, &inc FCONE);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&J, &w[i], &N, tmp_JD, &inc);
	  // Rprintf("logPostCand: %f\n", logPostCand); 
          if (corName == "matern"){
            logPostCand += log(nuCand - nuA[i]) + log(nuB[i] - nuCand); 
          }
	  if (sigmaSqIG == 0) {
            logPostCand += log(sigmaSqCand - sigmaSqA[i]) + log(sigmaSqB[i] - sigmaSqCand);
	  }

          /********************************
           * Current
           *******************************/
	  // Get logdetCov
          detCurr = 0.0;
	  // tmp_JJ has the current covariance matrix
	  F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
	  // Get log of the determinant of the covariance matrix. 
	  for (k = 0; k < J; k++) {
	    detCurr += 2.0 * log(tmp_JJ[k*J+k]);
	  } // k
	  // Inverse of current covariance matrix is already in C. 
          logPostCurr = 0.0; 
	  // Jacobian and Uniform prior. 
	  logPostCurr += log(phi[i] - phiA[i]) + log(phiB[i] - phi[i]); 
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &J, &one, C, &J, &w[i], &N, &zero, tmp_JD, &inc FCONE);
	  logPostCurr += -0.5*detCurr-0.5*F77_NAME(ddot)(&J, &w[i], &N, tmp_JD, &inc);
	  // Rprintf("logPostCurr: %f\n", logPostCurr); 
          if (corName == "matern"){
            logPostCurr += log(nu[i] - nuA[i]) + log(nuB[i] - nu[i]); 
          }
	  if (sigmaSqIG == 0) {
            logPostCurr += log(sigmaSq[i] - sigmaSqA[i]) + log(sigmaSqB[i] - sigmaSq[i]);
	  }

	  // MH Accept/Reject
	  logMHRatio = logPostCand - logPostCurr; 
	  if (runif(0.0, 1.0) <= exp(logMHRatio)) {
            theta[phiIndx * N + i] = phiCand; 
	    currTheta[phiIndx] = phiCand;
	    accept[phiIndx * N + i]++; 
            if (corName == "matern") {
              nu[i] = nuCand; 
	      currTheta[nuIndx] = nuCand; 
	      theta[nuIndx * N + i] = nu[i]; 
              accept[nuIndx * N + i]++; 
            }
	    if (sigmaSqIG == 0) {
              theta[sigmaSqIndx * N + i] = sigmaSqCand;
	      currTheta[sigmaSqIndx] = sigmaSqCand;
	      accept[sigmaSqIndx * N + i]++;
	    }
	    F77_NAME(dcopy)(&JJ, CCand, &inc, C, &inc); 
          }

          /********************************************************************
           *Update w (spatial random effects)
           *******************************************************************/
          /********************************
           * Compute b.w
           *******************************/
          for(j = 0; j < J; j++){
            tmp_JD[j] = kappaOcc[j] - omegaOcc[j] * (F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + betaStarSites[i * J + j]);
          }
          /********************************
           * Compute A.w
           *******************************/
	  // Copy inverse covariance matrix into tmp_JJ
	  F77_NAME(dcopy)(&JJ, C, &inc, tmp_JJ, &inc); 
	  for (k = 0; k < J; k++) {
	    tmp_JJ[k * J + k] += omegaOcc[k]; 
	  } // k
          F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf on A.w failed\n");}
          F77_NAME(dpotri)(lower, &J, tmp_JJ, &J, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri on A.w failed\n");}
          F77_NAME(dsymv)(lower, &J, &one, tmp_JJ, &J, tmp_JD, &inc, &zero, tmp_JD2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: dpotrf on A.w failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_J1, tmp_JD2, tmp_JJ, J);
          for (j = 0; j < J; j++) {
            w[j * N + i] = tmp_J1[j]; 
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
                  tmp_02 += betaStar[i * nOccRE + betaStarLongIndx[ll * J + j]];
	        } 
                  tmp_one[0] += kappaOcc[j] - (F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + tmp_02 - betaStar[i * nOccRE + l] + w[j * N + i]) * omegaOcc[j];
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
                betaStarSites[i * J + j] += betaStar[i * nOccRE + betaStarLongIndx[l * J + j]];
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
                  tmp_02 = 0.0;
                  for (ll = 0; ll < pDetRE; ll++) {
                    tmp_02 += alphaStar[i * nDetRE + alphaStarLongIndx[ll * nObs + r]];
	          } 
                  tmp_one[0] += kappaDet[r] - (F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + tmp_02 - alphaStar[i * nDetRE + l]) * omegaDet[r];
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
                alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + alphaStarLongIndx[l * nObs + r]]; 
              }
            }
          }

          /********************************************************************
           *Update Latent Occupancy
           *******************************************************************/
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

          // Compute occupancy probability 
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
	if (a >= nBurn) {
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
      } // b (end batch)

      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (i = 0; i < N; i++) {
        for (k = 0; k < nTheta; k++) {
          REAL(acceptSamples_r)[s * nThetaN + k * N + i] = accept[k * N + i]/batchLength; 
          REAL(tuningSamples_r)[s * nThetaN + k * N + i] = tuning[k * N + i]; 
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

    }
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

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

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

    // Rf_mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 10, Rf_mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 11, Rf_mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 12, Rf_mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 13, Rf_mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 14, Rf_mkChar("alpha.star.samples")); 
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


