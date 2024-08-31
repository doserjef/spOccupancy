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
  SEXP tMsPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
                SEXP consts_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
                SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
                SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
                SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
                SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
                SEXP betaStarStarting_r, SEXP alphaStarStarting_r, SEXP zLongIndx_r,
                SEXP zYearIndx_r, SEXP zDatIndx_r, SEXP zSiteIndx_r,
                SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
                SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
                SEXP muBetaComm_r, SEXP muAlphaComm_r, 
                SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
                SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r,  
                SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
                SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
                SEXP ar1_r, SEXP ar1Vals_r,
                SEXP tuning_r, SEXP nBatch_r, SEXP batchLength_r, 
                SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
                SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, g, t, h, r, l, ss, info, nProtect=0, rr;    

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
    int *XpRE = INTEGER(XpRE_r); 
    int *XRE = INTEGER(XRE_r);
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
    int nYearsMax = INTEGER(consts_r)[9];
    int indBetas = INTEGER(consts_r)[10];
    int indAlphas = INTEGER(consts_r)[11];
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
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *tuning = REAL(tuning_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *zYearIndx = INTEGER(zYearIndx_r); 
    int *zDatIndx = INTEGER(zDatIndx_r); 
    int *zSiteIndx = INTEGER(zSiteIndx_r);
    int yearIndx, siteIndx; 
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
    int sPost = 0; 
    /**********************************
     * AR1 Parameters 
     * *******************************/
    int ar1 = INTEGER(ar1_r)[0];
    double *rhoA = REAL(ar1Vals_r);
    double *rhoB = &REAL(ar1Vals_r)[N];
    double *sigmaSqTA = &REAL(ar1Vals_r)[2*N];
    double *sigmaSqTB = &REAL(ar1Vals_r)[3*N];

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
        Rprintf("Multi-season Multi-species Occupancy Model with Polya-Gamma latent\nvariables with %i sites, %i species, and %i primary time periods.\n\n", J, N, nYearsMax);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
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
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pOccN = pOcc * N; 
    int pDetN = pDet * N; 
    int nObsN = nObs * N; 
    int nDetREN = nDetRE * N; 
    int nOccREN = nOccRE * N; 
    int JN = J * N;
    int NnYears = nYearsMax * N;
    int JNnYears = JN * nYearsMax;
    int JnYears = J * nYearsMax;
    int JnYearspOcc = JnYears * pOcc;
    int JnYearspOccRE = JnYears * pOccRE;
    int nObspDet = nObs * pDet;
    int nObspDetRE = nObs * pDetRE;
    int nnYears = nYearsMax * nYearsMax;
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
    int *tmp_JnYearsInt = (int *) R_alloc(JnYears, sizeof(int));
    for (j = 0; j < JnYears; j++) {
      tmp_JnYearsInt[j] = 0; 
    }
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J, J);
    double *tmp_JnYears = (double *) R_alloc(JnYears, sizeof(double));
    zeros(tmp_JnYears, JnYears);
    double *tmp_JnYearspOcc = (double *) R_alloc(JnYears * pOcc, sizeof(double));
    zeros(tmp_JnYearspOcc, JnYears * pOcc);
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_NnYears = (double *) R_alloc(NnYears, sizeof(double));
    zeros(tmp_NnYears, NnYears);
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 

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
    // Temporal variance parameter
    double *sigmaSqT = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, &REAL(ar1Vals_r)[5 * N], &inc, sigmaSqT, &inc);
    // Temporal correlation parameter
    double *rho = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, &REAL(ar1Vals_r)[4 * N], &inc, rho, &inc);
    // Latent Occurrence
    double *z = (double *) R_alloc(JNnYears, sizeof(double));   
    F77_NAME(dcopy)(&JNnYears, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
    // Note, you can just write over the values for the detection
    // parameters. 
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(JNnYears, sizeof(double)); 
    zeros(omegaOcc, JNnYears);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(JNnYears, sizeof(double)); 
    zeros(kappaOcc, JNnYears);
    // Need this for all species
    double *zStar = (double *) R_alloc(JNnYears, sizeof(double));
    zeros(zStar, JNnYears);

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
    SEXP etaSamples_r; 
    SEXP thetaSamples_r;
    int nTheta = 2, sigmaSqTIndx = 0, rhoIndx = 1;
    int nThetaN = nTheta * N; 
    if (ar1) {
      PROTECT(etaSamples_r = Rf_allocMatrix(REALSXP, NnYears, nPost)); nProtect++; 
      zeros(REAL(etaSamples_r), NnYears * nPost);
      PROTECT(thetaSamples_r = Rf_allocMatrix(REALSXP, nThetaN, nPost)); nProtect++; 
      zeros(REAL(thetaSamples_r), nThetaN * nPost);
    }
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOccN, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOccN * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, pDetN, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDetN * nPost);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, JNnYears, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), JNnYears * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, JNnYears, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), JNnYears * nPost);
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
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, JNnYears, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), JNnYears * nPost);
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObsN, sizeof(double)); 
    double *yWAIC = (double *) R_alloc(JNnYears, sizeof(double)); ones(yWAIC, JNnYears);
    double *psi = (double *) R_alloc(JNnYears, sizeof(double)); 
    zeros(psi, JNnYears); 
    double *piProd = (double *) R_alloc(JNnYears, sizeof(double)); 
    ones(piProd, JNnYears); 
    double *piProdWAIC = (double *) R_alloc(JNnYears, sizeof(double)); 
    ones(piProdWAIC, JNnYears); 
    double *ySum = (double *) R_alloc(JNnYears, sizeof(double)); zeros(ySum, JNnYears); 

    // For normal community-level priors
    // Occurrence coefficients
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
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
    // Put community level occurrence variances in a pOcc x pOcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level detection variances in a pDet x pDet matrix. 
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauSqAlpha[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri TauAlphaInv failed\n");}

    /**********************************************************************
     * Prep for random effects (if they exist)
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JNnYears, sizeof(double)); 
    zeros(betaStarSites, JNnYears); 
    int *betaStarLongIndx = (int *) R_alloc(JnYearspOccRE, sizeof(int));
    // Initial sums (initiate with the first species)
    for (t = 0; t < nYearsMax; t++) {
      for (j = 0; j < J; j++) {
        for (l = 0; l < pOccRE; l++) {
          betaStarLongIndx[l * JnYears + t * J + j] = which(XRE[l * JnYears + t * J + j], 
        		                                    betaLevelIndx, nOccRE);
          // Rprintf("betaStarLongIndx[%i]: %f\n", l * JnYears + t * J + j, betaStarLongIndx[l * JnYears + t * J + j]);
          for (i = 0; i < N; i++) {
            betaStarSites[t * JN + j * N + i] += betaStar[i * nOccRE + betaStarLongIndx[l * JnYears + t * J + j]];
          }
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
     * Set up AR1 stuff
     *********************************************************************/
    double *theta = (double *) R_alloc(nThetaN, sizeof(double));
    for (i = 0; i < N; i++) {
      theta[sigmaSqTIndx * N + i] = sigmaSqT[i];
      theta[rhoIndx * N + i] = rho[i];
    }
    // Note you are overwriting these for each species each time. 
    double *SigmaEta = (double *) R_alloc(nnYears, sizeof(double));
    double *SigmaEtaCand = (double *) R_alloc(nnYears, sizeof(double));
    double *tmp_nYearsMax = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nYearsMax2 = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nnYears = (double *) R_alloc(nnYears, sizeof(double));
    // Initiate AR(1) matrix with first species. 
    if (ar1) {
      AR1(nYearsMax, theta[rhoIndx * N], theta[sigmaSqTIndx * N], SigmaEta);
      clearUT(SigmaEta, nYearsMax);
      F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){Rf_error("c++ error: Cholesky failed in initial time covariance matrix\n");}
      F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){Rf_error("c++ error: Cholesky inverse failed in initial time covariance matrix\n");}
    }
    // Initiate temporal random effects to zero.
    double *eta = (double *) R_alloc(NnYears, sizeof(double)); zeros(eta, NnYears);
    double *etaTRInv = (double *) R_alloc(nYearsMax, sizeof(double));
    double aSigmaSqTPost = 0.0;
    double bSigmaSqTPost = 0.0;

    /**********************************************************************
     Set up stuff for Adaptive MH and other misc
     * *******************************************************************/
    double logMHRatio, logPostCurr = 0.0, logPostCand = 0.0, detCand = 0.0, detCurr = 0.0, rhoCand = 0.0; 
    double *accept = (double *) R_alloc(nThetaN, sizeof(double)); zeros(accept, nThetaN); 
    // MCMC info if desired
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = Rf_allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nThetaN * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = Rf_allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nThetaN * nBatch);

    GetRNGstate();

    /**********************************************************************
     Start sampling
     * *******************************************************************/
    for (s = 0, g = 0; s < nBatch; s++) {
      for (ss = 0; ss < batchLength; ss++, g++) {

        /********************************************************************
         Update Community level Occupancy Coefficients
         *******************************************************************/
        if (indBetas == 0) {
          /********************************
           Compute b.beta.comm
           *******************************/
          zeros(tmp_pOcc, pOcc); 
          for (i = 0; i < N; i++) {
            F77_NAME(dgemv)(ytran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, &beta[i], &N, &one, tmp_pOcc, &inc FCONE); 
          } // i
          for (h = 0; h < pOcc; h++) {
            tmp_pOcc[h] += SigmaBetaCommInvMuBeta[h];  
          } // j

          /********************************
           Compute A.beta.comm
           *******************************/
          for (h = 0; h < ppOcc; h++) {
            tmp_ppOcc[h] = SigmaBetaCommInv[h] + N * TauBetaInv[h]; 
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
	}
        /********************************************************************
         Update Community level Detection Coefficients
         *******************************************************************/
	if (indAlphas == 0) {
          /********************************
           * Compute b.alpha.comm
           *******************************/
           zeros(tmp_pDet, pDet); 
           for (i = 0; i < N; i++) {
             F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, TauAlphaInv, &pDet, &alpha[i], &N, &one, tmp_pDet, &inc FCONE); 
           } // i
           for (h = 0; h < pDet; h++) {
             tmp_pDet[h] += SigmaAlphaCommInvMuAlpha[h];  
           } // j
          /********************************
           * Compute A.alpha.comm
           *******************************/
          for (h = 0; h < ppDet; h++) {
            tmp_ppDet[h] = SigmaAlphaCommInv[h] + N * TauAlphaInv[h]; 
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
	}

        /********************************************************************
         Update Community Occupancy Variance Parameter
        ********************************************************************/
	if (indBetas == 0) {
          for (h = 0; h < pOcc; h++) {
            tmp_0 = 0.0;  
            for (i = 0; i < N; i++) {
              tmp_0 += (beta[h * N + i] - betaComm[h]) * (beta[h * N + i] - betaComm[h]);
            } // i
            tmp_0 *= 0.5;
            tauSqBeta[h] = rigamma(tauSqBetaA[h] + N / 2.0, tauSqBetaB[h] + tmp_0); 
          } // h
          for (h = 0; h < pOcc; h++) {
            TauBetaInv[h * pOcc + h] = tauSqBeta[h]; 
          } // i
          F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
          F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}
	}
        /********************************************************************
         Update Community Detection Variance Parameter
        ********************************************************************/
	if (indAlphas == 0) {
          for (h = 0; h < pDet; h++) {
            tmp_0 = 0.0;  
            for (i = 0; i < N; i++) {
              tmp_0 += (alpha[h * N + i] - alphaComm[h]) * (alpha[h * N + i] - alphaComm[h]);
            } // i
            tmp_0 *= 0.5;
            tauSqAlpha[h] = rigamma(tauSqAlphaA[h] + N / 2.0, tauSqAlphaB[h] + tmp_0); 
          } // h
          for (h = 0; h < pDet; h++) {
            TauAlphaInv[h * pDet + h] = tauSqAlpha[h]; 
          } // i
          F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf TauAlphaInv failed\n");}
          F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri TauAlphaInv failed\n");}
	}

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
          for (t = 0; t < nYearsMax; t++) {
            for (j = 0; j < J; j++) {
              // Only calculate omegaOcc when there are observations at that 
              // site/year combination
              if (zDatIndx[t * J + j] == 1) {
                omegaOcc[t * JN + j * N + i] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, &beta[i], &N) + 
        			                        eta[i * nYearsMax + t] + betaStarSites[t * JN + j * N + i]);
                kappaOcc[t * JN + j * N + i] = z[t * JN + j * N + i] - 1.0 / 2.0; 
                zStar[t * JN + j * N + i] = kappaOcc[t * JN + j * N + i] / omegaOcc[t * JN + j * N + i];
              }
            } // j
          } // t
          /********************************************************************
           *Update Detection Auxiliary Variables 
           *******************************************************************/
           for (r = 0; r < nObs; r++) {
             yearIndx = zYearIndx[zLongIndx[r]]; 
             siteIndx = zSiteIndx[zLongIndx[r]];
             omegaDet[r] = 0.0; // reset
             // If the site is occupied and data were collected at that site/year combo.
             if ((z[yearIndx * JN + siteIndx * N + i] == 1.0) && (zDatIndx[zLongIndx[r]] == 1)) {
               omegaDet[r] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
             }
           } // r

          /********************************************************************
           *Update Occupancy Regression Coefficients
           *******************************************************************/
          zeros(tmp_JnYears, JnYears);
          for (t = 0; t < nYearsMax; t++) {
            for (j = 0; j < J; j++) {
              if (zDatIndx[t * J + j] == 1) {
                tmp_JnYears[t * J + j] = kappaOcc[t * JN + j * N + i] - omegaOcc[t * JN + j * N + i] * 
                          (eta[i * nYearsMax + t] + betaStarSites[t * JN + j * N + i]); 
              }
            } // j
          } // t
          /********************************
           * Compute b.beta
           *******************************/
          // t(X) * tmp_J1 + 0 * tmp_pOcc = tmp_pOcc. 
          F77_NAME(dgemv)(ytran, &JnYears, &pOcc, &one, X, &JnYears, 
                	  tmp_JnYears, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
          // TauBetaInv %*% betaComm + tmp_pOcc = tmp_pOcc
          F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, betaComm, &inc, &one, tmp_pOcc, &inc FCONE); 

          /********************************
           * Compute A.beta
           * *****************************/
          // t(X) %*% diag(omegaOcc)
          zeros(tmp_JnYearspOcc, JnYearspOcc);
          for (t = 0; t < nYearsMax; t++) {
            for(j = 0; j < J; j++){
              if (zDatIndx[t * J + j] == 1) {
                for(h = 0; h < pOcc; h++){
                  tmp_JnYearspOcc[h * JnYears + t * J + j] = X[h * JnYears + t * J + j] * omegaOcc[t * JN + j * N + i];
                }
              }
            }
          }
          // This finishes off A.beta
          // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
          F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &JnYears, &one, X, &JnYears, tmp_JnYearspOcc, &JnYears, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
          for (h = 0; h < ppOcc; h++) {
            tmp_ppOcc[h] += TauBetaInv[h]; 
          } // j
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf A.beta 2 failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_beta, tmp_pOcc2, tmp_ppOcc, pOcc);
          // Can eventually get rid of this and change order of beta. 
          for (h = 0; h < pOcc; h++) {
            beta[h * N + i] = tmp_beta[h]; 
          }
        
          /********************************************************************
           *Update Detection Regression Coefficients
           *******************************************************************/
          /********************************
           * Compute b.alpha
           *******************************/
          // Only calculate when z == 1 and the site/year combo has an observation. 
          for (r = 0; r < nObs; r++) {
            kappaDet[r] = 0.0;
            tmp_nObs[r] = 0.0;
            yearIndx = zYearIndx[zLongIndx[r]]; 
            siteIndx = zSiteIndx[zLongIndx[r]];
            // If the site is occupied and data were collected at that site/year combo:
            if ((z[yearIndx * JN + siteIndx * N + i] == 1.0) && (zDatIndx[zLongIndx[r]] == 1)) {
              kappaDet[r] = (y[r * N + i] - 1.0/2.0) * z[yearIndx * JN + siteIndx * N + i];
              tmp_nObs[r] = kappaDet[r] - omegaDet[r] * alphaStarObs[i * nObs + r]; 
              tmp_nObs[r] *= z[yearIndx * JN + siteIndx * N + i];
            }
          } // r
          
          F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
          F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, TauAlphaInv, &pDet, alphaComm, &inc, &one, tmp_pDet, &inc FCONE); 
          /********************************
           * Compute A.alpha
           * *****************************/
          for (r = 0; r < nObs; r++) {
            yearIndx = zYearIndx[zLongIndx[r]]; 
            siteIndx = zSiteIndx[zLongIndx[r]];
            for (h = 0; h < pDet; h++) {
              tmp_nObspDet[h*nObs + r] = Xp[h * nObs + r] * omegaDet[r] * z[yearIndx * JN + siteIndx * N + i];
            } // h
          } // r

          // This finishes off A.alpha
          // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
          F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

          for (h = 0; h < ppDet; h++) {
            tmp_ppDet[h] += TauAlphaInv[h]; 
          } // h
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
          if(info != 0){Rf_error("c++ error: dpotrf A.alpha 2 failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDet);
          for (h = 0; h < pDet; h++) {
            alpha[h * N + i] = tmp_alpha[h];
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
              tmp_one[0] = 0.0;
              tmp_0 = 0.0;	      
              // Only allow information to come from when XRE == betaLevelIndx[l]. 
              // aka information only comes from the sites with any given level 
              // of a random effect. 
              for (t = 0; t < nYearsMax; t++) {
                for (j = 0; j < J; j++) {
                  if (XRE[betaStarIndx[l] * JnYears + t * J + j] == betaLevelIndx[l]) {
                    if (zDatIndx[t * J + j] == 1) {
                      tmp_02 = 0.0;
        	      for (rr = 0; rr < pOccRE; rr++) {
                        tmp_02 += betaStar[i * nOccRE + betaStarLongIndx[rr * JnYears + t * J + j]]; 
        	      }
                      tmp_one[0] += kappaOcc[t * JN + j * N + i] - (F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, &beta[i], &N) + 
                		    tmp_02 - betaStar[i * nOccRE + l] + eta[i * nYearsMax + t]) * omegaOcc[t * JN + j * N + i];
                      tmp_0 += omegaOcc[t * JN + j * N + i];
                    }
                  }
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
            for (t = 0; t < nYearsMax; t++) {
              for (j = 0; j < J; j++) {
                betaStarSites[t * JN + j * N + i] = 0.0;
                for (l = 0; l < pOccRE; l++) {
                  betaStarSites[t * JN + j * N + i] += betaStar[i * nOccRE + betaStarLongIndx[l * JnYears + t * J + j]];
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
              for (r = 0; r < nObs; r++) {
                yearIndx = zYearIndx[zLongIndx[r]]; 
                siteIndx = zSiteIndx[zLongIndx[r]];
                if ((z[yearIndx * JN + siteIndx * N + i] == 1.0) && (XpRE[alphaStarIndx[l] * nObs + r] == alphaLevelIndx[l])) {
                  tmp_02 = 0.0;
        	  for (rr = 0; rr < pDetRE; rr++) {
                    tmp_02 += alphaStar[i * nDetRE + alphaStarLongIndx[rr * nObs + r]];
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
                alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE)]; 
              }
            }
          }

          if (ar1) {
            /********************************************************************
             *Update sigmaSqT
             *******************************************************************/
            // Form correlation matrix. 
            AR1(nYearsMax, theta[rhoIndx * N + i], 1.0, SigmaEta);
            clearUT(SigmaEta, nYearsMax);
            F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
            F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: Cholesky inverse failed in covariance matrix\n");}
            fillUTri(SigmaEta, nYearsMax);
            // Compute t(eta) %*% SigmaEta^-1 %*% eta
            for (t = 0; t < nYearsMax; t++) {
              etaTRInv[t] = F77_NAME(ddot)(&nYearsMax, &SigmaEta[t], &nYearsMax, 
            		               &eta[i * nYearsMax], &inc); 
            }
            bSigmaSqTPost = F77_NAME(ddot)(&nYearsMax, etaTRInv, &inc, &eta[i * nYearsMax], &inc);	
            bSigmaSqTPost /= 2.0;
            bSigmaSqTPost += sigmaSqTB[i];
	    aSigmaSqTPost = 0.5 * nYearsMax + sigmaSqTA[i];
            theta[sigmaSqTIndx * N + i] = rigamma(aSigmaSqTPost, bSigmaSqTPost);
            
            /********************************************************************
             *Update rho
             *******************************************************************/
            rho[i] = theta[rhoIndx * N + i];
            sigmaSqT[i] = theta[sigmaSqTIndx * N + i];
            rhoCand = logitInv(rnorm(logit(rho[i], rhoA[i], rhoB[i]), 
				     exp(tuning[rhoIndx * N + i])), rhoA[i], rhoB[i]); 
            theta[rhoIndx * N + i] = rhoCand; 

            // Construct proposal covariance matrix. 
            AR1(nYearsMax, theta[rhoIndx * N + i], theta[sigmaSqTIndx * N + i], SigmaEtaCand);
            clearUT(SigmaEtaCand, nYearsMax);

            /********************************
             * Proposal
             *******************************/
            // Invert SigmaEtaCand and log det cov. 
            detCand = 0.0;
            F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: Cholesky failed in proposal covariance matrix\n");}
            // Get log of the determinant of the covariance matrix. 
            for (k = 0; k < nYearsMax; k++) {
              detCand += 2.0 * log(SigmaEtaCand[k*nYearsMax+k]);
            } // k
            F77_NAME(dpotri)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: Cholesky inverse failed in proposal covariance matrix\n");}
            logPostCand = 0.0; 
            // Jacobian and Uniform prior. 
            logPostCand += log(rhoCand - rhoA[i]) + log(rhoB[i] - rhoCand); 
            F77_NAME(dsymv)(lower, &nYearsMax, &one,  SigmaEtaCand, &nYearsMax, 
			    &eta[i * nYearsMax], &inc, &zero, tmp_nYearsMax, &inc FCONE);
            logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&nYearsMax, &eta[i * nYearsMax], 
			                                   &inc, tmp_nYearsMax, &inc);
            /********************************
             * Current
             *******************************/
            theta[rhoIndx * N + i] = rho[i]; 
            AR1(nYearsMax, theta[rhoIndx * N + i], theta[sigmaSqTIndx * N + i], SigmaEta);
            clearUT(SigmaEta, nYearsMax);
            detCurr = 0.0;
            F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
            for (k = 0; k < nYearsMax; k++) {
              detCurr += 2.0 * log(SigmaEta[k*nYearsMax+k]);
            } // k
            F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: Cholesky inverse failed in covariance matrix\n");}
            logPostCurr = 0.0; 
            logPostCurr += log(rho[i] - rhoA[i]) + log(rhoB[i] - rho[i]); 
            // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
            F77_NAME(dsymv)(lower, &nYearsMax, &one, SigmaEta, &nYearsMax, 
			    &eta[i * nYearsMax], &inc, &zero, 
            		tmp_nYearsMax, &inc FCONE);
            logPostCurr += -0.5*detCurr-0.5*F77_NAME(ddot)(&nYearsMax, &eta[i * nYearsMax], 
			                                   &inc, tmp_nYearsMax, &inc);

            // MH Accept/Reject
            logMHRatio = logPostCand - logPostCurr; 
            if (runif(0.0, 1.0) <= exp(logMHRatio)) {
              theta[rhoIndx * N + i] = rhoCand;
              accept[rhoIndx * N + i]++;
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
                  tmp_nYearsMax[t] += kappaOcc[t * JN + j * N + i] - omegaOcc[t * JN + j * N + i] * 
			              (F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, &beta[i], &N) + 
				       betaStarSites[t * JN + j * N + i]);
                }
              }
            }
            /********************************
             * Compute A.w
             *******************************/
            // Copy inverse covariance matrix into tmp_nnYears
            F77_NAME(dcopy)(&nnYears, SigmaEta, &inc, tmp_nnYears, &inc); 
            for (j = 0; j < J; j++) {
              for (t = 0; t < nYearsMax; t++) {
                if (zDatIndx[t * J + j] == 1) {
                  tmp_nnYears[t * nYearsMax + t] += omegaOcc[t * JN + j * N + i];
                }
              } // t
            } // j

            // Cholesky of A.eta
            F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: dpotrf on A.eta failed\n");}
            // Inverse of A.eta
            F77_NAME(dpotri)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: dpotri on A.eta failed\n");}
            // A.eta.inv %*% b.eta. Stored in tmp_
            F77_NAME(dsymv)(lower, &nYearsMax, &one, tmp_nnYears, &nYearsMax, 
            		tmp_nYearsMax, &inc, &zero, tmp_nYearsMax2, &inc FCONE);
            F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
            if(info != 0){Rf_error("c++ error: dpotrf on A.eta failed\n");}
            // Args: destination, mu, cholesky of the covariance matrix, dimension
            mvrnorm(&eta[i * nYearsMax], tmp_nYearsMax2, tmp_nnYears, nYearsMax);
          }

          /********************************************************************
           *Update Latent Occupancy
           *******************************************************************/
          // Compute detection probability 
          for (r = 0; r < nObs; r++) {
            yearIndx = zYearIndx[zLongIndx[r]];
            siteIndx = zSiteIndx[zLongIndx[r]];
            detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
            if (tmp_JnYearsInt[zLongIndx[r]] == 0) {
              psi[yearIndx * JN + siteIndx * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[yearIndx * J + siteIndx], &JnYears, &beta[i], &N) + 
                	                                               eta[i * nYearsMax + yearIndx] + 
                					               betaStarSites[yearIndx * JN + siteIndx * N + i], zero, one);
            }
            piProd[yearIndx * JN + siteIndx * N + i] *= (1.0 - detProb[i * nObs + r]);
            piProdWAIC[yearIndx * JN + siteIndx * N + i] *= pow(detProb[i * nObs + r], y[r * N + i]);
            piProdWAIC[yearIndx * JN + siteIndx * N + i] *= pow(1.0 - detProb[i * nObs + r], 
                	                                                1.0 - y[r * N + i]);
            ySum[yearIndx * JN + siteIndx * N + i] += y[r * N + i]; 
            tmp_JnYearsInt[zLongIndx[r]]++;
          } // r
          // Compute occupancy probability and the integrated likelihood for WAIC
          for (t = 0; t < nYearsMax; t++ ) {
            for (j = 0; j < J; j++) {
              psiNum = psi[t * JN + j * N + i] * piProd[t * JN + j * N + i]; 
              if (zDatIndx[t * J + j] == 1) {
                if (ySum[t * JN + j * N + i] == zero) {
                  z[t * JN + j * N + i] = rbinom(one, psiNum / (psiNum + (1.0 - psi[t * JN + j * N + i])));
                  yWAIC[t * JN + j * N + i] = (1.0 - psi[t * JN + j * N + i]) + psi[t * JN + j * N + i] * piProdWAIC[t * JN + j * N + i];
                } else {
                  z[t * JN + j * N + i] = one; 
                  yWAIC[t * JN + j * N + i] = psi[t * JN + j * N + i] * piProdWAIC[t * JN + j * N + i]; 
                }
              } else {
                psi[t * JN + j * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, &beta[i], &N) + 
                 	                           eta[i * nYearsMax + t] + 
                 			           betaStarSites[t * JN + j * N + i], zero, one);
                z[t * JN + j * N + i] = rbinom(one, psi[t * JN + j * N + i]);		
              }
              // Reset variables
              piProd[t * JN + j * N + i] = one;
              piProdWAIC[t * JN + j * N + i] = one;
              ySum[t * JN + j * N + i] = zero; 
              tmp_JnYearsInt[t * J + j] = 0; 
            } // j
          } // t
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
            F77_NAME(dcopy)(&JNnYears, psi, &inc, &REAL(psiSamples_r)[sPost*JNnYears], &inc); 
            F77_NAME(dcopy)(&JNnYears, z, &inc, &REAL(zSamples_r)[sPost*JNnYears], &inc); 
            if (ar1) {
              F77_NAME(dcopy)(&nThetaN, theta, &inc, 
        		      &REAL(thetaSamples_r)[sPost*nThetaN], &inc);
              F77_NAME(dcopy)(&NnYears, eta, &inc, &REAL(etaSamples_r)[sPost*NnYears], &inc);
            }
            if (pDetRE > 0) {
              F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
              F77_NAME(dcopy)(&nDetREN, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetREN], &inc);
            }
            if (pOccRE > 0) {
              F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
              F77_NAME(dcopy)(&nOccREN, betaStar, &inc, &REAL(betaStarSamples_r)[sPost*nOccREN], &inc);
            }
            F77_NAME(dcopy)(&JNnYears, yWAIC, &inc, &REAL(likeSamples_r)[sPost*JNnYears], &inc); 
            sPost++; 
            thinIndx = 0; 
          }
        }
        R_CheckUserInterrupt();
      } // ss (end batch)

      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      if (ar1) {
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
      }
        
      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          if (ar1) {
            Rprintf("\tSpecies\t\tParameter\tAcceptance\tTuning\n");	  
            for (i = 0; i < N; i++) {
                Rprintf("\t%i\t\trho\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nThetaN + rhoIndx * N + i], exp(tuning[rhoIndx * N + i]));
            } // i
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

    // make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 11;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }
    if (ar1) {
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
    SET_VECTOR_ELT(result_r, 8, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 9, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 10, likeSamples_r); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 11, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 12, alphaStarSamples_r);
    }
    if (pOccRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 13;
      } else {
        tmp_0 = 11;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    int ar1Ind = 0;
    if (ar1) {
      if (pOccRE > 0) {
        ar1Ind = tmp_0 + 2;
      } else if (pDetRE > 0) {
        ar1Ind = 13; 
      } else {
        ar1Ind = 11;
      }
      SET_VECTOR_ELT(result_r, ar1Ind, etaSamples_r);
      SET_VECTOR_ELT(result_r, ar1Ind + 1, thetaSamples_r);
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
    SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 10, Rf_mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 11, Rf_mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 12, Rf_mkChar("alpha.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, Rf_mkChar("beta.star.samples")); 
    }
    if (ar1) {
      SET_VECTOR_ELT(resultName_r, ar1Ind, Rf_mkChar("eta.samples")); 
      SET_VECTOR_ELT(resultName_r, ar1Ind + 1, Rf_mkChar("theta.samples")); 
    }
   
    // Set the names of the output list.  
    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}




