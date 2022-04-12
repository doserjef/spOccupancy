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

extern "C" {
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
		 SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, g, t, h, r, l, info, nProtect=0, ii, ll;    

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
    // Xp is sorted by parameter, then by visit, then by site 
    double *Xp = REAL(Xp_r);
    int *XRE = INTEGER(XRE_r);
    int *XpRE = INTEGER(XpRE_r); 
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
    int q = INTEGER(consts_r)[9]; 
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
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    double *K = REAL(K_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nSamples = INTEGER(nSamples_r)[0];
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 
    int thinIndx = 0; 
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
        Rprintf("Latent Factor Multispecies Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        Rprintf("Using %i latent factors.\n\n", q);
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
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
    int Jq = J * q;
    int qq = q * q;
    int JN = J * N;
    int Nq = N * q;
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    int JJ = J * J;
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
    int *tmp_JInt = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_JInt[j] = 0; 
    }
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J, J);
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_qq = (double *) R_alloc(qq, sizeof(double));
    double *tmp_q = (double *) R_alloc(q, sizeof(double));
    double *tmp_q2 = (double *) R_alloc(q, sizeof(double));
    double *tmp_qq2 = (double *) R_alloc(qq, sizeof(double));
    double *tmp_Jq = (double *) R_alloc(Jq, sizeof(double));
    double *tmp_Nq = (double *) R_alloc(Nq, sizeof(double));
    double *tmp_N = (double *) R_alloc(N, sizeof(double));
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    int currDim = 0;

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
    // Latent effects
    double *w = (double *) R_alloc(Jq, sizeof(double)); zeros(w, Jq);
    // Latent spatial factors
    double *lambda = (double *) R_alloc(Nq, sizeof(double));
    F77_NAME(dcopy)(&Nq, REAL(lambdaStarting_r), &inc, lambda, &inc);
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREN, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREN, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Latent occurrence random effects
    double *betaStar = (double *) R_alloc(nOccREN, sizeof(double)); 
    F77_NAME(dcopy)(&nOccREN, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
    // Note, you can just write over the values for the detection
    // parameters. 
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(JN, sizeof(double)); zeros(omegaOcc, JN);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(JN, sizeof(double)); zeros(kappaOcc, JN);
    // Need this for all species
    double *zStar = (double *) R_alloc(JN, sizeof(double));

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
    // Factor model parameters
    SEXP lambdaSamples_r; 
    PROTECT(lambdaSamples_r = allocMatrix(REALSXP, Nq, nPost)); nProtect++;
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, Jq, nPost)); nProtect++; 
    // Detection random effects
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

    // For normal community-level priors
    // Occurrence coefficients
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    // dgemv computes linear combinations of different variables. 
    F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, SigmaBetaCommInvMuBeta, &inc FCONE); 	  
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, SigmaAlphaCommInvMuAlpha, &inc FCONE); 
    // Put community level occurrence variances in a pOcc x pOcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level detection variances in a pDet x pDet matrix. 
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
     Set up factor model parameters
     * *******************************************************************/
    // Species-level latent effects
    double *wStar = (double *) R_alloc(JN, sizeof(double)); zeros(wStar, JN);
    // Multiply Lambda %*% w[j] to get wStar. 
    for (j = 0; j < J; j++) {
      F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q], &inc, &zero, &wStar[j * N], &inc FCONE);
    }
    double *mu = (double *) R_alloc(q, sizeof(double));
    double *var = (double *) R_alloc(qq, sizeof(double));

    GetRNGstate();

    /**********************************************************************
     Start sampling
     * *******************************************************************/
    for (s = 0; s < nSamples; s++) {

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
      for (h = 0; h < pOcc; h++) {
        tmp_0 = 0.0;  
        for (i = 0; i < N; i++) {
          tmp_0 += (beta[h * N + i] - betaComm[h]) * (beta[h * N + i] - betaComm[h]);
        } // i
        tmp_0 *= 0.5;
        tauSqBeta[h] = rigamma(tauSqBetaA[h] + N / 2.0, tauSqBetaB[h] + tmp_0); 
      } // h
      // This is correct, nothing wrong here. 
      for (h = 0; h < pOcc; h++) {
        TauBetaInv[h * pOcc + h] = tauSqBeta[h]; 
      } // i
      F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
      if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
      if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
      /********************************************************************
       Update Community Detection Variance Parameter
      ********************************************************************/
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
          omegaOcc[j * N + i] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + wStar[j * N + i] + betaStarSites[i * J + j]);
        } // j
        /********************************************************************
         *Update Detection Auxiliary Variables 
         *******************************************************************/
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
          kappaOcc[j * N + i] = z[j * N + i] - 1.0 / 2.0; 
          tmp_J1[j] = kappaOcc[j * N + i] - omegaOcc[j * N + i] * 
		      (wStar[j * N + i] + betaStarSites[i * J + j]); 
	  // For later
	  zStar[j * N + i] = kappaOcc[j * N + i] / omegaOcc[j * N + i];
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
          for(h = 0; h < pOcc; h++){
            tmp_JpOcc[h*J+j] = X[h*J+j]*omegaOcc[j * N + i];
          }
        }
        // This finishes off A.beta
        // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
        for (h = 0; h < ppOcc; h++) {
          tmp_ppOcc[h] += TauBetaInv[h]; 
        } // j
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotri ABeta failed\n");}
        // A.beta.inv %*% b.beta
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
	if(info != 0){error("c++ error: dpotrf A.beta 2 failed\n");}
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
          for (h = 0; h < pDet; h++) {
            tmp_nObspDet[h*nObs + r] = Xp[h * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i];
          } // i
        } // j

        // This finishes off A.alpha
        // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
        F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

        for (h = 0; h < ppDet; h++) {
          tmp_ppDet[h] += TauAlphaInv[h]; 
        } // h
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
        if(info != 0){error("c++ error: dpotrf A.alpha 2 failed\n");}
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
            zeros(tmp_one, inc);
            tmp_0 = 0.0;	      
	    // Only allow information to come from when XRE == betaLevelIndx[l]. 
	    // aka information only comes from the sites with any given level 
	    // of a random effect. 
            for (j = 0; j < J; j++) {
              if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l]) {
                tmp_one[0] += kappaOcc[j * N + i] - (F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + betaStarSites[i * J + j] - betaStar[i * nOccRE + l] + wStar[j * N + i]) * omegaOcc[j * N + i];
	        tmp_0 += omegaOcc[j * N + i];
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
      }

      /********************************************************************
       *Update latent effects (w) 
       *******************************************************************/
      for (ii = 0; ii < J; ii++) {
        // tmp_qq = lambda' S_beta lambda 
        for (i = 0; i < N; i++) {
          for (ll = 0; ll < q; ll++) {
            tmp_Nq[ll * N + i] = lambda[ll * N + i] * omegaOcc[ii * N + i];
          } // ll
        } // i
	F77_NAME(dgemm)(ytran, ntran, &q, &q, &N, &one, tmp_Nq, &N, lambda, &N, &zero, var, &q FCONE FCONE);

	// var
	for (k = 0; k < q; k++) {
          var[k * q + k] += 1.0; 
        } // k
	F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE);
        if (info != 0){error("c++ error: dpotrf var failed\n");}
	F77_NAME(dpotri)(lower, &q, var, &q, &info FCONE);
        if (info != 0){error("c++ error: dpotri var failed\n");}

	// mu
	for (k = 0; k < N; k++) {
          tmp_N[k] = (zStar[ii * N + k] - F77_NAME(ddot)(&pOcc, &X[ii], &J, &beta[k], &N) - betaStarSites[k * J + ii]) * omegaOcc[ii * N + k];
        } // k

	F77_NAME(dgemv)(ytran, &N, &q, &one, lambda, &N, tmp_N, &inc, &zero, mu, &inc FCONE);

        F77_NAME(dsymv)(lower, &q, &one, var, &q, mu, &inc, &zero, tmp_N, &inc FCONE);

	F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf var 2 failed\n");}

	mvrnorm(&w[ii * q], tmp_N, var, q);

      } // ii
      /********************************************************************
       *Update spatial factors (lambda)
       *******************************************************************/
      for (i = 1; i < N; i++) {
        zeros(tmp_qq, qq);
        zeros(tmp_q, q);
        zeros(tmp_qq2, qq);
        // W' %*% S_beta %*% W
        for (k = 0; k < q; k++) {
          for (l = 0; l < q; l++) {
            for (j = 0; j < J; j++) {
              tmp_qq[k * q + l] += w[j * q + k] * w[j * q + l] * omegaOcc[j * N + i];
            } // j
          } // l
        } // k


        // currDim gives the mean dimension of mu and var. 
        if (i < q) {
          currDim = i;  
        } else {
          currDim = q;
        }
        /*****************************
         *mu
         *****************************/
        // zStar - X %*% beta
        for (j = 0; j < J; j++) {
          tmp_J[j] = zStar[j * N + i] - F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) -
		     betaStarSites[i * J + j];

          if (i < q) {
            tmp_J[j] -= w[j * q + i];
          }
        } // j

        // S_beta %*% W' = tmp_Jq
        // aka multiply W[j, ] by omegaOcc[j] of the current species you're on. 
        for (j = 0, l = 0; j < J; j++) {
          for (ll = 0; ll < q; ll++, l++) {
            tmp_Jq[l] = omegaOcc[j * N + i] * w[j * q + ll];  
          }
        }

        // tmp_Jq %*% tmp_J
        for (k = 0; k < currDim; k++) {
          for (j = 0; j < J; j++) {
            tmp_q[k] += tmp_Jq[j * q + k] * tmp_J[j];
          } // j
        } // k

        /*****************************
         *var
         *****************************/
        // Only get relevant columns of t(W) %*% W
        for (k = 0, l = 0; k < currDim; k++) {
          for (j = 0; j < currDim; j++, l++) {
            tmp_qq2[l] = tmp_qq[k * q + j];
          } // j
        } // k

        // Add 1
        for (j = 0; j < currDim; j++) {
          tmp_qq2[j * currDim + j] += 1.0;  
        } // j

        F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf for spatial factors failed\n");}
        F77_NAME(dpotri)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
        if(info != 0){error("c++ error: dpotri for spatial factors failed\n");}

        F77_NAME(dsymv)(lower, &currDim, &one, tmp_qq2, &currDim, tmp_q, &inc, &zero, tmp_q2, &inc FCONE);

        F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf for spatial factors 2 failed\n");}
        
        mvrnorm(tmp_q, tmp_q2, tmp_qq2, currDim);
        F77_NAME(dcopy)(&currDim, tmp_q, &inc, &lambda[i], &N);
      } // i

      // Multiply Lambda %*% w[j] to get wStar. 
      for (j = 0; j < J; j++) {
        F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q], &inc, &zero, &wStar[j * N], &inc FCONE);
      } // j

      /********************************************************************
       *Update Latent Occupancy
       *******************************************************************/
      for (i = 0; i < N; i++) {
        // Compute detection probability 
        if (nObs == J) {
          for (r = 0; r < nObs; r++) {
            detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
            psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + wStar[zLongIndx[r] * N + i] + betaStarSites[i * J + zLongIndx[r]], zero, one); 
            piProd[zLongIndx[r] * N + i] = pow(1.0 - detProb[i * nObs + r], K[r]);
	    piProdWAIC[zLongIndx[r] * N + i] *= pow(detProb[i * nObs + r], y[r * N + i]);
	    piProdWAIC[zLongIndx[r] * N + i] *= pow(1.0 - detProb[i * nObs + r], K[r] - y[r * N + i]);
            ySum[zLongIndx[r] * N + i] = y[r * N + i]; 
          } // r
        } else {
          for (r = 0; r < nObs; r++) {
            detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
            if (tmp_JInt[zLongIndx[r]] == 0) {
              psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + wStar[zLongIndx[r] * N + i] + betaStarSites[i * J + zLongIndx[r]], zero, one); 
            }
            piProd[zLongIndx[r] * N + i] *= (1.0 - detProb[i * nObs + r]);
	    piProdWAIC[zLongIndx[r] * N + i] *= pow(detProb[i * nObs + r], y[r * N + i]);
	    piProdWAIC[zLongIndx[r] * N + i] *= pow(1.0 - detProb[i * nObs + r], 
			                            1.0 - y[r * N + i]);
            ySum[zLongIndx[r] * N + i] += y[r * N + i]; 
            tmp_JInt[zLongIndx[r]]++;
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
          tmp_JInt[j] = 0; 
        } // j
      } // i

      /********************************************************************
       *Save samples
       *******************************************************************/
      if (s >= nBurn) {
        thinIndx++;
        if (thinIndx == nThin) {
          F77_NAME(dcopy)(&pOcc, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, alphaComm, &inc, &REAL(alphaCommSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&pOcc, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, tauSqAlpha, &inc, &REAL(tauSqAlphaSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&pOccN, beta, &inc, &REAL(betaSamples_r)[sPost*pOccN], &inc); 
          F77_NAME(dcopy)(&pDetN, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDetN], &inc); 
          F77_NAME(dcopy)(&Nq, lambda, &inc, &REAL(lambdaSamples_r)[sPost*Nq], &inc); 
          F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[sPost*JN], &inc); 
          F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[sPost*JN], &inc); 
          F77_NAME(dcopy)(&Jq, w, &inc, &REAL(wSamples_r)[sPost*Jq], &inc); 
	  if (pDetRE > 0) {
            F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
            F77_NAME(dcopy)(&nDetREN, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetREN], &inc);
	  }
	  if (pOccRE > 0) {
            F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
            F77_NAME(dcopy)(&nOccREN, betaStar, &inc, &REAL(betaStarSamples_r)[sPost*nOccREN], &inc);
	  }
          F77_NAME(dcopy)(&JN, yWAIC, &inc, &REAL(likeSamples_r)[sPost*JN], &inc); 
          sPost++; 
          thinIndx = 0; 
        }
      }

      /********************************************************************
       * Report
       *******************************************************************/
      if (status == nReport){
        if(verbose){
          Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
      	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
      	  R_FlushConsole();
          #endif
      	}
        status = 0;
      }
      
      status++;

      R_CheckUserInterrupt();

    }

    if (verbose) {
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
    }

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
    SET_VECTOR_ELT(result_r, 8, lambdaSamples_r);
    SET_VECTOR_ELT(result_r, 9, wSamples_r); 
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


    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("lambda.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 10, mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 11, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 12, mkChar("alpha.star.samples")); 
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


