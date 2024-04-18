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
	       SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, l, ll, q, r, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    // Sorted by data set, then visit, then site, then species. 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    // Sorted by parameter, then data set, site, visit.
    double *Xp = REAL(Xp_r);
    int *XRE = INTEGER(XRE_r); 
    // Load constants
    int N = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int nObs = INTEGER(consts_r)[2]; 
    int pOcc = INTEGER(consts_r)[3];
    int pOccRE = INTEGER(consts_r)[4];
    int nOccRE = INTEGER(consts_r)[5];
    int pDet = INTEGER(consts_r)[6];
    int nData = INTEGER(consts_r)[7];
    int nObsFull = INTEGER(consts_r)[8];
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *muAlphaComm = REAL(muAlphaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *sigmaAlphaComm = REAL(sigmaAlphaComm_r);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *tauSqAlphaA = REAL(tauSqAlphaA_r); 
    double *tauSqAlphaB = REAL(tauSqAlphaB_r); 
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    int *pDetLong = INTEGER(pDetLong_r);
    int *NLong = INTEGER(NLong_r);
    int *nObsLong = INTEGER(nObsLong_r);
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaCommIndx = INTEGER(alphaCommIndx_r);
    int *dataIndx = INTEGER(dataIndx_r);
    int *alphaSpIndx = INTEGER(alphaSpIndx_r);
    int *alphaDatIndx = INTEGER(alphaDatIndx_r);
    int *alphaPDetIndx = INTEGER(alphaPDetIndx_r);
    int *spLongIndx = INTEGER(spLongIndx_r);
    int *spDatLongIndx = INTEGER(spDatLongIndx_r);
    int *obsLongIndx = INTEGER(obsLongIndx_r);
    int *spSiteIndx = INTEGER(spSiteIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r);
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int status = 0; 
    int thinIndx = 0;
    int sPost = 0;  

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
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
        Rprintf("Integrated Multispecies Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Integrating %i occupancy data sets.\n\n", nData); 
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
	Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n", nPost * nChain); 
#ifdef _OPENMP
        Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
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
    int nOccREN = nOccRE * N; 
    int nAlpha = 0;
    for (i = 0; i < nData; i++) {
      nAlpha += NLong[i] * pDetLong[i];
    }
    int JN = J * N;
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    int JpOccRE = J * pOccRE; 
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double)); zeros(tmp_ppDet, ppDet);
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
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    // For conecting to different data sets
    int stAlpha = 0;
    int stAlphaComm = 0;

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *betaComm = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    double *alphaComm = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaCommStarting_r), &inc, alphaComm, &inc);
    double *tauSqAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(tauSqAlphaStarting_r), &inc, tauSqAlpha, &inc);
    double *beta = (double *) R_alloc(pOccN, sizeof(double));   
    F77_NAME(dcopy)(&pOccN, REAL(betaStarting_r), &inc, beta, &inc);
    double *alpha = (double *) R_alloc(nAlpha, sizeof(double));   
    F77_NAME(dcopy)(&nAlpha, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Occupancy random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nOccREN, sizeof(double)); 
    F77_NAME(dcopy)(&nOccREN, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
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
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++; 
    SEXP tauSqAlphaSamples_r; 
    PROTECT(tauSqAlphaSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++; 
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOccN, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, nAlpha, nPost)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = Rf_allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nOccREN, nPost)); nProtect++;
    }
    // TODO: eventually make this separate across different data sets.
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++;
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObsFull, sizeof(double)); 
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
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    // Get size of vector
    int currSize = 0;
    // Index of starting prior values
    int *alphaCommSigmaIndx = (int *) R_alloc(nData, sizeof(int));
    int *alphaCommMuIndx = (int *) R_alloc(nData, sizeof(int));
    int tmp0 = 0;
    int tmp02 = 0;
    for (q = 0; q < nData; q++) {
      currSize += pDetLong[q] * pDetLong[q];  
      alphaCommSigmaIndx[q] = tmp0; 
      tmp0 += pDetLong[q] * pDetLong[q]; 
      alphaCommMuIndx[q] = tmp02; 
      tmp02 += pDetLong[q]; 
    }
    double *SigmaAlphaCommInv = (double *) R_alloc(currSize, sizeof(double)); zeros(SigmaAlphaCommInv, currSize); 
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    // Fill SigmaAlphaComm
    for (q = 0, j = 0; q < nData; q++) {
      for (i = 0; i < pDetLong[q]; j++, i++) {
        SigmaAlphaCommInv[alphaCommSigmaIndx[q] + i * pDetLong[q] + i] = sigmaAlphaComm[j]; 
      } // i
      F77_NAME(dpotrf)(lower, &pDetLong[q], &SigmaAlphaCommInv[alphaCommSigmaIndx[q]], 
		       &pDetLong[q], &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
      F77_NAME(dpotri)(lower, &pDetLong[q], &SigmaAlphaCommInv[alphaCommSigmaIndx[q]], &pDetLong[q], &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
      F77_NAME(dsymv)(lower, &pDetLong[q], &one, &SigmaAlphaCommInv[alphaCommSigmaIndx[q]], &pDetLong[q], &muAlphaComm[alphaCommMuIndx[q]], &inc, &zero, 
                     &SigmaAlphaCommInvMuAlpha[alphaCommMuIndx[q]], &inc FCONE);
    } // q


    // Put community level variances in a pOcc x POcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = 1.0 / tauSqBeta[i]; 
    } // i
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = 1.0 / tauSqAlpha[i]; 
    } // i

    // Get starting location for each data set in Xp
    int *startXP = (int *) R_alloc(nData, sizeof(int));
    startXP[0] = 0;
    for (q = 1; q < nData; q++) {
      startXP[q] = startXP[q - 1] + nObsLong[q - 1];
    }

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
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pOccRE, sizeof(int)); 
    for (l = 0; l < pOccRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nOccRE); 
    }
    
    GetRNGstate();
    
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
      for (q = 0; q < pOcc; q++) {
        tmp_pOcc[q] += SigmaBetaCommInvMuBeta[q];  
      } // j

      /********************************
       Compute A.beta.comm
       *******************************/
      // Note that you really only need to do this for the off-diagonal elements, 
      // because the off-diagonal values just go to zero. This is what you do for the 
      // detection values. 
      for (q = 0; q < ppOcc; q++) {
        tmp_ppOcc[q] = SigmaBetaCommInv[q] + N * TauBetaInv[q]; 
      }
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri ABetaComm failed\n");}
      F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
      mvrnorm(betaComm, tmp_pOcc2, tmp_ppOcc, pOcc);

      /********************************************************************
       Update Community level Detection Coefficients
       *******************************************************************/
      /********************************
       * Compute b.alpha.comm
       *******************************/
       zeros(tmp_pDet, pDet); 
       for (i = 0; i < N; i++) {
         for (q = 0; q < nAlpha; q++) {
           if (alphaSpIndx[q] == i) {
             tmp_pDet[alphaPDetIndx[q]] += alpha[q] / tauSqAlpha[alphaPDetIndx[q]];   
	   }
	 }
       } // i
      /********************************
       * Compute A.alpha.comm
       *******************************/
      // NOTE: this differs from previously implementations because now the 
      //       number of species (aka the number of random effects) that informs
      //       each community-level effect depends on which data set the value 
      //       comes from. 
      zeros(tmp_ppDet, ppDet);
      for (q = 0; q < pDet; q++) {
        tmp_ppDet[q * pDet + q] = (1.0 / sigmaAlphaComm[q]) + 
                                  NLong[alphaCommIndx[q]] * TauAlphaInv[q * pDet + q]; 
      }
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf AAlphaComm1 failed\n");}
      F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri AAlphaComm failed\n");}
      F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf AAlphaComm2 failed\n");}
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
      } // q
      F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}

      /********************************************************************
       Update Community Detection Variance Parameter
      ********************************************************************/
      for (q = 0; q < pDet; q++) {
        tmp_0 = 0.0;  
	for (l = 0; l < nAlpha; l++) {
          if (alphaPDetIndx[l] == q) {
            tmp_0 += (alpha[l] - alphaComm[q]) * (alpha[l] - alphaComm[q]);
	  }
	}
        tmp_0 *= 0.5;
        tauSqAlpha[q] = rigamma(tauSqAlphaA[q] + NLong[alphaCommIndx[q]] / 2.0, tauSqAlphaB[q] + tmp_0); 
      } // q
      for (q = 0; q < pDet; q++) {
        TauAlphaInv[q * pDet + q] = tauSqAlpha[q]; 
      } // q
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
       * Species-specific occurrence parameters
       *******************************************************************/
      for (i = 0; i < N; i++) {  
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
        for (j = 0; j < J; j++) {
          omegaOcc[j] = 0.0;
          if (spSiteIndx[j * N + i] == 1) {
            omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + betaStarSites[i * J + j]);
	  }
        } // j
        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          kappaOcc[j] = z[j * N + i] - 1.0 / 2.0; 
          tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * betaStarSites[i * J + j]; 
	  if (spSiteIndx[j * N + i] == 0) {
            tmp_J1[j] = 0.0;
	  }
        } // j
        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
        // TauBetaInv %*% betaComm + tmp_pOcc = tmp_pOcc
        F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, betaComm, &inc, &one, tmp_pOcc, &inc FCONE); 

        /********************************
         * Compute A.beta
         * *****************************/
	// This is fine because it will go to zero if spSiteIndx == 0
        for(j = 0; j < J; j++){
          for(q = 0; q < pOcc; q++){
            tmp_JpOcc[q*J+j] = X[q*J+j]*omegaOcc[j];
          }
        }
        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
        for (q = 0; q < ppOcc; q++) {
          tmp_ppOcc[q] += TauBetaInv[q]; 
        } // q
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf ABeta failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri ABeta failed\n");}
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        mvrnorm(tmp_beta, tmp_pOcc2, tmp_ppOcc, pOcc);
        for (q = 0; q < pOcc; q++) {
          beta[q * N + i] = tmp_beta[q]; 
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
              if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l] && (spSiteIndx[j * N + i] == 1)) {
                tmp_02 = 0.0;
                for (ll = 0; ll < pOccRE; ll++) {
                  tmp_02 += betaStar[i * nOccRE + betaStarLongIndx[ll * J + j]];
                } 
                tmp_one[0] += kappaOcc[j] - (F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + tmp_02 - betaStar[i * nOccRE + l]) * omegaOcc[j];
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
         *Update Detection Auxiliary Variables and regression coefficients
         *******************************************************************/
        for (q = 0; q < nData; q++) {
          zeros(tmp_nObs, nObs);
	  zeros(omegaDet, nObs);
	  zeros(tmp_pDet, pDet);
          currSize = 0;
	  stAlphaComm = which(q, alphaCommIndx, pDet);
	  stAlpha = which(q, alphaDatIndx, nAlpha);
          for (r = 0; r < nObsFull; r++) {
            if ((spLongIndx[r] == i) && (dataIndx[r] == q)) {
              // First multiply kappDet * the current occupied values, such that values go 
              // to 0 if z == 0 and values go to kappaDet if z == 1
              kappaDet[obsLongIndx[r]] = (y[r] - 1.0/2.0) * z[zLongIndx[r] * N + i];
              tmp_nObs[currSize] = kappaDet[obsLongIndx[r]];
	      if (z[zLongIndx[r] * N + i] == 1.0) {
                omegaDet[currSize] = rpg(1.0, F77_NAME(ddot)(&pDetLong[q], &Xp[obsLongIndx[r]],
				       	                     &nObs, &alpha[stAlpha + spDatLongIndx[r]], 
							     &NLong[q]));

	      }
              currSize++;
            }
          } // r
          /********************************
           * Compute b.alpha
           *******************************/
	  // If the current species is sampled in the current data set. 
          if (currSize > 0) { 
            F77_NAME(dgemv)(ytran, &nObsLong[q], &pDetLong[q], &one, &Xp[startXP[q]], 
			    &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
	    for (r = 0; r < pDetLong[q]; r++) {
              tmp_pDet[r] += alphaComm[stAlphaComm + r] / tauSqAlpha[stAlphaComm + r];
	    }
            /********************************
             * Compute A.alpha
             * *****************************/
            zeros(tmp_nObspDet, nObspDet);
	    zeros(tmp_ppDet, ppDet);
            for (r = 0; r < nObsLong[q]; r++) {
              for (l = 0; l < pDetLong[q]; l++) {
                // Note that omegaDet[r] = 0.0 if z != 1, so you don't need to explicitly multiply 
	        // by z here. 
                tmp_nObspDet[l*nObsLong[q] + r] = Xp[startXP[q] + l * nObs + r] * omegaDet[r];
              } // l
            } // r

            F77_NAME(dgemm)(ytran, ntran, &pDetLong[q], &pDetLong[q], &nObsLong[q], 
	         	   &one, &Xp[startXP[q]], &nObs, tmp_nObspDet, &nObsLong[q], &zero, 
	         	   tmp_ppDet, &pDetLong[q] FCONE FCONE);

            for (r = 0; r < pDetLong[q]; r++) {
              tmp_ppDet[r * pDetLong[q] + r] += 1.0 / tauSqAlpha[stAlphaComm + r]; 
            } // r
            F77_NAME(dpotrf)(lower, &pDetLong[q], tmp_ppDet, &pDetLong[q], &info FCONE); 
            if(info != 0){Rf_error("c++ error: dpotrf A.alpha failed\n");}
            F77_NAME(dpotri)(lower, &pDetLong[q], tmp_ppDet, &pDetLong[q], &info FCONE); 
            if(info != 0){Rf_error("c++ error: dpotri A.alpha failed\n");}
            F77_NAME(dsymv)(lower, &pDetLong[q], &one, tmp_ppDet, &pDetLong[q], tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
            F77_NAME(dpotrf)(lower, &pDetLong[q], tmp_ppDet, &pDetLong[q], &info FCONE); 
            if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
            mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDetLong[q]);
	    currSize = 0; 
	    for (r = 0; r < nAlpha; r++) {
              if ((alphaSpIndx[r] == i) && alphaDatIndx[r] == q) {
                alpha[r] = tmp_alpha[currSize]; 
	        currSize++;
	      } 
	    }
          } 
        } 
      } // i (species)

      /********************************************************************
       *Update Latent Occupancy
       *******************************************************************/
      // Compute detection probability 
      for (r = 0; r < nObsFull; r++) {
        stAlpha = which(dataIndx[r], alphaDatIndx, nAlpha);
        detProb[r] = logitInv(F77_NAME(ddot)(&pDetLong[dataIndx[r]], &Xp[obsLongIndx[r]],
        			             &nObs, &alpha[stAlpha + spDatLongIndx[r]], 
        				     &NLong[dataIndx[r]]), zero, one);
        // if (tmp_J[zLongIndx[obsLongIndx[r]]] == 0) {
          psi[zLongIndx[r] * N + spLongIndx[r]] = 
            logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], 
        			    &J, &beta[spLongIndx[r]], &N) + 
        		            betaStarSites[spLongIndx[r] * J + zLongIndx[r]], zero, one); 
        // }
        piProd[zLongIndx[r] * N + spLongIndx[r]] *= (1.0 - detProb[r]);
        piProdWAIC[zLongIndx[r] * N + spLongIndx[r]] *= pow(detProb[r], y[r]);
        piProdWAIC[zLongIndx[r] * N + spLongIndx[r]] *= pow(1.0 - detProb[r], 
        		                                                 1.0 - y[r]);
        ySum[zLongIndx[r] * N + spLongIndx[r]] += y[r]; 
        // tmp_J[zLongIndx[obsLongIndx[r]]]++;
      } // r
      // Compute occupancy probability and the integrated likelihood for WAIC
      for (i = 0; i < N; i++) {
        for (j = 0; j < J; j++) {
          psiNum = psi[j * N + i] * piProd[j * N + i]; 
          if (ySum[j * N + i] == zero) {
            if (spSiteIndx[j * N + i] == 1) {
              z[j * N + i] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j * N + i])));
              yWAIC[j * N + i] = (1.0 - psi[j * N + i]) + psi[j * N + i] * piProdWAIC[j * N + i];
            } else {
              psi[j * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + 
                                                       betaStarSites[i * J + j], zero, one);
              z[j * N + i] = rbinom(one, psi[j * N + i]);
	      yWAIC[j * N + i] = NA_REAL;
            }
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
      }

      /********************************************************************
         *Save samples
      ********************************************************************/
      if (s >= nBurn) {
        thinIndx++; 
        if (thinIndx == nThin) {
          F77_NAME(dcopy)(&pOcc, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, alphaComm, &inc, &REAL(alphaCommSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&pOcc, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, tauSqAlpha, &inc, &REAL(tauSqAlphaSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&pOccN, beta, &inc, &REAL(betaSamples_r)[sPost*pOccN], &inc); 
          F77_NAME(dcopy)(&nAlpha, alpha, &inc, &REAL(alphaSamples_r)[sPost*nAlpha], &inc); 
          F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[sPost*JN], &inc); 
          F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[sPost*JN], &inc); 
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

    SEXP result_r, resultName_r;
    int nResultListObjs = 9;
    // if (pDetRE > 0) {
    //   nResultListObjs += 2; 
    // }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauSqAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, zSamples_r);
    SET_VECTOR_ELT(result_r, 7, psiSamples_r);
    SET_VECTOR_ELT(result_r, 8, likeSamples_r);
    // if (pDetRE > 0) {
    //   SET_VECTOR_ELT(result_r, 9, sigmaSqPSamples_r);
    //   SET_VECTOR_ELT(result_r, 10, alphaStarSamples_r);
    // }
    if (pOccRE > 0) {
      // if (pDetRE > 0) {
      //   tmp_0 = 11;
      // } else {
      //   tmp_0 = 9;
      // }
      tmp_0 = 9;
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("like.samples")); 
    // if (pDetRE > 0) {
    //   SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("sigma.sq.p.samples")); 
    //   SET_VECTOR_ELT(resultName_r, 10, Rf_mkChar("alpha.star.samples")); 
    // }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, Rf_mkChar("beta.star.samples")); 
    }
   
    Rf_namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


