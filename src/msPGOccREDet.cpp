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

extern "C" {
  SEXP msPGOccREDet(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XpRE_r, 
	            SEXP lambdaP_r, SEXP pocc_r, SEXP pdet_r, 
	            SEXP pDetRE_r, SEXP J_r, SEXP nObs_r, SEXP K_r, SEXP N_r, 
	            SEXP nDetRE_r, SEXP nDetRELong_r, 
	            SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
	            SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
	            SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
	            SEXP sigmaSqPStarting_r, SEXP alphaStarStarting_r, 
	            SEXP zLongIndx_r, SEXP alphaStarIndx_r, 
	            SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	            SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	            SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, 
	            SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	            SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
	            SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, SEXP currChain_r, 
		    SEXP nChain_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, l, q, r, info, nProtect=0;
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
    int *XpRE = INTEGER(XpRE_r); 
    double *lambdaP = REAL(lambdaP_r); 
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
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
    int pDetRE = INTEGER(pDetRE_r)[0]; 
    int nDetRE = INTEGER(nDetRE_r)[0]; 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int J = INTEGER(J_r)[0];
    double *K = REAL(K_r); 
    int nObs = INTEGER(nObs_r)[0];
    int N = INTEGER(N_r)[0]; 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int nThin = INTEGER(nThin_r)[0]; 
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    int currChain = INTEGER(currChain_r)[0];
    int nChain = INTEGER(nChain_r)[0];
    int status = 0; 
    int thinIndx = 0;
    int sPost = 0;  

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
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
        Rprintf("Multispecies Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
	Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
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
    }

    /**********************************************************************
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pOccN = pOcc * N; 
    int pDetN = pDet * N; 
    int nObsN = nObs * N; 
    int nDetREN = nDetRE * N; 
    int JN = J * N;
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
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
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));

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
    double *alpha = (double *) R_alloc(pDetN, sizeof(double));   
    F77_NAME(dcopy)(&pDetN, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREN, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREN, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double));
    double *omegaOcc = (double *) R_alloc(J, sizeof(double));
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); 
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); 

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
    SEXP sigmaSqPSamples_r; 
    PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
    SEXP alphaStarSamples_r; 
    PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetREN, nPost)); nProtect++;
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
   
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObsN, sizeof(double)); 
    double *psi = (double *) R_alloc(JN, sizeof(double)); 
    zeros(psi, JN); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J);

    // For normal priors
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, SigmaBetaCommInvMuBeta, &inc); 	  
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, SigmaAlphaCommInvMuAlpha, &inc); 
    // Put community level variances in a pOcc x POcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level variances in a pDet x pDet matrix. 
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauSqAlpha[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObsN, sizeof(double)); 
    zeros(alphaStarObs, nObsN); 
    // Get sums of the current REs for each site/visit combo for all species
    for (i = 0; i < N; i++) {
      for (r = 0; r < nObs; r++) {
        for (l = 0; l < pDetRE; l++) {
          alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + XpRE[l * nObs + r]];
        }
      }
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int)); 
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE); 
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
        F77_NAME(dgemv)(ytran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, &beta[i], &N, &one, tmp_pOcc, &inc); 
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
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotri ABetaComm failed\n");}
      F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
      mvrnorm(betaComm, tmp_pOcc2, tmp_ppOcc, pOcc);
      /********************************************************************
       Update Community level Detection Coefficients
       *******************************************************************/
      /********************************
       * Compute b.alpha.comm
       *******************************/
       zeros(tmp_pDet, pDet); 
       for (i = 0; i < N; i++) {
         F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, TauAlphaInv, &pDet, &alpha[i], &N, &one, tmp_pDet, &inc); 
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
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
      F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotri AAlphaComm failed\n");}
      F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc);
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
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
      F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
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
      } // q
      F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
      F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info); 
      if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

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

       
      for (i = 0; i < N; i++) {  
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
        for (j = 0; j < J; j++) {
          omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N));
        } // j
        /********************************************************************
         *Update Detection Auxiliary Variables 
         *******************************************************************/
        // Note that all of the variables are sampled, but only those at 
        // locations with z[j] == 1 actually effect the results. 
        if (nObs == J) {
          for (r = 0; r < nObs; r++) {
            omegaDet[r] = rpg(K[r], F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
          } // r
        } else {
          for (r = 0; r < nObs; r++) {
            omegaDet[r] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
          } // r
        }
           
        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          kappaOcc[j] = z[j * N + i] - 1.0 / 2.0; 
        } // j
        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, kappaOcc, &inc, &zero, tmp_pOcc, &inc); 	 
        // TauBetaInv %*% betaComm + tmp_pOcc = tmp_pOcc
        F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, betaComm, &inc, &one, tmp_pOcc, &inc); 

        /********************************
         * Compute A.beta
         * *****************************/
        for(j = 0; j < J; j++){
          for(q = 0; q < pOcc; q++){
            tmp_JpOcc[q*J+j] = X[q*J+j]*omegaOcc[j];
          }
        }
        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc);
        for (q = 0; q < ppOcc; q++) {
          tmp_ppOcc[q] += TauBetaInv[q]; 
        } // q
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotri ABeta failed\n");}
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        mvrnorm(tmp_beta, tmp_pOcc2, tmp_ppOcc, pOcc);
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
        // to 0 if z == 0 and values go to kappaDet if z == 1
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
        
        F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc); 	  
        F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, TauAlphaInv, &pDet, alphaComm, &inc, &one, tmp_pDet, &inc); 
        /********************************
         * Compute A.alpha
         * *****************************/
        for (r = 0; r < nObs; r++) {
          for (q = 0; q < pDet; q++) {
            tmp_nObspDet[q*nObs + r] = Xp[q * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i];
          } // i
        } // j

        F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet);

        for (q = 0; q < ppDet; q++) {
          tmp_ppDet[q] += TauAlphaInv[q]; 
        } // q
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
        if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info); 
        if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc);
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDet);
        for (q = 0; q < pDet; q++) {
          alpha[q * N + i] = tmp_alpha[q];
        }

        /********************************************************************
         *Update Detection random effects
         *******************************************************************/
        // Update each individual random effect one by one. 
        for (l = 0; l < nDetRE; l++) {
          /********************************
           * Compute b.alpha.star
           *******************************/
          for (r = 0; r < nObs; r++) {
            tmp_nObs[r] = kappaDet[r] - (F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r] - alphaStar[i * nDetRE + l]) * omegaDet[r];
            // Only allow information to come from when z == 1.
            tmp_nObs[r] *= z[zLongIndx[r] * N + i]; 
          }
          F77_NAME(dgemv)(ytran, &nObs, &inc, &one, &lambdaP[l * nObs], &nObs, tmp_nObs, &inc, &zero, tmp_one, &inc); 
          /********************************
           * Compute A.alpha.star
           *******************************/
          for (r = 0; r < nObs; r++) {
            tmp_nObs[r] = lambdaP[l * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i]; 
          }
          tmp_0 = F77_NAME(ddot)(&nObs, tmp_nObs, &inc, &lambdaP[l * nObs], &inc); 
          tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
          tmp_0 = 1.0 / tmp_0; 
          alphaStar[i * nDetRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
        }
	zeros(&alphaStarObs[i * nObs], nObs); 
        // Update the RE sums for the current species
        for (r = 0; r < nObs; r++) {
          for (l = 0; l < pDetRE; l++) {
            alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + XpRE[l * nObs + r]]; 
          }
        }

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
        if (nObs == J) {
          for (r = 0; r < nObs; r++) {
            detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
            psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N), zero, one); 
            piProd[zLongIndx[r]] = pow(1.0 - detProb[i * nObs + r], K[r]);
            ySum[zLongIndx[r]] = y[r * N + i]; 
          } // r
        } else {
          for (r = 0; r < nObs; r++) {
            detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
            if (tmp_J[zLongIndx[r]] == 0) {
              psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N), zero, one); 
            }
            piProd[zLongIndx[r]] *= (1.0 - detProb[i * nObs + r]);
            ySum[zLongIndx[r]] += y[r * N + i]; 	
            tmp_J[zLongIndx[r]]++;
          } // r
        }

        // Compute occupancy probability 
        for (j = 0; j < J; j++) {
          psiNum = psi[j * N + i] * piProd[j]; 
          if (ySum[j] == zero) {
            z[j * N + i] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j * N + i])));           
          } else {
            z[j * N + i] = one; 
          }
          // Reset variables
          piProd[j] = one;
          ySum[j] = zero; 
          tmp_J[j] = 0; 
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
          F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[sPost*JN], &inc); 
          F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[sPost*JN], &inc); 
          F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
          F77_NAME(dcopy)(&nDetREN, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetREN], &inc);
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
    int nResultListObjs = 10;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauSqAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, zSamples_r);
    SET_VECTOR_ELT(result_r, 7, psiSamples_r);
    SET_VECTOR_ELT(result_r, 8, sigmaSqPSamples_r);
    SET_VECTOR_ELT(result_r, 9, alphaStarSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("sigma.sq.p.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, mkChar("alpha.star.samples")); 
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


