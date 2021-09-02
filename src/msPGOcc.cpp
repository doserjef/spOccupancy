#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
  SEXP msPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pocc_r, SEXP pdet_r, 
	       SEXP J_r, SEXP K_r, SEXP N_r, SEXP betaStarting_r, 
	       SEXP alphaStarting_r, SEXP zStarting_r, SEXP betaCommStarting_r, 
	       SEXP alphaCommStarting_r, SEXP tauBetaStarting_r, SEXP tauAlphaStarting_r, 
	       SEXP zLongIndx_r, SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	       SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauBetaA_r, 
	       SEXP tauBetaB_r, SEXP tauAlphaA_r, SEXP tauAlphaB_r, 
	       SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, q, r, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *upper = "U";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *rside = "R";
    char const *lside = "L";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *Xp = REAL(Xp_r);
    double *muBetaComm = REAL(muBetaComm_r); 
    double *muAlphaComm = REAL(muAlphaComm_r); 
    double *SigmaBetaCommInv = REAL(SigmaBetaComm_r); 
    double *SigmaAlphaCommInv = REAL(SigmaAlphaComm_r); 
    double *tauBetaA = REAL(tauBetaA_r); 
    double *tauBetaB = REAL(tauBetaB_r); 
    double *tauAlphaA = REAL(tauAlphaA_r); 
    double *tauAlphaB = REAL(tauAlphaB_r); 
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
    int J = INTEGER(J_r)[0];
    int *K = INTEGER(K_r); 
    int N = INTEGER(N_r)[0]; 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int nObs = 0;
    for (j = 0; j < J; j++) {
      nObs += K[j]; 
    } // j
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 
    double *z = REAL(zStarting_r); 

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
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Multi-species Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      Rprintf("Sampling ... \n");
    }

    /**********************************************************************
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pOccN = pOcc * N; 
    int pDetN = pDet * N; 
    int nObsN = nObs * N; 
    int JN = J * N;
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    double tmp_0; 
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
    double *tauBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(tauBetaStarting_r), &inc, tauBeta, &inc);
    double *alphaComm = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaCommStarting_r), &inc, alphaComm, &inc);
    double *tauAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(tauAlphaStarting_r), &inc, tauAlpha, &inc);
    double *beta = (double *) R_alloc(pOccN, sizeof(double));   
    F77_NAME(dcopy)(&pOccN, REAL(betaStarting_r), &inc, beta, &inc);
    double *alpha = (double *) R_alloc(pDetN, sizeof(double));   
    F77_NAME(dcopy)(&pDetN, REAL(alphaStarting_r), &inc, alpha, &inc);
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
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pOcc, nSamples)); nProtect++;
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = allocMatrix(REALSXP, pDet, nSamples)); nProtect++;
    SEXP tauBetaSamples_r; 
    PROTECT(tauBetaSamples_r = allocMatrix(REALSXP, pOcc, nSamples)); nProtect++; 
    SEXP tauAlphaSamples_r; 
    PROTECT(tauAlphaSamples_r = allocMatrix(REALSXP, pDet, nSamples)); nProtect++; 
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOccN, nSamples)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDetN, nSamples)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, JN, nSamples)); nProtect++; 
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, JN, nSamples)); nProtect++; 
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(INTSXP, nObsN, nSamples)); nProtect++; 
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
   
    // For latent occupancy
    double psiNum; 
    double psiNew; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
    double *psi = (double *) R_alloc(JN, sizeof(double)); 
    zeros(psi, JN); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); 
    int *yRep = (int *) R_alloc(nObsN, sizeof(int)); 

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
      TauBetaInv[i * pOcc + i] = tauBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level variances in a pDet x pDet matrix. 
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauAlpha[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

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
        tauBeta[q] = rigamma(tauBetaA[q] + N / 2.0, tauBetaB[q] + tmp_0); 
      } // q
      for (q = 0; q < pOcc; q++) {
        TauBetaInv[q * pOcc + q] = tauBeta[q]; 
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
        tauAlpha[q] = rigamma(tauAlphaA[q] + N / 2.0, tauAlphaB[q] + tmp_0); 
      } // q
      for (q = 0; q < pDet; q++) {
        TauAlphaInv[q * pDet + q] = tauAlpha[q]; 
      } // q
      F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
      F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info); 
      if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}
       
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
        for (r = 0; r < nObs; r++) {
          omegaDet[r] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N));
        } // i
           
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
        for (r = 0; r < nObs; r++) {
          kappaDet[r] = (y[r * N + i] - 1.0/2.0) * z[zLongIndx[r] * N + i];
        } // r
        
        F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, kappaDet, &inc, &zero, tmp_pDet, &inc); 	  
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
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
        for (r = 0; r < nObs; r++) {
          detProb[r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N), zero, one);
          if (tmp_J[zLongIndx[r]] == 0) {
            psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N), zero, one); 
          }
          piProd[zLongIndx[r]] *= (1.0 - detProb[r]);
          ySum[zLongIndx[r]] += y[r * N + i]; 	
          tmp_J[zLongIndx[r]]++;
        } // r
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

        /********************************************************************
         *Replicate data set for GoF
         *******************************************************************/
        for (r = 0; r < nObs; r++) {
          yRep[r * N + i] = rbinom(one, detProb[i] * z[zLongIndx[r] * N + i]);
          INTEGER(yRepSamples_r)[s * nObsN + r * N + i] = yRep[r * N + i]; 
        } // r
    } // i


     /********************************************************************
      *Save samples
      *******************************************************************/
      F77_NAME(dcopy)(&pOcc, betaComm, &inc, &REAL(betaCommSamples_r)[s*pOcc], &inc);
      F77_NAME(dcopy)(&pDet, alphaComm, &inc, &REAL(alphaCommSamples_r)[s*pDet], &inc);
      F77_NAME(dcopy)(&pOcc, tauBeta, &inc, &REAL(tauBetaSamples_r)[s*pOcc], &inc);
      F77_NAME(dcopy)(&pDet, tauAlpha, &inc, &REAL(tauAlphaSamples_r)[s*pDet], &inc);
      F77_NAME(dcopy)(&pOccN, beta, &inc, &REAL(betaSamples_r)[s*pOccN], &inc); 
      F77_NAME(dcopy)(&pDetN, alpha, &inc, &REAL(alphaSamples_r)[s*pDetN], &inc); 
      F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[s*JN], &inc); 
      F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[s*JN], &inc); 

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
    PutRNGstate();

    SEXP result_r, resultName_r;
    int nResultListObjs = 9;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, zSamples_r);
    SET_VECTOR_ELT(result_r, 7, psiSamples_r);
    SET_VECTOR_ELT(result_r, 8, yRepSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("y.rep.samples")); 
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


