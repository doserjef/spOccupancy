#define USE_FC_LEN_T
#include <string>
#include "util.h"

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
  SEXP postHocLM(SEXP y_r, SEXP X_r, SEXP XRE_r, 
                 SEXP consts_r, SEXP nRELong_r, 
                 SEXP betaStarting_r, SEXP tauSqStarting_r,  
                 SEXP sigmaSqStarting_r, SEXP betaStarStarting_r, 
                 SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	         SEXP muBeta_r, SEXP SigmaBeta_r, 
		 SEXP tauSqA_r, SEXP tauSqB_r, 
                 SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nSamples_r, 
                 SEXP verbose_r, SEXP nReport_r, 
	         SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, r, l, ll, info, nProtect=0;
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
    int *XRE = INTEGER(XRE_r); 
    // Load constants
    int N = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1]; 
    int p = INTEGER(consts_r)[2];
    int pRE = INTEGER(consts_r)[3];
    int nRE = INTEGER(consts_r)[4];
    int pp = p * p; 
    double *muBeta = (double *) R_alloc(p, sizeof(double));   
    F77_NAME(dcopy)(&p, REAL(muBeta_r), &inc, muBeta, &inc);
    double *SigmaBetaInv = (double *) R_alloc(pp, sizeof(double));   
    F77_NAME(dcopy)(&pp, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *sigmaSqA = REAL(sigmaSqA_r); 
    double *sigmaSqB = REAL(sigmaSqB_r); 
    double *tauSqA = REAL(tauSqA_r); 
    double *tauSqB = REAL(tauSqB_r); 
    int *nRELong = INTEGER(nRELong_r); 
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nSamples = INTEGER(nSamples_r)[0];
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 

    /**********************************************************************
     * Print Information 
     * *******************************************************************/
    if(verbose){
      if (currChain == 1) {
        Rprintf("----------------------------------------\n");
        Rprintf("\tModel description\n");
        Rprintf("----------------------------------------\n");
        Rprintf("Post-hoc linear model fit with %i observations.\n\n", N);
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nSamples * nChain); 
      }
      Rprintf("----------------------------------------\n");
      Rprintf("\tChain %i\n", currChain);
      Rprintf("----------------------------------------\n");
      Rprintf("Sampling ... \n");
    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Occupancy covariates
    double *beta = (double *) R_alloc(p, sizeof(double));   
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSq = (double *) R_alloc(pRE, sizeof(double)); 
    F77_NAME(dcopy)(&pRE, REAL(sigmaSqStarting_r), &inc, sigmaSq, &inc); 
    // Residual variance
    double *tauSq = (double *) R_alloc(1, sizeof(double)); 
    F77_NAME(dcopy)(&inc, REAL(tauSqStarting_r), &inc, tauSq, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nRE, sizeof(double)); 
    F77_NAME(dcopy)(&nRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    SEXP tauSqSamples_r;
    PROTECT(tauSqSamples_r = allocMatrix(REALSXP, 1, nSamples)); nProtect++;
    SEXP yHatSamples_r; 
    PROTECT(yHatSamples_r = allocMatrix(REALSXP, N, nSamples)); nProtect++; 
    // Occurrence random effects
    SEXP sigmaSqSamples_r; 
    SEXP betaStarSamples_r; 
    if (pRE > 0) {
      PROTECT(sigmaSqSamples_r = allocMatrix(REALSXP, pRE, nSamples)); nProtect++;
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nRE, nSamples)); nProtect++;
    }
    
    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int Np = N * p; 
    int NpRE = N * pRE;
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); 
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    int *tmp_N = (int *) R_alloc(N, sizeof(int));
    for (j = 0; j < N; j++) {
      tmp_N[j] = 0; 
    }
    double *tmp_Np = (double *) R_alloc(Np, sizeof(double));
    double *tmp_N1 = (double *) R_alloc(N, sizeof(double));
   
    double *yHat = (double *) R_alloc(N, sizeof(double));
    zeros(yHat, N); 
    double *yCurr = (double *) R_alloc(N, sizeof(double));
    zeros(yCurr, N);

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &p, SigmaBetaInv, &p, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &p, SigmaBetaInv, &p, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(p, sizeof(double)); 
    F77_NAME(dsymv)(lower, &p, &one, SigmaBetaInv, &p, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);
    double *XtX = (double *) R_alloc(pp, sizeof(double));
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &N, &one, X, &N, 
		    X, &N, &zero, XtX, &p FCONE FCONE);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the random effects
    double *betaStarSites = (double *) R_alloc(N, sizeof(double)); 
    zeros(betaStarSites, N); 
    int *betaStarLongIndx = (int *) R_alloc(NpRE, sizeof(int));
    // Initial sums
    for (j = 0; j < N; j++) {
      for (l = 0; l < pRE; l++) {
        betaStarLongIndx[l * N + j] = which(XRE[l * N + j], betaLevelIndx, nRE);
        betaStarSites[j] += betaStar[betaStarLongIndx[l * N + j]];
      }
    }

    int *betaStarStart = (int *) R_alloc(pRE, sizeof(int)); 
    for (l = 0; l < pRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nRE); 
    }

    GetRNGstate(); 

    for (s = 0; s < nSamples; s++) {
      /********************************************************************
       *Extract current y values
       *******************************************************************/
      for (j = 0; j < N; j++) {
        yCurr[j] = y[s * N + j];
      } // j

      /********************************************************************
       *Update Occupancy Regression Coefficients
       *******************************************************************/
      for (j = 0; j < N; j++) {
        tmp_N1[j] = (yCurr[j] - betaStarSites[j]) / tauSq[0]; 
      } // j

      /********************************
       * Compute b.beta
       *******************************/
      F77_NAME(dgemv)(ytran, &N, &p, &one, X, &N, tmp_N1, &inc, &zero, tmp_p, &inc FCONE); 	 
      for (j = 0; j < p; j++) {
        tmp_p[j] += SigmaBetaInvMuBeta[j]; 
      } // j 

      /********************************
       * Compute A.beta
       * *****************************/
      for (i = 0; i < pp; i++) {
         tmp_pp[i] = XtX[i] / tauSq[0];
         tmp_pp[i] += SigmaBetaInv[i];
      }

      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); 
      if(info != 0){error("c++ error: dpotri here failed\n");}
      F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      mvrnorm(beta, tmp_p2, tmp_pp, p);
      
      /********************************************************************
       *Update Residual Variance
       *******************************************************************/
      for (j = 0; j < N; j++) {
        tmp_N1[j] = yCurr[j] - betaStarSites[j] - F77_NAME(ddot)(&p, &X[j], &N, beta, &inc);
      }
      tauSq[0] = rigamma(tauSqA[0] + N / 2.0, tauSqB[0] + 0.5 * 
			                      F77_NAME(ddot)(&N, tmp_N1, &inc, tmp_N1, 
							     &inc));
      
      /********************************************************************
       *Update Occupancy random effects variance
       *******************************************************************/
      for (l = 0; l < pRE; l++) {
        tmp_0 = F77_NAME(ddot)(&nRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
        tmp_0 *= 0.5; 
        sigmaSq[l] = rigamma(sigmaSqA[l] + nRELong[l] / 2.0, sigmaSqB[l] + tmp_0); 
      }

      /********************************************************************
       *Update Occupancy random effects
       *******************************************************************/
      if (pRE > 0) {
        // Update each individual random effect one by one. 
        for (l = 0; l < nRE; l++) {
          /********************************
           * Compute b.beta.star
           *******************************/
          zeros(tmp_one, inc);
          tmp_0 = 0.0;	      
          // Only allow information to come from when XRE == betaLevelIndx[l]. 
          // aka information only comes from the sites with any given level 
          // of a random effect. 
          for (j = 0; j < N; j++) {
            if (XRE[betaStarIndx[l] * N + j] == betaLevelIndx[l]) {
              tmp_02 = 0.0;     
	      for (ll = 0; ll < pRE; ll++) {
                tmp_02 += betaStar[betaStarLongIndx[ll * N + j]];
	      }
              tmp_one[0] += (yCurr[j] - F77_NAME(ddot)(&p, &X[j], &N, beta, &inc) - 
        		    tmp_02 + betaStar[l]) / tauSq[0];
              tmp_0 += 1.0/tauSq[0];
            }
          }
          /********************************
           * Compute A.beta.star
           *******************************/
          tmp_0 += 1.0 / sigmaSq[betaStarIndx[l]]; 
          tmp_0 = 1.0 / tmp_0; 
          betaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
        }
      
        // Update the RE sums for the current species
        zeros(betaStarSites, N);
        for (j = 0; j < N; j++) {
          for (l = 0; l < pRE; l++) {
            betaStarSites[j] += betaStar[betaStarLongIndx[l * N + j]];
          }
        }
      }

      /********************************************************************
       *Get predicted mean values for Bayes R2. 
       *******************************************************************/
      for (j = 0; j < N; j++) {
        yHat[j] = F77_NAME(ddot)(&p, &X[j], &N, beta, &inc) + betaStarSites[j];
      }

      /********************************************************************
       *Save samples
       *******************************************************************/
      F77_NAME(dcopy)(&p, beta, &inc, &REAL(betaSamples_r)[s*p], &inc);
      F77_NAME(dcopy)(&inc, tauSq, &inc, &REAL(tauSqSamples_r)[s], &inc);
      F77_NAME(dcopy)(&N, yHat, &inc, &REAL(yHatSamples_r)[s*N], &inc); 
      if (pRE > 0) {
        F77_NAME(dcopy)(&pRE, sigmaSq, &inc, 
            	    &REAL(sigmaSqSamples_r)[s*pRE], &inc);
        F77_NAME(dcopy)(&nRE, betaStar, &inc, 
            	    &REAL(betaStarSamples_r)[s*nRE], &inc);
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
    int nResultListObjs = 3;
    if (pRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, tauSqSamples_r);
    SET_VECTOR_ELT(result_r, 2, yHatSamples_r); 
    if (pRE > 0) {
      SET_VECTOR_ELT(result_r, 3, sigmaSqSamples_r);
      SET_VECTOR_ELT(result_r, 4, betaStarSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("tau.sq.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("y.hat.samples")); 
    if (pRE > 0) {
      SET_VECTOR_ELT(resultName_r, 3, mkChar("sigma.sq.samples")); 
      SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.star.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

