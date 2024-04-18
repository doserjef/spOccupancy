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
  SEXP lfJSDM(SEXP y_r, SEXP X_r, SEXP XRE_r, SEXP consts_r, SEXP nOccRELong_r, 
		 SEXP betaStarting_r, SEXP betaCommStarting_r, 
		 SEXP tauSqBetaStarting_r, SEXP lambdaStarting_r, 
		 SEXP sigmaSqPsiStarting_r, SEXP betaStarStarting_r, 
		 SEXP betaStarIndx_r, SEXP betaLevelIndx_r,  
		 SEXP muBetaComm_r, SEXP SigmaBetaComm_r, 
		 SEXP tauSqBetaA_r, SEXP tauSqBetaB_r, 
		 SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r,
		 SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		 SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, h, l, info, nProtect=0, ii, ll;    

    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    // Sorted by then by site, then by species within site. 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int *XRE = INTEGER(XRE_r);
    // Load constants
    int N = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int pOcc = INTEGER(consts_r)[2];
    int pOccRE = INTEGER(consts_r)[3];
    int nOccRE = INTEGER(consts_r)[4];
    int q = INTEGER(consts_r)[5]; 
    int ppOcc = pOcc * pOcc; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
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
        Rprintf("Latent Factor JSDM with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
	if (q > 0) {
          Rprintf("Using %i latent factors.\n\n", q);
	} else {
          Rprintf("Assuming no residual species correlations.\n\n");
	}
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
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    /**********************************************************************
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pOccN = pOcc * N; 
    int nOccREN = nOccRE * N; 
    int Jq = J * q;
    int qq = q * q;
    int JN = J * N;
    int Nq = N * q;
    int JpOcc = J * pOcc; 
    int JpOccRE = J * pOccRE; 
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_beta = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    int *tmp_JInt = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_JInt[j] = 0; 
    }
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J, J);
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_qq = (double *) R_alloc(qq, sizeof(double));
    double *tmp_q = (double *) R_alloc(q, sizeof(double));
    double *tmp_q2 = (double *) R_alloc(q, sizeof(double));
    double *tmp_qq2 = (double *) R_alloc(qq, sizeof(double));
    double *tmp_Jq = (double *) R_alloc(Jq, sizeof(double));
    double *tmp_Nq = (double *) R_alloc(Nq, sizeof(double));
    double *tmp_N = (double *) R_alloc(N, sizeof(double));
    int currDim = 0;

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Community level
    double *betaComm = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    // Species level
    double *beta = (double *) R_alloc(pOccN, sizeof(double));   
    F77_NAME(dcopy)(&pOccN, REAL(betaStarting_r), &inc, beta, &inc);
    // Occurrence random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent effects
    double *w = (double *) R_alloc(Jq, sizeof(double)); zeros(w, Jq);
    // Latent spatial factors
    double *lambda = (double *) R_alloc(Nq, sizeof(double));
    F77_NAME(dcopy)(&Nq, REAL(lambdaStarting_r), &inc, lambda, &inc);
    // Latent occurrence random effects
    double *betaStar = (double *) R_alloc(nOccREN, sizeof(double)); 
    F77_NAME(dcopy)(&nOccREN, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Auxiliary variables
    double *omegaOcc = (double *) R_alloc(JN, sizeof(double)); zeros(omegaOcc, JN);
    double *kappaOcc = (double *) R_alloc(JN, sizeof(double)); zeros(kappaOcc, JN);
    // Need this for all species
    double *yStar = (double *) R_alloc(JN, sizeof(double));

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), pOcc * nPost);
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++; 
    zeros(REAL(tauSqBetaSamples_r), pOcc * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOccN, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOccN * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    // Factor model parameters
    SEXP lambdaSamples_r; 
    PROTECT(lambdaSamples_r = Rf_allocMatrix(REALSXP, Nq, nPost)); nProtect++;
    zeros(REAL(lambdaSamples_r), Nq * nPost);
    SEXP wSamples_r; 
    PROTECT(wSamples_r = Rf_allocMatrix(REALSXP, Jq, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), Jq * nPost);
    // Fitted values
    double *z = (double *) R_alloc(JN, sizeof(double)); zeros(z, JN);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), JN * nPost);
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
    double *psi = (double *) R_alloc(JN, sizeof(double)); 
    zeros(psi, JN); 
    double *yWAIC = (double *) R_alloc(JN, sizeof(double)); zeros(yWAIC, JN);


    // For normal community-level priors
    // Occurrence coefficients
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    // dgemv computes linear combinations of different variables. 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Put community level occurrence variances in a pOcc x pOcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}

    /**********************************************************************
     * Prep for random effects (if they exist)
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
    /**********************************************************************
     Set up factor model parameters
     * *******************************************************************/
    // Species-level latent effects
    double *wStar = (double *) R_alloc(JN, sizeof(double)); zeros(wStar, JN);
    // Multiply Lambda %*% w[j] to get wStar. 
    if (q > 0) {
      for (j = 0; j < J; j++) {
        F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q], &inc, &zero, &wStar[j * N], &inc FCONE);
      }
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
      if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}

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
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          kappaOcc[j * N + i] = y[j * N + i] - 1.0 / 2.0; 
          tmp_J1[j] = kappaOcc[j * N + i] - omegaOcc[j * N + i] * 
		      (wStar[j * N + i] + betaStarSites[i * J + j]); 
	  // For later
	  yStar[j * N + i] = kappaOcc[j * N + i] / omegaOcc[j * N + i];
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
                tmp_one[0] += kappaOcc[j * N + i] - (F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + tmp_02 - betaStar[i * nOccRE + l] + wStar[j * N + i]) * omegaOcc[j * N + i];
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
              betaStarSites[i * J + j] += betaStar[i * nOccRE + betaStarLongIndx[l * J + j]];
            } //l
          } // j
	}
      } // i
      /********************************************************************
       *Update latent effects (w) 
       *******************************************************************/
      if (q > 0) {
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
          if (info != 0){Rf_error("c++ error: dpotrf var failed\n");}
          F77_NAME(dpotri)(lower, &q, var, &q, &info FCONE);
          if (info != 0){Rf_error("c++ error: dpotri var failed\n");}

          // mu
          for (k = 0; k < N; k++) {
            tmp_N[k] = (yStar[ii * N + k] - F77_NAME(ddot)(&pOcc, &X[ii], &J, &beta[k], &N) - betaStarSites[k * J + ii]) * omegaOcc[ii * N + k];
          } // k

          F77_NAME(dgemv)(ytran, &N, &q, &one, lambda, &N, tmp_N, &inc, &zero, mu, &inc FCONE);

          F77_NAME(dsymv)(lower, &q, &one, var, &q, mu, &inc, &zero, tmp_N, &inc FCONE);

          F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf var 2 failed\n");}

          mvrnorm(&w[ii * q], tmp_N, var, q);

        } // ii
      }
      /********************************************************************
       *Update spatial factors (lambda)
       *******************************************************************/
      if (q > 0) {
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
            tmp_J[j] = yStar[j * N + i] - F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) -
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
          if(info != 0){Rf_error("c++ error: dpotrf for spatial factors failed\n");}
          F77_NAME(dpotri)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri for spatial factors failed\n");}

          F77_NAME(dsymv)(lower, &currDim, &one, tmp_qq2, &currDim, tmp_q, &inc, &zero, tmp_q2, &inc FCONE);

          F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf for spatial factors 2 failed\n");}
          
          mvrnorm(tmp_q, tmp_q2, tmp_qq2, currDim);
          F77_NAME(dcopy)(&currDim, tmp_q, &inc, &lambda[i], &N);
        } // i

        // Multiply Lambda %*% w[j] to get wStar. 
        for (j = 0; j < J; j++) {
          F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q], &inc, &zero, &wStar[j * N], &inc FCONE);
        } // j
      }
      /********************************************************************
       *Get fitted values and occurrence probability
       *******************************************************************/
      for (i = 0; i < N; i++) {
        for (j = 0; j < J; j++) {
          psi[j * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + wStar[j * N + i] + betaStarSites[i * J + j], zero, one); 
          z[j * N + i] = rbinom(one, psi[j * N + i]);           
	  if (y[j * N + i] == 1) {
            yWAIC[j * N + i] = psi[j * N + i];
	  } else {
            yWAIC[j * N + i] = 1.0 - psi[j * N + i];
	  }
        } // j
      } // i

      /********************************************************************
       *Save samples
       *******************************************************************/
      if (s >= nBurn) {
        thinIndx++;
        if (thinIndx == nThin) {
          F77_NAME(dcopy)(&pOcc, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pOcc, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pOccN, beta, &inc, &REAL(betaSamples_r)[sPost*pOccN], &inc); 
          F77_NAME(dcopy)(&Nq, lambda, &inc, &REAL(lambdaSamples_r)[sPost*Nq], &inc); 
          F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[sPost*JN], &inc); 
          F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[sPost*JN], &inc); 
          F77_NAME(dcopy)(&Jq, w, &inc, &REAL(wSamples_r)[sPost*Jq], &inc); 
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
    int nResultListObjs = 8;
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 2, betaSamples_r);
    SET_VECTOR_ELT(result_r, 3, zSamples_r);
    SET_VECTOR_ELT(result_r, 4, psiSamples_r);
    SET_VECTOR_ELT(result_r, 5, lambdaSamples_r);
    SET_VECTOR_ELT(result_r, 6, wSamples_r); 
    SET_VECTOR_ELT(result_r, 7, likeSamples_r); 
    if (pOccRE > 0) {
      SET_VECTOR_ELT(result_r, 8, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, 9, betaStarSamples_r);
    }


    // Rf_mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("lambda.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("like.samples")); 
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("beta.star.samples")); 
    }
   
    // Set the names of the output list.  
    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


