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
  SEXP spMsPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coordsD_r, SEXP pocc_r, SEXP pdet_r, 
	         SEXP J_r, SEXP nObs_r, SEXP K_r, SEXP N_r, SEXP betaStarting_r, 
	         SEXP alphaStarting_r, SEXP zStarting_r, SEXP betaCommStarting_r, 
	         SEXP alphaCommStarting_r, SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		 SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		 SEXP nuStarting_r, SEXP zLongIndx_r, SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	         SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	         SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP phiA_r, 
		 SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, 
		 SEXP nuA_r, SEXP nuB_r, SEXP tuning_r, 
		 SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, 
	         SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, SEXP nBurn_r, 
		 SEXP nThin_r, SEXP nPost_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, a, b, q, r, info, nProtect=0;
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
    // Sorted by visit, then by site, then by species. 
    // (e.g., visit 1, site 1, sp 1, v1, s1, sp2, 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *coordsD = REAL(coordsD_r); 
    // Xp is sorted by parameter then site, then visit (parameter 1 site 1, p1 site 2, etc) 
    double *Xp = REAL(Xp_r);
    double *muBetaComm = REAL(muBetaComm_r); 
    double *muAlphaComm = REAL(muAlphaComm_r); 
    double *SigmaBetaCommInv = REAL(SigmaBetaComm_r); 
    double *SigmaAlphaCommInv = REAL(SigmaAlphaComm_r); 
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
    double *tuning = REAL(tuning_r); 
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
    int J = INTEGER(J_r)[0];
    double *K = REAL(K_r); 
    int nObs = INTEGER(nObs_r)[0];
    int N = INTEGER(N_r)[0]; 
    double *phiAccept = (double *) R_alloc(N, sizeof(double)); 
    zeros(phiAccept, N); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nThin = INTEGER(nThin_r)[0];
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    double acceptRate = REAL(acceptRate_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 
    // Sorted by site, then species (e.g., site 1 sp 1, site 1 sp2, ...)
    double *z = REAL(zStarting_r); 
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
      // Rprintf allows you to print messages and value on the R console screen. 
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Spatial Multispecies Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
      Rprintf("Number of MCMC samples: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
      Rprintf("Burn-in: %i \n", nBurn); 
      Rprintf("Thinning Rate: %i \n", nThin); 
      Rprintf("Total Posterior Samples: %i \n\n", nPost); 
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      Rprintf("Adaptive Metropolis with target acceptance rate: %.1f\n", 100*acceptRate);
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
    int JJ = J * J; 
    double tmp_0; 
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
    // Detection covariates
    double *alpha = (double *) R_alloc(pDetN, sizeof(double));   
    F77_NAME(dcopy)(&pDetN, REAL(alphaStarting_r), &inc, alpha, &inc);
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
    // Auxiliary variables
    // Only need to set aside J locations since you don't save these 
    // for each species
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
    // Spatial parameters
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 

    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
   
    // For latent occupancy
    double psiNum; 
    double psiNew; 
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
    // dgemv computes linear combinations of different variables. 
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
    spCorLT(coordsD, J, currTheta, corName, R); 
    // Get spatial covariance matrix 
    spCovLT(coordsD, J, currTheta, corName, C); 
    F77_NAME(dpotrf)(lower, &J, C, &J, &info); 
    if(info != 0){error("c++ error: Cholesky failed in initial covariance matrix\n");}
    F77_NAME(dpotri)(lower, &J, C, &J, &info); 
    if(info != 0){error("c++ error: Cholesky inverse failed in initial covariance matrix\n");}
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nThetaN, nPost)); nProtect++; 

    /**********************************************************************
     Set up stuff for Adaptive MH and other misc
     * *******************************************************************/
    double logMHRatio, logPostCurr = 0.0, logPostCand = 0.0, detCand = 0.0, detCurr = 0.0; 
    logPostCurr = R_NegInf; 
    double *accept = (double *) R_alloc(nThetaN, sizeof(double)); zeros(accept, nThetaN); 
    double phiCand = 0.0, nuCand = 0.0; 
    double logDet; 
    // For sigmaSq Sampler
    double aSigmaSqPost = 0.0; 
    double bSigmaSqPost = 0.0; 
    // MCMC info if desired
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
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
        // A.beta.inv %*% b.beta
        // 1 * tmp_ppOcc * tmp_pOcc + 0 * tmp_pOcc2  = tmp_pOcc2
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
        // Computes cholesky of tmp_pp again stored back in tmp_ppOcc. This chol(A.beta.inv)
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
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
        // A.alpha.inv %*% b.alpha
        // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2  = tmp_pDet2
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc);
        // Computes cholesky of tmp_pp again stored back in tmp_ppDet. This chol(A.alpha.inv)
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
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
          // Rprintf("tauSqBetaB[%i]: %f\n", q, tauSqBetaB[q]); 
          // Rprintf("tmp_0[%i]: %f\n", q, tmp_0 + tauSqBetaB[q]); 
        } // q
        // This is correct, nothing wrong here. 
        for (q = 0; q < pOcc; q++) {
          TauBetaInv[q * pOcc + q] = tauSqBeta[q]; 
          // Rprintf("TauBetaInv[%i]: %f\n", q * pOcc + q, TauBetaInv[q * pOcc + q]); 
        } // i
        // for (q = 0; q < ppOcc; q++) {
        //   Rprintf("TauBetaInv[%i]: %f\n", q, TauBetaInv[q]); 
        // }
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
          // Rprintf("TauAlphaInv[%i]: %f\n", q * pDet + q, TauAlphaInv[q * pDet + q]); 
        } // i
        F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info); 
        if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
        F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info); 
        if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}
         
        for (i = 0; i < N; i++) {  
          /********************************************************************
           *Update Occupancy Auxiliary Variables 
           *******************************************************************/
          for (j = 0; j < J; j++) {
            omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + w[j * N + i]);
          } // j
          /********************************************************************
           *Update Detection Auxiliary Variables 
           *******************************************************************/
          // Note that all of the variables are sampled, but only those at 
          // locations with z[j] == 1 actually effect the results. 
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              omegaDet[r] = rpg(K[r], F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N));
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              omegaDet[r] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N));
            } // r
          }
             
          /********************************************************************
           *Update Occupancy Regression Coefficients
           *******************************************************************/
          for (j = 0; j < J; j++) {
            kappaOcc[j] = z[j * N + i] - 1.0 / 2.0; 
	    tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * w[j * N + i]; 
          } // j
          /********************************
           * Compute b.beta
           *******************************/
          // t(X) * tmp_J1 + 0 * tmp_pOcc = tmp_pOcc. 
          // dgemv computes linear combinations of different variables. 
          F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc); 	 
          // TauBetaInv %*% betaComm + tmp_pOcc = tmp_pOcc
          F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, betaComm, &inc, &one, tmp_pOcc, &inc); 

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
          F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc);
          for (q = 0; q < ppOcc; q++) {
            tmp_ppOcc[q] += TauBetaInv[q]; 
          } // j
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
          if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
          if(info != 0){error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); if(info != 0){error("c++ error: dpotrf here failed\n");}
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
            } // r
          } else { 
            for (r = 0; r < nObs; r++) {
              kappaDet[r] = (y[r * N + i] - 1.0/2.0) * z[zLongIndx[r] * N + i];
            } // r
          }
          
          F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, kappaDet, &inc, &zero, tmp_pDet, &inc); 	  
          F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, TauAlphaInv, &pDet, alphaComm, &inc, &one, tmp_pDet, &inc); 
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
          F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet);

          for (q = 0; q < ppDet; q++) {
            tmp_ppDet[q] += TauAlphaInv[q]; 
            // Rprintf("TauAlphaInv: %f\n", TauAlphaInv[q]); 
          } // q
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
          if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
          F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info); 
          if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
          // A.alpha.inv %*% b.alpha
          // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2 
          // (which is currently nothing) = tmp_pDet2
          F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc);
          // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
          if(info != 0){error("c++ error: dpotrf here failed\n");}
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
	  // Get inverse correlation matrix
          spCorLT(coordsD, J, currTheta, corName, R); 
	  fillUTri(R, J); 
	  // Save R in tmp_JJ
          F77_NAME(dcopy)(&JJ, R, &inc, tmp_JJ, &inc); 
          F77_NAME(dpotrf)(lower, &J, R, &J, &info); 
          if(info != 0){error("c++ error: Cholesky failed in correlation matrix\n");}
          F77_NAME(dpotri)(lower, &J, R, &J, &info); 
          if(info != 0){error("c++ error: Cholesky inverse failed in correlation matrix\n");}
	  fillUTri(R, J); 

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

	  // Construct proposal covariance matrix (stored in CCand). 
	  spCovLT(coordsD, J, currTheta, corName, CCand); 

          /********************************
           * Proposal
           *******************************/
	  // Invert CCand and log det cov. 
          detCand = 0.0;
	  F77_NAME(dpotrf)(lower, &J, CCand, &J, &info); 
	  if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	  // Get log of the determinant of the covariance matrix. 
	  for (k = 0; k < J; k++) {
            // Rprintf("detCand Value: %f\n", 2.0 * log(CCand[k * J + k])); 
	    detCand += 2.0 * log(CCand[k*J+k]);
	  } // k
	  F77_NAME(dpotri)(lower, &J, CCand, &J, &info); 
	  if(info != 0){error("c++ error: Cholesky inverse failed in covariance matrix\n");}
          logPostCand = 0.0; 
	  // Jacobian and Uniform prior. 
	  logPostCand += log(phiCand - phiA[i]) + log(phiB[i] - phiCand); 
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &J, &one,  CCand, &J, &w[i], &N, &zero, tmp_JD, &inc);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&J, &w[i], &N, tmp_JD, &inc);
	  // Rprintf("logPostCand: %f\n", logPostCand); 
          if (corName == "matern"){
            logPostCand += log(nuCand - nuA[i]) + log(nuB[i] - nuCand); 
          }

          /********************************
           * Current
           *******************************/
	  // Get logdetCov
          detCurr = 0.0;
	  // tmp_JJ has the current covariance matrix
	  F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info); 
	  if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	  // Get log of the determinant of the covariance matrix. 
	  for (k = 0; k < J; k++) {
	    detCurr += 2.0 * log(tmp_JJ[k*J+k]);
	  } // k
	  // Inverse of current covariance matrix is already in C. 
          logPostCurr = 0.0; 
	  // Jacobian and Uniform prior. 
	  logPostCurr += log(phi[i] - phiA[i]) + log(phiB[i] - phi[i]); 
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &J, &one, C, &J, &w[i], &N, &zero, tmp_JD, &inc);
	  logPostCurr += -0.5*detCurr-0.5*F77_NAME(ddot)(&J, &w[i], &N, tmp_JD, &inc);
	  // Rprintf("logPostCurr: %f\n", logPostCurr); 
          if (corName == "matern"){
            logPostCurr += log(nu[i] - nuA[i]) + log(nuB[i] - nu[i]); 
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
	    F77_NAME(dcopy)(&JJ, CCand, &inc, C, &inc); 
          }

          /********************************************************************
           *Update w (spatial random effects)
           *******************************************************************/
          /********************************
           * Compute b.w
           *******************************/
          for(j = 0; j < J; j++){
            tmp_JD[j] = kappaOcc[j] - F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) * omegaOcc[j];
          }
          /********************************
           * Compute A.w
           *******************************/
	  // Copy inverse covariance matrix into tmp_JJ
	  F77_NAME(dcopy)(&JJ, C, &inc, tmp_JJ, &inc); 
	  for (k = 0; k < J; k++) {
	    tmp_JJ[k * J + k] += omegaOcc[k]; 
	  } // k
          F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info); 
          if(info != 0){error("c++ error: dpotrf on A.w failed\n");}
          F77_NAME(dpotri)(lower, &J, tmp_JJ, &J, &info); 
          if(info != 0){error("c++ error: dpotri on A.w failed\n");}
          F77_NAME(dsymv)(lower, &J, &one, tmp_JJ, &J, tmp_JD, &inc, &zero, tmp_JD2, &inc);
          F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info); 
	  if(info != 0){error("c++ error: dpotrf on A.w failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_J1, tmp_JD2, tmp_JJ, J);
          for (j = 0; j < J; j++) {
            w[j * N + i] = tmp_J1[j]; 
          }

          /********************************************************************
           *Update Latent Occupancy
           *******************************************************************/
          // Compute detection probability 
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N), zero, one);
              psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + w[zLongIndx[r] * N + i], zero, one); 
              piProd[zLongIndx[r]] = pow(1.0 - detProb[i * nObs + r], K[r]);
              ySum[zLongIndx[r]] = y[r * N + i]; 
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N), zero, one);
              if (tmp_J[zLongIndx[r]] == 0) {
                psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + w[zLongIndx[r] * N + i], zero, one); 
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
            piProd[j] = one;
            ySum[j] = zero; 
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
	  Rprintf("\tSpecies\t\tAcceptance\tTuning\n");	  
	  for (i = 0; i < N; i++) {
	    Rprintf("\t%i\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nThetaN + phiIndx * N + i], exp(tuning[phiIndx * N + i]));
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
    int nResultListObjs = 12;

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
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


