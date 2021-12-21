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
  SEXP spIntPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coordsD_r, 
		  SEXP pOcc_r, SEXP pDet_r, SEXP pDetLong_r, 
	          SEXP J_r, SEXP JLong_r, SEXP K_r, SEXP nObs_r, SEXP nObsLong_r, SEXP nData_r, 
		  SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		  SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		  SEXP nuStarting_r, SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
		  SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP sigmaAlpha_r, 
		  SEXP phiA_r, SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, 
		  SEXP nuA_r, SEXP nuB_r, SEXP tuning_r, 
		  SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		  SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		  SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, SEXP currChain_r, SEXP nChain_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, t, s, r, q, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    // Sorted by data set, then visit, then by site. 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    // Sorted by parameter, then data set, site, visit
    double *Xp = REAL(Xp_r);
    double *coordsD = REAL(coordsD_r); 
    int pOcc = INTEGER(pOcc_r)[0];
    int pDet = INTEGER(pDet_r)[0];
    int nData = INTEGER(nData_r)[0]; 
    int *pDetLong = INTEGER(pDetLong_r); 
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    // Covariance model
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    // Priors for regression coefficients
    double *muBeta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(muBeta_r), &inc, muBeta, &inc);
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *SigmaBetaInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double sigmaAlpha = REAL(sigmaAlpha_r)[0]; 
    double phiA = REAL(phiA_r)[0];
    double phiB = REAL(phiB_r)[0]; 
    double sigmaSqA = REAL(sigmaSqA_r)[0]; 
    double sigmaSqB = REAL(sigmaSqB_r)[0]; 
    double nuA = REAL(nuA_r)[0]; 
    double nuB = REAL(nuB_r)[0]; 
    double *tuning = REAL(tuning_r); 
    int J = INTEGER(J_r)[0];
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int nObs = INTEGER(nObs_r)[0]; 
    // Rprintf("nObs: %i\n", nObs); 
    int *nObsLong = INTEGER(nObsLong_r); 
    int *dataIndx = INTEGER(dataIndx_r); 
    int *alphaIndx = INTEGER(alphaIndx_r); 
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
    int status = 0; 
    // For looping through data sets
    int stNObs = 0; 
    int stAlpha = 0; 
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
        Rprintf("Spatial Integrated Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
        Rprintf("Integrating %i occupancy data sets.\n\n", nData); 
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
    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Occupancy covariates
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    // Detection covariates
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Spatial random effects
    double *w = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(wStarting_r), &inc, w, &inc);
    // Latent Occurrence
    double *z = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(zStarting_r), &inc, z, &inc);
    // Spatial smooth parameter for matern. 
    double nu = REAL(nuStarting_r)[0]; 
    // Auxiliary variables
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double));
    double *omegaOcc = (double *) R_alloc(J, sizeof(double));
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); 
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); 

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(INTSXP, nObs, nPost)); nProtect++; 
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    int JJ = J * J; 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    double * tmp_JJ = (double *) R_alloc(JJ, sizeof(double)); 
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = zero; 
    }
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
   
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J);
    int *yRep = (int *) R_alloc(nObs, sizeof(int)); 

    // For normal priors
    // Occupancy regression coefficient priors. 
    // Compute cholesky
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    // Compute inverse
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    // dgemv computes linear combinations of different variables. 
    F77_NAME(dgemv)(ytran, &pOcc, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, SigmaBetaInvMuBeta, &inc); 	  
    // Detection regression coefficient priors. 
    // Have "separate" multivariate normal priors for the different sets of coefficients
    // that vary across the data sets. 
    // Get size of vector
    int currSize = 0; 
    // Index of starting prior values. 
    int *alphaSigmaIndx = (int *) R_alloc(nData, sizeof(int)); 
    int *alphaMuIndx = (int *) R_alloc(nData, sizeof(int)); 
    int tmp0 = 0; 
    int tmp02 = 0; 
    for (q = 0; q < nData; q++) {
      currSize += pDetLong[q] * pDetLong[q];  
      alphaSigmaIndx[q] = tmp0; 
      tmp0 += pDetLong[q] * pDetLong[q]; 
      alphaMuIndx[q] = tmp02; 
      tmp02 += pDetLong[q]; 
    } // q
    double *SigmaAlphaInv = (double *) R_alloc(currSize, sizeof(double)); zeros(SigmaAlphaInv, currSize); 
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    // Fill SigmaAlpha
    for (q = 0; q < nData; q++) {
      for (i = 0; i < pDetLong[q]; i++) {
        SigmaAlphaInv[alphaSigmaIndx[q] + i * pDetLong[q] + i] = sigmaAlpha; 
	// Rprintf("Index: %i\n", alphaSigmaIndx[q] + i * pDetLong[q] + i); 
      } // i
      F77_NAME(dpotrf)(lower, &pDetLong[q], &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &info); 
      if(info != 0){error("c++ error: dpotrf SigmaAlphaInv failed\n");}
      F77_NAME(dpotri)(lower, &pDetLong[q], &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &info); 
      if(info != 0){error("c++ error: dpotri SigmaAlphaInv failed\n");}
      F77_NAME(dgemv)(ytran, &pDetLong[q], &pDetLong[q], &one, &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &muAlpha[alphaMuIndx[q]], &inc, &zero, &SigmaAlphaInvMuAlpha[alphaMuIndx[q]], &inc); 	  
    } // q

    /**********************************************************************
     * Set up spatial stuff and MH stuff
     * *******************************************************************/
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }  
    double *accept = (double *) R_alloc(nTheta, sizeof(double)); zeros(accept, nTheta); 
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    double logMHRatio, logPostCurr = 0.0, logPostCand = 0.0, detCand = 0.0, detCurr = 0.0;
    double phiCand = 0.0, nuCand = 0.0;  
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nTheta, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nTheta, nBatch)); nProtect++; 
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nPost)); nProtect++; 
    // Initiate spatial values
    theta[sigmaSqIndx] = REAL(sigmaSqStarting_r)[0]; 
    double phi = REAL(phiStarting_r)[0]; 
    theta[phiIndx] = phi; 
    if (corName == "matern") {
      theta[nuIndx] = nu; 
    }
    double *C = (double *) R_alloc(JJ, sizeof(double));
    double *CCand = (double *) R_alloc(JJ, sizeof(double));
    double *tmp_JD = (double *) R_alloc(J, sizeof(double));
    double *tmp_JD2 = (double *) R_alloc(J, sizeof(double));
    double *R = (double *) R_alloc(JJ, sizeof(double)); 
    // Get spatial correlation matrix
    spCorLT(coordsD, J, theta, corName, R); 
    logPostCurr = R_NegInf; 
    // Get spatial covariance matrix 
    spCovLT(coordsD, J, theta, corName, C); 
    // Get cholesky of C
    F77_NAME(dpotrf)(lower, &J, C, &J, &info); 
    if(info != 0){error("c++ error: Cholesky failed in initial covariance matrix\n");}
    // Get inverse Cholesky of C. 
    F77_NAME(dpotri)(lower, &J, C, &J, &info); 
    if(info != 0){error("c++ error: Cholesky inverse failed in initial covariance matrix\n");}
    // C now contains the inverse of the covariance matrix. 
    // For sigmaSq sampler
    double aSigmaSqPost = 0.5 * J + sigmaSqA; 
    double bSigmaSqPost = 0.0; 
    double *wTRInv = (double *) R_alloc(J, sizeof(double)); 


    GetRNGstate();
   
    for (s = 0, t = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, t++) {
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
        for (j = 0; j < J; j++) {
          omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + w[j]);
        } // j
        /********************************************************************
         *Update Detection Auxiliary Variables 
         *******************************************************************/
        // Note that all of the variables are sampled, but only those at 
        // locations with z[j] == 1 actually effect the results. 
        for (i = 0; i < nObs; i++) {
          stAlpha = which(dataIndx[i], alphaIndx, pDet); 
          omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, &alpha[stAlpha], &inc));
          // Rprintf("omegaDet[%i]: %f\n", i, omegaDet[i]); 
        } // i
           
        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          kappaOcc[j] = z[j] - 1.0 / 2.0; 
	  tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * w[j]; 
        } // j
        /********************************
         * Compute b.beta
         *******************************/
        // X * tmp_J1 + 0 * tmp_p. Output is stored in tmp_p
        F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc); 	 
        for (j = 0; j < pOcc; j++) {
          tmp_pOcc[j] += SigmaBetaInvMuBeta[j]; 
        } // j 
        /********************************
         * Compute A.beta
         * *****************************/
        // tmp_JpOcc is X %*% omegaOcc. 
        for(j = 0; j < J; j++){
          for(i = 0; i < pOcc; i++){
            tmp_JpOcc[i*J+j] = X[i*J+j]*omegaOcc[j];
          }
        }
        // This finishes off A.beta
        // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc);
        for (j = 0; j < ppOcc; j++) {
          tmp_ppOcc[j] += SigmaBetaInv[j]; 
        } // j
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotri here failed\n");}
        // 1 * tmp_ppOcc * tmp_pOcc + 0 * tmp_pOcc2 
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); if(info != 0){error("c++ error: dpotrf here failed\n");}
        // Args: destination, mu, cholesky of the covariance matrix, dimension
        mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);

      
        /********************************************************************
         *Update Detection Regression Coefficients
         *******************************************************************/
        for (q = 0; q < nData; q++) {
          // Rprintf("q: %i\n", q); 
          // Starting locations
          stNObs = which(q, dataIndx, nObs); 
          stAlpha = which(q, alphaIndx, pDet); 
          // Rprintf("nObsLong[%i]: %i\n", q, nObsLong[q]); 
          /********************************
           * Compute b.alpha
           *******************************/
          // First multiply kappaDet * the current occupied values, such that values go 
          // to 0 if z == 0 and values go to kappaDet if z == 1
          for (i = 0; i < nObsLong[q]; i++) {
            // 1.0 is currently hardcoded in for occupancy data
            kappaDet[stNObs + i] = (y[stNObs + i] - 1.0/2.0) * z[zLongIndx[stNObs + i]];
          } // i
          // Xp * kappaDet + 0 * tmp_pDet. Output is stored in tmp_pDet
          F77_NAME(dgemv)(ytran, &nObsLong[q], &pDetLong[q], &one, &Xp[stNObs], &nObs, &kappaDet[stNObs], &inc, &zero, &tmp_pDet[stAlpha], &inc); 	  
          for (j = 0; j < pDetLong[q]; j++) {
            tmp_pDet[stAlpha + j] += SigmaAlphaInvMuAlpha[stAlpha + j]; 
            // Rprintf("tmp_pDet: %f\n", tmp_pDet[stAlpha + j]); 
          } // j


          /********************************
           * Compute A.alpha
           * *****************************/
          for (j = 0; j < nObsLong[q]; j++) {
            for (i = 0; i < pDetLong[q]; i++) {
              tmp_nObspDet[stNObs + i*nObs + j] = Xp[stNObs + i * nObs + j] * omegaDet[stNObs + j] * z[zLongIndx[stNObs + j]];
              // Rprintf("tmp_nObspDet: %f\n", tmp_nObspDet[stNObs + i*nObs + j]);  
              // Rprintf("omegaDet[%i]: %f\n", stNObs + j, omegaDet[stNObs + j]); 
            } // i
          } // j

          // This finishes off A.alpha
          // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
          F77_NAME(dgemm)(ytran, ntran, &pDetLong[q], &pDetLong[q], &nObsLong[q], &one, &Xp[stNObs], &nObs, &tmp_nObspDet[stNObs], &nObs, &zero, &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q]);

          for (j = 0; j < pDetLong[q] * pDetLong[q]; j++) {
            tmp_ppDet[alphaSigmaIndx[q] + j] += SigmaAlphaInv[alphaSigmaIndx[q] + j]; 
            // Rprintf("tmp_ppDet: %f\n", tmp_ppDet[alphaSigmaIndx[q] + j]); 
          } // j

          // This gives the Cholesky of A.alpha
          // Computes cholesky of tmp_ppDet. Output stored in tmp_ppOcc
          F77_NAME(dpotrf)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info); 
          if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
          // Computes the inverse tmp_ppOcc. Stored in tmp_ppOcc. This is A.beta.inv. 
          F77_NAME(dpotri)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info); 
          if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
          // A.alpha.inv %*% b.alpha
          // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2 
          // (which is currently nothing) = tmp_pDet2
          F77_NAME(dsymv)(lower, &pDetLong[q], &one, &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &tmp_pDet[stAlpha], &inc, &zero, &tmp_pDet2[stAlpha], &inc);
          // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
          F77_NAME(dpotrf)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info); 
          if(info != 0){error("c++ error: dpotrf here failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(&alpha[stAlpha], &tmp_pDet2[stAlpha], &tmp_ppDet[alphaSigmaIndx[q]], pDetLong[q]);
          // for (j = 0; j < pDetLong[q]; j++) {
          //   Rprintf("alpha[%i]: %f\n", stAlpha + j, alpha[stAlpha + j]); 
          // } // j
        } // q

	/********************************************************************
         *Update sigmaSq
         *******************************************************************/
	// Get inverse correlation matrix in reverse from inverse covariance matrix
	// Remember: C currently contains the inverse of covariance matrix. 
	fillUTri(C, J); 
	for (j = 0; j < JJ; j++) {
          R[j] = theta[sigmaSqIndx] * C[j]; 
	} // j
	// Compute t(w) %*% R^-1 %*% w / 
	// t(w) %*% R^-1
	// Probably a better way to do this operation. 
	for (j = 0; j < J; j++) {
          wTRInv[j] = F77_NAME(ddot)(&J, &R[j], &J, w, &inc);  
	  // Rprintf("wTRInv[j]: %f\n", wTRInv[j]); 
        } // j
	// // // wTRInv %*% w
	bSigmaSqPost = F77_NAME(ddot)(&J, wTRInv, &inc, w, &inc); 
	bSigmaSqPost /= 2.0; 
	bSigmaSqPost += sigmaSqB; 
	// Rprintf("bSigmaSqPost: %f\n", bSigmaSqPost); 
	theta[sigmaSqIndx] = rigamma(aSigmaSqPost, bSigmaSqPost); 

        /********************************************************************
         *Update phi (and nu if matern)
         *******************************************************************/
	if (corName == "matern") {
          nu = theta[nuIndx]; 
	  nuCand = logitInv(rnorm(logit(theta[nuIndx], nuA, nuB), exp(tuning[nuIndx])), nuA, nuB); 
        }
	phi = theta[phiIndx]; 
	phiCand = logitInv(rnorm(logit(phi, phiA, phiB), exp(tuning[phiIndx])), phiA, phiB); 
	theta[phiIndx] = phiCand; 
	theta[nuIndx] = nuCand; 

	// Construct covariance matrix (stored in C). 
	spCovLT(coordsD, J, theta, corName, CCand); 

        /********************************
         * Proposal
         *******************************/
	// Invert CCand and log det cov. 
        detCand = 0.0;
	F77_NAME(dpotrf)(lower, &J, CCand, &J, &info); 
	if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	// Get log of the determinant of the covariance matrix. 
	for (k = 0; k < J; k++) {
	  detCand += 2.0 * log(CCand[k*J+k]);
	} // k
	F77_NAME(dpotri)(lower, &J, CCand, &J, &info); 
	if(info != 0){error("c++ error: Cholesky inverse failed in covariance matrix\n");}
        logPostCand = 0.0; 
	// Jacobian and Uniform prior. 
	logPostCand += log(phiCand - phiA) + log(phiB - phiCand); 
	// (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	F77_NAME(dsymv)(lower, &J, &one,  CCand, &J, w, &inc, &zero, tmp_JD, &inc);
	logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&J, w, &inc, tmp_JD, &inc);
        if (corName == "matern"){
          logPostCand += log(nuCand - nuA) + log(nuB - nuCand); 
        }

        /********************************
         * Current
         *******************************/
	theta[nuIndx] = nu; 
	theta[phiIndx] = phi; 
	// Construct covariance matrix (stored in C). 
	spCovLT(coordsD, J, theta, corName, C); 
        detCurr = 0.0;
	F77_NAME(dpotrf)(lower, &J, C, &J, &info); 
	if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	for (k = 0; k < J; k++) {
	  detCurr += 2.0 * log(C[k*J+k]);
	} // k
	F77_NAME(dpotri)(lower, &J, C, &J, &info); 
	if(info != 0){error("c++ error: Cholesky inverse failed in covariance matrix\n");}
        logPostCurr = 0.0; 
	logPostCurr += log(phi - phiA) + log(phiB - phi); 
	// (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	F77_NAME(dsymv)(lower, &J, &one, C, &J, w, &inc, &zero, tmp_JD, &inc);
	logPostCurr += -0.5*detCurr-0.5*F77_NAME(ddot)(&J, w, &inc, tmp_JD, &inc);
        if (corName == "matern"){
          logPostCurr += log(nu - nuA) + log(nuB - nu); 
        }

	// MH Accept/Reject
	logMHRatio = logPostCand - logPostCurr; 
	if (runif(0.0, 1.0) <= exp(logMHRatio)) {
          theta[phiIndx] = phiCand;
          accept[phiIndx]++;
          if (corName == "matern") {
            theta[nuIndx] = nuCand; 
            accept[nuIndx]++; 
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
          tmp_JD[j] = kappaOcc[j] - F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) * omegaOcc[j];
        }
        /********************************
         * Compute A.w
         *******************************/
	// Copy inverse covariance matrix into tmp_JJ
	F77_NAME(dcopy)(&JJ, C, &inc, tmp_JJ, &inc); 
	for (k = 0; k < J; k++) {
	  tmp_JJ[k * J + k] += omegaOcc[k]; 
	} // k

        // Cholesky of A.w
        F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info); 
        if(info != 0){error("c++ error: dpotrf on A.w failed\n");}
	// Inverse of A.w
        F77_NAME(dpotri)(lower, &J, tmp_JJ, &J, &info); 
        if(info != 0){error("c++ error: dpotri on A.w failed\n");}
        // A.w.inv %*% b.w. Stored in tmp_JD2
        F77_NAME(dsymv)(lower, &J, &one, tmp_JJ, &J, tmp_JD, &inc, &zero, tmp_JD2, &inc);
        // Computes cholesky of tmp_JJ again stored back in tmp_JJ. This chol(A.beta.inv)
        F77_NAME(dpotrf)(lower, &J, tmp_JJ, &J, &info); 
	if(info != 0){error("c++ error: dpotrf on A.w failed\n");}
        // Args: destination, mu, cholesky of the covariance matrix, dimension
        mvrnorm(w, tmp_JD2, tmp_JJ, J);

	// At end of each MCMC, the variable C contains the inverse of the 
	// current covariance matrix. 
     
        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
        for (i = 0; i < nObs; i++) {
          stAlpha = which(dataIndx[i], alphaIndx, pDet); 
          detProb[i] = logitInv(F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, &alpha[stAlpha], &inc), zero, one);
          if (tmp_J[zLongIndx[i]] == 0) {
            psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + w[zLongIndx[i]], zero, one); 
          }
          piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
          ySum[zLongIndx[i]] += y[i]; 	
          tmp_J[zLongIndx[i]]++;
        } // i
        // Compute occupancy probability 
        for (j = 0; j < J; j++) {
          psiNum = psi[j] * piProd[j]; 
          if (ySum[j] == zero) {
            z[j] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j])));           
          } else {
            z[j] = one; 
          }
          piProd[j] = one;
          ySum[j] = zero; 
          tmp_J[j] = 0; 
        } // j


        /********************************************************************
         *Save samples
         *******************************************************************/
	if (t >= nBurn) {
          thinIndx++; 
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
            F77_NAME(dcopy)(&J, w, &inc, &REAL(wSamples_r)[sPost*J], &inc); 
	    F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[sPost*nTheta], &inc); 
	    F77_NAME(dcopy)(&J, z, &inc, &REAL(zSamples_r)[sPost*J], &inc); 
	    // Replicate data set for GoF
            for (i = 0; i < nObs; i++) {
              yRep[i] = rbinom(one, detProb[i] * z[zLongIndx[i]]);
              INTEGER(yRepSamples_r)[sPost * nObs + i] = yRep[i]; 
            } // i
	    sPost++; 
	    thinIndx = 0; 
	  }
	}

        R_CheckUserInterrupt();
      } // end batch

      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (j = 0; j < nTheta; j++) {
        REAL(acceptSamples_r)[s * nTheta + j] = accept[j]/batchLength; 
        REAL(tuningSamples_r)[s * nTheta + j] = tuning[j]; 
        if (accept[j] / batchLength > acceptRate) {
          tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        } else{
            tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          }
        accept[j] = 0;
      }
      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
	  Rprintf("\tAcceptance\tTuning\n");	  
	  Rprintf("\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nTheta + phiIndx], exp(tuning[phiIndx]));
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;        
    } // s (sample loop)
    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }
    // This is necessary when generating random numbers in C.     
    PutRNGstate();

    //make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 9;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, yRepSamples_r);
    SET_VECTOR_ELT(result_r, 5, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 6, wSamples_r); 
    SET_VECTOR_ELT(result_r, 7, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 8, acceptSamples_r); 
    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("y.rep.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("phi.tune")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("phi.accept")); 
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

