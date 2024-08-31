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
  SEXP intPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
                SEXP consts_r, SEXP pDetLong_r, 
                SEXP JLong_r, SEXP nObsLong_r, SEXP K_r, SEXP nOccRELong_r,
                SEXP nDetRELong_r, 
                SEXP betaStarting_r, SEXP alphaStarting_r, SEXP sigmaSqPsiStarting_r, 
                SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, 
                SEXP alphaStarStarting_r, SEXP zStarting_r, 
                SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
                SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
                SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, SEXP alphaNREIndx_r,
                SEXP alphaColIndx_r, SEXP waicJIndx_r, SEXP waicNObsIndx_r, SEXP WAIC_r,
                SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP sigmaAlpha_r, 
                SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
                SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
                SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
                SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, l, ll, q, info, nProtect=0;
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
    int *XRE = INTEGER(XRE_r);
    // Sorted by parameter, then data set, site, visit
    double *Xp = REAL(Xp_r);
    int *XpRE = INTEGER(XpRE_r);
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1];
    int pOcc = INTEGER(consts_r)[2];
    int pOccRE = INTEGER(consts_r)[3];
    int nOccRE = INTEGER(consts_r)[4];
    int pDet = INTEGER(consts_r)[5];
    int pDetRE = INTEGER(consts_r)[6];
    int nDetRE = INTEGER(consts_r)[7];
    int nData = INTEGER(consts_r)[8]; 
    int *pDetLong = INTEGER(pDetLong_r); 
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    // Priors for regression coefficients
    double *muBeta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(muBeta_r), &inc, muBeta, &inc);
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *SigmaBetaInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *sigmaAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(sigmaAlpha_r), &inc, sigmaAlpha, &inc);
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    int *JLong = INTEGER(JLong_r);
    // Total number of sites across all data sets, including
    // sites that are sampled by multiple data sources
    int JSum = 0; 
    for (i = 0; i < nData; i++) {
      JSum += JLong[i];
    }
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *nObsLong = INTEGER(nObsLong_r); 
    int *dataIndx = INTEGER(dataIndx_r); 
    int *alphaIndx = INTEGER(alphaIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int *alphaNREIndx = INTEGER(alphaNREIndx_r);
    int *alphaColIndx = INTEGER(alphaColIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int *waicJIndx = INTEGER(waicJIndx_r);
    int waic = INTEGER(WAIC_r)[0];
    int *waicNObsIndx = INTEGER(waicNObsIndx_r);
    int *nOccRELong = INTEGER(nOccRELong_r); 
    int *nDetRELong = INTEGER(nDetRELong_r);
    int nSamples = INTEGER(nSamples_r)[0];
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int nThreads = INTEGER(nThreads_r)[0];
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
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
        Rprintf("Integrated Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
        Rprintf("Integrating %i occupancy data sets.\n\n", nData); 
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
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
     * Parameters
     * *******************************************************************/
    // Occupancy covariates
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&nOccRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Detection covariates
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(J, sizeof(double)); zeros(omegaOcc, J);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); zeros(kappaOcc, J);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOcc * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDet * nPost);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), J * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), J * nPost);
    SEXP pSamples_r; 
    PROTECT(pSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
    zeros(REAL(pSamples_r), nObs * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = Rf_allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = Rf_allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetRE * nPost);
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = Rf_allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPsiSamples_r), pOccRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nOccRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nOccRE * nPost);
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, JSum, nPost)); nProtect++;
    
    /**********************************************************************
     * Some constants and temporary variables to be used later
     * *******************************************************************/
    int JpOcc = J * pOcc; 
    int JpOccRE = J * pOccRE;
    int nObspDet = nObs * pDet;
    int nObspDetRE = nObs * pDetRE;
    double tmp_0 = 0.0; 
    double tmp_02 = 0.0; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = zero; 
    }
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(tmp_nObs, nObs);
   
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); zeros(detProb, nObs);
    double *yWAIC = (double *) R_alloc(JSum, sizeof(double)); zeros(yWAIC, JSum);
    double *piProdWAIC = (double *) R_alloc(JSum, sizeof(double)); 
    ones(piProdWAIC, JSum); 
    double *yWAICSum = (double *) R_alloc(JSum, sizeof(double)); zeros(yWAICSum, JSum);
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J); 
    // For normal priors
    // Occupancy regression coefficient priors. 
    // Compute cholesky
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaInv failed\n");}
    // Compute inverse
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);
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
    for (q = 0, j = 0; q < nData; q++) {
      for (i = 0; i < pDetLong[q]; i++, j++) {
        SigmaAlphaInv[alphaSigmaIndx[q] + i * pDetLong[q] + i] = sigmaAlpha[j]; 
	// Rprintf("Index: %i\n", alphaSigmaIndx[q] + i * pDetLong[q] + i); 
      } // i
      F77_NAME(dpotrf)(lower, &pDetLong[q], &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf SigmaAlphaInv failed\n");}
      F77_NAME(dpotri)(lower, &pDetLong[q], &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri SigmaAlphaInv failed\n");}
      F77_NAME(dsymv)(lower, &pDetLong[q], &one, &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &muAlpha[alphaMuIndx[q]], &inc, &zero, 
                     &SigmaAlphaInvMuAlpha[alphaMuIndx[q]], &inc FCONE);
    } // q

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(J, sizeof(double)); 
    zeros(betaStarSites, J); 
    int *betaStarLongIndx = (int *) R_alloc(JpOccRE, sizeof(int));
    // Initial sums
    for (j = 0; j < J; j++) {
      for (l = 0; l < pOccRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nOccRE);
        betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]];
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(alphaStarObs, nObs); 
    int *alphaStarLongIndx = (int *) R_alloc(nObspDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < alphaNREIndx[i]; l++) {
        alphaStarLongIndx[l * nObs + i] = which(XpRE[l * nObs + i], alphaLevelIndx, nDetRE);
        alphaStarObs[i] += alphaStar[alphaStarLongIndx[l * nObs + i]];
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

    GetRNGstate();
    
    for (s = 0; s < nSamples; s++) {
      /********************************************************************
       *Update Occupancy Auxiliary Variables 
       *******************************************************************/
      for (j = 0; j < J; j++) {
        omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + betaStarSites[j]);
      } // j
      /********************************************************************
       *Update Detection Auxiliary Variables 
       *******************************************************************/
      // Note that all of the variables are sampled, but only those at 
      // locations with z[j] == 1 actually effect the results. 
      for (i = 0; i < nObs; i++) {
        stAlpha = which(dataIndx[i], alphaIndx, pDet); 
	if (z[zLongIndx[i]] == 1.0) {
          omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, 
				                &alpha[stAlpha], &inc) + alphaStarObs[i]);
	}
        // Rprintf("omegaDet[%i]: %f\n", i, omegaDet[i]); 
      } // i
           
      /********************************************************************
       *Update Occupancy Regression Coefficients
       *******************************************************************/
      for (j = 0; j < J; j++) {
        kappaOcc[j] = z[j] - 1.0 / 2.0; 
        tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * betaStarSites[j]; 
      } // j
      /********************************
       * Compute b.beta
       *******************************/
      // X * kappaOcc + 0 * tmp_p. Output is stored in tmp_p
      F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
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
      F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
      for (j = 0; j < ppOcc; j++) {
        tmp_ppOcc[j] += SigmaBetaInv[j]; 
      } // j
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri here failed\n");}
      // A.beta.inv %*% b.beta
      F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
      // Args: destination, mu, cholesky of the covariance matrix, dimension
      mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
      
      /********************************************************************
       *Update Detection Regression Coefficients
       *******************************************************************/
      for (q = 0; q < nData; q++) {
        zeros(tmp_nObs, nObs);
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
	  tmp_nObs[stNObs + i] = kappaDet[stNObs + i] - 
		                 omegaDet[stNObs + i] * alphaStarObs[stNObs + i];
	  tmp_nObs[stNObs + i] *= z[zLongIndx[stNObs + i]];
        } // i
        // Xp * kappaDet + 0 * tmp_pDet. Output is stored in tmp_pDet
        F77_NAME(dgemv)(ytran, &nObsLong[q], &pDetLong[q], &one, &Xp[stNObs], &nObs, &tmp_nObs[stNObs], &inc, &zero, &tmp_pDet[stAlpha], &inc FCONE); 	  
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
        F77_NAME(dgemm)(ytran, ntran, &pDetLong[q], &pDetLong[q], &nObsLong[q], &one, &Xp[stNObs], 
			&nObs, &tmp_nObspDet[stNObs], &nObs, &zero, &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q] FCONE FCONE);

        for (j = 0; j < pDetLong[q] * pDetLong[q]; j++) {
          tmp_ppDet[alphaSigmaIndx[q] + j] += SigmaAlphaInv[alphaSigmaIndx[q] + j]; 
          // Rprintf("tmp_ppDet: %f\n", tmp_ppDet[alphaSigmaIndx[q] + j]); 
        } // j

        // This gives the Cholesky of A.alpha
        // Computes cholesky of tmp_ppDet. Output stored in tmp_ppOcc
        F77_NAME(dpotrf)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf A.alpha failed\n");}
        // Computes the inverse tmp_ppOcc. Stored in tmp_ppOcc. This is A.beta.inv. 
        F77_NAME(dpotri)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri A.alpha failed\n");}
        // A.alpha.inv %*% b.alpha
        // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2 
        // (which is currently nothing) = tmp_pDet2
        F77_NAME(dsymv)(lower, &pDetLong[q], &one, &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &tmp_pDet[stAlpha], &inc, &zero, &tmp_pDet2[stAlpha], &inc FCONE);
        // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
        F77_NAME(dpotrf)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        // Args: destination, mu, cholesky of the covariance matrix, dimension
        mvrnorm(&alpha[stAlpha], &tmp_pDet2[stAlpha], &tmp_ppDet[alphaSigmaIndx[q]], pDetLong[q]);
      } // q

      /********************************************************************
       *Update Occupancy random effects variance
       *******************************************************************/
      for (l = 0; l < pOccRE; l++) {
        tmp_0 = F77_NAME(ddot)(&nOccRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
        tmp_0 *= 0.5; 
        sigmaSqPsi[l] = rigamma(sigmaSqPsiA[l] + nOccRELong[l] / 2.0, sigmaSqPsiB[l] + tmp_0); 
      }

      /********************************************************************
       *Update Detection random effects variance
       *******************************************************************/
      for (l = 0; l < pDetRE; l++) {
        tmp_0 = F77_NAME(ddot)(&nDetRELong[l], &alphaStar[alphaStarStart[l]], &inc, &alphaStar[alphaStarStart[l]], &inc); 
        tmp_0 *= 0.5; 
        sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] / 2.0, sigmaSqPB[l] + tmp_0); 
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
                tmp_02 += betaStar[betaStarLongIndx[ll * J + j]];
              } 
              tmp_one[0] += kappaOcc[j] - (F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc) + 
        		    tmp_02 - betaStar[l]) * omegaOcc[j];
              tmp_0 += omegaOcc[j];
            }
          }
          /********************************
           * Compute A.beta.star
           *******************************/
          tmp_0 += 1.0 / sigmaSqPsi[betaStarIndx[l]]; 
          tmp_0 = 1.0 / tmp_0; 
          betaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
        }
      
        // Update the RE sums for the current species
        zeros(betaStarSites, J);
        for (j = 0; j < J; j++) {
          for (l = 0; l < pOccRE; l++) {
            betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]];
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
          for (i = 0; i < nObs; i++) {
            stAlpha = which(dataIndx[i], alphaIndx, pDet); 
            if ((z[zLongIndx[i]] == 1.0) && (XpRE[alphaColIndx[l] * nObs + i] == alphaLevelIndx[l])) {
              tmp_02 = 0.0;
      	for (ll = 0; ll < alphaNREIndx[i]; ll++) {
                tmp_02 += alphaStar[alphaStarLongIndx[ll * nObs + i]];
      	}
              tmp_one[0] += kappaDet[i] - (F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, &alpha[stAlpha], &inc) + tmp_02 - alphaStar[l]) * omegaDet[i];
      	tmp_0 += omegaDet[i];
            }
          }
          /********************************
           * Compute A.alpha.star
           *******************************/
          tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
          tmp_0 = 1.0 / tmp_0; 
          alphaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
        }
        zeros(alphaStarObs, nObs); 
        // Update the RE sums
        for (i = 0; i < nObs; i++) {
          for (l = 0; l < alphaNREIndx[i]; l++) {
            alphaStarObs[i] += alphaStar[alphaStarLongIndx[l * nObs + i]]; 
          }
        }
      }

      /********************************************************************
       *Update Latent Occupancy
       *******************************************************************/
      // Compute detection probability 
      for (i = 0; i < nObs; i++) {
        stAlpha = which(dataIndx[i], alphaIndx, pDet); 
        detProb[i] = logitInv(F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, &alpha[stAlpha], &inc) + alphaStarObs[i], zero, one);
        if (tmp_J[zLongIndx[i]] == 0) {
          psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + 
			               betaStarSites[zLongIndx[i]], zero, one); 
        }
        piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
	if (waic > 0) {
	  piProdWAIC[waicNObsIndx[i]] *= pow(detProb[i], y[i]);
	  piProdWAIC[waicNObsIndx[i]] *= pow(1.0 - detProb[i], 1 - y[i]);
	  yWAICSum[waicNObsIndx[i]] += y[i];
	}
        ySum[zLongIndx[i]] += y[i]; 	
        tmp_J[zLongIndx[i]]++;
      } // i
      // Compute occupancy probability and the integrated likelihood for WAIC
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
      // Calculating WAIC
      if (waic > 0) {
        for (j = 0; j < JSum; j++) {
          if (yWAICSum[j] == zero) {
            yWAIC[j] = (1.0 - psi[waicJIndx[j]]) + psi[waicJIndx[j]] * piProdWAIC[j];
          } else {
            yWAIC[j] = psi[waicJIndx[j]] * piProdWAIC[j];
          }
          piProdWAIC[j] = one;
          yWAICSum[j] = zero;
        } // j
      }

     /********************************************************************
      *Save samples
      *******************************************************************/
      if (s >= nBurn) {
        thinIndx++; 
	if (thinIndx == nThin) {
          F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
          F77_NAME(dcopy)(&J, z, &inc, &REAL(zSamples_r)[sPost*J], &inc); 
          if (pOccRE > 0) {
            F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, 
        		    &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
            F77_NAME(dcopy)(&nOccRE, betaStar, &inc, 
        		    &REAL(betaStarSamples_r)[sPost*nOccRE], &inc);
          }
          if (pDetRE > 0) {
            F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, 
        		    &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
            F77_NAME(dcopy)(&nDetRE, alphaStar, &inc, 
        		    &REAL(alphaStarSamples_r)[sPost*nDetRE], &inc);
          }
          for (i = 0; i < nObs; i++) {
            REAL(pSamples_r)[sPost * nObs + i] = detProb[i]; 
          } // i
	  if (waic > 0) {
            F77_NAME(dcopy)(&JSum, yWAIC, &inc, 
        		    &REAL(likeSamples_r)[sPost*JSum], &inc);
	  }
          sPost++; 
	  thinIndx = 0; 
	}
      }

      R_CheckUserInterrupt();

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


    }
    if (verbose) {
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
    }
    // This is necessary when generating random numbers in C.     
    PutRNGstate();

    //make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 6;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, pSamples_r);
    SET_VECTOR_ELT(result_r, 5, likeSamples_r);
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 6, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 7, alphaStarSamples_r);
    }
    if (pOccRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 8;
      } else {
        tmp_0 = 6;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    // Rf_mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("p.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("alpha.star.samples")); 
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



