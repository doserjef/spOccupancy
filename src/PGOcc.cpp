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
  SEXP PGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
             SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
             SEXP betaStarting_r, SEXP alphaStarting_r, 
             SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
             SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
             SEXP zStarting_r, SEXP zLongIndx_r, 
             SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
             SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
             SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
             SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, SEXP sigmaSqPA_r, 
             SEXP sigmaSqPB_r, SEXP nSamples_r, SEXP nThreads_r, 
             SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
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
    double *Xp = REAL(Xp_r);
    int *XRE = INTEGER(XRE_r); 
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
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    double *muBeta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(muBeta_r), &inc, muBeta, &inc);
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *SigmaBetaInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *SigmaAlphaInv = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlpha_r), &inc, SigmaAlphaInv, &inc);
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
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
        Rprintf("Occupancy model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
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
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), J * nPost);
    
    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int JpOcc = J * pOcc; 
    int JpOccRE = J * pOccRE; 
    int nObspDet = nObs * pDet;
    int nObspDetRE = nObs * pDetRE;
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
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
   
    // For latent occupancy and WAIC
    double psiNum; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); zeros(detProb, nObs);
    double *yWAIC = (double *) R_alloc(J, sizeof(double)); zeros(yWAIC, J);
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *piProdWAIC = (double *) R_alloc(J, sizeof(double)); 
    ones(piProdWAIC, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J);

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, 
                   SigmaAlphaInvMuAlpha, &inc FCONE);

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
      for (l = 0; l < pDetRE; l++) {
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
      if (nObs == J) {
        for (i = 0; i < nObs; i++) {
          if (z[zLongIndx[i]] == 1.0) {
            omegaDet[i] = rpg(K[i], F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
        		      alphaStarObs[i]);
          }
        } // i
      } else {
        for (i = 0; i < nObs; i++) {
          if (z[zLongIndx[i]] == 1.0) {
            omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
        		      alphaStarObs[i]);
          }
        } // i
      }

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
      F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
      mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
      
      /********************************************************************
       *Update Detection Regression Coefficients
       *******************************************************************/
      /********************************
       * Compute b.alpha
       *******************************/
      if (nObs == J) {
        for (i = 0; i < nObs; i++) {
          kappaDet[i] = (y[i] - K[i]/2.0) * z[zLongIndx[i]];
          tmp_nObs[i] = kappaDet[i] - omegaDet[i] * alphaStarObs[i]; 
          tmp_nObs[i] *= z[zLongIndx[i]]; 
        } // i
      } else {
        for (i = 0; i < nObs; i++) {
          kappaDet[i] = (y[i] - 1.0/2.0) * z[zLongIndx[i]];
          tmp_nObs[i] = kappaDet[i] - omegaDet[i] * alphaStarObs[i]; 
          tmp_nObs[i] *= z[zLongIndx[i]]; 
        } // i
      }
      
      F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
      for (j = 0; j < pDet; j++) {
        tmp_pDet[j] += SigmaAlphaInvMuAlpha[j]; 
      } // j

      /********************************
       * Compute A.alpha
       * *****************************/
      for (j = 0; j < nObs; j++) {
        for (i = 0; i < pDet; i++) {
          tmp_nObspDet[i*nObs + j] = Xp[i * nObs + j] * omegaDet[j] * z[zLongIndx[j]];
        } // i
      } // j

      F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

      for (j = 0; j < ppDet; j++) {
        tmp_ppDet[j] += SigmaAlphaInv[j]; 
      } // j

      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf A.alpha failed\n");}
      F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri A.alpha failed\n");}
      F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
      // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
      mvrnorm(alpha, tmp_pDet2, tmp_ppDet, pDet);

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
          for (r = 0; r < nObs; r++) {
            if ((z[zLongIndx[r]] == 1.0) && (XpRE[alphaStarIndx[l] * nObs + r] == alphaLevelIndx[l])) {
              tmp_02 = 0.0;
              for (ll = 0; ll < pDetRE; ll++) {
                tmp_02 += alphaStar[alphaStarLongIndx[ll * nObs + r]];
	      } 
              tmp_one[0] += kappaDet[r] - (F77_NAME(ddot)(&pDet, &Xp[r], &nObs, alpha, &inc) + tmp_02 - alphaStar[l]) * omegaDet[r];
      	      tmp_0 += omegaDet[r];
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
        // Update the RE sums for the current species
        for (r = 0; r < nObs; r++) {
          for (l = 0; l < pDetRE; l++) {
            alphaStarObs[r] += alphaStar[alphaStarLongIndx[l * nObs + r]]; 
          }
        }
      }
      /********************************************************************
       *Update Latent Occupancy
       *******************************************************************/
      // Compute detection probability 
      if (nObs == J) {
        for (i = 0; i < nObs; i++) {
          detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
          psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + betaStarSites[zLongIndx[i]], zero, one); 
          piProd[zLongIndx[i]] = pow(1.0 - detProb[i], K[i]);
	  piProdWAIC[zLongIndx[i]] *= pow(detProb[i], y[i]);
	  piProdWAIC[zLongIndx[i]] *= pow(1.0 - detProb[i], K[i] - y[i]);
          ySum[zLongIndx[i]] = y[i]; 	
        } // i
      } else {
        for (i = 0; i < nObs; i++) {
          detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
          if (tmp_J[zLongIndx[i]] == 0) {
            psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + betaStarSites[zLongIndx[i]], zero, one); 
          }
          piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
	  piProdWAIC[zLongIndx[i]] *= pow(detProb[i], y[i]);
	  piProdWAIC[zLongIndx[i]] *= pow(1.0 - detProb[i], 1 - y[i]);
          ySum[zLongIndx[i]] += y[i]; 	
          tmp_J[zLongIndx[i]]++;
        } // i
      }

      // Compute occupancy probability and the integrated likelihood for WAIC
      for (j = 0; j < J; j++) {
        psiNum = psi[j] * piProd[j]; 
        if (ySum[j] == zero) {
          z[j] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j])));          
	  yWAIC[j] = (1.0 - psi[j]) + psi[j] * piProdWAIC[j]; 
        } else {
          z[j] = one; 
	  yWAIC[j] = psi[j] * piProdWAIC[j];
        }
        piProd[j] = one;
	piProdWAIC[j] = one;
        ySum[j] = zero; 
        tmp_J[j] = 0; 
      } // j

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
          F77_NAME(dcopy)(&J, yWAIC, &inc, 
        		  &REAL(likeSamples_r)[sPost*J], &inc);
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
    int nResultListObjs = 5;
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
    SET_VECTOR_ELT(result_r, 4, likeSamples_r);
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 5, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 6, alphaStarSamples_r);
    }
    if (pOccRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 7;
      } else {
        tmp_0 = 5;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("like.samples"));
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("alpha.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, Rf_mkChar("beta.star.samples")); 
    }
   
    Rf_namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}



