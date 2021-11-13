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
  SEXP PGOccREDet(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XpRE_r, SEXP lambdaP_r, 
		  SEXP pocc_r, SEXP pdet_r, SEXP pDetRE_r,  
                  SEXP J_r, SEXP nObs_r, SEXP K_r, SEXP nDetRE_r, SEXP nDetRELong_r, 
		  SEXP betaStarting_r, SEXP alphaStarting_r, 
	          SEXP sigmaSqPStarting_r, SEXP alphaStarStarting_r, 
	          SEXP zStarting_r, SEXP zLongIndx_r, 
		  SEXP alphaStarIndx_r, SEXP muBeta_r, SEXP muAlpha_r, 
		  SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP sigmaSqPA_r, 
	          SEXP sigmaSqPB_r, SEXP nSamples_r, SEXP nThreads_r, 
		  SEXP verbose_r, SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, 
		  SEXP nPost_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, l, ll, info, nProtect=0;
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
    int *XpRE = INTEGER(XpRE_r); 
    double *lambdaP = REAL(lambdaP_r); 
    double *muBeta = REAL(muBeta_r); 
    double *muAlpha = REAL(muAlpha_r); 
    double *SigmaBetaInv = REAL(SigmaBeta_r); 
    double *SigmaAlphaInv = REAL(SigmaAlpha_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
    int pDetRE = INTEGER(pDetRE_r)[0]; 
    int nDetRE = INTEGER(nDetRE_r)[0]; 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int J = INTEGER(J_r)[0];
    double *K = REAL(K_r); 
    int nObs = INTEGER(nObs_r)[0];
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int nThin = INTEGER(nThin_r)[0]; 
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    int status = 0; 
    double *z = REAL(zStarting_r); 
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
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Occupancy model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
      Rprintf("Number of MCMC samples: %i \n", nSamples);
      Rprintf("Burn-in: %i \n", nBurn); 
      Rprintf("Thinning Rate: %i \n", nThin); 
      Rprintf("Total Posterior Samples: %i \n\n", nPost); 
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
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
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
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
    SEXP sigmaSqPSamples_r; 
    PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
    SEXP alphaStarSamples_r; 
    PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    int nnDetRE = nDetRE * nDetRE; 
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    double tmp_0; 
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
   
    // For latent occupancy
    double psiNum; 
    double psiNew; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); zeros(detProb, nObs);
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J);

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, SigmaBetaInvMuBeta, &inc); 	  
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, SigmaAlphaInvMuAlpha, &inc); 	  


    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(alphaStarObs, nObs); 
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarObs[i] += alphaStar[XpRE[l * nObs + i]];
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
       *Update Occupancy Auxiliary Variables 
       *******************************************************************/
      for (j = 0; j < J; j++) {
        omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc));
      } // j

      /********************************************************************
       *Update Detection Auxiliary Variables 
       *******************************************************************/
      if (nObs == J) {
        for (i = 0; i < nObs; i++) {
          omegaDet[i] = rpg(K[i], F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i]);
        } // i
      } else {
        for (i = 0; i < nObs; i++) {
          omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i]);
        } // i
      }

      /********************************************************************
       *Update Occupancy Regression Coefficients
       *******************************************************************/
      for (j = 0; j < J; j++) {
	kappaOcc[j] = z[j] - 1.0 / 2.0; 
      } // j

      /********************************
       * Compute b.beta
       *******************************/
      F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, kappaOcc, &inc, &zero, tmp_pOcc, &inc); 	 
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
      F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);

      /********************************************************************
       *Update Detection Regression Coefficients
       *******************************************************************/
      // /********************************
      //  * Compute b.alpha
      //  *******************************/
      // First multiply kappDet * the current occupied values, such that values go 
      // to 0 if z == 0 and values go to kappaDet if z == 1
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
      
      F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc); 	  
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

      F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet);

      for (j = 0; j < ppDet; j++) {
        tmp_ppDet[j] += SigmaAlphaInv[j]; 
      } // j

      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
      F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
      F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc);
      // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      mvrnorm(alpha, tmp_pDet2, tmp_ppDet, pDet);

      /********************************************************************
       *Update Detection random effects variance
       *******************************************************************/
      for (l = 0; l < pDetRE; l++) {
        tmp_0 = F77_NAME(ddot)(&nDetRELong[l], &alphaStar[alphaStarStart[l]], &inc, &alphaStar[alphaStarStart[l]], &inc); 
        tmp_0 *= 0.5; 
        sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] / 2.0, sigmaSqPB[l] + tmp_0); 
      }

      /********************************************************************
       *Update Detection random effects
       *******************************************************************/
      // Update each individual random effect one by one. 
      for (l = 0; l < nDetRE; l++) {
        /********************************
         * Compute b.alpha.star
         *******************************/
        for (i = 0; i < nObs; i++) {
          tmp_nObs[i] = kappaDet[i] - (F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i] - alphaStar[l]) * omegaDet[i];
          // Only allow information to come from when z == 1.
          tmp_nObs[i] *= z[zLongIndx[i]]; 
        }
        F77_NAME(dgemv)(ytran, &nObs, &inc, &one, &lambdaP[l * nObs], &nObs, tmp_nObs, &inc, &zero, tmp_one, &inc); 
        /********************************
         * Compute A.alpha.star
         *******************************/
        for (i = 0; i < nObs; i++) {
          tmp_nObs[i] = lambdaP[l * nObs + i] * omegaDet[i] * z[zLongIndx[i]]; 
        }
        tmp_0 = F77_NAME(ddot)(&nObs, tmp_nObs, &inc, &lambdaP[l * nObs], &inc); 
        tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
        tmp_0 = 1.0 / tmp_0; 
        alphaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
      }
      
      zeros(alphaStarObs, nObs); 
      // Get sums of the current REs for each site/visit combo
      for (i = 0; i < nObs; i++) {
        for (l = 0; l < pDetRE; l++) {
          alphaStarObs[i] += alphaStar[XpRE[l * nObs + i]];
        }
      }


      /********************************************************************
       *Update Latent Occupancy
       *******************************************************************/
      // Compute detection probability 
      if (nObs == J) {
        for (i = 0; i < nObs; i++) {
          detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
          psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc), zero, one); 
          piProd[zLongIndx[i]] *= pow(1.0 - detProb[i], K[i]);
          ySum[zLongIndx[i]] = y[i]; 	
        } // i
      } else {
        for (i = 0; i < nObs; i++) {
          detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
          if (tmp_J[zLongIndx[i]] == 0) {
            psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc), zero, one); 
          }
          piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
          ySum[zLongIndx[i]] += y[i]; 	
          tmp_J[zLongIndx[i]]++;
        } // i
      }
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
      if (s >= nBurn) {
        thinIndx++; 
	if (thinIndx == nThin) {
          F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
          F77_NAME(dcopy)(&J, z, &inc, &REAL(zSamples_r)[sPost*J], &inc); 
          F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
          F77_NAME(dcopy)(&nDetRE, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetRE], &inc);
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
    int nResultListObjs = 6;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, sigmaSqPSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaStarSamples_r);

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("sigma.sq.p.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.star.samples")); 
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

