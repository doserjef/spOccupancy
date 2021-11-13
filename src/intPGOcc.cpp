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
  SEXP intPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pOcc_r, SEXP pDet_r, SEXP pDetLong_r, 
	        SEXP J_r, SEXP JLong_r, SEXP K_r, SEXP nObs_r, SEXP nObsLong_r, SEXP nData_r, 
		SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
		SEXP muBeta_r, SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP sigmaAlpha_r, 
		SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, q, info, nProtect=0;
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
    // Sorted by data set, then visit, then by site. 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    // Sorted by parameter, then data set, site, visit
    double *Xp = REAL(Xp_r);
    // Priors for regression coefficients
    double *muBeta = REAL(muBeta_r); 
    double *muAlpha = REAL(muAlpha_r); 
    double *SigmaBetaInv = REAL(SigmaBeta_r); 
    double sigmaAlpha = REAL(sigmaAlpha_r)[0]; 
    int pOcc = INTEGER(pOcc_r)[0];
    int pDet = INTEGER(pDet_r)[0];
    int nData = INTEGER(nData_r)[0]; 
    int *pDetLong = INTEGER(pDetLong_r); 
    int J = INTEGER(J_r)[0];
    int *JLong = INTEGER(JLong_r); 
    int *K = INTEGER(K_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int nObs = INTEGER(nObs_r)[0]; 
    // Rprintf("nObs: %i\n", nObs); 
    int *nObsLong = INTEGER(nObsLong_r); 
    int *dataIndx = INTEGER(dataIndx_r); 
    int *alphaIndx = INTEGER(alphaIndx_r); 
    int nSamples = INTEGER(nSamples_r)[0];
    int nThin = INTEGER(nThin_r)[0]; 
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0]; 
    int status = 0; 
    // z starting values 
    double *z = REAL(zStarting_r); 
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
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Integrated Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
      Rprintf("Integrating %i occupancy data sets.\n\n", nData); 
      Rprintf("Number of MCMC samples: %i\n", nSamples);
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
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
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
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
   
    // For latent occupancy
    double psiNum; 
    double psiNew; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); zeros(detProb, nObs);
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


    GetRNGstate();
    
    for (s = 0; s < nSamples; s++) {
    // for (s = 0; s < 1; s++) {
    // for (s = 0; s < 1; s++) {
      /********************************************************************
       *Update Occupancy Auxiliary Variables 
       *******************************************************************/
      for (j = 0; j < J; j++) {
        omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc));
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
      } // j
      /********************************
       * Compute b.beta
       *******************************/
      // X * kappaOcc + 0 * tmp_p. Output is stored in tmp_p
      F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, kappaOcc, &inc, &zero, tmp_pOcc, &inc); 	 
      for (j = 0; j < pOcc; j++) {
        // Rprintf("Value of SigmaBetaInv is %f \n", SigmaBetaInv[j]);
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
      // A.beta.inv %*% b.beta
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
        // TODO: check this if a bug arises. 
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
        // TODO: check this if need be. 
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
       *Update Latent Occupancy
       *******************************************************************/
      // Compute detection probability 
      for (i = 0; i < nObs; i++) {
        stAlpha = which(dataIndx[i], alphaIndx, pDet); 
        detProb[i] = logitInv(F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, &alpha[stAlpha], &inc), zero, one);
        if (tmp_J[zLongIndx[i]] == 0) {
          psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc), zero, one); 
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
      if (s >= nBurn) {
        thinIndx++; 
	if (thinIndx == nThin) {
          F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
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

      // This allows for users to interrupt the C code when running, 
      // which is not normally allowed. 
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
    int nResultListObjs = 5;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, yRepSamples_r);
    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("y.rep.samples")); 
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

