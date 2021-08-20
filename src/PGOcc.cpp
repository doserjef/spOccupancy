#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#include "rpg.h"

// Optionally include OPENMP for parallelization if it exists. 
#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
  SEXP PGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pocc_r, SEXP pdet_r, 
             SEXP J_r, SEXP K_r, SEXP betaStarting_r, SEXP alphaStarting_r, 
	     SEXP zStarting_r, SEXP zLongIndx_r, SEXP muBeta_r, 
	     SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    // These are specified as pointers since functions require the addresses.
    char const *lower = "L";
    char const *upper = "U";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *rside = "R";
    char const *lside = "L";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    // The REAL or INTEGER are helper functions that allow you to access
    // the C array inside the R objects that are read in as inputs or 
    // created in the function. 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    // Xp is sorted by parameter then visit (parameter 1 v1, p1v2, etc..) 
    double *Xp = REAL(Xp_r);
    // Priors for regression coefficients
    double *muBeta = REAL(muBeta_r); 
    double *muAlpha = REAL(muAlpha_r); 
    double *SigmaBetaInv = REAL(SigmaBeta_r); 
    double *SigmaAlphaInv = REAL(SigmaAlpha_r); 
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
    int J = INTEGER(J_r)[0];
    int *K = INTEGER(K_r); 
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
    // z starting values 
    double *z = REAL(zStarting_r); 

    // Rprintf("Value of nObs is %i\n\n", nObs);

// For parallelization.  
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
      Rprintf("Occupancy model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
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
    // This copies the starting values provided as user input into beta.  
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
    // Create an R-level matrix. The PROTECT is necessary to ensure that 
    // the R objects you want for output are not deleted even if the garbage
    // collector is activated. 
    // The nProtect is used to track the number of protected objects, which 
    // is added to as additional objects are protected. 
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOcc, nSamples)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDet, nSamples)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, J, nSamples)); nProtect++; 
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, J, nSamples)); nProtect++; 
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(INTSXP, nObs, nSamples)); nProtect++; 
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    // R_alloc is used to allocate memory. 
    // The memory allocated with R_alloc is automatically released when 
    // R returns from .Call. 
    // R_alloc is used when one wants to rrepresent native c data types
    // rather than R objects. 
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
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    int *ySum = (int *) R_alloc(J, sizeof(int)); 
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
    // Compute cholesky
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    // Compute inverse
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, SigmaAlphaInvMuAlpha, &inc); 	  


    // This is necessary for generating random numbers in C 
    GetRNGstate();
    
    for (s = 0; s < nSamples; s++) {
    // for (s = 0; s < 10; s++) {
    // for (s = 0; s < 1; s++) {
      /********************************************************************
       *Update Occupancy Auxiliary Variables 
       *******************************************************************/
      for (j = 0; j < J; j++) {
        // ddot forms the dot product of two vectors. Note how the third argument
        // of ddot that is the storage spacing. So the elements grabbed from X are
        // X[j], X[j + J],  which is a row of the design matrix. 
        omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc));
      } // j
      /********************************************************************
       *Update Detection Auxiliary Variables 
       *******************************************************************/
      // Note that all of the variables are sampled, but only those at 
      // locations with z[j] == 1 actually effect the results. 
      for (i = 0; i < nObs; i++) {
        omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc));
	// Rprintf("Value of omega currently %f \n", omegaDet[i]);
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
      // dgemv computes linear combinations of different variables. 
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

      // This gives the Cholesky of A.beta
      // Computes cholesky of tmp_ppOcc. Output stored in tmp_ppOcc
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      // Computes the inverse tmp_ppOcc. Stored in tmp_ppOcc. This is A.beta.inv. 
      F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotri here failed\n");}
      // A.beta.inv %*% b.beta
      // 1 * tmp_ppOcc * tmp_pOcc + 0 * tmp_pOcc2 
      // (which is currently nothing) = tmp_pOcc2
      F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
      // Computes cholesky of tmp_pp again stored back in tmp_ppOcc. This chol(A.beta.inv)
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); if(info != 0){error("c++ error: dpotrf here failed\n");}
      // Args: destination, mu, cholesky of the covariance matrix, dimension
      mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
      // for (j = 0; j < pOcc; j++) {
      //   Rprintf("Beta: %f\n", beta[j]); 
      // } // j
      
      /********************************************************************
       *Update Detection Regression Coefficients
       *******************************************************************/
      // /********************************
      //  * Compute b.alpha
      //  *******************************/
      // First multiply kappDet * the current occupied values, such that values go 
      // to 0 if they z == 0 and values go to kappaDet if z == 1
      for (i = 0; i < nObs; i++) {
        // 1.0 is currently hardcoded in for occupancy data
        kappaDet[i] = (y[i] - 1.0/2.0) * z[zLongIndx[i]];
        // Rprintf("kappa value %i is %f \n", i, kappaDet[i]);
        // Rprintf("zLongIndx value %i is %i \n", i, zLongIndx[i]);
      } // i
      
      // Xp * kappaDet + 0 * tmp_pDet. Output is stored in tmp_pDet
      // dgemv computes linear combinations of different variables. 
      F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, kappaDet, &inc, &zero, tmp_pDet, &inc); 	  
      for (j = 0; j < pDet; j++) {
        tmp_pDet[j] += SigmaAlphaInvMuAlpha[j]; 
	// Rprintf("b.alpha: %f\n", tmp_pDet[j]); 
      } // j

      /********************************
       * Compute A.alpha
       * *****************************/
      for (j = 0; j < nObs; j++) {
        for (i = 0; i < pDet; i++) {
          tmp_nObspDet[i*nObs + j] = Xp[i * nObs + j] * omegaDet[j] * z[zLongIndx[j]];
	  //Rprintf("tmp_nObspDet: %f\n", tmp_nObspDet[i*nObs + j]);  
        } // i
      } // j

      // This finishes off A.alpha
      // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
      F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet);

      for (j = 0; j < ppDet; j++) {
        tmp_ppDet[j] += SigmaAlphaInv[j]; 
	// Rprintf("tmp_ppDet: %f\n", tmp_ppDet[j]); 
      } // j

      // This gives the Cholesky of A.alpha
      // Computes cholesky of tmp_ppDet. Output stored in tmp_ppOcc
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
      // Computes the inverse tmp_ppOcc. Stored in tmp_ppOcc. This is A.beta.inv. 
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
      mvrnorm(alpha, tmp_pDet2, tmp_ppDet, pDet);
      // for (j = 0; j < pDet; j++) {
      //   Rprintf("alpha: %f\n", alpha[j]); 
      // } // j

     
      /********************************************************************
       *Update Latent Occupancy
       *******************************************************************/
      // Compute detection probability 
      for (i = 0; i < nObs; i++) {
        detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc), zero, one);
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
        // Save z samples along the way. 
        REAL(zSamples_r)[s * J + j] = z[j]; 
        piProd[j] = one;
        ySum[j] = zero; 
	tmp_J[j] = 0; 
      } // j

      /********************************************************************
       *Replicate data set for GoF
       *******************************************************************/
      for (i = 0; i < nObs; i++) {
        yRep[i] = rbinom(one, detProb[i] * z[zLongIndx[i]]);
	INTEGER(yRepSamples_r)[s * nObs + i] = yRep[i]; 
      } // i


      /********************************************************************
       *Save samples
       *******************************************************************/
      F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[s*pOcc], &inc);
      F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[s*pDet], &inc);
      F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[s*J], &inc); 

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

