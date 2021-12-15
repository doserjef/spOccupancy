#include <string>
#include "util.h"

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

  SEXP spPGOccPredict(SEXP J_r, SEXP pOcc_r, SEXP X0_r, SEXP q_r, 
		      SEXP obsD_r, SEXP obsPredD_r, SEXP betaSamples_r, 
		      SEXP thetaSamples_r, SEXP wSamples_r, 
		      SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, 
		      SEXP verbose_r, SEXP nReport_r){
    

    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, b, s, ii, jj, info, nProtect= 0;
    const char *lower = "L";
    const char *upper = "U";
    const char *ntran = "N";
    const char *ytran = "T";
    const char *rside = "R";
    const char *lside = "L";
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    const int inc = 1;

    /*****************************************
                     Set-up
    *****************************************/

    int J = INTEGER(J_r)[0];
    int pOcc = INTEGER(pOcc_r)[0];
    double *X0 = REAL(X0_r);
    int q = INTEGER(q_r)[0];

    double *obsD = REAL(obsD_r);
    double *obsPredD = REAL(obsPredD_r);
    double *betaSamples = REAL(betaSamples_r);
    double *thetaSamples = REAL(thetaSamples_r);
    double *wSamples = REAL(wSamples_r);
    
    int nSamples = INTEGER(nSamples_r)[0];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int nThreads = INTEGER(nThreads_r)[0]; 
    int verbose = INTEGER(verbose_r)[0]; 
    int nReport = INTEGER(nReport_r)[0];

    /*****************************************
                     Display
    *****************************************/
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tPrediction description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Spatial Occupancy model with Polya-Gamma latent\nvariable fit with %i observations.\n\n", J);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", pOcc);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Predicting at %i non-sampled locations.\n\n", q);  
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
 
    /*****************************************
         Set-up sample matrices etc.
    *****************************************/
    // parameters
    int nTheta, sigmaSqIndx,  phiIndx, nuIndx;

    if (corName != "matern") {
      nTheta = 2; //sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
      } else{
	nTheta = 3; //sigma^2, phi, nu
	sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
      }
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    double JJ = J * J; 
    double qJ = q * J; 
    
    SEXP w0_r, psi0_r, z0_r;

    PROTECT(w0_r = allocMatrix(REALSXP, q, nSamples)); nProtect++; 
    double *w0 = REAL(w0_r);

    PROTECT(psi0_r = allocMatrix(REALSXP, q, nSamples)); nProtect++; 
    double *psi0 = REAL(psi0_r);
    PROTECT(z0_r = allocMatrix(REALSXP, q, nSamples)); nProtect++; 
    double *z0 = REAL(z0_r);
    
    double *S_obs = (double *) R_alloc(JJ, sizeof(double));
    double *S_obsPred = (double *) R_alloc(qJ, sizeof(double));

    double *beta = (double *) R_alloc(pOcc, sizeof(double));
    double phi, nu, sigmaSq; 
   
    double *tmp_J = (double *) R_alloc(J, sizeof(double));  
    double *tmp_q = (double *) R_alloc(q, sizeof(double));
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_one2 = (double *) R_alloc(inc, sizeof(double)); 
    
    int status = 0;
    
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){
      
      F77_NAME(dcopy)(&pOcc, &betaSamples[s*pOcc], &inc, beta, &inc);
      phi = thetaSamples[s * nTheta + phiIndx]; 
      if (corName == "matern") {
        nu = thetaSamples[s * nTheta + nuIndx]; 
      }
      sigmaSq = thetaSamples[s * nTheta + sigmaSqIndx]; 
      theta[sigmaSqIndx] = sigmaSq; 
      theta[phiIndx] = phi; 
      theta[nuIndx] = nu; 

      // Get covariance matrices
      spCov(obsD, JJ, theta, corName, S_obs); 
      spCov(obsPredD, qJ, theta, corName, S_obsPred); 
      F77_NAME(dpotrf)(lower, &J, S_obs, &J, &info); 
      if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &J, S_obs, &J, &info); 
      if(info != 0){error("c++ error: dpotri failed\n");}	 

      F77_NAME(dgemv)(ntran, &q, &pOcc, &one, X0, &q, beta, &inc, &zero, tmp_q, &inc);
   
      // Predicting each element one at a time, instead of doing joint prediction.  
      for(j = 0; j < q; j++){

	//get Mu
	F77_NAME(dsymm)(lside, lower, &J, &inc, &one, S_obs, &J, &S_obsPred[j*J], &J, &zero, tmp_J, &J);
	F77_NAME(dgemv)(ytran, &J, &inc, &one, tmp_J, &J, &wSamples[s*J], &inc, &zero, tmp_one, &inc);
	
	//get Sigma
	F77_NAME(dgemm)(ytran, ntran, &inc, &inc, &J, &one, tmp_J, &J, &S_obsPred[j*J], &J, &zero, tmp_one2, &inc);
        tmp_one2[0] = sigmaSq - tmp_one2[0];
	w0[s * q + j] = rnorm(tmp_one[0], sqrt(tmp_one2[0])); 
	psi0[s * q + j] = logitInv(tmp_q[j] + w0[s * q + j], zero, one); 
	z0[s * q + j] = rbinom(one, psi0[s * q + j]);
      }

      //report
      if(verbose){
	if(status == nReport){
	  Rprintf("Samples: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      
      R_CheckUserInterrupt();
      
     } //end sample loop

     if(verbose){
       Rprintf("Samples: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
       #ifdef Win32
       R_FlushConsole();
       #endif
     }
     
     PutRNGstate();

     //make return object
     SEXP result, resultNames;
     int nResultListObjs = 0;
     nResultListObjs = 3;
     
     PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
     PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
          
     SET_VECTOR_ELT(result, 0, w0_r);
     SET_VECTOR_ELT(resultNames, 0, mkChar("w.0.samples"));

     SET_VECTOR_ELT(result, 1, psi0_r);
     SET_VECTOR_ELT(resultNames, 1, mkChar("psi.0.samples"));

     SET_VECTOR_ELT(result, 2, z0_r);
     SET_VECTOR_ELT(resultNames, 2, mkChar("z.0.samples"));
     
     namesgets(result, resultNames);
     
     //unprotect
     UNPROTECT(nProtect);
     
     return(result);

  }
}

