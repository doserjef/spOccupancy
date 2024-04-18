#define USE_FC_LEN_T
#include <string>
#include "util.h"

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

  SEXP svcPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP pOcc_r, SEXP pTilde_r, 
		           SEXP m_r, SEXP X0_r, SEXP Xw0_r, SEXP coords0_r, 
			   SEXP weights0_r,
			   SEXP JStr_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			   SEXP thetaSamples_r, SEXP wSamples_r, 
			   SEXP betaStarSiteSamples_r, SEXP nSamples_r, 
			   SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			   SEXP nReport_r, SEXP Jw0_r, SEXP Jw_r, SEXP gridIndx0_r, 
			   SEXP sitesLink_r, SEXP sites0Sampled_r){

    int i, j, ll, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    
    double *coords = REAL(coords_r);
    int J = INTEGER(J_r)[0];
    int pOcc = INTEGER(pOcc_r)[0];
    int pTilde = INTEGER(pTilde_r)[0];
    int Jw0 = INTEGER(Jw0_r)[0];
    int Jw = INTEGER(Jw_r)[0];

    double *X0 = REAL(X0_r);
    double *Xw0 = REAL(Xw0_r);
    double *coords0 = REAL(coords0_r);
    double *weights0 = REAL(weights0_r);
    int JStr = INTEGER(JStr_r)[0];
    int m = INTEGER(m_r)[0]; 
    int mm = m * m; 
    int Jw0pTilde = Jw0 * pTilde; 
    int JwpTilde = Jw * pTilde;

    int *nnIndx0 = INTEGER(nnIndx0_r);        
    double *beta = REAL(betaSamples_r);
    double *theta = REAL(thetaSamples_r);
    double *w = REAL(wSamples_r);
    double *betaStarSite = REAL(betaStarSiteSamples_r);
    int *gridIndx0 = INTEGER(gridIndx0_r);
    int *sitesLink = INTEGER(sitesLink_r);
    int *sites0Sampled = INTEGER(sites0Sampled_r);
    
    int nSamples = INTEGER(nSamples_r)[0];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tPrediction description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("NNGP Occupancy model with Polya-Gamma latent\nvariable fit with %i observations.\n\n", J);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", pOcc);
      Rprintf("Number of spatially-varying covariates %i (including intercept if specified).\n\n", pTilde);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Predicting at %i non-sampled locations.\n\n", JStr);  
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    // parameters
    int nTheta, sigmaSqIndx,  phiIndx, nuIndx;

    if (corName != "matern") {
      nTheta = 2; //sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
      } else{
	nTheta = 3; //sigma^2, phi, nu
	sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
      }
    int nThetapTilde = nTheta * pTilde;
    
    // get max nu
    double *nuMax = (double *) R_alloc(pTilde, sizeof(double));
    int *nb = (int *) R_alloc(pTilde, sizeof(nb)); 
    // Fill in with zeros
    for (ll = 0; ll < pTilde; ll++) {
      nuMax[ll] = 0.0; 
      nb[ll] = 0; 
    }
    
    if(corName == "matern"){
      for (ll = 0; ll < pTilde; ll++) {
        for(s = 0; s < nSamples; s++){
          if(theta[s*nThetapTilde + nuIndx * pTilde + ll] > nuMax[ll]){
            nuMax[ll] = theta[s*nThetapTilde + nuIndx * pTilde + ll];
          }
        }
        nb[ll] = 1+static_cast<int>(floor(nuMax[ll]));
      }
    }

    int nbMax = 0; 
    for (ll = 0; ll < pTilde; ll++) {
      if (nb[ll] > nbMax) {
        nbMax = nb[ll]; 
      }
    }

    double *bk = (double *) R_alloc(nThreads*nbMax, sizeof(double));
   
    double *C = (double *) R_alloc(nThreads*mm, sizeof(double)); zeros(C, nThreads*mm);
    double *c = (double *) R_alloc(nThreads*m, sizeof(double)); zeros(c, nThreads*m);
    double *tmp_m  = (double *) R_alloc(nThreads*m, sizeof(double));
    double phi = 0, nu = 0, sigmaSq = 0, d;
    int threadID = 0, status = 0;

    SEXP z0_r, w0_r, psi0_r;
    PROTECT(z0_r = Rf_allocMatrix(REALSXP, JStr, nSamples)); nProtect++; 
    PROTECT(psi0_r = Rf_allocMatrix(REALSXP, JStr, nSamples)); nProtect++; 
    PROTECT(w0_r = Rf_allocMatrix(REALSXP, Jw0pTilde, nSamples)); nProtect++;
    double *z0 = REAL(z0_r);
    double *psi0 = REAL(psi0_r); 
    double *w0 = REAL(w0_r);
    double wSites; 
 
    if (verbose) {
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tPredicting\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    int vIndx = -1;
    double *wV = (double *) R_alloc(Jw0pTilde*nSamples, sizeof(double));

    GetRNGstate();
    
    for(i = 0; i < Jw0pTilde*nSamples; i++){
      wV[i] = rnorm(0.0,1.0);
    }
    
    for(j = 0; j < Jw0; j++){
      for (ll = 0; ll < pTilde; ll++) {
#ifdef _OPENMP
#pragma omp parallel for private(threadID, phi, nu, sigmaSq, k, l, d, info)
#endif     
        for(s = 0; s < nSamples; s++){
#ifdef _OPENMP
	  threadID = omp_get_thread_num();
#endif 	
	  if (sites0Sampled[j] == 1) {
	    w0[s * Jw0pTilde + j * pTilde + ll] = w[s * JwpTilde + sitesLink[j] * pTilde + ll];
	  } else {
	    phi = theta[s * nThetapTilde + phiIndx * pTilde + ll];
	    if(corName == "matern"){
	      nu = theta[s * nThetapTilde + nuIndx * pTilde + ll];
	    }
	    sigmaSq = theta[s * nThetapTilde + sigmaSqIndx * pTilde + ll];

	    for(k = 0; k < m; k++){
	      d = dist2(coords[nnIndx0[j+Jw0*k]], coords[Jw+nnIndx0[j+Jw0*k]], coords0[j], coords0[Jw0+j]);
	      c[threadID*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb[ll]]);
	      for(l = 0; l < m; l++){
	        d = dist2(coords[nnIndx0[j+Jw0*k]], coords[Jw+nnIndx0[j+Jw0*k]], coords[nnIndx0[j+Jw0*l]], coords[Jw+nnIndx0[j+Jw0*l]]);
	        C[threadID*mm+l*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb[ll]]);
	      }
	    }

	    F77_NAME(dpotrf)(lower, &m, &C[threadID*mm], &m, &info FCONE); 
	    if(info != 0){Rf_error("c++ error: dpotrf failed\n");}
	    F77_NAME(dpotri)(lower, &m, &C[threadID*mm], &m, &info FCONE); 
	    if(info != 0){Rf_error("c++ error: dpotri failed\n");}

	    F77_NAME(dsymv)(lower, &m, &one, &C[threadID*mm], &m, &c[threadID*m], &inc, &zero, &tmp_m[threadID*m], &inc FCONE);

	    d = 0;
	    for(k = 0; k < m; k++){
	      d += tmp_m[threadID*m+k]*w[s*JwpTilde+nnIndx0[j+Jw0*k] * pTilde + ll];
	    }

	    #ifdef _OPENMP
            #pragma omp atomic
            #endif   
	    vIndx++;

	    w0[s * Jw0pTilde + j * pTilde + ll] = sqrt(sigmaSq - F77_NAME(ddot)(&m, &tmp_m[threadID*m], &inc, &c[threadID*m], &inc))*wV[vIndx] + d;
	  }

        } // sample
      } // covariate
      
      if(verbose){
	if(status == nReport){
	  Rprintf("Location: %i of %i, %3.2f%%\n", j, Jw0, 100.0*j/Jw0);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      R_CheckUserInterrupt();
    } // location
    
    if(verbose){
      Rprintf("Location: %i of %i, %3.2f%%\n", j, Jw0, 100.0*j/Jw0);
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    // Generate latent occurrence state after the fact.
    if (verbose) {
      Rprintf("Generating latent occupancy state\n");
    }
    for(j = 0; j < JStr; j++){
      for(s = 0; s < nSamples; s++){
        wSites = F77_NAME(ddot)(&pTilde, &Xw0[j], &JStr, 
			        &w0[s * Jw0pTilde + gridIndx0[j] * pTilde], &inc);
	psi0[s * JStr + j] = logitInv(F77_NAME(ddot)(&pOcc, &X0[j], &JStr, 
				                     &beta[s*pOcc], &inc) + 
			              wSites + betaStarSite[s * JStr + j], zero, one);
        z0[s * JStr + j] = rbinom(weights0[j], psi0[s * JStr + j]);
      } // s
    } // i

    PutRNGstate();
    

    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 3;

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, z0_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("z.0.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, w0_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("w.0.samples"));

    SET_VECTOR_ELT(result_r, 2, psi0_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("psi.0.samples")); 

    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  
  }
}
