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

  SEXP sfMsPGOccNNGPPredict(SEXP coords_r, SEXP J_r, SEXP N_r, SEXP q_r,
                            SEXP pOcc_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, 
                            SEXP JStr_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
                            SEXP thetaSamples_r, SEXP lambdaSamples_r, 
                            SEXP wSamples_r, SEXP betaStarSiteSamples_r, 
                            SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
                            SEXP nReport_r, SEXP sharedSpatial_r, SEXP sitesLink_r, 
                            SEXP sites0Sampled_r, SEXP Jw0_r, SEXP Jw_r, SEXP gridIndx0_r){

    int i, j, k, l, s, info, nProtect=0, ll;
    const int inc = 1;
    const double one = 1.0;
    char const *ntran = "N";
    const double zero = 0.0;
    char const *lower = "L";
    
    double *coords = REAL(coords_r);
    int J = INTEGER(J_r)[0];
    int N = INTEGER(N_r)[0];
    int q = INTEGER(q_r)[0]; 
    int pOcc = INTEGER(pOcc_r)[0];
    int Jw0 = INTEGER(Jw0_r)[0];
    int Jw = INTEGER(Jw_r)[0];
    int pOccN = pOcc * N;
    int Nq = N * q;

    double *X0 = REAL(X0_r);
    double *coords0 = REAL(coords0_r);
    int JStr = INTEGER(JStr_r)[0];
    int JStrN = JStr * N;
    int Jw0N = Jw0 * N;
    int m = INTEGER(m_r)[0]; 
    int mm = m * m; 
    int *gridIndx0 = INTEGER(gridIndx0_r);
    int JStrpOcc = JStr * pOcc;

    int *nnIndx0 = INTEGER(nnIndx0_r);        
    double *beta = REAL(betaSamples_r);
    double *theta = REAL(thetaSamples_r);
    double *lambda = REAL(lambdaSamples_r);
    double *w = REAL(wSamples_r);
    double *betaStarSite = REAL(betaStarSiteSamples_r);
    int *sitesLink = INTEGER(sitesLink_r);
    int *sites0Sampled = INTEGER(sites0Sampled_r);
    int Jwq = Jw * q;
    int Jw0q = Jw0 * q;
    
    int nSamples = INTEGER(nSamples_r)[0];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int sharedSpatial = INTEGER(sharedSpatial_r)[0];
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tPrediction description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Spatial Factor NNGP Multispecies Occupancy model with Polya-Gamma latent\nvariable fit with %i observations.\n\n", J);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", pOcc);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n", m);
      if (sharedSpatial == 0) {
        Rprintf("Using %i latent spatial factors.\n\n", q);
      } else {
          Rprintf("Using a common spatial process across species.\n\n");
      }
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Predicting at %i non-sampled locations.\n\n", JStr);  
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    // parameters
    int nTheta, phiIndx, nuIndx, sigmaSqIndx;

    // NOTE: this differs from other non-factor modeling functions
    if (sharedSpatial == 0) {
    if (corName != "matern") {
      nTheta = 1; //phi
      phiIndx = 0;
    } else{
	nTheta = 2; //phi, nu
	phiIndx = 0; nuIndx = 1;
      }
    } else {
      if (corName != "matern") {
        nTheta = 2; //sigmaSq, phi
        sigmaSqIndx = 0; phiIndx = 1;
      } else{
          nTheta = 3; //sigmaSq, phi, nu
          sigmaSqIndx = 0; phiIndx = 2; nuIndx = 2;
        }

    }
   
    int nThetaq = nTheta * q; 
    // get max nu
    double *nuMax = (double *) R_alloc(q, sizeof(double));
    int *nb = (int *) R_alloc(q, sizeof(nb)); 
    // Fill in with zeros
    for (ll = 0; ll < q; ll++) {
      nuMax[ll] = 0.0; 
      nb[ll] = 0; 
    }
    
    if(corName == "matern"){
      for (ll = 0; ll < q; ll++) {
        for(s = 0; s < nSamples; s++){
          if(theta[s*nThetaq + nuIndx * q + ll] > nuMax[ll]){
            nuMax[ll] = theta[s*nThetaq + nuIndx * q + ll];
          }
        }
        nb[ll] = 1+static_cast<int>(floor(nuMax[ll]));
      }
    }

    int nbMax = 0; 
    for (ll = 0; ll < q; ll++) {
      if (nb[ll] > nbMax) {
        nbMax = nb[ll]; 
      }
    }

    double *bk = (double *) R_alloc(nThreads*nbMax, sizeof(double));
    double *C = (double *) R_alloc(nThreads*mm, sizeof(double)); zeros(C, nThreads*mm);
    double *c = (double *) R_alloc(nThreads*m, sizeof(double)); zeros(c, nThreads*m);
    double *tmp_m  = (double *) R_alloc(nThreads*m, sizeof(double));
    double phi = 0, nu = 0, sigmaSq = 1.0, d;
    int threadID = 0, status = 0;

    SEXP z0_r, w0_r, psi0_r;
    PROTECT(z0_r = Rf_allocMatrix(REALSXP, JStrN, nSamples)); nProtect++; 
    PROTECT(psi0_r = Rf_allocMatrix(REALSXP, JStrN, nSamples)); nProtect++; 
    PROTECT(w0_r = Rf_allocMatrix(REALSXP, Jw0q, nSamples)); nProtect++;
    double *z0 = REAL(z0_r);
    double *psi0 = REAL(psi0_r); 
    double *w0 = REAL(w0_r);
    double *w0Star = (double *) R_alloc(nSamples * Jw0N, sizeof(double));
    if (verbose) {
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tPredicting\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    int vIndx = -1;
    double *wV = (double *) R_alloc(Jw0q*nSamples, sizeof(double));

    GetRNGstate();
    
    for(i = 0; i < Jw0q*nSamples; i++){
      wV[i] = rnorm(0.0,1.0);
    }
    
    for(j = 0; j < Jw0; j++){
      for (ll = 0; ll < q; ll++) {
#ifdef _OPENMP
#pragma omp parallel for private(threadID, phi, nu, sigmaSq, k, l, d, info)
#endif     
        for(s = 0; s < nSamples; s++){
#ifdef _OPENMP
	  threadID = omp_get_thread_num();
#endif 	
	  if (sites0Sampled[j] == 1) {
            w0[s * Jw0q + j * q + ll] = w[s * Jwq + sitesLink[j] * q + ll];
	  } else {
	    phi = theta[s * nThetaq + phiIndx * q + ll];
	    if(corName == "matern"){
	      nu = theta[s * nThetaq + nuIndx * q + ll];
	    }
	    if (sharedSpatial == 0) {
	      sigmaSq = 1.0;
	    } else {
              sigmaSq = theta[s * nThetaq + sigmaSqIndx * q + ll];
	    }

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
	      d += tmp_m[threadID*m+k]*w[s*Jwq+nnIndx0[j+Jw0*k] * q + ll];
	    }

	    #ifdef _OPENMP
            #pragma omp atomic
            #endif   
	    vIndx++;

	    w0[s * Jw0q + j * q + ll] = sqrt(sigmaSq - F77_NAME(ddot)(&m, &tmp_m[threadID*m], &inc, &c[threadID*m], &inc))*wV[vIndx] + d;
	  }

        } // sample
      } // factor

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
    // Temporary fix. Will embed this in the above loop at some point.
    if (verbose) {
      Rprintf("Generating latent occupancy state\n");
    }
    for (j = 0; j < JStr; j++) {
      for (s = 0; s < nSamples; s++) {
          F77_NAME(dgemv)(ntran, &N, &q, &one, &lambda[s * Nq], &N, &w0[s * Jw0q + gridIndx0[j]*q], &inc, &zero, &w0Star[s * Jw0N + gridIndx0[j] * N], &inc FCONE);
	  for (i = 0; i < N; i++) {
	    psi0[s * JStrN + j * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X0[i * JStrpOcc + j], &JStr, &beta[s*pOccN + i], &N) + w0Star[s * Jw0N + gridIndx0[j] * N + i] + betaStarSite[s*JStrN + j * N + i], zero, one);
	    z0[s * JStrN + j * N + i] = rbinom(one, psi0[s * JStrN + j * N + i]);
	  } // i
      } // s
      R_CheckUserInterrupt();
    } // j

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

    
