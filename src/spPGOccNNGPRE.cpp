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

void updateBF1RE(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';

  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuUnifb));
  int threadID = 0;
  double e;
  int mm = m*m;

#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
    for(i = 0; i < n; i++){
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
      if(i > 0){
	for(k = 0; k < nnIndxLU[n+i]; k++){
	  e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
	  c[m*threadID+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
	    C[mm*threadID+l*nnIndxLU[n+i]+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  }
	}
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
  SEXP spPGOccNNGPRE(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, SEXP XpRE_r, 
		     SEXP lambdaP_r, SEXP pocc_r, SEXP pdet_r, SEXP pDetRE_r,
	             SEXP J_r, SEXP nObs_r, SEXP K_r, SEXP nDetRE_r, SEXP nDetRELong_r, 
		     SEXP m_r, SEXP nnIndx_r, 
		     SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		     SEXP betaStarting_r, SEXP alphaStarting_r, 
		     SEXP sigmaSqPStarting_r, SEXP alphaStarStarting_r, 
	             SEXP zStarting_r, SEXP phiStarting_r,
	             SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
	             SEXP zLongIndx_r, SEXP alphaStarIndx_r, 
		     SEXP muBeta_r, SEXP muAlpha_r, 
	             SEXP SigmaBeta_r, SEXP SigmaAlpha_r, SEXP phiA_r, SEXP phiB_r, 
	             SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		     SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		     SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	             SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	             SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, l, k, s, r, q, info, nProtect=0;
    int status = 0; // For AMCMC. 
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
    int m = INTEGER(m_r)[0]; 
    // Priors
    double *muBeta = REAL(muBeta_r); 
    double *muAlpha = REAL(muAlpha_r); 
    double *SigmaBetaInv = REAL(SigmaBeta_r); 
    double *SigmaAlphaInv = REAL(SigmaAlpha_r); 
    double phiA = REAL(phiA_r)[0];
    double phiB = REAL(phiB_r)[0]; 
    double nuA = REAL(nuA_r)[0]; 
    double nuB = REAL(nuB_r)[0]; 
    double sigmaSqA = REAL(sigmaSqA_r)[0]; 
    double sigmaSqB = REAL(sigmaSqB_r)[0]; 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double *tuning = REAL(tuning_r); 
    double *coords = REAL(coords_r);
    int pDetRE = INTEGER(pDetRE_r)[0]; 
    int nDetRE = INTEGER(nDetRE_r)[0]; 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
    int J = INTEGER(J_r)[0];
    double *K = REAL(K_r); 
    int nObs = INTEGER(nObs_r)[0];
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nThin = INTEGER(nThin_r)[0]; 
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    double acceptRate = REAL(acceptRate_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    double *z = REAL(zStarting_r); 
    int nReport = INTEGER(nReport_r)[0];
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
      Rprintf("NNGP Occupancy model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
      Rprintf("Number of MCMC samples: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
      Rprintf("Burn-in: %i \n", nBurn); 
      Rprintf("Thinning Rate: %i \n", nThin); 
      Rprintf("Total Posterior Samples: %i \n\n", nPost); 
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
#ifdef _OPENMP
      Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      Rprintf("Adaptive Metropolis with target acceptance rate: %.1f\n", 100*acceptRate);
      Rprintf("Sampling ... \n");
      #ifdef Win32
        R_FlushConsole();
      #endif

    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    double *w = (double *) R_alloc(J, sizeof(double)); zeros(w, J);   
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    double nu = REAL(nuStarting_r)[0]; 
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
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
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
    int JpOcc = J * pOcc; 
    int JJ = J * J; 
    int nObspDet = nObs * pDet;
    int jj, kk;
    double tmp_0; 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_one = (double *) R_alloc(1, sizeof(double)); 
    double * tmp_JJ = (double *) R_alloc(JJ, sizeof(double)); 
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = zero; 
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
   
    // For latent occupancy
    double psiNum; 
    double psiNew; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
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
    F77_NAME(dgemv)(ytran, &pOcc, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, SigmaBetaInvMuBeta, &inc); 	  
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, SigmaAlphaInvMuAlpha, &inc); 	  

    /**********************************************************************
     * Prep for detection random effects
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
    double logMHRatio, logPostCurrent = 0.0, logPostCand = 0.0, detCand = 0.0, detCurr = 0.0;
    double logDet;  
    double phiCand = 0.0, nuCand = 0.0;  
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nTheta, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nTheta, nBatch)); nProtect++; 
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nPost)); nProtect++; 
    double a, v, b, e, mu, var, aij; 
    // Initiate spatial values
    theta[sigmaSqIndx] = REAL(sigmaSqStarting_r)[0]; 
    theta[phiIndx] = REAL(phiStarting_r)[0]; 
    if (corName == "matern") {
      theta[nuIndx] = nu; 
    } 
    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(J, sizeof(double));
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads, sizeof(double));

    double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuB))), sizeof(double));

    if (corName == "matern") {
      nu = theta[nuIndx];
    }
    updateBF1RE(B, F, c, C, coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx], theta[phiIndx], nu, covModel, bk, nuB);


    GetRNGstate();
   
    /**********************************************************************
     * Begin Sampler 
     * *******************************************************************/
    for (s = 0, q = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, q++) {
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
          tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * w[j]; 
        } // j
        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc); 	 
        for (j = 0; j < pOcc; j++) {
          tmp_pOcc[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
        for(j = 0; j < J; j++){
          for(i = 0; i < pOcc; i++){
            tmp_JpOcc[i*J+j] = X[i*J+j]*omegaOcc[j];
          }
        }

        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc);
        for (j = 0; j < ppOcc; j++) {
          tmp_ppOcc[j] += SigmaBetaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotri here failed\n");}
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); if(info != 0){error("c++ error: dpotrf here failed\n");}
        mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
        
        /********************************************************************
         *Update Detection Regression Coefficients
         *******************************************************************/
        // /********************************
        //  * Compute b.alpha
        //  *******************************/
        // First multiply kappDet * the current occupied values, such that values go 
        // to 0 if they z == 0 and values go to kappaDet if z == 1
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
         *Update w (spatial random effects)
         *******************************************************************/
	for (i = 0; i < J; i++ ) {
          a = 0;
	  v = 0;
	  if (uIndxLU[J + i] > 0){ // is i a neighbor for anybody
	    for (j = 0; j < uIndxLU[J+i]; j++){ // how many locations have i as a neighbor
	      b = 0;
	      // now the neighbors for the jth location who has i as a neighbor
	      jj = uIndx[uIndxLU[i]+j]; // jj is the index of the jth location who has i as a neighbor
	      for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	        kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	        if(kk != i){ //if the neighbor of jj is not i
	  	b += B[nnIndxLU[jj]+k]*w[kk]; //covariance between jj and kk and the random effect of kk
	        }
	      }
	      aij = w[jj] - b;
	      a += B[nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]]*aij/F[jj];
	      v += pow(B[nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]],2)/F[jj];
	    }
	  }
	  
	  e = 0;
	  for(j = 0; j < nnIndxLU[J+i]; j++){
	    e += B[nnIndxLU[i]+j]*w[nnIndx[nnIndxLU[i]+j]];
	  }
	  
	  mu = (kappaOcc[i] / omegaOcc[i] - F77_NAME(ddot)(&pOcc, &X[i], &J, beta, &inc))*omegaOcc[i] + e/F[i] + a;
	  
	  var = 1.0/(omegaOcc[i] + 1.0/F[i] + v);
	  
	  w[i] = rnorm(mu*var, sqrt(var));

        } // i 

        /********************************************************************
         *Update sigmaSq
         *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
         for (j = 0; j < J; j++){
           if(nnIndxLU[J+j] > 0){
             e = 0;
             for(i = 0; i < nnIndxLU[J+j]; i++){
               e += B[nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i]];
             }
             b = w[j] - e;
           }else{
             b = w[j];
           }	
           a += b*b/F[j];
         }

	 theta[sigmaSqIndx] = rigamma(sigmaSqA + J / 2.0, sigmaSqB + 0.5 * a * theta[sigmaSqIndx]); 

        /********************************************************************
         *Update phi (and nu if matern)
         *******************************************************************/
        // Current
        if (corName == "matern"){ nu = theta[nuIndx]; }
        updateBF1RE(B, F, c, C, coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx], theta[phiIndx], nu, covModel, bk, nuB);
        
        a = 0;
        logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
        for (j = 0; j < J; j++){
          if (nnIndxLU[J+j] > 0){
            e = 0;
            for (i = 0; i < nnIndxLU[J+j]; i++){
              e += B[nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i]];
            }
            b = w[j] - e;
          } else{
            b = w[j];
          }	
          a += b*b/F[j];
          logDet += log(F[j]);
        }
      
        logPostCurrent = -0.5*logDet - 0.5*a;
        logPostCurrent += log(theta[phiIndx] - phiA) + log(phiB - theta[phiIndx]); 
        if(corName == "matern"){
        	logPostCurrent += log(theta[nuIndx] - nuA) + log(nuB - theta[nuIndx]); 
        }
        
        // Candidate
        phiCand = logitInv(rnorm(logit(theta[phiIndx], phiA, phiB), exp(tuning[phiIndx])), phiA, phiB);
        if (corName == "matern"){
      	  nuCand = logitInv(rnorm(logit(theta[nuIndx], nuA, nuB), exp(tuning[nuIndx])), nuA, nuB);
        }
      
        updateBF1RE(BCand, FCand, c, C, coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx], phiCand, nuCand, covModel, bk, nuB);
      
        a = 0;
        logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
        for (j = 0; j < J; j++){
          if (nnIndxLU[J+j] > 0){
            e = 0;
            for (i = 0; i < nnIndxLU[J+j]; i++){
              e += BCand[nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i]];
            }
            b = w[j] - e;
          } else{
            b = w[j];
            }	
            a += b*b/FCand[j];
            logDet += log(FCand[j]);
        }
        
        logPostCand = -0.5*logDet - 0.5*a;      
        logPostCand += log(phiCand - phiA) + log(phiB - phiCand); 
        if (corName == "matern"){
          logPostCand += log(nuCand - nuA) + log(nuB - nuCand); 
        }

        if (runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)) {

          std::swap(BCand, B);
          std::swap(FCand, F);
          
          theta[phiIndx] = phiCand;
          accept[phiIndx]++;
          if(corName == "matern"){
            theta[nuIndx] = nuCand; 
            accept[nuIndx]++; 
          }
        }

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
        if (nObs == J) {
          for (i = 0; i < nObs; i++) {
            detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
            psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + w[zLongIndx[i]], zero, one); 
            piProd[zLongIndx[i]] *= pow(1.0 - detProb[i], K[i]);
            ySum[zLongIndx[i]] = y[i]; 	
          } // i
        } else {
          for (i = 0; i < nObs; i++) {
            detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + alphaStarObs[i], zero, one);
            if (tmp_J[zLongIndx[i]] == 0) {
              psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc) + w[zLongIndx[i]], zero, one); 
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
	if (q >= nBurn) {
          thinIndx++; 
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
            F77_NAME(dcopy)(&J, w, &inc, &REAL(wSamples_r)[sPost*J], &inc); 
	    F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[sPost*nTheta], &inc); 
	    F77_NAME(dcopy)(&J, z, &inc, &REAL(zSamples_r)[sPost*J], &inc); 
            F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
            F77_NAME(dcopy)(&nDetRE, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetRE], &inc);
	    // Replicate data set for GoF
	    sPost++; 
	    thinIndx = 0; 
	  }
	}

        R_CheckUserInterrupt();
      } // r (end batch)


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
    int nResultListObjs = 11;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 5, wSamples_r); 
    SET_VECTOR_ELT(result_r, 6, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 7, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 8, sigmaSqPSamples_r);
    SET_VECTOR_ELT(result_r, 9, alphaStarSamples_r);

    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("sigma.sq.p.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, mkChar("alpha.star.samples")); 
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

