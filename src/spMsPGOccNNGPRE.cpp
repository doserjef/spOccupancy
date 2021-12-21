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

void updateBF1MsRE(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){
    
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
  SEXP spMsPGOccNNGPRE(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, 
		       SEXP XpRE_r, SEXP lambdaP_r, SEXP pocc_r, SEXP pdet_r, 
		       SEXP pDetRE_r, SEXP J_r, SEXP nObs_r, SEXP K_r, SEXP N_r, 
		       SEXP nDetRE_r, SEXP nDetRELong_r, SEXP m_r, SEXP nnIndx_r, 
		       SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
		       SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		       SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		       SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		       SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
		       SEXP nuStarting_r, SEXP sigmaSqPStarting_r, SEXP alphaStarStarting_r, 
		       SEXP zLongIndx_r, SEXP alphaStarIndx_r, SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	               SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	               SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP phiA_r, 
		       SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
		       SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		       SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		       SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		       SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, SEXP currChain_r, 
		       SEXP nChain_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, l, s, g, t, q, r, info, nProtect=0, ii;    

    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    // Sorted by visit, then by site, then by species. 
    // (e.g., visit 1, site 1, sp 1, v1, s1, sp2, 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *coords = REAL(coords_r); 
    int *XpRE = INTEGER(XpRE_r); 
    double *lambdaP = REAL(lambdaP_r); 
    int m = INTEGER(m_r)[0]; 
    // Xp is sorted by parameter then site, then visit (parameter 1 site 1, p1 site 2, etc) 
    double *Xp = REAL(Xp_r);
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *muAlphaComm = REAL(muAlphaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *SigmaAlphaCommInv = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlphaComm_r), &inc, SigmaAlphaCommInv, &inc);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *tauSqAlphaA = REAL(tauSqAlphaA_r); 
    double *tauSqAlphaB = REAL(tauSqAlphaB_r); 
    double *phiA = REAL(phiA_r); 
    double *phiB = REAL(phiB_r); 
    double *nuA = REAL(nuA_r); 
    double *nuB = REAL(nuB_r); 
    double *sigmaSqA = REAL(sigmaSqA_r); 
    double *sigmaSqB = REAL(sigmaSqB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double *tuning = REAL(tuning_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int pDetRE = INTEGER(pDetRE_r)[0]; 
    int nDetRE = INTEGER(nDetRE_r)[0]; 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int J = INTEGER(J_r)[0];
    double *K = REAL(K_r); 
    int nObs = INTEGER(nObs_r)[0];
    int N = INTEGER(N_r)[0]; 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
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
        Rprintf("NNGP Multispecies Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
        Rprintf("Using %i nearest neighbors.\n\n", m);
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
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
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pOccN = pOcc * N; 
    int pDetN = pDet * N; 
    int nObsN = nObs * N; 
    int nDetREN = nDetRE * N; 
    int JN = J * N;
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    int JJ = J * J;
    int jj, kk;
    double tmp_0; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_beta = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_alpha = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = 0; 
    }
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Community level
    double *betaComm = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dcopy)(&pOcc, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    double *alphaComm = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaCommStarting_r), &inc, alphaComm, &inc);
    double *tauSqAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(tauSqAlphaStarting_r), &inc, tauSqAlpha, &inc);
    // Species level
    double *beta = (double *) R_alloc(pOccN, sizeof(double));   
    F77_NAME(dcopy)(&pOccN, REAL(betaStarting_r), &inc, beta, &inc);
    // Detection covariates
    double *alpha = (double *) R_alloc(pDetN, sizeof(double));   
    F77_NAME(dcopy)(&pDetN, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREN, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREN, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Spatial random effects
    double *w = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(wStarting_r), &inc, w, &inc);
    // Spatial variance
    double *sigmaSq = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, REAL(sigmaSqStarting_r), &inc, sigmaSq, &inc); 
    // Spatial range parameter
    double *phi = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, REAL(phiStarting_r), &inc, phi, &inc); 
    // Spatial smoothing parameter for Matern
    double *nu = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, REAL(nuStarting_r), &inc, nu, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
    // Only need to set aside J locations since you don't save these 
    // for each species
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double));
    double *omegaOcc = (double *) R_alloc(J, sizeof(double));
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); 
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); 

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++; 
    SEXP tauSqAlphaSamples_r; 
    PROTECT(tauSqAlphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++; 
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOccN, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDetN, nPost)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    SEXP sigmaSqPSamples_r; 
    PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
    SEXP alphaStarSamples_r; 
    PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetREN, nPost)); nProtect++;
    // Spatial random effects
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 

    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
   
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObsN, sizeof(double)); 
    double *psi = (double *) R_alloc(JN, sizeof(double)); 
    zeros(psi, JN); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J);

    // For normal priors
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    // dgemv computes linear combinations of different variables. 
    F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, SigmaBetaCommInvMuBeta, &inc); 	  
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, SigmaAlphaCommInvMuAlpha, &inc); 
    // Put community level variances in a pOcc x POcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level variances in a pDet x pDet matrix. 
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauSqAlpha[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObsN, sizeof(double)); 
    zeros(alphaStarObs, nObsN); 
    // Get sums of the current REs for each site/visit combo for all species
    for (i = 0; i < N; i++) {
      for (r = 0; r < nObs; r++) {
        for (l = 0; l < pDetRE; l++) {
          alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + XpRE[l * nObs + r]];
        }
      }
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int)); 
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE); 
    }

    /**********************************************************************
     Set up spatial stuff
     * *******************************************************************/
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }
    int nThetaN = nTheta * N; 
    double *theta = (double *) R_alloc(nThetaN, sizeof(double));
    for (i = 0; i < N; i++) {
      theta[sigmaSqIndx * N + i] = sigmaSq[i]; 
      theta[phiIndx * N + i] = phi[i]; 
      if (corName == "matern") {
        theta[nuIndx * N + i] = nu[i]; 
      } 
    } // i
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nThetaN, nPost)); nProtect++; 
    // For NNGP
    double a, v, b, e, mu, var, aij; 

    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP. Create a copy of these for each species. Increases storage 
    // space that is needed, but reduces amount of computations. 
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * N, sizeof(double));
    double *F = (double *) R_alloc(J * N, sizeof(double));
    // Only need one of these. 
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads*N, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads*N, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(N*sizeBK, sizeof(double));

    // Initiate B and F for each species
    for (i = 0; i < N; i++) {
    updateBF1MsRE(&B[i * nIndx], &F[i*J], &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * N + i], theta[phiIndx * N + i], nu[i], covModel, &bk[i * sizeBK], nuB[0]);
    }

    /**********************************************************************
     Set up stuff for Adaptive MH and other misc
     * *******************************************************************/
    double logPostCurr = 0.0, logPostCand = 0.0; 
    logPostCurr = R_NegInf; 
    double *accept = (double *) R_alloc(nThetaN, sizeof(double)); zeros(accept, nThetaN); 
    double phiCand = 0.0, nuCand = 0.0; 
    double logDet; 
    // MCMC info if desired
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nThetaN, nBatch)); nProtect++; 

    GetRNGstate();

    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {
    // for (s = 0, g = 0; s < 1; s++) {
    //   for (t = 0; t < 10; t++, g++) {

        /********************************************************************
         Update Community level Occupancy Coefficients
         *******************************************************************/
        /********************************
         Compute b.beta.comm
         *******************************/
        zeros(tmp_pOcc, pOcc); 
        for (i = 0; i < N; i++) {
          F77_NAME(dgemv)(ytran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, &beta[i], &N, &one, tmp_pOcc, &inc); 
        } // i
        for (q = 0; q < pOcc; q++) {
          tmp_pOcc[q] += SigmaBetaCommInvMuBeta[q];  
        } // j

        /********************************
         Compute A.beta.comm
         *******************************/
        for (q = 0; q < ppOcc; q++) {
          tmp_ppOcc[q] = SigmaBetaCommInv[q] + N * TauBetaInv[q]; 
        }
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotri ABetaComm failed\n");}
        // A.beta.inv %*% b.beta
        // 1 * tmp_ppOcc * tmp_pOcc + 0 * tmp_pOcc2  = tmp_pOcc2
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
        // Computes cholesky of tmp_pp again stored back in tmp_ppOcc. This chol(A.beta.inv)
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
        // Args: destination, mu, cholesky of the inverse covariance matrix, dimension
        mvrnorm(betaComm, tmp_pOcc2, tmp_ppOcc, pOcc);
        /********************************************************************
         Update Community level Detection Coefficients
         *******************************************************************/
        /********************************
         * Compute b.alpha.comm
         *******************************/
         zeros(tmp_pDet, pDet); 
         for (i = 0; i < N; i++) {
           F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, TauAlphaInv, &pDet, &alpha[i], &N, &one, tmp_pDet, &inc); 
         } // i
         for (q = 0; q < pDet; q++) {
           tmp_pDet[q] += SigmaAlphaCommInvMuAlpha[q];  
         } // j
        /********************************
         * Compute A.alpha.comm
         *******************************/
        for (q = 0; q < ppDet; q++) {
          tmp_ppDet[q] = SigmaAlphaCommInv[q] + N * TauAlphaInv[q]; 
        }
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info); 
        if(info != 0){error("c++ error: dpotri AAlphaComm failed\n");}
        // A.alpha.inv %*% b.alpha
        // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2  = tmp_pDet2
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc);
        // Computes cholesky of tmp_pp again stored back in tmp_ppDet. This chol(A.alpha.inv)
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        // Args: destination, mu, cholesky of the inverse covariance matrix, dimension
        mvrnorm(alphaComm, tmp_pDet2, tmp_ppDet, pDet);

        /********************************************************************
         Update Community Occupancy Variance Parameter
        ********************************************************************/
        for (q = 0; q < pOcc; q++) {
          tmp_0 = 0.0;  
          for (i = 0; i < N; i++) {
            tmp_0 += (beta[q * N + i] - betaComm[q]) * (beta[q * N + i] - betaComm[q]);
          } // i
          tmp_0 *= 0.5;
          tauSqBeta[q] = rigamma(tauSqBetaA[q] + N / 2.0, tauSqBetaB[q] + tmp_0); 
        } // q
        // This is correct, nothing wrong here. 
        for (q = 0; q < pOcc; q++) {
          TauBetaInv[q * pOcc + q] = tauSqBeta[q]; 
        } // i
        F77_NAME(dpotrf)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, TauBetaInv, &pOcc, &info); 
        if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
        /********************************************************************
         Update Community Detection Variance Parameter
        ********************************************************************/
        for (q = 0; q < pDet; q++) {
          tmp_0 = 0.0;  
          for (i = 0; i < N; i++) {
            tmp_0 += (alpha[q * N + i] - alphaComm[q]) * (alpha[q * N + i] - alphaComm[q]);
          } // i
          tmp_0 *= 0.5;
          tauSqAlpha[q] = rigamma(tauSqAlphaA[q] + N / 2.0, tauSqAlphaB[q] + tmp_0); 
        } // q
        for (q = 0; q < pDet; q++) {
          TauAlphaInv[q * pDet + q] = tauSqAlpha[q]; 
        } // i
        F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info); 
        if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
        F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info); 
        if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

        /********************************************************************
         *Update Detection random effects variance
         *******************************************************************/
        for (l = 0; l < pDetRE; l++) {
          tmp_0 = 0.0; 
          for (i = 0; i < N; i++) {
            tmp_0 += F77_NAME(ddot)(&nDetRELong[l], &alphaStar[i*nDetRE + alphaStarStart[l]], &inc, &alphaStar[i*nDetRE + alphaStarStart[l]], &inc); 
          }
          tmp_0 *= 0.5; 
          sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] * N / 2.0, sigmaSqPB[l] + tmp_0);
        }
         
        for (i = 0; i < N; i++) {  
          /********************************************************************
           *Update Occupancy Auxiliary Variables 
           *******************************************************************/
          for (j = 0; j < J; j++) {
            omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, &beta[i], &N) + w[j * N + i]);
          } // j
          /********************************************************************
           *Update Detection Auxiliary Variables 
           *******************************************************************/
          // Note that all of the variables are sampled, but only those at 
          // locations with z[j] == 1 actually effect the results. 
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              omegaDet[r] = rpg(K[r], F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              omegaDet[r] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
            } // r
          }
             
          /********************************************************************
           *Update Occupancy Regression Coefficients
           *******************************************************************/
          for (j = 0; j < J; j++) {
            kappaOcc[j] = z[j * N + i] - 1.0 / 2.0; 
	    tmp_J1[j] = kappaOcc[j] - omegaOcc[j] * w[j * N + i]; 
          } // j
          /********************************
           * Compute b.beta
           *******************************/
          // t(X) * tmp_J1 + 0 * tmp_pOcc = tmp_pOcc. 
          // dgemv computes linear combinations of different variables. 
          F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, tmp_J1, &inc, &zero, tmp_pOcc, &inc); 	 
          // TauBetaInv %*% betaComm + tmp_pOcc = tmp_pOcc
          F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, betaComm, &inc, &one, tmp_pOcc, &inc); 

          /********************************
           * Compute A.beta
           * *****************************/
          // t(X) %*% diag(omegaOcc)
          for(j = 0; j < J; j++){
            for(q = 0; q < pOcc; q++){
              tmp_JpOcc[q*J+j] = X[q*J+j]*omegaOcc[j];
            }
          }
          // This finishes off A.beta
          // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
          F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc);
          for (q = 0; q < ppOcc; q++) {
            tmp_ppOcc[q] += TauBetaInv[q]; 
          } // j
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
          if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
          if(info != 0){error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); if(info != 0){error("c++ error: dpotrf here failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_beta, tmp_pOcc2, tmp_ppOcc, pOcc);
          // Can eventually get rid of this and change order of beta. 
          for (q = 0; q < pOcc; q++) {
            beta[q * N + i] = tmp_beta[q]; 
          }
        
          /********************************************************************
           *Update Detection Regression Coefficients
           *******************************************************************/
          /********************************
           * Compute b.alpha
           *******************************/
          // First multiply kappDet * the current occupied values, such that values go 
          // to 0 if they z == 0 and values go to kappaDet if z == 1
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              kappaDet[r] = (y[r * N + i] - K[r]/2.0) * z[zLongIndx[r] * N + i];
              tmp_nObs[r] = kappaDet[r] - omegaDet[r] * alphaStarObs[i * nObs + r]; 
              tmp_nObs[r] *= z[zLongIndx[r] * N + i]; 
            } // r
          } else { 
            for (r = 0; r < nObs; r++) {
              kappaDet[r] = (y[r * N + i] - 1.0/2.0) * z[zLongIndx[r] * N + i];
              tmp_nObs[r] = kappaDet[r] - omegaDet[r] * alphaStarObs[i * nObs + r]; 
              tmp_nObs[r] *= z[zLongIndx[r] * N + i]; 
            } // r
          }
          
          F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc);  
          F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, TauAlphaInv, &pDet, alphaComm, &inc, &one, tmp_pDet, &inc); 

          /********************************
           * Compute A.alpha
           * *****************************/
          for (r = 0; r < nObs; r++) {
            for (q = 0; q < pDet; q++) {
              tmp_nObspDet[q*nObs + r] = Xp[q * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i];
            } // i
          } // j

          // This finishes off A.alpha
          // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
          F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet);

          for (q = 0; q < ppDet; q++) {
            tmp_ppDet[q] += TauAlphaInv[q]; 
          } // q
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
          if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
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
          mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDet);
          for (q = 0; q < pDet; q++) {
            alpha[q * N + i] = tmp_alpha[q];
          }

          /********************************************************************
           *Update w (spatial random effects)
           *******************************************************************/
	  for (ii = 0; ii < J; ii++ ) {
            a = 0;
	    v = 0;
	    if (uIndxLU[J + ii] > 0){ // is ii a neighbor for anybody
	      for (j = 0; j < uIndxLU[J+ii]; j++){ // how many locations have ii as a neighbor
	        b = 0;
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[ii]+j]; // jj is the index of the jth location who has ii as a neighbor
	        for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	          kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	          if(kk != ii){ //if the neighbor of jj is not ii
	    	    b += B[i*nIndx + nnIndxLU[jj]+k]*w[kk * N + i]; //covariance between jj and kk and the random effect of kk
	          }
	        }
	        aij = w[jj * N + i] - b;
	        a += B[i*nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]]*aij/F[i*J + jj];
	        v += pow(B[i * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]],2)/F[i * J + jj];
	      }
	    }
	    
	    e = 0;
	    for(j = 0; j < nnIndxLU[J+ii]; j++){
	      e += B[i * nIndx + nnIndxLU[ii]+j]*w[nnIndx[nnIndxLU[ii]+j] * N + i];
	    }
	    
	    mu = (kappaOcc[ii] / omegaOcc[ii] - F77_NAME(ddot)(&pOcc, &X[ii], &J, &beta[i], &N))*omegaOcc[ii] + e/F[i*J + ii] + a;
	    
	    var = 1.0/(omegaOcc[ii] + 1.0/F[i * J + ii] + v);
	    
	    w[ii * N + i] = rnorm(mu*var, sqrt(var));

          } // ii

          /********************************************************************
           *Update sigmaSq
           *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:a, logDet)
#endif

          for (j = 0; j < J; j++){
            if(nnIndxLU[J+j] > 0){
              e = 0;
              for(ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += B[i * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * N + i];
              }
              b = w[j * N + i] - e;
            }else{
              b = w[j * N + i];
            }	
            a += b*b/F[i * J + j];
          }

	  theta[sigmaSqIndx * N + i] = rigamma(sigmaSqA[i] + J / 2.0, sigmaSqB[i] + 0.5 * a * theta[sigmaSqIndx * N + i]); 

          /********************************************************************
           *Update phi (and nu if matern)
           *******************************************************************/
          // Current
          if (corName == "matern"){ 
	    nu[i] = theta[nuIndx * N + i];
       	  }
          updateBF1MsRE(&B[i * nIndx], &F[i*J], &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * N + i], theta[phiIndx * N + i], nu[i], covModel, &bk[i * sizeBK], nuB[i]);
          a = 0;
          logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:a, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += B[i * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * N + i];
              }
              b = w[j * N + i] - e;
            } else{
              b = w[j * N + i];
            }	
            a += b*b/F[i * J + j];
            logDet += log(F[i * J + j]);
          }
      
          logPostCurr = -0.5 * logDet - 0.5 * a;
          logPostCurr += log(theta[phiIndx * N + i] - phiA[i]) + log(phiB[i] - theta[phiIndx * N + i]); 
          if(corName == "matern"){
       	    logPostCurr += log(theta[nuIndx * N + i] - nuA[i]) + log(nuB[i] - theta[nuIndx * N + i]); 
          }
          
          // Candidate
          phiCand = logitInv(rnorm(logit(theta[phiIndx * N + i], phiA[i], phiB[i]), exp(tuning[phiIndx * N + i])), phiA[i], phiB[i]);
          if (corName == "matern"){
      	    nuCand = logitInv(rnorm(logit(theta[nuIndx * N + i], nuA[i], nuB[i]), exp(tuning[nuIndx * N + i])), nuA[i], nuB[i]);
          }
      
          updateBF1MsRE(BCand, FCand, &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * N + i], phiCand, nuCand, covModel, &bk[i * sizeBK], nuB[i]);
      
          a = 0;
          logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:a, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += BCand[nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * N + i];
              }
              b = w[j * N + i] - e;
            } else{
              b = w[j * N + i];
              }	
              a += b*b/FCand[j];
              logDet += log(FCand[j]);
          }
          
          logPostCand = -0.5*logDet - 0.5*a;      
          logPostCand += log(phiCand - phiA[i]) + log(phiB[i] - phiCand); 
          if (corName == "matern"){
            logPostCand += log(nuCand - nuA[i]) + log(nuB[i] - nuCand); 
          }

          if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {

            F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[i * nIndx], &inc);
            F77_NAME(dcopy)(&J, FCand, &inc, &F[i * J], &inc);
            
	    theta[phiIndx * N + i] = phiCand;
            accept[phiIndx * N + i]++;
            if (corName == "matern") {
              nu[i] = nuCand; 
	      theta[nuIndx * N + i] = nu[i]; 
              accept[nuIndx * N + i]++; 
            }
          }

          /********************************************************************
           *Update Detection random effects
           *******************************************************************/
          // Update each individual random effect one by one. 
          for (l = 0; l < nDetRE; l++) {
            /********************************
             * Compute b.alpha.star
             *******************************/
            for (r = 0; r < nObs; r++) {
              tmp_nObs[r] = kappaDet[r] - (F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r] - alphaStar[i * nDetRE + l]) * omegaDet[r];
              // Only allow information to come from when z == 1.
              tmp_nObs[r] *= z[zLongIndx[r] * N + i]; 
            }
            F77_NAME(dgemv)(ytran, &nObs, &inc, &one, &lambdaP[l * nObs], &nObs, tmp_nObs, &inc, &zero, tmp_one, &inc); 
            /********************************
             * Compute A.alpha.star
             *******************************/
            for (r = 0; r < nObs; r++) {
              tmp_nObs[r] = lambdaP[l * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i]; 
            }
            tmp_0 = F77_NAME(ddot)(&nObs, tmp_nObs, &inc, &lambdaP[l * nObs], &inc); 
            tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
            tmp_0 = 1.0 / tmp_0; 
            alphaStar[i * nDetRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
          }
	  zeros(&alphaStarObs[i * nObs], nObs); 
          // Update the RE sums for the current species
          for (r = 0; r < nObs; r++) {
            for (l = 0; l < pDetRE; l++) {
              alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + XpRE[l * nObs + r]]; 
            }
          }

          /********************************************************************
           *Update Latent Occupancy
           *******************************************************************/
          // Compute detection probability 
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
              psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + w[zLongIndx[r] * N + i], zero, one); 
              piProd[zLongIndx[r]] = pow(1.0 - detProb[i * nObs + r], K[r]);
              ySum[zLongIndx[r]] = y[r * N + i]; 
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
              if (tmp_J[zLongIndx[r]] == 0) {
                psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[r]], &J, &beta[i], &N) + w[zLongIndx[r] * N + i], zero, one); 
              }
              piProd[zLongIndx[r]] *= (1.0 - detProb[i * nObs + r]);
              ySum[zLongIndx[r]] += y[r * N + i]; 	
              tmp_J[zLongIndx[r]]++;
            } // r
          }
          // Compute occupancy probability 
          for (j = 0; j < J; j++) {
            psiNum = psi[j * N + i] * piProd[j]; 
            if (ySum[j] == zero) {
              z[j * N + i] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j * N + i])));           
            } else {
              z[j * N + i] = one; 
            }
            piProd[j] = one;
            ySum[j] = zero; 
            tmp_J[j] = 0; 
          } // j
        } // i

        /********************************************************************
         *Save samples
         *******************************************************************/
	if (g >= nBurn) {
	  thinIndx++;
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, alphaComm, &inc, &REAL(alphaCommSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&pOcc, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, tauSqAlpha, &inc, &REAL(tauSqAlphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&pOccN, beta, &inc, &REAL(betaSamples_r)[sPost*pOccN], &inc); 
            F77_NAME(dcopy)(&pDetN, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDetN], &inc); 
            F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&JN, w, &inc, &REAL(wSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&nThetaN, theta, &inc, &REAL(thetaSamples_r)[sPost*nThetaN], &inc); 
            F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
            F77_NAME(dcopy)(&nDetREN, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetREN], &inc);
	    sPost++; 
	    thinIndx = 0; 
	  }
	}
        R_CheckUserInterrupt();
      } // t (end batch)

      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (i = 0; i < N; i++) {
        for (k = 0; k < nTheta; k++) {
          REAL(acceptSamples_r)[s * nThetaN + k * N + i] = accept[k * N + i]/batchLength; 
          REAL(tuningSamples_r)[s * nThetaN + k * N + i] = tuning[k * N + i]; 
	  // Rprintf("accept[k * N + i]: %f\n", accept[k * N + i]); 
          if (accept[k * N + i] / batchLength > acceptRate) {
            tuning[k * N + i] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          } else{
              tuning[k * N + i] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
            }
          accept[k * N + i] = 0.0;
        } // k
      } // i

      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
	  Rprintf("\tSpecies\t\tAcceptance\tTuning\n");	  
	  for (i = 0; i < N; i++) {
	    Rprintf("\t%i\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nThetaN + phiIndx * N + i], exp(tuning[phiIndx * N + i]));
	  } // i
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      } 
      status++;        

    } // all batches
    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }
  
    // This is necessary when generating random numbers in C.     
    PutRNGstate();

    // make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 14;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauSqAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, zSamples_r);
    SET_VECTOR_ELT(result_r, 7, psiSamples_r);
    SET_VECTOR_ELT(result_r, 8, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 9, wSamples_r); 
    SET_VECTOR_ELT(result_r, 10, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 11, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 12, sigmaSqPSamples_r);
    SET_VECTOR_ELT(result_r, 13, alphaStarSamples_r);

    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 10, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 11, mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 12, mkChar("sigma.sq.p.samples")); 
    SET_VECTOR_ELT(resultName_r, 13, mkChar("alpha.star.samples")); 
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


