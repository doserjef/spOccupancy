#define USE_FC_LEN_T
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
#ifndef FCONE
# define FCONE
#endif

void updateBFsvcMs(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){
    
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
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
  SEXP svcMsPGOccNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r, SEXP XpRE_r, 
		      SEXP rangeInd_r, SEXP consts_r, SEXP K_r, SEXP nOccRELong_r, SEXP nDetRELong_r, 
		      SEXP m_r, SEXP nnIndx_r, 
		      SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
		      SEXP betaStarting_r, SEXP alphaStarting_r, SEXP zStarting_r, 
		      SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
		      SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
		      SEXP phiStarting_r, SEXP lambdaStarting_r, SEXP wStarting_r, 
		      SEXP nuStarting_r, 
		      SEXP sigmaSqPsiStarting_r, SEXP sigmaSqPStarting_r, 
		      SEXP betaStarStarting_r, SEXP alphaStarStarting_r, SEXP zLongIndx_r,
		      SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		      SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
		      SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	              SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP tauSqBetaA_r, 
	              SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP phiA_r, 
		      SEXP phiB_r, SEXP nuA_r, SEXP nuB_r, 
		      SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		      SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
		      SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, 
		      SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
		      SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    // Indexing key (this is general, might deviate in a couple spots): 
    // i = species 
    // j = site
    // r = observation
    // ll = latent factor
    // h = coefficients (occurrence or detection)
    // rr = spatially-varying coefficient. 
    int i, j, k, s, g, t, h, r, l, info, nProtect=0, ii, ll, rr, rrr;    

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
    // Order: covariate, site
    double *Xw = REAL(Xw_r);
    double *coords = REAL(coords_r); 
    double *rangeInd = REAL(rangeInd_r);
    int *XpRE = INTEGER(XpRE_r); 
    int *XRE = INTEGER(XRE_r);
    int m = INTEGER(m_r)[0]; 
    // Xp is sorted by parameter, then by visit, then by site 
    double *Xp = REAL(Xp_r);
    // Load constants
    int N = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int nObs = INTEGER(consts_r)[2]; 
    int pOcc = INTEGER(consts_r)[3];
    int pOccRE = INTEGER(consts_r)[4];
    int nOccRE = INTEGER(consts_r)[5];
    int pDet = INTEGER(consts_r)[6];
    int pDetRE = INTEGER(consts_r)[7];
    int nDetRE = INTEGER(consts_r)[8];
    int q = INTEGER(consts_r)[9]; 
    int pTilde = INTEGER(consts_r)[10];
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
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *tuning = REAL(tuning_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    double *K = REAL(K_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
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
        Rprintf("Spatial Factor NNGP Multispecies Occupancy Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, N);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
	Rprintf("Number of spatially-varying coefficients: %i \n", pTilde);
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
        Rprintf("Using %i latent spatial factors.\n", q);
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
    int nOccREN = nOccRE * N; 
    int Jq = J * q;
    int qq = q * q;
    int JN = J * N;
    int Nq = N * q;
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    int JpOccRE = J * pOccRE; 
    int nObspDetRE = nObs * pDetRE;
    int JpTilde = J * pTilde;
    int qpTilde = q * pTilde;
    int JqpTilde = J * q * pTilde;
    int JNpTilde = J * N * pTilde;
    int NqpTilde = N * q * pTilde;
    int jj, kk;
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_beta = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_alpha = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    int *tmp_JInt = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_JInt[j] = 0; 
    }
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J, J);
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_qq = (double *) R_alloc(qq, sizeof(double));
    double *tmp_q = (double *) R_alloc(q, sizeof(double));
    double *tmp_q2 = (double *) R_alloc(q, sizeof(double));
    double *tmp_qq2 = (double *) R_alloc(qq, sizeof(double));
    double *tmp_Jq = (double *) R_alloc(Jq, sizeof(double));
    double *tmp_Jq2 = (double *) R_alloc(Jq, sizeof(double));
    double *tmp_Nq = (double *) R_alloc(Nq, sizeof(double));
    double *tmp_Nq2 = (double *) R_alloc(Nq, sizeof(double));
    double *tmp_N = (double *) R_alloc(N, sizeof(double));
    double *tmp_N2 = (double *) R_alloc(N, sizeof(double));
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    int currDim = 0;

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
    // Occurrence random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Detection covariates
    double *alpha = (double *) R_alloc(pDetN, sizeof(double));   
    F77_NAME(dcopy)(&pDetN, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Spatial random effects
    double *w = (double *) R_alloc(JqpTilde, sizeof(double));
    F77_NAME(dcopy)(&JqpTilde, REAL(wStarting_r), &inc, w, &inc);
    // Latent spatial factor loadings 
    double *lambda = (double *) R_alloc(NqpTilde, sizeof(double));
    F77_NAME(dcopy)(&NqpTilde, REAL(lambdaStarting_r), &inc, lambda, &inc);
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREN, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREN, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Latent occurrence random effects
    double *betaStar = (double *) R_alloc(nOccREN, sizeof(double)); 
    F77_NAME(dcopy)(&nOccREN, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Spatial range parameter
    double *phi = (double *) R_alloc(qpTilde, sizeof(double)); 
    F77_NAME(dcopy)(&qpTilde, REAL(phiStarting_r), &inc, phi, &inc); 
    // Spatial smoothing parameter for Matern
    double *nu = (double *) R_alloc(qpTilde, sizeof(double)); 
    F77_NAME(dcopy)(&qpTilde, REAL(nuStarting_r), &inc, nu, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(JN, sizeof(double));   
    F77_NAME(dcopy)(&JN, REAL(zStarting_r), &inc, z, &inc);
    // Auxiliary variables
    // Note, you can just write over the values for the detection
    // parameters. 
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double));zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(JN, sizeof(double)); zeros(omegaOcc, JN);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(JN, sizeof(double)); zeros(kappaOcc, JN);
    // Need this for all species
    double *zStar = (double *) R_alloc(JN, sizeof(double));

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), pOcc * nPost);
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaCommSamples_r), pDet * nPost);
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++; 
    zeros(REAL(tauSqBetaSamples_r), pOcc * nPost);
    SEXP tauSqAlphaSamples_r; 
    PROTECT(tauSqAlphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++; 
    zeros(REAL(tauSqAlphaSamples_r), pDet * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOccN, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOccN * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDetN, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDetN * nPost);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), JN * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), JN * nPost);
    // Spatial parameters
    SEXP lambdaSamples_r; 
    PROTECT(lambdaSamples_r = allocMatrix(REALSXP, NqpTilde, nPost)); nProtect++;
    zeros(REAL(lambdaSamples_r), NqpTilde * nPost);
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, JqpTilde, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), JqpTilde * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetREN, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetREN * nPost);
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPsiSamples_r), pOccRE * nPost);
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nOccREN, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nOccREN * nPost);
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), JN * nPost);
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObsN, sizeof(double)); 
    double *yWAIC = (double *) R_alloc(JN, sizeof(double)); ones(yWAIC, JN);
    double *psi = (double *) R_alloc(JN, sizeof(double)); 
    zeros(psi, JN); 
    double *piProd = (double *) R_alloc(JN, sizeof(double)); 
    ones(piProd, JN); 
    double *piProdWAIC = (double *) R_alloc(JN, sizeof(double)); 
    ones(piProdWAIC, JN); 
    double *ySum = (double *) R_alloc(JN, sizeof(double)); zeros(ySum, JN); 

    // For normal community-level priors
    // Occurrence coefficients
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaCommInv, &pOcc, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaCommInv, &pOcc, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, 
                   SigmaAlphaCommInvMuAlpha, &inc FCONE);
    // Put community level occurrence variances in a pOcc x pOcc matrix.
    double *TauBetaInv = (double *) R_alloc(ppOcc, sizeof(double)); zeros(TauBetaInv, ppOcc); 
    for (i = 0; i < pOcc; i++) {
      TauBetaInv[i * pOcc + i] = 1.0 / tauSqBeta[i]; 
    } // i
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = 1.0 / tauSqAlpha[i]; 
    } // i

    /**********************************************************************
     * Prep for random effects (if they exist)
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JN, sizeof(double)); 
    zeros(betaStarSites, JN); 
    int *betaStarLongIndx = (int *) R_alloc(JpOccRE, sizeof(int));
    // Initial sums (initiate with the first species)
    for (j = 0; j < J; j++) {
      for (l = 0; l < pOccRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nOccRE);
        for (i = 0; i < N; i++) {
          betaStarSites[i * J + j] += betaStar[i * nOccRE + betaStarLongIndx[l * J + j]];
        }
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObsN, sizeof(double)); 
    zeros(alphaStarObs, nObsN); 
    int *alphaStarLongIndx = (int *) R_alloc(nObspDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo for all species
    for (r = 0; r < nObs; r++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarLongIndx[l * nObs + r] = which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE);
        for (i = 0; i < N; i++) {
          alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + alphaStarLongIndx[l * nObs + r]];
        }
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

    /**********************************************************************
     Set up spatial stuff
     * *******************************************************************/
    // Note that even though sigmaSq is fixed at 1 for spatial factor models, 
    // still maintaining the same approach to leverage the same correlation
    // functions underneath. 
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }
    int nThetaq = nTheta * q; 
    int nThetaqpTilde = nThetaq * pTilde;
    int nThetaqpTildeSave = (nTheta - 1) * qpTilde;
    double *theta = (double *) R_alloc(nThetaqpTilde, sizeof(double));
    for (rr = 0; rr < pTilde; rr++) {
      for (ll = 0; ll < q; ll++) {
        theta[phiIndx * qpTilde + rr * q + ll] = phi[rr * q + ll];
        // sigmaSq by default is 1 for spatial factor models. 
        theta[sigmaSqIndx * qpTilde + rr * q + ll] = 1.0;
        if (corName == "matern") {
          theta[nuIndx * qpTilde + rr * q + ll] = nu[rr * q + ll]; 
        } 
      } // ll
    } // rr 
    SEXP thetaSamples_r; 
    // Note you don't save any of the sigma.sqs which are fixed at 1. 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nThetaqpTildeSave, nPost)); nProtect++; 
    zeros(REAL(thetaSamples_r), nThetaqpTildeSave * nPost);
    // Species-level spatial random effects for each SVC. 
    double *wStar = (double *) R_alloc(JNpTilde, sizeof(double)); zeros(wStar, JNpTilde);
    // Multiply Lambda %*% w[j] to get wStar for each SVC. 
    for (rr = 0; rr < pTilde; rr++) {
      for (j = 0; j < J; j++) {
        F77_NAME(dgemv)(ntran, &N, &q, &one, &lambda[rr * Nq], &N, &w[rr * Jq + j * q], 
			&inc, &zero, &wStar[rr * JN + j * N], &inc FCONE);
      } // j
    } // rr
    // For NNGP.
    double b, e, aij, aa; 
    double *a = (double *) R_alloc(q, sizeof(double));
    double *v = (double *) R_alloc(q, sizeof(double));
    double *mu = (double *) R_alloc(q, sizeof(double));
    double *var = (double *) R_alloc(qq, sizeof(double));
    double *ff = (double *) R_alloc(q, sizeof(double));
    double *gg = (double *) R_alloc(q, sizeof(double));

    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP. Create a copy of these for each species. Increases storage 
    // space that is needed, but reduces amount of computations. 
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * q, sizeof(double)); 
    double *F = (double *) R_alloc(J * q, sizeof(double));
    // Only need one of these. 
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads*q, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads*q, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(q*sizeBK, sizeof(double));

    // wStar is JN x pTilde
    // Spatial process sums for each site and species
    double *wSites = (double *) R_alloc(JN, sizeof(double));
    // For each species and each location, multiply w x Xw
    for (i = 0; i < N; i++) {
      for (j = 0; j < J; j++) {
        wSites[j * N + i] = 0.0;
        for (rr = 0; rr < pTilde; rr++) {
          wSites[j * N + i] += wStar[rr * JN + j * N + i] * Xw[i * JpTilde + rr * J + j];
        } // rr
      } // j
    } // i

    /**********************************************************************
     Set up stuff for Adaptive MH and other misc
     * *******************************************************************/
    double logPostCurr = 0.0, logPostCand = 0.0; 
    logPostCurr = R_NegInf; 
    double *accept = (double *) R_alloc(nThetaqpTilde, sizeof(double)); zeros(accept, nThetaqpTilde); 
    double phiCand = 0.0, nuCand = 0.0; 
    double logDet; 
    // MCMC info if desired
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nThetaqpTilde, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nThetaqpTilde * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nThetaqpTilde, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nThetaqpTilde * nBatch);

    GetRNGstate();

    /**********************************************************************
     Start sampling
     * *******************************************************************/
    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {

        /********************************************************************
         Update Community level Occupancy Coefficients
         *******************************************************************/
        /********************************
         Compute b.beta.comm
         *******************************/
        zeros(tmp_pOcc, pOcc); 
        for (i = 0; i < N; i++) {
          F77_NAME(dgemv)(ytran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, &beta[i], &N, &one, tmp_pOcc, &inc FCONE); 
        } // i
        for (h = 0; h < pOcc; h++) {
          tmp_pOcc[h] += SigmaBetaCommInvMuBeta[h];  
        } // j

        /********************************
         Compute A.beta.comm
         *******************************/
        for (h = 0; h < ppOcc; h++) {
          tmp_ppOcc[h] = SigmaBetaCommInv[h] + N * TauBetaInv[h]; 
        }
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){error("c++ error: dpotri ABetaComm failed\n");}
        // A.beta.inv %*% b.beta
        // 1 * tmp_ppOcc * tmp_pOcc + 0 * tmp_pOcc2  = tmp_pOcc2
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        // Computes cholesky of tmp_pp again stored back in tmp_ppOcc. This chol(A.beta.inv)
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
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
           F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, TauAlphaInv, &pDet, &alpha[i], &N, &one, tmp_pDet, &inc FCONE); 
         } // i
         for (h = 0; h < pDet; h++) {
           tmp_pDet[h] += SigmaAlphaCommInvMuAlpha[h];  
         } // j
        /********************************
         * Compute A.alpha.comm
         *******************************/
        for (h = 0; h < ppDet; h++) {
          tmp_ppDet[h] = SigmaAlphaCommInv[h] + N * TauAlphaInv[h]; 
        }
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotri AAlphaComm failed\n");}
        // A.alpha.inv %*% b.alpha
        // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2  = tmp_pDet2
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
        // Computes cholesky of tmp_pp again stored back in tmp_ppDet. This chol(A.alpha.inv)
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        // Args: destination, mu, cholesky of the inverse covariance matrix, dimension
        mvrnorm(alphaComm, tmp_pDet2, tmp_ppDet, pDet);

        /********************************************************************
         Update Community Occupancy Variance Parameter
        ********************************************************************/
        for (h = 0; h < pOcc; h++) {
          tmp_0 = 0.0;  
          for (i = 0; i < N; i++) {
            tmp_0 += (beta[h * N + i] - betaComm[h]) * (beta[h * N + i] - betaComm[h]);
          } // i
          tmp_0 *= 0.5;
          tauSqBeta[h] = rigamma(tauSqBetaA[h] + N / 2.0, tauSqBetaB[h] + tmp_0); 
        } // h
        for (h = 0; h < pOcc; h++) {
          TauBetaInv[h * pOcc + h] = 1.0 / tauSqBeta[h]; 
        } // i
        /********************************************************************
         Update Community Detection Variance Parameter
        ********************************************************************/
        for (h = 0; h < pDet; h++) {
          tmp_0 = 0.0;  
          for (i = 0; i < N; i++) {
            tmp_0 += (alpha[h * N + i] - alphaComm[h]) * (alpha[h * N + i] - alphaComm[h]);
          } // i
          tmp_0 *= 0.5;
          tauSqAlpha[h] = rigamma(tauSqAlphaA[h] + N / 2.0, tauSqAlphaB[h] + tmp_0); 
        } // h
        for (h = 0; h < pDet; h++) {
          TauAlphaInv[h * pDet + h] = 1.0 / tauSqAlpha[h]; 
        } // i

        /********************************************************************
         *Update Occupancy random effects variance
         *******************************************************************/
        for (l = 0; l < pOccRE; l++) {
          tmp_0 = 0.0; 
          for (i = 0; i < N; i++) {
            tmp_0 += F77_NAME(ddot)(&nOccRELong[l], &betaStar[i*nOccRE + betaStarStart[l]], &inc, &betaStar[i*nOccRE + betaStarStart[l]], &inc); 
          }
          tmp_0 *= 0.5; 
          sigmaSqPsi[l] = rigamma(sigmaSqPsiA[l] + nOccRELong[l] * N / 2.0, sigmaSqPsiB[l] + tmp_0);
        }

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

        /********************************************************************
         *Update Species-Specific Regression Parameters
         *******************************************************************/
        for (i = 0; i < N; i++) {  
          /********************************************************************
           *Update Occupancy Auxiliary Variables 
           *******************************************************************/
          for (j = 0; j < J; j++) {
            // Only calculate if site is within the range. Note this will not change
	    // by iteration.  
            if (rangeInd[j * N + i] == 1.0) {
              omegaOcc[j * N + i] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[i * JpOcc + j], &J, &beta[i], &N) + wSites[j * N + i] + betaStarSites[i * J + j]);
	    }
          } // j
          /********************************************************************
           *Update Detection Auxiliary Variables 
           *******************************************************************/
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              if (z[zLongIndx[r] * N + i] == 1.0) {
                omegaDet[r] = rpg(K[r], F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
	      }
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              if (z[zLongIndx[r] * N + i] == 1.0) {
                omegaDet[r] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r]);
	      }
            } // r
          }

          /********************************************************************
           *Update Occupancy Regression Coefficients
           *******************************************************************/
	  zeros(tmp_J1, J);
          for (j = 0; j < J; j++) {
            if (rangeInd[j * N + i] == 1.0) {
              kappaOcc[j * N + i] = z[j * N + i] - 1.0 / 2.0; 
              tmp_J1[j] = kappaOcc[j * N + i] - omegaOcc[j * N + i] * 
	                  (wSites[j * N + i] + betaStarSites[i * J + j]); 
	      // For later
	      zStar[j * N + i] = kappaOcc[j * N + i] / omegaOcc[j * N + i];
	    }
          } // j
          /********************************
           * Compute b.beta
           *******************************/
          // t(X) * tmp_J1 + 0 * tmp_pOcc = tmp_pOcc. 
          F77_NAME(dgemv)(ytran, &J, &pOcc, &one, &X[i * JpOcc], &J, tmp_J1, 
			  &inc, &zero, tmp_pOcc, &inc FCONE); 	 
          // TauBetaInv %*% betaComm + tmp_pOcc = tmp_pOcc
          F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, TauBetaInv, &pOcc, 
			  betaComm, &inc, &one, tmp_pOcc, &inc FCONE); 

          /********************************
           * Compute A.beta
           * *****************************/
          // t(X) %*% diag(omegaOcc)
	  zeros(tmp_JpOcc, JpOcc);
          for(j = 0; j < J; j++){
            if (rangeInd[j * N + i] == 1.0) {
              for(h = 0; h < pOcc; h++){
                tmp_JpOcc[h*J+j] = X[i * JpOcc + h * J + j] * omegaOcc[j * N + i];
              }
	    }
          }
          // This finishes off A.beta
          // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
          F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, &X[i * JpOcc], &J, 
			  tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
          for (h = 0; h < ppOcc; h++) {
            tmp_ppOcc[h] += TauBetaInv[h]; 
          } // j
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
          if(info != 0){error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, 
			  &inc, &zero, tmp_pOcc2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
	  if(info != 0){error("c++ error: dpotrf A.beta 2 failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_beta, tmp_pOcc2, tmp_ppOcc, pOcc);
          // Can eventually get rid of this and change order of beta. 
          for (h = 0; h < pOcc; h++) {
            beta[h * N + i] = tmp_beta[h]; 
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
          
          F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
          F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, TauAlphaInv, &pDet, alphaComm, &inc, &one, tmp_pDet, &inc FCONE); 
          /********************************
           * Compute A.alpha
           * *****************************/
          for (r = 0; r < nObs; r++) {
            for (h = 0; h < pDet; h++) {
              tmp_nObspDet[h*nObs + r] = Xp[h * nObs + r] * omegaDet[r] * z[zLongIndx[r] * N + i];
            } // h
          } // j

          // This finishes off A.alpha
          // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
          F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

          for (h = 0; h < ppDet; h++) {
            tmp_ppDet[h] += TauAlphaInv[h]; 
          } // h
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
          F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
          // A.alpha.inv %*% b.alpha
          // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2 
          // (which is currently nothing) = tmp_pDet2
          F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
          // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf A.alpha 2 failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDet);
          for (h = 0; h < pDet; h++) {
            alpha[h * N + i] = tmp_alpha[h];
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
                if (rangeInd[j * N + i] == 1.0) { 
                  if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l]) {
                    tmp_02 = 0.0;
                    for (ll = 0; ll < pOccRE; ll++) {
                      tmp_02 += betaStar[i * nOccRE + betaStarLongIndx[ll * J + j]];
	            } 
                    tmp_one[0] += kappaOcc[j * N + i] - (F77_NAME(ddot)(&pOcc, &X[i * JpOcc + j], &J, &beta[i], &N) + tmp_02 - betaStar[i * nOccRE + l] + wSites[j * N + i]) * omegaOcc[j * N + i];
	            tmp_0 += omegaOcc[j * N + i];
	          }
		}
              }
              /********************************
               * Compute A.beta.star
               *******************************/
              tmp_0 += 1.0 / sigmaSqPsi[betaStarIndx[l]]; 
              tmp_0 = 1.0 / tmp_0; 
              betaStar[i * nOccRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
            }

            // Update the RE sums for the current species
            zeros(&betaStarSites[i * J], J);
            for (j = 0; j < J; j++) {
              if (rangeInd[j * N + i] == 1.0) {
                for (l = 0; l < pOccRE; l++) {
                  betaStarSites[i * J + j] += betaStar[i * nOccRE + betaStarLongIndx[l * J + j]];
                }
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
                if ((z[zLongIndx[r] * N + i] == 1.0) && (XpRE[alphaStarIndx[l] * nObs + r] == alphaLevelIndx[l])) {
                  tmp_02 = 0.0;
                  for (ll = 0; ll < pDetRE; ll++) {
                    tmp_02 += alphaStar[i * nDetRE + alphaStarLongIndx[ll * nObs + r]];
	          } 
                  tmp_one[0] += kappaDet[r] - (F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + tmp_02 - alphaStar[i * nDetRE + l]) * omegaDet[r];
	  	  tmp_0 += omegaDet[r];
	        }
	      }
              /********************************
               * Compute A.alpha.star
               *******************************/
              tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
              tmp_0 = 1.0 / tmp_0; 
              alphaStar[i * nDetRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
            }
            zeros(&alphaStarObs[i * nObs], nObs); 
            // Update the RE sums for the current species
            for (r = 0; r < nObs; r++) {
              for (l = 0; l < pDetRE; l++) {
                alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + alphaStarLongIndx[l * nObs + r]]; 
              }
            }
          }
	} // i (species)

        /********************************************************************
         *Update Spatial Random Effects (w)
         *******************************************************************/

        for (rr = 0; rr < pTilde; rr++) { // svc

          // Update B and F for current SVC
	  for (ll = 0; ll < q; ll++) {
            // Current
            if (corName == "matern"){ 
	      nu[rr * q + ll] = theta[nuIndx * qpTilde + rr * q + ll];
       	    }
            updateBFsvcMs(&B[ll * nIndx], &F[ll * J], &c[ll * m*nThreads], 
	  		&C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, 
	  		theta[sigmaSqIndx * qpTilde + rr * q + ll], 
	  		theta[phiIndx * qpTilde + rr * q + ll], nu[rr * q + ll], 
	  	        covModel, &bk[ll * sizeBK], nuB[rr * q + ll]);
	  } // ll

	  for (ii = 0; ii < J; ii++) { // site
	    for (ll = 0; ll < q; ll++) { // factor

              a[ll] = 0; 
	      v[ll] = 0; 

	      if (uIndxLU[J + ii] > 0){ // is ii a neighbor for anybody
	        for (j = 0; j < uIndxLU[J+ii]; j++){ // how many locations have ii as a neighbor
	          b = 0;
	          // now the neighbors for the jth location who has ii as a neighbor
	          jj = uIndx[uIndxLU[ii]+j]; // jj is the index of the jth location who has ii as a neighbor
	          for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	            kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	            if(kk != ii){ //if the neighbor of jj is not ii
	      	    b += B[ll * nIndx + nnIndxLU[jj]+k] * w[rr * Jq + kk * q + ll]; //covariance between jj and kk and the random effect of kk
	            }
	          } // k
	          aij = w[rr * Jq + jj * q + ll] - b;
	          a[ll] += B[ll * nIndx + nnIndxLU[jj] + uiIndx[uIndxLU[ii] + j]] * 
			   aij / F[ll * J + jj];
	          v[ll] += pow(B[ll * nIndx + nnIndxLU[jj] + uiIndx[uIndxLU[ii] + j]], 2) / 
			   F[ll * J + jj];
	        } // j
	      }
	      
	      e = 0;
	      for(j = 0; j < nnIndxLU[J+ii]; j++){
	        e += B[ll * nIndx + nnIndxLU[ii] + j] * 
                     w[rr * Jq + nnIndx[nnIndxLU[ii]+j] * q + ll];
	      }

	      ff[ll] = 1.0 / F[ll * J + ii];
	      gg[ll] = e / F[ll * J + ii];
	    } // ll (factor)

	  // var
          // tmp_qq = (lambda * Xw)' S_beta (lambda * Xw)
	  zeros(tmp_Nq2, Nq);
	  zeros(tmp_Nq, Nq);
	  for (i = 0; i < N; i++) {
            if (rangeInd[ii * N + i] == 1.0) {
              for (ll = 0; ll < q; ll++) {
                tmp_Nq2[ll * N + i] = lambda[rr * Nq + ll * N + i] * Xw[i * JpTilde + rr * J + ii];
                tmp_Nq[ll * N + i] = tmp_Nq2[ll * N + i] * omegaOcc[ii * N + i];
              } // ll
	    }
          } // i
	  F77_NAME(dgemm)(ytran, ntran, &q, &q, &N, &one, tmp_Nq, &N, 
			  tmp_Nq2, &N, &zero, tmp_qq, &q FCONE FCONE);
	  F77_NAME(dcopy)(&qq, tmp_qq, &inc, var, &inc);
	  for (ll = 0; ll < q; ll++) {
            var[ll * q + ll] += ff[ll] + v[ll]; 
          } // ll
	  F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE);
          if(info != 0){error("c++ error: dpotrf var failed\n");}
	  F77_NAME(dpotri)(lower, &q, var, &q, &info FCONE);
          if(info != 0){error("c++ error: dpotri var failed\n");}

	  // mu
          // //Get Z*Lambda*w for all columns except for the kth column for the current location.
          zeros(tmp_N2, N);
          for (rrr = 0; rrr < pTilde; rrr++) {
            if (rrr != rr) { // if not on the current svc
              for (i = 0; i < N; i++) {
                if (rangeInd[ii * N + i] == 1.0) {
                  tmp_N2[i] += wStar[rrr * JN + ii * N + i] * Xw[i * JpTilde + rrr * J + ii];
		}
	      }
	    }
	  } // rrr svc)

	  zeros(tmp_N, N);
	  for (i = 0; i < N; i++) {
            if (rangeInd[ii * N + i] == 1.0) {
              tmp_N[i] = (zStar[ii * N + i] - F77_NAME(ddot)(&pOcc, &X[i * JpOcc + ii], 
				                             &J, &beta[i], &N) - 
			                    betaStarSites[i * J + ii] - tmp_N2[i]) * 
		       omegaOcc[ii * N + i];
	    }
          } // i

	  F77_NAME(dgemv)(ytran, &N, &q, &one, tmp_Nq2, &N, tmp_N, &inc, &zero, mu, &inc FCONE);

	  for (ll = 0; ll < q; ll++) {
            mu[ll] += gg[ll] + a[ll];
	  } // ll

	  F77_NAME(dsymv)(lower, &q, &one, var, &q, mu, &inc, &zero, tmp_q, &inc FCONE);

	  F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf var 2 failed\n");}

	  mvrnorm(&w[rr * Jq + ii * q], tmp_q, var, q);
	  } // ii (site)

  	  // Update wStar to make sure updated w values are used in lambda update.  
          for (j = 0; j < J; j++) {
            F77_NAME(dgemv)(ntran, &N, &q, &one, &lambda[rr * Nq], &N, 
			    &w[rr * Jq + j * q], &inc, &zero, &wStar[rr * JN + j * N], &inc FCONE);
          } // j

          /********************************************************************
           *Update spatial factors (lambda)
           *******************************************************************/
          for (i = 1; i < N; i++) {
            zeros(tmp_qq, qq);
            zeros(tmp_q, q);
	    zeros(tmp_qq2, qq);
	    zeros(tmp_Jq, Jq);
	    // Compute Xtilde %*% Z* %*% W for the current svc. 
	    for (ll = 0; ll < q; ll++) {
	      for (j = 0; j < J; j++) {
                if (rangeInd[j * N + i] == 1.0) {
                  tmp_Jq[j * q + ll] = w[rr * Jq + j * q + ll] * Xw[i * JpTilde + rr * J + j];
		}
	      } // j (site)
	    } // ll (factor)
	    // tmp_Jq' %*% S_beta %*% tmp_Jq 
            for (ll = 0; ll < q; ll++) {
              for (l = 0; l < q; l++) {
                for (j = 0; j < J; j++) {
                  if (rangeInd[j * N + i] == 1.0) {
                    tmp_qq[ll * q + l] += tmp_Jq[j * q + ll] * tmp_Jq[j * q + l] * omegaOcc[j * N + i];
		  }
                } // j
              } // l
            } // ll


            // currDim gives the mean dimension of mu and var. 
	    if (i < q) {
              currDim = i;  
            } else {
              currDim = q;
            }
            /*****************************
	     *mu
             *****************************/

	    // zStar - X %*% beta
	    for (j = 0; j < J; j++) {
              tmp_J[j] = 0.0;
              // //Get X*Lambda*w for all columns except for the kth column for the current location.
	      tmp_0 = 0.0;
	      if (rangeInd[j * N + i] == 1.0) {
                for (rrr = 0; rrr < pTilde; rrr++) {
                  if (rrr != rr) { // if not on the current svc
                    tmp_0 += wStar[rrr * JN + j * N + i] * Xw[i * JpTilde + rrr * J + j];
	          }
	        } // rrr (svc)
                tmp_J[j] = zStar[j * N + i] - F77_NAME(ddot)(&pOcc, &X[i * JpOcc + j], 
				                             &J, &beta[i], &N) - 
	                 betaStarSites[i * J + j] - tmp_0;

	        if (i < q) {
                  tmp_J[j] -= w[rr * Jq + j * q + i] * Xw[i * JpTilde + rr * J + j];
                }
	      }
            } // j

	    // S_beta %*% (W * Xw)' = tmp_Jq2
	    for (j = 0; j < J; j++) {
              if (rangeInd[j * N + i] == 1.0) {
                for (ll = 0; ll < q; ll++) {
                  tmp_Jq2[j * q + ll] = omegaOcc[j * N + i] * tmp_Jq[j * q + ll];  
                }
	      }
            }

	    // tmp_Jq2 %*% tmp_J
	    for (k = 0; k < currDim; k++) {
              for (j = 0; j < J; j++) {
                if (rangeInd[j * N + i] == 1.0) {
                  tmp_q[k] += tmp_Jq2[j * q + k] * tmp_J[j];
		}
              } // j
            } // k

            /*****************************
	     *var
             *****************************/
	    // Only get relevant columns of t(W) %*% W
	    for (k = 0, l = 0; k < currDim; k++) {
              for (ll = 0; ll < currDim; ll++, l++) {
                tmp_qq2[l] = tmp_qq[k * q + ll];
	      } // j
            } // k

	    // Add 1
	    for (ll = 0; ll < currDim; ll++) {
              tmp_qq2[ll * currDim + ll] += 1.0;  
            } // ll

            F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
	    if(info != 0){error("c++ error: dpotrf for spatial factors failed\n");}
            F77_NAME(dpotri)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
	    if(info != 0){error("c++ error: dpotri for spatial factors failed\n");}

            F77_NAME(dsymv)(lower, &currDim, &one, tmp_qq2, &currDim, tmp_q, &inc, &zero, tmp_q2, &inc FCONE);

            F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
	    if(info != 0){error("c++ error: dpotrf for spatial factors 2 failed\n");}
            
            mvrnorm(tmp_q, tmp_q2, tmp_qq2, currDim);
            F77_NAME(dcopy)(&currDim, tmp_q, &inc, &lambda[rr * Nq + i], &N);
          } // i

          // // Multiply Lambda %*% w[j] to get wStar.
          for (j = 0; j < J; j++) {
            F77_NAME(dgemv)(ntran, &N, &q, &one, &lambda[rr * Nq], &N, &w[rr * Jq + j * q], &inc, &zero, &wStar[rr * JN + j * N], &inc FCONE);
          } // j

          /********************************************************************
           *Update phi (and nu if matern)
           *******************************************************************/
	  for (ll = 0; ll < q; ll++) {
            aa = 0;
            logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
            for (j = 0; j < J; j++){
              if (nnIndxLU[J+j] > 0){
                e = 0;
                for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                  e += B[ll * nIndx + nnIndxLU[j]+ii]*w[rr * Jq + nnIndx[nnIndxLU[j]+ii] * q + ll];
                }
                b = w[rr * Jq + j * q + ll] - e;
              } else{
                b = w[rr * Jq + j * q + ll];
              }	
              aa += b*b/F[ll * J + j];
              logDet += log(F[ll * J + j]);
            }
      
            logPostCurr = -0.5 * logDet - 0.5 * aa;
            logPostCurr += log(theta[phiIndx * qpTilde + rr * q + ll] - 
	          	 phiA[rr * q + ll]) + log(phiB[rr * q + ll] - 
	          	 theta[phiIndx * qpTilde + rr * q + ll]); 
            if(corName == "matern"){
       	      logPostCurr += log(theta[nuIndx * qpTilde + rr * q + ll] - nuA[rr * q + ll]) + 
	                         log(nuB[rr * q + ll] - theta[nuIndx * qpTilde + rr * q + ll]); 
            }
            
            // Candidate
            phiCand = logitInv(rnorm(logit(theta[phiIndx * qpTilde + rr * q + ll], phiA[rr * q + ll], phiB[rr * q + ll]), exp(tuning[phiIndx * qpTilde + rr * q + ll])), phiA[rr * q + ll], phiB[rr * q + ll]);
            if (corName == "matern"){
      	      nuCand = logitInv(rnorm(logit(theta[nuIndx * qpTilde + rr * q + ll], nuA[rr * q + ll], nuB[rr * q + ll]), exp(tuning[nuIndx * qpTilde + rr * q + ll])), nuA[rr * q + ll], nuB[rr * q + ll]);
            }
      
            updateBFsvcMs(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * qpTilde + rr * q + ll], phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[rr * q + ll]);
      
            aa = 0;
            logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
            for (j = 0; j < J; j++){
              if (nnIndxLU[J+j] > 0){
                e = 0;
                for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                  e += BCand[nnIndxLU[j]+ii]*w[rr * Jq + nnIndx[nnIndxLU[j]+ii] * q + ll];
                }
                b = w[rr * Jq + j * q + ll] - e;
              } else{
                b = w[rr * Jq + j * q + ll];
                }	
                aa += b*b/FCand[j];
                logDet += log(FCand[j]);
            }
            
            logPostCand = -0.5*logDet - 0.5*aa;      
            logPostCand += log(phiCand - phiA[rr * q + ll]) + log(phiB[rr * q + ll] - phiCand); 
            if (corName == "matern"){
              logPostCand += log(nuCand - nuA[rr * q + ll]) + log(nuB[rr * q + ll] - nuCand); 
            }

            if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {

              F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[ll * nIndx], &inc);
              F77_NAME(dcopy)(&J, FCand, &inc, &F[ll * J], &inc);
              
	      theta[phiIndx * qpTilde + rr * q + ll] = phiCand;
              accept[phiIndx * qpTilde + rr * q + ll]++;
              if (corName == "matern") {
                nu[rr * q + ll] = nuCand; 
	        theta[nuIndx * qpTilde + rr * q + ll] = nu[rr * q + ll]; 
                accept[nuIndx * qpTilde + rr * q + ll]++; 
              }
            }
	  } // ll (factor)
        } // rr (svc)

        // Spatial process sums for each site and species
        for (i = 0; i < N; i++) {
          for (j = 0; j < J; j++) {
            wSites[j * N + i] = 0.0;
            if (rangeInd[j * N + i] == 1.0) {
              for (rr = 0; rr < pTilde; rr++) {
                wSites[j * N + i] += wStar[rr * JN + j * N + i] * Xw[i * JpTilde + rr * J + j];
              } // rr
	    }
          } // j
        } // i

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        for (i = 0; i < N; i++) {
          // Compute detection probability 
          if (nObs == J) {
            for (r = 0; r < nObs; r++) {
              if (rangeInd[zLongIndx[r] * N + i] == 1.0) {
                detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
                psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[i * JpOcc + zLongIndx[r]], &J, &beta[i], &N) + wSites[zLongIndx[r] * N + i] + betaStarSites[i * J + zLongIndx[r]], zero, one); 
                piProd[zLongIndx[r] * N + i] = pow(1.0 - detProb[i * nObs + r], K[r]);
	        piProdWAIC[zLongIndx[r] * N + i] *= pow(detProb[i * nObs + r], y[r * N + i]);
	        piProdWAIC[zLongIndx[r] * N + i] *= pow(1.0 - detProb[i * nObs + r], K[r] - y[r * N + i]);
                ySum[zLongIndx[r] * N + i] = y[r * N + i]; 
	      }
            } // r
          } else {
            for (r = 0; r < nObs; r++) {
              if (rangeInd[zLongIndx[r] * N + i]) {
                detProb[i * nObs + r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &N) + alphaStarObs[i * nObs + r], zero, one);
                if (tmp_JInt[zLongIndx[r]] == 0) {
                  psi[zLongIndx[r] * N + i] = logitInv(F77_NAME(ddot)(&pOcc, &X[i * JpOcc + zLongIndx[r]], &J, &beta[i], &N) + wSites[zLongIndx[r] * N + i] + betaStarSites[i * J + zLongIndx[r]], zero, one); 
                }
                piProd[zLongIndx[r] * N + i] *= (1.0 - detProb[i * nObs + r]);
	        piProdWAIC[zLongIndx[r] * N + i] *= pow(detProb[i * nObs + r], y[r * N + i]);
	        piProdWAIC[zLongIndx[r] * N + i] *= pow(1.0 - detProb[i * nObs + r], 
	            	                            1.0 - y[r * N + i]);
                ySum[zLongIndx[r] * N + i] += y[r * N + i]; 
                tmp_JInt[zLongIndx[r]]++;
	      }
            } // r
          }
          // Compute occupancy probability and the integrated likelihood for WAIC
          for (j = 0; j < J; j++) {
            if (rangeInd[j * N + i] == 1.0) {
              psiNum = psi[j * N + i] * piProd[j * N + i]; 
              if (ySum[j * N + i] == zero) {
                z[j * N + i] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j * N + i])));
                yWAIC[j * N + i] = (1.0 - psi[j * N + i]) + psi[j * N + i] * piProdWAIC[j * N + i];
              } else {
                z[j * N + i] = one; 
                yWAIC[j * N + i] = psi[j * N + i] * piProdWAIC[j * N + i]; 
              }
	    } else {
              z[j * N + i] = 0.0;
	      psi[j * N + i] = 0.0;
	    }
            // Reset variables
            piProd[j * N + i] = one;
	    piProdWAIC[j * N + i] = one;
            ySum[j * N + i] = zero; 
            tmp_JInt[j] = 0; 
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
            F77_NAME(dcopy)(&NqpTilde, lambda, &inc, &REAL(lambdaSamples_r)[sPost*NqpTilde], &inc); 
            F77_NAME(dcopy)(&JN, psi, &inc, &REAL(psiSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&JN, z, &inc, &REAL(zSamples_r)[sPost*JN], &inc); 
            F77_NAME(dcopy)(&JqpTilde, w, &inc, &REAL(wSamples_r)[sPost*JqpTilde], &inc); 
            F77_NAME(dcopy)(&nThetaqpTildeSave, &theta[phiIndx * qpTilde], &inc, &REAL(thetaSamples_r)[sPost*nThetaqpTildeSave], &inc); 
	    if (pDetRE > 0) {
              F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
              F77_NAME(dcopy)(&nDetREN, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetREN], &inc);
	    }
	    if (pOccRE > 0) {
              F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
              F77_NAME(dcopy)(&nOccREN, betaStar, &inc, &REAL(betaStarSamples_r)[sPost*nOccREN], &inc);
	    }
            F77_NAME(dcopy)(&JN, yWAIC, &inc, &REAL(likeSamples_r)[sPost*JN], &inc); 
            sPost++; 
            thinIndx = 0; 
          }
        }
        R_CheckUserInterrupt();
      } // t (end batch)

      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (rr = 0; rr < pTilde; rr++) {
        for (ll = 0; ll < q; ll++) {
          for (k = 0; k < nTheta; k++) {
            REAL(acceptSamples_r)[s * nThetaqpTilde + k * qpTilde + rr * q + ll] = accept[k * qpTilde + rr * q + ll] / batchLength; 
            REAL(tuningSamples_r)[s * nThetaqpTilde + k * qpTilde + rr * q + ll] = tuning[k * qpTilde + rr * q + ll]; 
            if (accept[k * qpTilde + rr * q + ll] / batchLength > acceptRate) {
              tuning[k * qpTilde + rr * q + ll] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
            } else{
                tuning[k * qpTilde + rr * q + ll] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
              }
            accept[k * qpTilde + rr * q + ll] = 0.0;
          } // k
        } // ll
      } // rr

      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tCoefficient\tLatent Factor\tAcceptance\tTuning\n");	  
	  for (rr = 0; rr < pTilde; rr++) {
            for (ll = 0; ll < q; ll++) {
              Rprintf("\t%i\t\t%i\t\t%3.1f\t\t%1.5f\n", rr + 1, ll + 1, 100.0*REAL(acceptSamples_r)[s * nThetaqpTilde + phiIndx * qpTilde + rr * q + ll], exp(tuning[phiIndx * qpTilde + rr * q + ll]));
            } // ll
	  } // rr
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
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }

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
    SET_VECTOR_ELT(result_r, 8, lambdaSamples_r);
    SET_VECTOR_ELT(result_r, 9, wSamples_r); 
    SET_VECTOR_ELT(result_r, 10, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 11, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 12, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 13, likeSamples_r); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 14, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 15, alphaStarSamples_r);
    }
    if (pOccRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 16;
      } else {
        tmp_0 = 14;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }

    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("psi.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("lambda.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 10, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 11, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 12, mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 13, mkChar("like.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 14, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 15, mkChar("alpha.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, mkChar("beta.star.samples")); 
    }
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


