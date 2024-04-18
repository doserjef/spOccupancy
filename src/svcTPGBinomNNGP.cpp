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

void updateBFSVCTBin(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

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
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){Rf_error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){Rf_error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
  SEXP svcTPGBinomNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP coords_r, SEXP XRE_r, 
		       SEXP consts_r, SEXP weights_r, SEXP nRELong_r, 
		       SEXP m_r, SEXP nnIndx_r, SEXP nnIndxLU_r, 
		       SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
		       SEXP betaStarting_r, SEXP sigmaSqPsiStarting_r, SEXP betaStarStarting_r, 
	               SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
		       SEXP wStarting_r, SEXP zYearIndx_r, SEXP zDatIndx_r, 
		       SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
		       SEXP muBeta_r, SEXP SigmaBeta_r, 
		       SEXP phiA_r, SEXP phiB_r, SEXP sigmaSqA_r, SEXP sigmaSqB_r,
		       SEXP nuA_r, SEXP nuB_r, SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
		       SEXP ar1_r, SEXP ar1Vals_r,
		       SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
	               SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r, 
	               SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, 
		       SEXP currChain_r, SEXP nChain_r, SEXP sigmaSqIG_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, r, ll, rr, ii, t, info, nProtect=0, l;
    int status = 0; // For AMCMC. 
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
    // Order: covariate, site
    double *Xw = REAL(Xw_r);
    int *XRE = INTEGER(XRE_r); 
    double *weights = REAL(weights_r);
    int m = INTEGER(m_r)[0]; 
    // Load constants
    int J = INTEGER(consts_r)[0];
    int p = INTEGER(consts_r)[1];
    int pRE = INTEGER(consts_r)[2];
    int nRE = INTEGER(consts_r)[3];
    int nYearsMax = INTEGER(consts_r)[4];
    int pTilde = INTEGER(consts_r)[5];
    int pp = p * p;
    int ppTilde = pTilde * pTilde;
    int JpTilde = J * pTilde;
    /**********************************
     * Priors
     * *******************************/
    double *muBeta = (double *) R_alloc(p, sizeof(double));   
    F77_NAME(dcopy)(&p, REAL(muBeta_r), &inc, muBeta, &inc);
    double *SigmaBetaInv = (double *) R_alloc(pp, sizeof(double));   
    F77_NAME(dcopy)(&pp, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *phiA = REAL(phiA_r);
    double *phiB = REAL(phiB_r); 
    double *nuA = REAL(nuA_r); 
    double *nuB = REAL(nuB_r); 
    double *sigmaSqA = REAL(sigmaSqA_r); 
    double *sigmaSqB = REAL(sigmaSqB_r); 
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *tuning = REAL(tuning_r); 
    double *coords = REAL(coords_r);
    int *nRELong = INTEGER(nRELong_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int *zDatIndx = INTEGER(zDatIndx_r); 
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    /**********************************
     * AR1 Parameters 
     * *******************************/
    int ar1 = INTEGER(ar1_r)[0];
    double rhoA = REAL(ar1Vals_r)[0];
    double rhoB = REAL(ar1Vals_r)[1];
    double sigmaSqTA = REAL(ar1Vals_r)[2];
    double sigmaSqTB = REAL(ar1Vals_r)[3];
    int sigmaSqIG = INTEGER(sigmaSqIG_r)[0];
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
    int thinIndx = 0; 
    int sPost = 0; 

    // Some constants
    int JnYears = J * nYearsMax;

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
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
        Rprintf("Spatial NNGP Multi-season Binomial Model with Polya-Gamma latent\nvariable fit with %i sites and %i years.\n\n", J, nYearsMax);
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
	Rprintf("Number of spatially-varying coefficients: %i \n", pTilde);
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
        Rprintf("Using %i nearest neighbors.\n\n", m);
	if (ar1) {
          Rprintf("Using an AR(1) temporal autocorrelation matrix.\n\n");
	}
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
     #ifdef Win32
       R_FlushConsole();
     #endif

    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Occurrence fixed effects
    double *beta = (double *) R_alloc(p, sizeof(double));   
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pRE, sizeof(double)); 
    F77_NAME(dcopy)(&pRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nRE, sizeof(double)); 
    F77_NAME(dcopy)(&nRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Spatial processes
    double *w = (double *) R_alloc(JpTilde, sizeof(double));
    F77_NAME(dcopy)(&JpTilde, REAL(wStarting_r), &inc, w, &inc); 
    // Spatial variance
    double *sigmaSq = (double *) R_alloc(pTilde, sizeof(double)); 
    F77_NAME(dcopy)(&pTilde, REAL(sigmaSqStarting_r), &inc, sigmaSq, &inc); 
    // Spatial range parameter
    double *phi = (double *) R_alloc(pTilde, sizeof(double)); 
    F77_NAME(dcopy)(&pTilde, REAL(phiStarting_r), &inc, phi, &inc); 
    // Spatial smoothing parameter for Matern
    double *nu = (double *) R_alloc(pTilde, sizeof(double)); 
    F77_NAME(dcopy)(&pTilde, REAL(nuStarting_r), &inc, nu, &inc); 
    // PG auxiliary variables 
    double *omega = (double *) R_alloc(JnYears, sizeof(double)); zeros(omega, JnYears);
    double *kappa = (double *) R_alloc(JnYears, sizeof(double)); zeros(kappa, JnYears);
    double *yStar = (double *) R_alloc(JnYears, sizeof(double)); zeros(yStar, JnYears);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, p, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), p * nPost);
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = Rf_allocMatrix(REALSXP, JnYears, nPost)); nProtect++; 
    zeros(REAL(yRepSamples_r), JnYears * nPost);
    SEXP wSamples_r; 
    PROTECT(wSamples_r = Rf_allocMatrix(REALSXP, JpTilde, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), JpTilde * nPost);
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, JnYears, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), JnYears * nPost);
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = Rf_allocMatrix(REALSXP, pRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPsiSamples_r), pRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nRE * nPost);
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, JnYears, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), JnYears * nPost);
    // AR(1) stuff. 
    SEXP etaSamples_r; 
    if (ar1) {
      PROTECT(etaSamples_r = Rf_allocMatrix(REALSXP, nYearsMax, nPost)); nProtect++; 
      zeros(REAL(etaSamples_r), nYearsMax * nPost);
    }
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int JnYearspRE = J * nYearsMax * pRE; 
    int jj, kk;
    double tmp_0 = 0.0;
    double tmp_02;
    int nnYears = nYearsMax * nYearsMax;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); 
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double));
    zeros(tmp_pp2, pp);
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_one = (double *) R_alloc(1, sizeof(double)); 
    int *tmp_JnYearsInt = (int *) R_alloc(JnYears, sizeof(int));
    for (j = 0; j < JnYears; j++) {
      tmp_JnYearsInt[j] = zero; 
    }
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J1, J);
    double *tmp_JnYears = (double *) R_alloc(JnYears, sizeof(double));
    zeros(tmp_JnYears, JnYears);
    double *tmp_JnYearsp = (double *) R_alloc(JnYears * p, sizeof(double));
    zeros(tmp_JnYearsp, JnYears * p);
    double *tmp_pTilde = (double *) R_alloc(pTilde, sizeof(double));
    zeros(tmp_pTilde, pTilde);
    double * tmp_ppTilde = (double *) R_alloc(ppTilde, sizeof(double));
    zeros(tmp_ppTilde, ppTilde);

    // For latent occupancy
    double *psi = (double *) R_alloc(JnYears, sizeof(double)); 
    zeros(psi, JnYears); 
    double *like = (double *) R_alloc(JnYears, sizeof(double)); ones(like, JnYears);
    double *yRep = (double *) R_alloc(JnYears, sizeof(double)); zeros(yRep, JnYears);

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &p, SigmaBetaInv, &p, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &p, SigmaBetaInv, &p, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(p, sizeof(double)); 
    F77_NAME(dsymv)(lower, &p, &one, SigmaBetaInv, &p, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site/year-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JnYears, sizeof(double)); 
    zeros(betaStarSites, JnYears); 
    int *betaStarLongIndx = (int *) R_alloc(JnYearspRE, sizeof(int));
    // Initial sums
    for (t = 0; t < nYearsMax; t++) {
      for (j = 0; j < J; j++) {
        for (l = 0; l < pRE; l++) {
          betaStarLongIndx[l * JnYears + t * J + j] = which(XRE[l * JnYears + t * J + j], 
                                                            betaLevelIndx, nRE);
          betaStarSites[t * J + j] += betaStar[betaStarLongIndx[l * JnYears + t * J + j]];
        } // l
      } // j
    } // t
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pRE, sizeof(int)); 
    for (l = 0; l < pRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nRE); 
    }

    /**********************************************************************
     * Set up spatial stuff, MH stuff, and AR1 stuff
     *********************************************************************/
    int nTheta, sigmaSqIndx, phiIndx, nuIndx, sigmaSqTIndx, rhoIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }  
    int nThetapTilde = nTheta * pTilde;
    int nThetaAll = nThetapTilde;
    sigmaSqTIndx = nThetapTilde;
    rhoIndx = sigmaSqTIndx + 1;
    if (ar1 == 1) {
      nThetaAll += 2;
    }
    double *accept = (double *) R_alloc(nThetaAll, sizeof(double)); zeros(accept, nThetaAll); 
    double *theta = (double *) R_alloc(nThetaAll, sizeof(double));
    double logPostCurr = 0.0, logPostCand = 0.0;
    double logDet;  
    double phiCand = 0.0, nuCand = 0.0, rhoCand = 0.0, sigmaSqCand = 0.0;  
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = Rf_allocMatrix(REALSXP, nThetaAll, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nThetaAll * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = Rf_allocMatrix(REALSXP, nThetaAll, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nThetaAll * nBatch);
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = Rf_allocMatrix(REALSXP, nThetaAll, nPost)); nProtect++; 
    zeros(REAL(thetaSamples_r), nThetaAll * nPost);
    double b, e, aij, aa; 
    double *a = (double *) R_alloc(pTilde, sizeof(double));
    double *v = (double *) R_alloc(pTilde, sizeof(double));
    double *mu = (double *) R_alloc(pTilde, sizeof(double));
    double *var = (double *) R_alloc(ppTilde, sizeof(double)); zeros(var, ppTilde);
    double *ff = (double *) R_alloc(pTilde, sizeof(double));
    double *gg = (double *) R_alloc(pTilde, sizeof(double));
    // Initiate spatial values
    for (i = 0; i < pTilde; i++) {
      theta[sigmaSqIndx * pTilde + i] = sigmaSq[i]; 
      theta[phiIndx * pTilde + i] = phi[i]; 
      if (corName == "matern") {
        theta[nuIndx * pTilde + i] = nu[i]; 
      } 
    } // i
    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * pTilde, sizeof(double));
    double *F = (double *) R_alloc(J * pTilde, sizeof(double));
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads * pTilde, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads * pTilde, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(pTilde*sizeBK, sizeof(double));

    // Initiate B and F for each SVC
    for (i = 0; i < pTilde; i++) {
    updateBFSVCTBin(&B[i * nIndx], &F[i*J], &c[i * m*nThreads], &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * pTilde + i], theta[phiIndx * pTilde + i], nu[i], covModel, &bk[i * sizeBK], nuB[i]);
    }
    // Spatial process sums for each site/year. 
    double *wSites = (double *) R_alloc(JnYears, sizeof(double)); zeros(wSites, JnYears);
    // For each location, multiply w x Xw
    for (t = 0; t < nYearsMax; t++) {
      wSites[t * J + j] = 0.0;
      for (j = 0; j < J; j++) {
        for (ll = 0; ll < pTilde; ll++) {
          wSites[t * J + j] += w[j * pTilde + ll] * Xw[ll * JnYears + t * J + j];
        }
      }
    }

    /**********************************************************************
     * Set up AR1 stuff
     *********************************************************************/
    double rho; 
    if (ar1) {
      theta[sigmaSqTIndx] = REAL(ar1Vals_r)[5];
      theta[rhoIndx] = REAL(ar1Vals_r)[4];
      rho = theta[rhoIndx];
    }
    double *SigmaEta = (double *) R_alloc(nnYears, sizeof(double));
    double *SigmaEtaCand = (double *) R_alloc(nnYears, sizeof(double));
    double *tmp_nYearsMax = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nYearsMax2 = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nnYears = (double *) R_alloc(nnYears, sizeof(double));
    if (ar1) {
      AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
      clearUT(SigmaEta, nYearsMax);
      F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){Rf_error("c++ error: Cholesky failed in initial AR(1) covariance matrix\n");}
      F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){Rf_error("c++ error: Cholesky inverse failed in initial AR(1) covariance matrix\n");}
    }
    double *eta = (double *) R_alloc(nYearsMax, sizeof(double)); zeros(eta, nYearsMax);
    // For sigmaSqT sampler
    double aSigmaSqTPost = 0.5 * nYearsMax + sigmaSqTA;
    double bSigmaSqTPost = 0.0;
    double *etaTRInv = (double *) R_alloc(nYearsMax, sizeof(double));

    GetRNGstate();

    /**********************************************************************
     * Begin Sampler 
     * *******************************************************************/
    for (s = 0, rr = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, rr++) {
        
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
	zeros(tmp_JnYears, JnYears);
	for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            // Only calculate omega when there are observations at that 
	    // site/year combo. Otherwise, you're just wasting time. 
	    if (zDatIndx[t * J + j] == 1) { 
              tmp_JnYears[t * J + j] = F77_NAME(ddot)(&p, &X[t * J + j], &JnYears, 
	         	                      beta, &inc);
	      omega[t * J + j] = rpg(weights[t * J + j], tmp_JnYears[t * J + j] + 
			                                 wSites[t * J + j] + 
							 betaStarSites[t * J + j] + 
							 eta[t]);
	      // Update kappa values along the way. 
              kappa[t * J + j] = y[t * J + j] - weights[t * J + j] / 2.0; 
	      yStar[t * J + j] = kappa[t * J + j] / omega[t * J + j];
	    }
	  } // j 
	} // t
        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        zeros(tmp_JnYears, JnYears);
	for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            if (zDatIndx[t * J + j] == 1) {
              tmp_JnYears[t * J + j] = kappa[t * J + j] - omega[t * J + j] * (wSites[t * J + j] + betaStarSites[t * J + j] + eta[t]); 
	    }
	  } // j
	} // t
        /********************************
         * Compute b.beta
         *******************************/
	// This is fine, because the elements in tmp_JnYears corresponding
	// to unobserve site/time locations is set to 0 and not changed. 
        F77_NAME(dgemv)(ytran, &JnYears, &p, &one, X, &JnYears, 
			tmp_JnYears, &inc, &zero, tmp_p, &inc FCONE); 	 
        for (j = 0; j < p; j++) {
          tmp_p[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
	// This is fine, because omega == 0 for the site/year combos
	// that don't have any observations at them, which will cause this
	// whole product to go to 0.  
	for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            if (zDatIndx[t * J + j] == 1) {
              for(i = 0; i < p; i++){
                tmp_JnYearsp[i * JnYears + t * J + j] = X[i * JnYears + t * J + j] * omega[t * J + j];
              } // i
	    }
	  } // j
	} // t

        F77_NAME(dgemm)(ytran, ntran, &p, &p, &JnYears, &one, X, &JnYears, tmp_JnYearsp, &JnYears, &zero, tmp_pp, &p FCONE FCONE);
        for (j = 0; j < pp; j++) {
          tmp_pp[j] += SigmaBetaInv[j]; 
        } // j


        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf A.beta failed\n");}
        F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri A.beta failed\n");}
        F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
	if(info != 0){Rf_error("c++ error: dpotrf A.beta2 failed\n");}
        mvrnorm(beta, tmp_p2, tmp_pp, p);

        /********************************************************************
         *Update Occupancy random effects variance
         *******************************************************************/
        for (l = 0; l < pRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqPsi[l] = rigamma(sigmaSqPsiA[l] + nRELong[l] / 2.0, sigmaSqPsiB[l] + tmp_0); 
        }

        /********************************************************************
         *Update Occupancy random effects
         *******************************************************************/
        if (pRE > 0) {
          // Update each individual random effect one by one. 
          for (l = 0; l < nRE; l++) {
            /********************************
             * Compute b.beta.star
             *******************************/
            zeros(tmp_one, inc);
            tmp_0 = 0.0;	      
            // Only allow information to come from when XRE == betaLevelIndx[l]. 
            // aka information only comes from the sites with any given level 
            // of a random effect. 
	    for (t = 0; t < nYearsMax; t++) {
              for (j = 0; j < J; j++) {
                if (XRE[betaStarIndx[l] * JnYears + t * J + j] == betaLevelIndx[l]) {
                   if (zDatIndx[t * J + j] == 1) {
                    tmp_02 = 0.0;
                    for (ll = 0; ll < pRE; ll++) {
                      tmp_02 += betaStar[betaStarLongIndx[ll * JnYears + t * J + j]];
	            } 
                    tmp_one[0] += kappa[t * J + j] - (F77_NAME(ddot)(&p, &X[t * J + j], &JnYears, beta, &inc) + 
                  	          tmp_02 - betaStar[l] + wSites[t * J + j] + eta[t]) * omega[t * J + j];
                    tmp_0 += omega[t * J + j];
		  }
                }
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
          zeros(betaStarSites, JnYears);
	  for (t = 0; t < nYearsMax; t++) {
            for (j = 0; j < J; j++) {
              for (l = 0; l < pRE; l++) {
                betaStarSites[t * J + j] += betaStar[betaStarLongIndx[l * JnYears + t * J + j]];
              }
            }
	  }
        }

        /********************************************************************
         *Update w (spatial random effects)
         *******************************************************************/
	for (ii = 0; ii < J; ii++ ) {
          zeros(tmp_ppTilde, ppTilde);
          for (ll = 0; ll < pTilde; ll++) {
            for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + ii] == 1) {
                tmp_0 = Xw[ll * JnYears + t * J + ii] * omega[t * J + ii];
	        for (k = 0; k < pTilde; k++) {
                  tmp_ppTilde[ll * pTilde + k] += tmp_0 * Xw[k * JnYears + t * J + ii];
	        }
	      }
	    } // t
            a[ll] = 0.0;
	    v[ll] = 0.0;
	    if (uIndxLU[J + ii] > 0){ // is i a neighbor for anybody
	      for (j = 0; j < uIndxLU[J+ii]; j++){ // how many locations have ii as a neighbor
	        b = 0.0;
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[ii]+j]; // jj is the index of the jth location who has ii as a neighbor
	        for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	          kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	          if(kk != ii){ //if the neighbor of jj is not ii
	    	    b += B[ll * nIndx + nnIndxLU[jj]+k]*w[kk * pTilde + ll]; //covariance between jj and kk and the random effect of kk
	          }
	        }
	        aij = w[jj * pTilde + ll] - b;
	        a[ll] += B[ll * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]]*aij/F[ll * J + jj];
	        v[ll] += pow(B[ll * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]],2)/F[ll * J + jj];
	      }
	    }
	    e = 0.0;
	    for(j = 0; j < nnIndxLU[J+ii]; j++){
	      e += B[ll * nIndx + nnIndxLU[ii]+j]*w[nnIndx[nnIndxLU[ii]+j] * pTilde + ll];
	    }
	    ff[ll] = 1.0 / F[ll * J + ii];
	    gg[ll] = e / F[ll * J + ii];
	  } // ll (svc)

	  // var
	  F77_NAME(dcopy)(&ppTilde, tmp_ppTilde, &inc, var, &inc);
	  for (k = 0; k < pTilde; k++) {
            var[k * pTilde + k] += ff[k] + v[k]; 
          } // k
	  F77_NAME(dpotrf)(lower, &pTilde, var, &pTilde, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotrf var failed\n");}
	  F77_NAME(dpotri)(lower, &pTilde, var, &pTilde, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri var failed\n");}

	  // mu
	  for (k = 0; k < pTilde; k++) {
            mu[k] = 0.0;
	    for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + ii] == 1) {
                mu[k] += (yStar[t * J + ii] - F77_NAME(ddot)(&p, &X[t * J + ii], &JnYears, beta, &inc) - betaStarSites[t * J + ii] - eta[t]) * omega[t * J + ii] * Xw[k * JnYears + t * J + ii];
	      }
	    } // t
	    mu[k] += gg[k] + a[k];
	  } // k
	  F77_NAME(dsymv)(lower, &pTilde, &one, var, &pTilde, mu, &inc, &zero, tmp_pTilde, &inc FCONE);

	  F77_NAME(dpotrf)(lower, &pTilde, var, &pTilde, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf var 2 failed\n");}

	  mvrnorm(&w[ii * pTilde], tmp_pTilde, var, pTilde);
        } // ii (site)

	// Compute wSites. 
        for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            wSites[t * J + j] = 0.0;
            for (ll = 0; ll < pTilde; ll++) {
              wSites[t * J + j] += w[j * pTilde + ll] * Xw[ll * JnYears + t * J + j];
	    }
          }
        }

        /********************************************************************
         *Update spatial covariance parameters
         *******************************************************************/
	for (ll = 0; ll < pTilde; ll++) {
          /******************************************************************
           *Update sigmaSq
           *****************************************************************/
          if (sigmaSqIG) {
            aa = 0;
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:aa, logDet)
#endif
            for (j = 0; j < J; j++){
              if(nnIndxLU[J+j] > 0){
                e = 0;
                for(i = 0; i < nnIndxLU[J+j]; i++){
                  e += B[ll * nIndx + nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i] * pTilde + ll];
                }
                b = w[j * pTilde + ll] - e;
              }else{
                b = w[j * pTilde + ll];
              }	
              aa += b*b/F[ll * J + j];
            }

	    theta[sigmaSqIndx * pTilde + ll] = rigamma(sigmaSqA[ll] + J / 2.0, sigmaSqB[ll] + 0.5 * aa * theta[sigmaSqIndx * pTilde + ll]); 
	  }

          /******************************************************************
           *Update phi (and nu if matern)
           *****************************************************************/
          // Current
          if (corName == "matern"){ 
	    nu[ll] = theta[nuIndx * pTilde + ll];
       	  }
          updateBFSVCTBin(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * pTilde + ll], theta[phiIndx * pTilde + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[ll]);
          aa = 0;
          logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += B[ll * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * pTilde + ll];
              }
              b = w[j * pTilde + ll] - e;
            } else{
              b = w[j * pTilde + ll];
            }	
            aa += b*b/F[ll * J + j];
            logDet += log(F[ll * J + j]);
          }
      
          logPostCurr = -0.5 * logDet - 0.5 * aa;
          logPostCurr += log(theta[phiIndx * pTilde + ll] - phiA[ll]) + log(phiB[ll] - theta[phiIndx * pTilde + ll]); 
          if(corName == "matern"){
       	    logPostCurr += log(theta[nuIndx * pTilde + ll] - nuA[ll]) + log(nuB[ll] - theta[nuIndx * pTilde + ll]); 
          }
	  if (sigmaSqIG == 0) {
            logPostCurr += log(theta[sigmaSqIndx * pTilde + ll] - sigmaSqA[ll]) +
		           log(sigmaSqB[ll] - theta[sigmaSqIndx * pTilde + ll]);
	  }
          
          // Candidate
          phiCand = logitInv(rnorm(logit(theta[phiIndx * pTilde + ll], phiA[ll], phiB[ll]), exp(tuning[phiIndx * pTilde + ll])), phiA[ll], phiB[ll]);
          if (corName == "matern"){
      	    nuCand = logitInv(rnorm(logit(theta[nuIndx * pTilde + ll], nuA[ll], nuB[ll]), exp(tuning[nuIndx * pTilde + ll])), nuA[ll], nuB[ll]);
          }
	  if (sigmaSqIG == 0) {
            sigmaSqCand = logitInv(rnorm(logit(theta[sigmaSqIndx * pTilde + ll], sigmaSqA[ll], sigmaSqB[ll]),
				         exp(tuning[sigmaSqIndx * pTilde + ll])), sigmaSqA[ll], sigmaSqB[ll]);
	  }
      
      	  if (sigmaSqIG) { 
            updateBFSVCTBin(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * pTilde + ll], phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
	  } else {
            updateBFSVCTBin(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, sigmaSqCand, phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
	  }

          aa = 0;
          logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += BCand[nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * pTilde + ll];
              }
              b = w[j * pTilde + ll] - e;
            } else{
              b = w[j * pTilde + ll];
              }	
              aa += b*b/FCand[j];
              logDet += log(FCand[j]);
          }
          
          logPostCand = -0.5*logDet - 0.5*aa;      
          logPostCand += log(phiCand - phiA[ll]) + log(phiB[ll] - phiCand); 
          if (corName == "matern"){
            logPostCand += log(nuCand - nuA[ll]) + log(nuB[ll] - nuCand); 
          }
	  if (sigmaSqIG == 0) {
            logPostCand += log(sigmaSqCand - sigmaSqA[ll]) + log(sigmaSqB[ll] - sigmaSqCand);
	  }

          if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {

            F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[ll * nIndx], &inc);
            F77_NAME(dcopy)(&J, FCand, &inc, &F[ll * J], &inc);
            
	    theta[phiIndx * pTilde + ll] = phiCand;
            accept[phiIndx * pTilde + ll]++;
            if (corName == "matern") {
              nu[ll] = nuCand; 
	      theta[nuIndx * pTilde + ll] = nu[ll]; 
              accept[nuIndx * pTilde + ll]++; 
            }
	    if (sigmaSqIG == 0) {
              theta[sigmaSqIndx * pTilde + ll] = sigmaSqCand;
	      accept[sigmaSqIndx * pTilde + ll]++;
	    }
          }
	} // ll

	if (ar1) {
          /********************************************************************
           *Update sigmaSqT
           *******************************************************************/
	  // Form correlation matrix. 
          AR1(nYearsMax, theta[rhoIndx], 1.0, SigmaEta);
	  clearUT(SigmaEta, nYearsMax);
	  F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
	  F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky inverse failed in covariance matrix\n");}
	  fillUTri(SigmaEta, nYearsMax);
	  // Compute t(eta) %*% SigmaEta^-1 %*% eta
          for (t = 0; t < nYearsMax; t++) {
            etaTRInv[t] = F77_NAME(ddot)(&nYearsMax, &SigmaEta[t], &nYearsMax, 
	  		               eta, &inc); 
	  }
          bSigmaSqTPost = F77_NAME(ddot)(&nYearsMax, etaTRInv, &inc, eta, &inc);	
	  bSigmaSqTPost /= 2.0;
	  bSigmaSqTPost += sigmaSqTB;
	  theta[sigmaSqTIndx] = rigamma(aSigmaSqTPost, bSigmaSqTPost);
          
          /********************************************************************
           *Update rho
           *******************************************************************/
	  rho = theta[rhoIndx];
	  rhoCand = logitInv(rnorm(logit(rho, rhoA, rhoB), exp(tuning[rhoIndx])), rhoA, rhoB); 
	  theta[rhoIndx] = rhoCand; 

	  // Construct proposal covariance matrix. 
          AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEtaCand);
	  clearUT(SigmaEtaCand, nYearsMax);

          /********************************
           * Proposal
           *******************************/
	  // Invert SigmaEtaCand and log det cov. 
          logPostCand = 0.0;
	  F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky failed in proposal covariance matrix\n");}
	  // Get log of the determinant of the covariance matrix. 
	  for (k = 0; k < nYearsMax; k++) {
	    logPostCand += 2.0 * log(SigmaEtaCand[k*nYearsMax+k]);
	  } // k
	  F77_NAME(dpotri)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky inverse failed in proposal covariance matrix\n");}
          logPostCand = 0.0; 
	  // Jacobian and Uniform prior. 
	  logPostCand += log(rhoCand - rhoA) + log(rhoB - rhoCand); 
	  F77_NAME(dsymv)(lower, &nYearsMax, &one,  SigmaEtaCand, &nYearsMax, eta, &inc, &zero, tmp_nYearsMax, &inc FCONE);
	  logPostCand += -0.5*logPostCand-0.5*F77_NAME(ddot)(&nYearsMax, eta, &inc, tmp_nYearsMax, &inc);
          /********************************
           * Current
           *******************************/
	  theta[rhoIndx] = rho; 
          AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
	  clearUT(SigmaEta, nYearsMax);
          logPostCurr = 0.0;
	  F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
	  for (k = 0; k < nYearsMax; k++) {
	    logPostCurr += 2.0 * log(SigmaEta[k*nYearsMax+k]);
	  } // k
	  F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: Cholesky inverse failed in covariance matrix\n");}
          logPostCurr = 0.0; 
	  logPostCurr += log(rho - rhoA) + log(rhoB - rho); 
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &nYearsMax, &one, SigmaEta, &nYearsMax, eta, &inc, &zero, 
	  		tmp_nYearsMax, &inc FCONE);
	  logPostCurr += -0.5*logPostCurr-0.5*F77_NAME(ddot)(&nYearsMax, eta, &inc, tmp_nYearsMax, &inc);

	  // MH Accept/Reject
	  if (runif(0.0, 1.0) <= exp(logPostCand - logPostCurr)) {
            theta[rhoIndx] = rhoCand;
            accept[rhoIndx]++;
	    F77_NAME(dcopy)(&nnYears, SigmaEtaCand, &inc, SigmaEta, &inc); 
          }
          /********************************************************************
           *Update eta 
           *******************************************************************/
          /********************************
           * Compute b.w
           *******************************/
	  zeros(tmp_nYearsMax, nYearsMax);
          for (j = 0; j < J; j++) {
            for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + j] == 1) {
                tmp_nYearsMax[t] += kappa[t * J + j] - omega[t * J + j] * (F77_NAME(ddot)(&p, &X[t * J + j], &JnYears, beta, &inc) + betaStarSites[t * J + j] + wSites[t * J + j]);
	      }
	    }
          }
          /********************************
           * Compute A.w
           *******************************/
	  // Copy inverse covariance matrix into tmp_JJ
	  F77_NAME(dcopy)(&nnYears, SigmaEta, &inc, tmp_nnYears, &inc); 
          for (j = 0; j < J; j++) {
	    for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + j] == 1) {
                tmp_nnYears[t * nYearsMax + t] += omega[t * J + j];
	      }
	    } // t
	  } // j

          // Cholesky of A.eta
          F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf on A.eta failed\n");}
	  // Inverse of A.eta
          F77_NAME(dpotri)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri on A.eta failed\n");}
          // A.eta.inv %*% b.eta. Stored in tmp_
          F77_NAME(dsymv)(lower, &nYearsMax, &one, tmp_nnYears, &nYearsMax, 
	  		tmp_nYearsMax, &inc, &zero, tmp_nYearsMax2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
	  if(info != 0){Rf_error("c++ error: dpotrf on A.eta failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(eta, tmp_nYearsMax2, tmp_nnYears, nYearsMax);
	}

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
	for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
              psi[t * J + j] = logitInv(F77_NAME(ddot)(&p, &X[t * J + j], 
					&JnYears, beta, &inc) + wSites[t * J + j] + 
			                betaStarSites[t * J + j] + eta[t], zero, one); 
              if (zDatIndx[t * J + j] == 1) {
                yRep[t * J + j] = rbinom(weights[t * J + j], psi[t * J + j]);
	        like[t * J + j] = dbinom(y[t * J + j], weights[t * J + j], psi[t * J + j], 0);
	      }
          } // j
	} // t


        /********************************************************************
         *Save samples
         *******************************************************************/
	if (rr >= nBurn) {
          thinIndx++; 
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&p, beta, &inc, &REAL(betaSamples_r)[sPost*p], &inc);
            F77_NAME(dcopy)(&JnYears, psi, &inc, &REAL(psiSamples_r)[sPost*JnYears], &inc); 
            F77_NAME(dcopy)(&JpTilde, w, &inc, &REAL(wSamples_r)[sPost*JpTilde], &inc); 
	    F77_NAME(dcopy)(&nThetaAll, theta, &inc, 
			    &REAL(thetaSamples_r)[sPost*nThetaAll], &inc); 
	    F77_NAME(dcopy)(&JnYears, yRep, &inc, &REAL(yRepSamples_r)[sPost*JnYears], &inc); 
            if (pRE > 0) {
              F77_NAME(dcopy)(&pRE, sigmaSqPsi, &inc, 
                  	    &REAL(sigmaSqPsiSamples_r)[sPost*pRE], &inc);
              F77_NAME(dcopy)(&nRE, betaStar, &inc, 
                  	    &REAL(betaStarSamples_r)[sPost*nRE], &inc);
            }
            F77_NAME(dcopy)(&JnYears, like, &inc, 
        		    &REAL(likeSamples_r)[sPost*JnYears], &inc);
	    if (ar1) {
	      F77_NAME(dcopy)(&nYearsMax, eta, &inc, &REAL(etaSamples_r)[sPost*nYearsMax], &inc);
	    }
	    sPost++; 
	    thinIndx = 0; 
	  }
	}

        R_CheckUserInterrupt();
      } // r (end batch)


      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (ll = 0; ll < nThetaAll; ll++) {
        REAL(acceptSamples_r)[s * nThetaAll + ll] = accept[ll]/batchLength; 
        REAL(tuningSamples_r)[s * nThetaAll + ll] = tuning[ll]; 
        if (accept[ll] / batchLength > acceptRate) {
          tuning[ll] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        } else{
            tuning[ll] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          }
        accept[ll] = 0;
      } // ll
      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tCoefficient\tParameter\tAcceptance\tTuning\n");	  
          for (ll = 0; ll < pTilde; ll++) {
            Rprintf("\t%i\t\tphi\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nThetaAll + phiIndx * pTilde + ll], exp(tuning[phiIndx * pTilde + ll]));
	    if (corName == "matern") {
              Rprintf("\t%i\t\tnu\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nThetaAll + nuIndx * pTilde + ll], exp(tuning[nuIndx * pTilde + ll]));
	    }
	    if (sigmaSqIG == 0) {
	      Rprintf("\t%i\t\tsigmaSq\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nThetapTilde + sigmaSqIndx * pTilde + ll], exp(tuning[sigmaSqIndx * pTilde + ll]));
	    }
          } // ll
	  if (ar1) {
	    Rprintf("\t%i\t\trho\t\t%3.1f\t\t%1.5f\n", 1, 100.0*REAL(acceptSamples_r)[s * nThetaAll + rhoIndx], exp(tuning[rhoIndx]));
	  }
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
 
    //make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 9;
    if (pRE > 0) {
      nResultListObjs += 2;
    }
    if (ar1) {
      nResultListObjs += 1;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, yRepSamples_r); 
    SET_VECTOR_ELT(result_r, 2, psiSamples_r);
    SET_VECTOR_ELT(result_r, 3, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 4, wSamples_r); 
    SET_VECTOR_ELT(result_r, 5, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 6, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 7, likeSamples_r); 
    if (pRE > 0) {
      SET_VECTOR_ELT(result_r, 8, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, 9, betaStarSamples_r);
    }
    if (ar1 > 0) {
      if (pRE > 0) {
        tmp_0 = 10;
      } else {
        tmp_0 = 8;
      }
      SET_VECTOR_ELT(result_r, tmp_0, etaSamples_r);
    }    

    // Rf_mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("y.rep.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("like.samples")); 
    if (pRE > 0) {
      SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("beta.star.samples")); 
    }
    if (ar1) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("eta.samples")); 
    }
    
    // Set the names of the output list.  
    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

