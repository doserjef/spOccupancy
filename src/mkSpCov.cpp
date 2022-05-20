#define USE_FC_LEN_T
#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#ifndef FCONE
# define FCONE
#endif

extern"C" {
  
  SEXP mkSpCov(SEXP coords_r, SEXP n_r, SEXP m_r, SEXP Psi_r, SEXP V_r, SEXP theta_r, SEXP covModel_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, ii, jj, info;
    char const *lower = "L";
    const int incOne = 1;
    
    double *coords = REAL(coords_r);
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *Psi = REAL(Psi_r);
    double *V = REAL(V_r);
    double *theta = REAL(theta_r);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    
    double *gamma = (double *) R_alloc(2, sizeof(double));
    
    int mm = m*m;
    int nn = n*n;
    int nm = n*m;
    
    double *D = (double *) R_alloc(nn, sizeof(double));
    for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
	D[n*j+i] = sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2));
      }
    }
    
    SEXP C;
    PROTECT(C = allocMatrix(REALSXP, nm, nm)); 
    
    //Get A
    double *A = (double *) R_alloc(mm, sizeof(double));
    F77_NAME(dcopy)(&mm, V, &incOne, A, &incOne);
    F77_NAME(dpotrf)(lower, &m, A, &m, &info FCONE); if(info != 0){error("Cholesky failed\n");}
    clearUT(A, m);    
    
    for(jj = 0; jj < n; jj++){
      for(ii = 0; ii < n; ii++){	
	for(k = 0; k < m; k++){
	  for(l = 0; l < m; l++){
	    REAL(C)[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
	    for(h = 0; h < m; h++){
	      gamma[0] = theta[h];
	      if(covModel == "matern"){
		gamma[1] = theta[m+h];
	      }
	      REAL(C)[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(D[jj*n+ii], gamma, covModel);
	    }
	  }
	}
      }
    }
    
    for(i = 0; i < n; i++){
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  REAL(C)[(i*m+l)*nm+(i*m+k)] += Psi[l*m+k];
	}
      }
    }
    
    UNPROTECT(1);
    
    return(C);
    
  } 
}

