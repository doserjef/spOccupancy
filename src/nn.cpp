#define USE_FC_LEN_T
#include <iomanip>
#include <string>
#include <limits>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "util.h"
#include "nn.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

///////////////////////////////////////////////////////////////////
//u index 
///////////////////////////////////////////////////////////////////
SEXP mkUIndx(SEXP n_r, SEXP m_r, SEXP nnIndx_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP nnIndxLU_r, SEXP searchType_r){

  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int *nnIndx = INTEGER(nnIndx_r);
  int *uIndx = INTEGER(uIndx_r);
  int *uIndxLU = INTEGER(uIndxLU_r);
  int *uiIndx = INTEGER(uiIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  int searchType = INTEGER(searchType_r)[0];
  int i, j, k;

  if(searchType == 0){
    mkUIndx0(n, m, nnIndx, uIndx, uIndxLU);
  }else if(searchType == 1){
    mkUIndx1(n, m, nnIndx, uIndx, uIndxLU);
  }else{
    mkUIndx2(n, m, nnIndx, nnIndxLU, uIndx, uIndxLU);
  }
  
  //u lists those locations that have the i-th location as a neighbor
  //then for each of those locations that have i as a neighbor, we need to know the index of i in each of their B vectors (i.e. where does i fall in their neighbor set)
  for(i = 0; i < n; i++){//for each i
    for(j = 0; j < uIndxLU[n+i]; j++){//for each location that has i as a neighbor
  	k = uIndx[uIndxLU[i]+j];//index of a location that has i as a neighbor
  	uiIndx[uIndxLU[i]+j] = which(i, &nnIndx[nnIndxLU[k]], nnIndxLU[n+k]);
    }
  }
  
  return R_NilValue;
}


///////////////////////////////////////////////////////////////////
//Brute force 
///////////////////////////////////////////////////////////////////

SEXP mkNNIndx(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP nnIndx_r, SEXP nnDist_r, SEXP nnIndxLU_r, SEXP nThreads_r){
  
  int i, j, iNNIndx, iNN;
  double d;
  
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  double *coords = REAL(coords_r);
  int *nnIndx = INTEGER(nnIndx_r);
  double *nnDist = REAL(nnDist_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  int nThreads = INTEGER(nThreads_r)[0];
    
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif

  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }
  
#ifdef _OPENMP
#pragma omp parallel for private(j, iNNIndx, iNN, d)
#endif
  for(i = 0; i < n; i++){
    getNNIndx(i, m, iNNIndx, iNN);
    nnIndxLU[i] = iNNIndx;
    nnIndxLU[n+i] = iNN;   
    if(i != 0){  
      for(j = 0; j < i; j++){	
	d = dist2(coords[i], coords[n+i], coords[j], coords[n+j]);	
	if(d < nnDist[iNNIndx+iNN-1]){	  
	  nnDist[iNNIndx+iNN-1] = d;
	  nnIndx[iNNIndx+iNN-1] = j;
	  rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	}	
      }
    }
  }
  
  return R_NilValue;
}


///////////////////////////////////////////////////////////////////
//code book
///////////////////////////////////////////////////////////////////

//Description: using the fast mean-distance-ordered nn search by Ra and Kim 1993
//Input:
//ui = is the index for which we need the m nearest neighbors
//m = number of nearest neighbors
//n = number of observations, i.e., length of u
//sIndx = the NNGP ordering index of length n that is pre-sorted by u
//u = x+y vector of coordinates assumed sorted on input
//rSIndx = vector or pointer to a vector to store the resulting nn sIndx (this is at most length m for ui >= m)
//rNNDist = vector or point to a vector to store the resulting nn Euclidean distance (this is at most length m for ui >= m)  

double dmi(double *x, double *c, int inc){
    return pow(x[0]+x[inc]-c[0]-c[inc], 2);
}

double dei(double *x, double *c, int inc){
  return pow(x[0]-c[0],2)+pow(x[inc]-c[inc],2);
}

void fastNN(int m, int n, double *coords, int ui, double *u, int *sIndx, int *rSIndx, double *rSNNDist){
  
  int i,j;
  bool up, down;
  double dm, de;
  
  //rSNNDist will hold de (i.e., squared Euclidean distance) initially.
  for(i = 0; i < m; i++){
    rSNNDist[i] = std::numeric_limits<double>::infinity();
  }
  
  i = j = ui;
  
  up = down = true;
  
  while(up || down){
    
    if(i == 0){
      down = false;
    }

    if(j == (n-1)){
      up = false;
    }

    if(down){
      
      i--;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[i]], n);
      
      if(dm > 2*rSNNDist[m-1]){
	down = false;
	
      }else{
	de = dei(&coords[sIndx[ui]], &coords[sIndx[i]], n);

	if(de < rSNNDist[m-1] && sIndx[i] < sIndx[ui]){
	  rSNNDist[m-1] = de;
	  rSIndx[m-1] = sIndx[i];
	  rsort_with_index(rSNNDist, rSIndx, m);
	}
	
      }
    }//end down
    
    if(up){
      
      j++;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[j]], n);
      
      if(dm > 2*rSNNDist[m-1]){
	up = false;
	
      }else{
	de = dei(&coords[sIndx[ui]], &coords[sIndx[j]], n);

	if(de < rSNNDist[m-1] && sIndx[j] < sIndx[ui]){
	  rSNNDist[m-1] = de;
	  rSIndx[m-1] = sIndx[j];
	  rsort_with_index(rSNNDist, rSIndx, m);
	}
	
      }
      
    }//end up
    
  }
  
  for(i = 0; i < m; i++){
    rSNNDist[i] = sqrt(rSNNDist[i]);
  }

  return;
}

extern "C" {
  SEXP mkNNIndxCB(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP nnIndx_r, SEXP nnDist_r, SEXP nnIndxLU_r, SEXP nThreads_r){
    
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);
    int *nnIndx = INTEGER(nnIndx_r);
    double *nnDist = REAL(nnDist_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int nThreads = INTEGER(nThreads_r)[0];
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    int i, iNNIndx, iNN;
    
    int *sIndx = new int[n];
    double *u = new double[n];
    
    for(i = 0; i < n; i++){
      sIndx[i] = i;
      u[i] = coords[i]+coords[n+i];
    }
    
    rsort_with_index(u, sIndx, n); 
    
    //make nnIndxLU and fill nnIndx and d
#ifdef _OPENMP
#pragma omp parallel for private(iNNIndx, iNN)
#endif  
    for(i = 0; i < n; i++){ //note this i indexes the u vector
      getNNIndx(sIndx[i], m, iNNIndx, iNN);
      nnIndxLU[sIndx[i]] = iNNIndx;
      nnIndxLU[n+sIndx[i]] = iNN;   
      fastNN(iNN, n, coords, i, u, sIndx, &nnIndx[iNNIndx], &nnDist[iNNIndx]);
    } 
    
    return R_NilValue;
  }
}
