#define USE_FC_LEN_T
#include <string>
#include <limits>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

  void zeros(double *a, int n){
    for(int i = 0; i < n; i++)
      a[i] = 0.0;
  }

  void ones(double *a, int n){
    for(int i = 0; i < n; i++)
      a[i] = 1.0;
  }
  
  void mvrnorm(double *des, double *mu, double *cholCov, int dim){
    
    int i;
    int inc = 1;
    double one = 1.0;
    
    for(i = 0; i < dim; i++){
      des[i] = rnorm(0, 1);
    }
   
    F77_NAME(dtrmv)("L", "N", "N", &dim, cholCov, &dim, des, &inc FCONE FCONE FCONE);
    F77_NAME(daxpy)(&dim, &one, mu, &inc, des, &inc);
  }
  
  double logit(double theta, double a, double b){
    return log((theta-a)/(b-theta));
  }
  
  double logitInv(double z, double a, double b){
    return b-(b-a)/(1+exp(z));
  }
  
  double dist2(double &a1, double &a2, double &b1, double &b2){
    return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
  }
  
  void getNNIndx(int i, int m, int &iNNIndx, int &iNN){
    
    if(i == 0){
      iNNIndx = 0;//this should never be accessed
      iNN = 0;
      return;
    }else if(i < m){
      iNNIndx = static_cast<int>(static_cast<double>(i)/2*(i-1));
      iNN = i;
      return;
    }else{
      iNNIndx = static_cast<int>(static_cast<double>(m)/2*(m-1)+(i-m)*m);
      iNN = m;
      return;
    } 
  }
  
  void mkUIndx0(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU){ 
    
    int iNNIndx, iNN, i, j, k, l, h;
    
    for(i = 0, l = 0; i < n; i++){    
      uIndxLU[i] = l; 
      for(j = 0, h = 0; j < n; j++){   
        getNNIndx(j, m, iNNIndx, iNN);  
        for(k = 0; k < iNN; k++){      	
  	if(nnIndx[iNNIndx+k] == i){
  	  uIndx[l+h] = j;
  	  h++;
  	}    
        }
      }
      l += h;
      uIndxLU[n+i] = h;
      R_CheckUserInterrupt();
    }
  }
  
  void mkUIndx1(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU){ 
    
    int iNNIndx, iNN, i, j, k, l, h;
    
    for(i = 0, l = 0; i < n; i++){    
      uIndxLU[i] = l; 
      for(j = n-1, h = 0; j > i; j--){   
        getNNIndx(j, m, iNNIndx, iNN);  
        for(k = 0; k < iNN; k++){      	
  	if(nnIndx[iNNIndx+k] == i){
  	  uIndx[l+h] = j;
  	  h++;
  	}    
        }
      }
      l += h;
      uIndxLU[n+i] = h;
      R_CheckUserInterrupt();
    }
  }
  
  
  void mkUIndx2(int n, int m, int* nnIndx, int *nnIndxLU, int* uIndx, int* uIndxLU){ 
  
    int i, k;
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    
    //int *j_A = new int[nIndx]; is nnIndx
    int *i_nnIndx = new int[n+1];
    //int *j_A_csc = new int[nIndx];//uIndx
    int *i_A_csc = new int[n+1];
  
    for(i = 0, k = 0; i < n; i++){
      if(nnIndxLU[n+i] == 0){//excludes rows with no elements, i.e., the first row because it is zero by design A[0,0] = 0
        i_nnIndx[0] = 0;
      }else{
        i_nnIndx[k] = i_nnIndx[k-1]+nnIndxLU[n+i-1];
      }
      k++;
    }
    i_nnIndx[n] = i_nnIndx[0]+nIndx;
      
    crs_csc(n, i_nnIndx, nnIndx, i_A_csc, uIndx);
    
    for(i = 0; i < n; i++){
      uIndxLU[i] = i_A_csc[i];
      uIndxLU[i+n] = i_A_csc[i+1]-i_A_csc[i];
    }
    
    delete[] i_nnIndx;
    delete[] i_A_csc;
    
  }
  
  
  void crs_csc(int n, int *i_A, int *j_A, int *i_B, int *j_B){
  
    int i, j, col, cumsum, temp, row, dest, last;
    
    int nnz = i_A[n];
  
    for(i = 0; i < n; i++){
      i_B[i] = 0;
    }
    
    for(i = 0; i < nnz; i++){            
      i_B[j_A[i]]++;
    }
    
    //cumsum the nnz per column to get i_B[]
    for(col = 0, cumsum = 0; col < n; col++){     
      temp  = i_B[col];
      i_B[col] = cumsum;
      cumsum += temp;
    }
    i_B[n] = nnz; 
    
    for(row = 0; row < n; row++){
      for(j = i_A[row]; j < i_A[row+1]; j++){
        col  = j_A[j];
        dest = i_B[col];
        
        j_B[dest] = row;
        i_B[col]++;
      }
    }  
    
    for(col = 0, last = 0; col <= n; col++){
      temp  = i_B[col];
      i_B[col] = last;
      last = temp;
    }
  } 
  
  
  
  
  
  
  std::string getCorName(int i){
  
    if(i == 0){
      return "exponential";
    }else if(i == 1){
      return "spherical";
    }else if(i == 2){
      return "matern";
    }else if(i == 3){
      return "gaussian";
    }else{
      error("c++ error: cov.model is not correctly specified");
    }
    
  }
  
  //which index of b equals a, where b is of length n
  int which(int a, int *b, int n){
    int i;
    for(i = 0; i < n; i++){
      if(a == b[i]){
        return(i);
      }
    }
  
    error("c++ error: which failed");
    return -9999;
  }
  
  //Description: computes the quadratic term.
  double Q(double *B, double *F, double *u, double *v, int n, int *nnIndx, int *nnIndxLU){
    
    double a, b, q = 0;
    int i, j;
    
  #ifdef _OPENMP
  #pragma omp parallel for private(a, b, j) reduction(+:q)
  #endif  
    for(i = 0; i < n; i++){
      a = 0;
      b = 0;
      for(j = 0; j < nnIndxLU[n+i]; j++){
        a += B[nnIndxLU[i]+j]*u[nnIndx[nnIndxLU[i]+j]];
        b += B[nnIndxLU[i]+j]*v[nnIndx[nnIndxLU[i]+j]];
      }
      q += (u[i] - a)*(v[i] - b)/F[i];
    }
    
    return(q);
  }
  
  
  void printMtrx(double *m, int nRow, int nCol){
  
    int i, j;
  
    for(i = 0; i < nRow; i++){
      Rprintf("\t");
      for(j = 0; j < nCol; j++){
        Rprintf("%.10f\t", m[j*nRow+i]);
      }
      Rprintf("\n");
    }
  }
  
  
  void printMtrxInt(int *m, int nRow, int nCol){
  
    int i, j;
  
    for(i = 0; i < nRow; i++){
      Rprintf("\t");
      for(j = 0; j < nCol; j++){
        Rprintf("%i\t", m[j*nRow+i]);
      }
      Rprintf("\n");
    }
  
  }

void spCovLT(double *D, int n, double *theta, std::string &covModel, double *C){
  int i,j;
  
  if(covModel == "exponential"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	C[i*n+j] = theta[0]*exp(-1.0*theta[1]*D[i*n+j]);
      }
    }
    
  }else if(covModel == "spherical"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	if(D[i*n+j] > 0 && D[i*n+j] <= 1.0/theta[1]){
	  C[i*n+j] = theta[0]*(1.0 - 1.5*theta[1]*D[i*n+j] + 0.5*pow(theta[1]*D[i*n+j],3));
	}else if(D[i*n+j] >= 1.0/theta[1]){
	  C[i*n+j] = 0.0;
	}else{
	  C[i*n+j] = theta[0];
	}
      }
    }
    
  }else if(covModel == "gaussian"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	C[i*n+j] = theta[0]*exp(-1.0*(pow(theta[1]*D[i*n+j],2)));
      }
    }

  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	if(D[i*n+j]*theta[1] > 0.0){
	  C[i*n+j] = theta[0]*pow(D[i*n+j]*theta[1], theta[2])/(pow(2, theta[2]-1)*gammafn(theta[2]))*bessel_k(D[i*n+j]*theta[1], theta[2], 1.0);
	}else{
	  C[i*n+j] = theta[0];
	}
      }
    }
    
 }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

void spCorLT(double *D, int n, double *theta, std::string &covModel, double *R){
  int i,j;
  
  if(covModel == "exponential"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	R[i*n+j] = exp(-1.0*theta[1]*D[i*n+j]);
      }
    }

  }else if(covModel == "spherical"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	if(D[i*n+j] > 0 && D[i*n+j] <= 1.0/theta[1]){
	  R[i*n+j] = (1.0 - 1.5*theta[1]*D[i*n+j] + 0.5*pow(theta[1]*D[i*n+j],3));
	}else if(D[i*n+j] >= 1.0/theta[1]){
	  R[i*n+j] = 0.0;
	}else{
	  R[i*n+j] = 1.0;
	}
      }
    }
    
  }else if(covModel == "gaussian"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	R[i*n+j] = exp(-1.0*(pow(theta[1]*D[i*n+j],2)));
      }
    }

  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	if(D[i*n+j]*theta[1] > 0.0){
	  R[i*n+j] = pow(D[i*n+j]*theta[1], theta[2])/(pow(2, theta[2]-1)*gammafn(theta[2]))*bessel_k(D[i*n+j]*theta[1], theta[2], 1.0);
	}else{
	  R[i*n+j] = 1.0;
	}
      }
    }
    
 }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

double rigamma(double a, double b) {
  return 1.0 / rgamma(a, 1.0 / b);
}

void spCov(double *D, int n, double *theta, std::string &covModel, double *C){
  int i;
  
  if(covModel == "exponential"){
    
    for(i = 0; i < n; i++){
      C[i] = theta[0]*exp(-1.0*theta[1]*D[i]);
    }
    
  }else if(covModel == "spherical"){
    
    for(i = 0; i < n; i++){
      if(D[i] > 0 && D[i] <= 1.0/theta[1]){
	C[i] = theta[0]*(1.0 - 1.5*theta[1]*D[i] + 0.5*pow(theta[1]*D[i],3));
      }else if(D[i] >= 1.0/theta[1]){
	C[i] = 0.0;
      }else{
	C[i] = theta[0];
      }
    }
    
  }else if(covModel == "gaussian"){
    
    for(i = 0; i < n; i++){
      C[i] = theta[0]*exp(-1.0*(pow(theta[1]*D[i],2)));
    }
    
  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    for(i = 0; i < n; i++){
      if(D[i]*theta[1] > 0.0){
	C[i] = theta[0]*pow(D[i]*theta[1], theta[2])/(pow(2, theta[2]-1)*gammafn(theta[2]))*bessel_k(D[i]*theta[1], theta[2], 1.0);
      }else{
	C[i] = theta[0];
      }
    }
    
 }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

void fillUTri(double *v, int m){
  int i, j;
  for(i = 0; i < m; i++){
    for(j = i; j < m; j++){
      v[j*m+i] = v[i*m+j];
    }
  }
}

double spCor(double D, double *theta, std::string &covModel){

  if(covModel == "exponential"){

    return exp(-1.0*theta[0]*D);

  }else if(covModel == "spherical"){

    if(D > 0 && D <= 1.0/theta[0]){
      return 1.0 - 1.5*theta[0]*D + 0.5*pow(theta[0]*D,3);
    }else if(D >= 1.0/theta[0]){
      return 0.0;
    }else{
      return 1.0;
    }

  }else if(covModel == "gaussian"){

    return exp(-1.0*(pow(theta[0]*D,2)));

  }else if(covModel == "matern"){

    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)

    if(D*theta[0] > 0.0){
      return pow(D*theta[0], theta[1])/(pow(2, theta[1]-1)*gammafn(theta[1]))*bessel_k(D*theta[0], theta[1], 1.0);
    }else{
      return 1.0;
    }

  }else{
    error("c++ error: cov.model is not correctly specified");
    return 0;
  }
}

double spCor(double &D, double &phi, double &nu, int &covModel, double *bk){

  //0 exponential
  //1 spherical
  //2 matern
  //3 gaussian

  if(covModel == 0){//exponential

    return exp(-phi*D);

  }else if(covModel == 1){//spherical

    if(D > 0 && D <= 1.0/phi){
      return 1.0 - 1.5*phi*D + 0.5*pow(phi*D,3);
    }else if(D >= 1.0/phi){
      return 0.0;
    }else{
      return 1.0;
    }
  }else if(covModel == 2){//matern

    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)

    if(D*phi > 0.0){
      return pow(D*phi, nu)/(pow(2, nu-1)*gammafn(nu))*bessel_k_ex(D*phi, nu, 1.0, bk);//thread safe bessel
    }else{
      return 1.0;
    }
  }else if(covModel == 3){//gaussian

    return exp(-1.0*(pow(phi*D,2)));

  }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

void clearUT(double *m, int n){
  for(int i = 1; i < n; i++){
    for(int j = 0; j < i; j++){
      m[i*n+j] = 0;
    }
  }
}

