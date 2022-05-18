#include <string>


  void zeros(double *a, int n);

  void ones(double *a, int n); 
  
  void mvrnorm(double *des, double *mu, double *cholCov, int dim);
  
  double logit(double theta, double a, double b);
  
  double logitInv(double z, double a, double b);
  
  double dist2(double &a1, double &a2, double &b1, double &b2);
  
  //Description: given a location's index i and number of neighbors m this function provides the index to i and number of neighbors in nnIndx
  void getNNIndx(int i, int m, int &iNNIndx, int &iNN);
  
  //Description: given the nnIndex this function fills uIndx for identifying those locations that have the i-th location as a neighbor.
  //Input:
  //n = number of locations
  //m = number of nearest neighbors
  //nnIndx = set of nearest neighbors for all n locations
  //Output:
  //uIndx = holds the indexes for locations that have each location as a neighbor
  //uIndxLU = nx2 look-up matrix with row values correspond to each location's index in uIndx and number of neighbors (columns 1 and 2, respectively)
  //Note: uIndx must be of length (1+m)/2*m+(n-m-1)*m on input. uINdxLU must also be allocated on input.
  void mkUIndx0(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU);
  void mkUIndx1(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU);
  void mkUIndx2(int n, int m, int* nnIndx, int *nnIndxLU, int* uIndx, int* uIndxLU);
  void crs_csc(int n, int *i_A, int *j_A, int *i_B, int *j_B);
  
  std::string getCorName(int i);
  
  int which(int a, int *b, int n);
  
  double Q(double *B, double *F, double *u, double *v, int n, int *nnIndx, int *nnIndxLU);
  
  void printMtrx(double *m, int nRow, int nCol);
  
  void printMtrxInt(int *m, int nRow, int nCol);
  void spCovLT(double *D, int n, double *theta, std::string &covModel, double *C); 
  void spCov(double *D, int n, double *theta, std::string &covModel, double *C); 
  void spCorLT(double *D, int n, double *theta, std::string &covModel, double *R); 

  double rigamma(double a, double b); 
  void fillUTri(double *v, int m);
  double spCor(double &D, double &phi, double &nu, int &covModel, double *bk); 
  double spCor(double D, double *theta, std::string &covModel);
  void clearUT(double *m, int n);

