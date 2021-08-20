#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

extern "C" {

  SEXP idist(SEXP coords1_r, SEXP n1_r, SEXP coords2_r, SEXP n2_r, SEXP p_r, SEXP D_r){
    
    int i, j, k;
    double dist = 0.0;

    for(i = 0; i < INTEGER(n1_r)[0]; i++){
      for(j = 0; j < INTEGER(n2_r)[0]; j++){
	dist = 0.0;
	for(k = 0; k < INTEGER(p_r)[0]; k++){
	  dist += pow(REAL(coords1_r)[k*INTEGER(n1_r)[0]+i]-REAL(coords2_r)[k*INTEGER(n2_r)[0]+j],2);
	}
	REAL(D_r)[INTEGER(n1_r)[0]*j+i] = sqrt(dist);
      }
    }

    return(R_NilValue);
  }  
}
