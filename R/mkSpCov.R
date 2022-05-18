# This function comes directly from the spBayes R package. 
# Author: Andrew O. Finley
mkSpCov <- function(coords, K, Psi, theta, cov.model){

  if(missing(coords)){stop("error: coords must be specified")}
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  n <- nrow(coords)
  if(ncol(coords) != 2 ){
    stop("error: coords have more than two columns")
  }
  
  if(!is.matrix(K)){stop("error: K must be a m-by-m cross-covaraince matrix")}
  if(!is.matrix(Psi)){stop("error: Psi must be a m-by-m cross-covaraince matrix")}

  m <- nrow(K)

  if(!all(c(dim(K), dim(Psi)) == m)){stop("error: K and Psi be m-by-m cross-covaraince matrices")}
  
  storage.mode(coords) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(Psi) <- "double"
  storage.mode(K) <- "double"
  storage.mode(theta) <- "double"
  
  .Call("mkSpCov", coords, n, m, Psi, K, theta, cov.model)
}

