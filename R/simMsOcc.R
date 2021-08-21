simMsOcc <- function(J.x, J.y, K, N, beta, alpha, sigma.sq = 2, phi = 3/0.5, 
		          sp = FALSE) {

  # Subroutines -----------------------------------------------------------
  # MVN 
  rmvn <- function(n, mu=0, V = matrix(1)) {
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }

  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}
  
  # Form occupancy covariates (if any) ------------------------------------
  J <- J.x * J.y
  p.occ <- ncol(beta)
  X.psi <- matrix(1, nrow = J, ncol = p.occ) 
  if (p.occ > 1) {
    for (i in 2:p.occ) {
      X.psi[, i] <- rnorm(J)
    } # i
  }

  # Form detection covariate (if any) -------------------------------------
  p.det <- ncol(alpha)
  X.p <- array(NA, dim = c(J, max(K), p.det))
  X.p[, , 1] <- 1
  if (p.det > 1) {
    for (i in 2:p.det) {
      for (j in 1:J) {
        X.p[j, 1:K[j], i] <- rnorm(K[j])
      } # j
    } # i
  }

  # Simulate spatial random effect for each species -----------------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- expand.grid(s.x, s.y)
  # Distance matrix
  D <- as.matrix(dist(coords))
  # Spatial random effects for each species (assuming exponential decay)
  w <- matrix(NA, nrow = N, ncol = J)
  for (i in 1:N) {
    w[i, ] <- rmvn(1, rep(0, J), sigma.sq[i] * exp(-phi[i] * D))
  }

  # Latent Occupancy Process ----------------------------------------------
  psi <- matrix(NA, nrow = N, ncol = J)
  z <- matrix(NA, nrow = N, ncol = J)
  for (i in 1:N) {
    if (sp) {
      psi[i, ] <- logit.inv(X.psi %*% as.matrix(beta[i, ]) + w[i, ])
    } else {
        psi[i, ] <- logit.inv(X.psi %*% as.matrix(beta[i, ]))
      }
    z[i, ] <- rbinom(J, 1, psi[i, ])
  }

  # Data Formation --------------------------------------------------------
  p <- array(NA, dim = c(N, J, max(K)))
  y <- array(NA, dim = c(N, J, max(K)))
  for (i in 1:N) {
    for (j in 1:J) {
        p[i, j, 1:K[j]] <- logit.inv(X.p[j, 1:K[j], ] %*% as.matrix(alpha[i, ]))
        y[i, j, 1:K[j]] <- rbinom(K[j], 1, p[i, j, 1:K[j]] * z[i, j]) 
    } # j
  } # i
  return(
    list(X.psi = X.psi, X.p = X.p, coords = coords, D = D,
	 w = w, psi = psi, z = z, p = p, y = y, s.x = s.x, s.y = s.y)
  )
}
