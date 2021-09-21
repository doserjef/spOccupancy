simMsOcc <- function(J.x, J.y, n.rep, N, beta, alpha, psi.RE = list(), 
		     p.RE = list(), sigma.sq = rep(2, N), 
		     phi = rep(3/0.5, N), sp = FALSE) {

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
  X <- matrix(1, nrow = J, ncol = p.occ) 
  if (p.occ > 1) {
    for (i in 2:p.occ) {
      X[, i] <- rnorm(J)
    } # i
  }

  # Form detection covariate (if any) -------------------------------------
  p.det <- ncol(alpha)
  X.p <- array(NA, dim = c(J, max(n.rep), p.det))
  X.p[, , 1] <- 1
  if (p.det > 1) {
    for (i in 2:p.det) {
      for (j in 1:J) {
        X.p[j, 1:n.rep[j], i] <- rnorm(n.rep[j])
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

  # Random effects --------------------------------------------------------
  if (length(psi.RE) > 0) {
    p.occ.re <- length(psi.RE$levels)
    sigma.sq.psi <- rep(NA, p.occ.re)
    n.occ.re.long <- psi.RE$levels
    n.occ.re <- sum(n.occ.re.long)
    beta.star.indx <- rep(1:p.occ.re, n.occ.re.long)
    beta.star <- matrix(0, N, n.occ.re)
    X.re <- matrix(NA, J, p.occ.re)
    for (l in 1:p.occ.re) {
      X.re[, l] <- sample(1:psi.RE$levels[l], J, replace = TRUE)         
      for (i in 1:N) {
        beta.star[i, which(beta.star.indx == l)] <- rnorm(psi.RE$levels[l], 0, 
							  sqrt(psi.RE$sigma.sq.psi[l]))
      }
    }
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1])
      }
    } 
    beta.star.sites <- matrix(NA, N, J)
    for (i in 1:N) {
      beta.star.sites[i, ] <- apply(X.re, 1, function(a) sum(beta.star[i, a]))
    }
  } else {
    X.re <- NA
    beta.star <- NA
  }
  if (length(p.RE) > 0) {
    p.det.re <- length(p.RE$levels)
    sigma.sq.p <- rep(NA, p.det.re)
    n.det.re.long <- p.RE$levels
    n.det.re <- sum(n.det.re.long)
    alpha.star.indx <- rep(1:p.det.re, n.det.re.long)
    alpha.star <- matrix(0, N, n.det.re)
    X.p.re <- array(NA, dim = c(J, max(n.rep), p.det.re))
    for (l in 1:p.det.re) {
      X.p.re[, , l] <- matrix(sample(1:p.RE$levels[l], J * max(n.rep), replace = TRUE), 
		              J, max(n.rep))	      
      for (i in 1:N) {
        alpha.star[i, which(alpha.star.indx == l)] <- rnorm(p.RE$levels[l], 0, sqrt(p.RE$sigma.sq.p[l]))
      }
    }
    for (j in 1:J) {
      X.p.re[j, -(1:n.rep[j]), ] <- NA
    }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        X.p.re[, , j] <- X.p.re[, , j] + max(X.p.re[, , j - 1]) 
      }
    }
    alpha.star.sites <- array(NA, c(N, J, max(n.rep)))
    for (i in 1:N) {
      alpha.star.sites[i, , ] <- apply(X.p.re, c(1, 2), function(a) sum(alpha.star[i, a]))
    }
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }

  # Latent Occupancy Process ----------------------------------------------
  psi <- matrix(NA, nrow = N, ncol = J)
  z <- matrix(NA, nrow = N, ncol = J)
  for (i in 1:N) {
    if (sp) {
      if (length(psi.RE) > 0) {
        psi[i, ] <- logit.inv(X %*% as.matrix(beta[i, ]) + w[i, ] + beta.star.sites[i, ])
      } else {
        psi[i, ] <- logit.inv(X %*% as.matrix(beta[i, ]) + w[i, ])
      }
    } else {
      if (length(psi.RE) > 0) {
        psi[i, ] <- logit.inv(X %*% as.matrix(beta[i, ]) + beta.star.sites[i, ])
      } else {
        psi[i, ] <- logit.inv(X %*% as.matrix(beta[i, ]))
      }
    }
    z[i, ] <- rbinom(J, 1, psi[i, ])
  }

  # Data Formation --------------------------------------------------------
  p <- array(NA, dim = c(N, J, max(n.rep)))
  y <- array(NA, dim = c(N, J, max(n.rep)))
  for (i in 1:N) {
    for (j in 1:J) {
      if (length(p.RE) > 0) {
        p[i, j, 1:n.rep[j]] <- logit.inv(X.p[j, 1:n.rep[j], ] %*% as.matrix(alpha[i, ]) + alpha.star.sites[i, j, 1:n.rep[j]])
      } else {
        p[i, j, 1:n.rep[j]] <- logit.inv(X.p[j, 1:n.rep[j], ] %*% as.matrix(alpha[i, ]))
      }
 
        y[i, j, 1:n.rep[j]] <- rbinom(n.rep[j], 1, p[i, j, 1:n.rep[j]] * z[i, j]) 
    } # j
  } # i
  return(
    list(X = X, X.p = X.p, coords = coords, D = D,
	 w = w, psi = psi, z = z, p = p, y = y, X.p.re = X.p.re, 
	 X.re = X.re, alpha.star = alpha.star, beta.star = beta.star)
  )
}
