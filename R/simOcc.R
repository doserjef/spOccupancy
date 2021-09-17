simOcc <- function(J.x, J.y, n.rep, beta, alpha, psi.RE = list(), p.RE = list(), 
		   sigma.sq = 2, phi = 3/0.5, sp = FALSE) {

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
  n.beta <- length(beta)
  X <- matrix(1, nrow = J, ncol = n.beta)
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      X[, i] <- rnorm(J)
    } # i
  }

  # Form detection covariate (if any) -------------------------------------
  n.alpha <- length(alpha)
  X.p <- array(NA, dim = c(J, max(n.rep), n.alpha))
  X.p[, , 1] <- 1
  if (n.alpha > 1) {
    for (i in 2:n.alpha) {
      for (j in 1:J) {
        X.p[j, 1:n.rep[j], i] <- rnorm(n.rep[j])
      } # j
    } # i
  }

  # Simulate spatial random effect ----------------------------------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- expand.grid(s.x, s.y)
  # Distance matrix
  D <- as.matrix(dist(coords))
  # Exponential correlation function
  R <- exp(-phi * D)
  # Random spatial process
  w <- rmvn(1, rep(0, J), sigma.sq * R)

  # Random effects --------------------------------------------------------
  if (length(psi.RE) > 0) {
    p.occ.re <- length(psi.RE$levels)
    sigma.sq.psi <- rep(NA, p.occ.re)
    n.occ.re.long <- psi.RE$levels
    n.occ.re <- sum(n.occ.re.long)
    beta.star.indx <- rep(1:p.occ.re, n.occ.re.long)
    beta.star <- rep(0, n.occ.re)
    X.re <- matrix(NA, J, p.occ.re)
    for (i in 1:p.occ.re) {
      X.re[, i] <- sample(1:psi.RE$levels[i], J, replace = TRUE)         
      beta.star[which(beta.star.indx == i)] <- rnorm(psi.RE$levels[i], 0, sqrt(psi.RE$sigma.sq.psi[i]))
    }
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1])
      }
    } 
    beta.star.sites <- apply(X.re, 1, function(a) sum(beta.star[a]))
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
    alpha.star <- rep(0, n.det.re)
    X.p.re <- array(NA, dim = c(J, max(n.rep), p.det.re))
    for (i in 1:p.det.re) {
      X.p.re[, , i] <- matrix(sample(1:p.RE$levels[i], J * max(n.rep), replace = TRUE), 
		              J, max(n.rep))	      
      alpha.star[which(alpha.star.indx == i)] <- rnorm(p.RE$levels[i], 0, sqrt(p.RE$sigma.sq.p[i]))
    }
    for (j in 1:J) {
      X.p.re[j, -(1:n.rep[j]), ] <- NA
    }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        X.p.re[, , j] <- X.p.re[, , j] + max(X.p.re[, , j - 1]) 
      }
    }
    alpha.star.sites <- apply(X.p.re, c(1, 2), function(a) sum(alpha.star[a]))
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }


  # Latent Occupancy Process ----------------------------------------------
  if (sp) {
    if (length(psi.RE) > 0) {
      psi <- logit.inv(X %*% as.matrix(beta) + w + beta.star.sites)
    } else {
      psi <- logit.inv(X %*% as.matrix(beta) + w)
    }
  } else {
    if (length(psi.RE) > 0) {
      psi <- logit.inv(X %*% as.matrix(beta) + beta.star.sites)
    } else {
      psi <- logit.inv(X %*% as.matrix(beta))
    }
  }
  z <- rbinom(J, 1, psi)

  # Data Formation --------------------------------------------------------
  p <- matrix(NA, nrow = J, ncol = max(n.rep))
  y <- matrix(NA, nrow = J, ncol = max(n.rep))
  for (j in 1:J) {
    if (length(p.RE) > 0) {
      p[j, 1:n.rep[j]] <- logit.inv(X.p[j, 1:n.rep[j], ] %*% as.matrix(alpha) + 
				    alpha.star.sites[j, 1:n.rep[j]])
    } else {
      p[j, 1:n.rep[j]] <- logit.inv(X.p[j, 1:n.rep[j], ] %*% as.matrix(alpha))
    }
    y[j, 1:n.rep[j]] <- rbinom(n.rep[j], 1, p[j, 1:n.rep[j]] * z[j])
  } # j

  return(
    list(X = X, X.p = X.p, coords = coords, D = D,
         w = w, psi = psi, z = z, p = p, y = y, 
	 X.p.re = X.p.re, X.re = X.re, 
	 alpha.star = alpha.star, beta.star = beta.star)
  )
}
