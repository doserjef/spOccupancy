simMsOcc <- function(J.x, J.y, n.rep, n.rep.max, N, beta, alpha, psi.RE = list(), 
                     p.RE = list(), sp = FALSE, svc.cols = 1, cov.model, 
                     sigma.sq, phi, nu, factor.model = FALSE, n.factors, 
                     range.probs, shared.spatial = FALSE, grid, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Subroutines -----------------------------------------------------------
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  # Check function inputs -------------------------------------------------
  # J.x -------------------------------
  if (missing(J.x)) {
    stop("error: J.x must be specified")
  }
  if (length(J.x) != 1) {
    stop("error: J.x must be a single numeric value.")
  }
  # J.y -------------------------------
  if (missing(J.y)) {
    stop("error: J.y must be specified")
  }
  if (length(J.y) != 1) {
    stop("error: J.y must be a single numeric value.")
  }
  J <- J.x * J.y
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (length(n.rep) != J) {
    stop(paste("error: n.rep must be a vector of length ", J, sep = ''))
  }
  if (missing(n.rep.max)) {
    n.rep.max <- max(n.rep)
  }
  # N ---------------------------------
  if (missing(N)) {
    stop("error: N must be specified")
  }
  if (length(N) != 1) {
    stop("error: N must be a single numeric value.")
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified")
  }
  if (!is.matrix(beta)) {
    stop(paste("error: beta must be a numeric matrix with ", N, " rows", sep = ''))
  }
  if (nrow(beta) != N) {
    stop(paste("error: beta must be a numeric matrix with ", N, " rows", sep = ''))
  }
  # alpha -----------------------------
  if (missing(alpha)) {
    stop("error: alpha must be specified.")
  }
  if (!is.matrix(alpha)) {
    stop(paste("error: alpha must be a numeric matrix with ", N, " rows", sep = ''))
  }
  if (nrow(alpha) != N) {
    stop(paste("error: alpha must be a numeric matrix with ", N, " rows", sep = ''))
  }
  # Check spatial stuff ---------------
  if (sp & !factor.model) {
    N.p.svc <- N * length(svc.cols)
    # sigma.sq --------------------------
    if (missing(sigma.sq)) {
      stop("error: sigma.sq must be specified when sp = TRUE")
    }
    if (length(sigma.sq) != N.p.svc) {
      stop(paste("error: sigma.sq must be a vector of length ", N.p.svc, sep = ''))
    }
    # phi -------------------------------
    if(missing(phi)) {
      stop("error: phi must be specified when sp = TRUE")
    }
    if (length(phi) != N.p.svc) {
      stop(paste("error: phi must be a vector of length ", N.p.svc, sep = ''))
    }
  }
  if (sp) {
    # Covariance model ----------------
    if(missing(cov.model)) {
      stop("error: cov.model must be specified when sp = TRUE")
    }
    cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
    if(! cov.model %in% cov.model.names){
      stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ", 
           paste(cov.model.names, collapse=", ", sep="") ,".")
    }
    if (cov.model == 'matern' & missing(nu)) {
      stop("error: nu must be specified when cov.model = 'matern'")
    }
  }
  p.svc <- length(svc.cols)
  if (factor.model) {
    # n.factors -----------------------
    if (missing(n.factors)) {
      stop("error: n.factors must be specified when factor.model = TRUE")
    }
    q.p.svc <- n.factors * length(svc.cols)
    if (sp) {
      if (!missing(sigma.sq) & !shared.spatial) {
        message("sigma.sq is specified but will be set to 1 for spatial latent factor model")
      }
      if (missing(sigma.sq) & shared.spatial) {
        stop("sigma.sq must be specified when shared.spatial = TRUE")
      }
      if(missing(phi)) {
        stop("error: phi must be specified when sp = TRUE")
      }
      if (length(phi) != q.p.svc) {
        stop(paste("error: phi must be a vector of length ", q.p.svc, sep = ''))
      }
    }
    if (!sp & length(svc.cols) > 1) {
      stop("error: length(svc.cols) > 1 when sp = FALSE. Set sp = TRUE to simulate data with spatially-varying coefficients")
    }
  }
  # psi.RE ----------------------------
  names(psi.RE) <- tolower(names(psi.RE))
  if (!is.list(psi.RE)) {
    stop("error: if specified, psi.RE must be a list with tags 'levels' and 'sigma.sq.psi'")
  }
  if (length(names(psi.RE)) > 0) {
    if (!'sigma.sq.psi' %in% names(psi.RE)) {
      stop("error: sigma.sq.psi must be a tag in psi.RE with values for the occurrence random effect variances")
    }
    if (!'levels' %in% names(psi.RE)) {
      stop("error: levels must be a tag in psi.RE with the number of random effect levels for each occurrence random intercept.")
    }
  }
  # p.RE ----------------------------
  names(p.RE) <- tolower(names(p.RE))
  if (!is.list(p.RE)) {
    stop("error: if specified, p.RE must be a list with tags 'levels' and 'sigma.sq.p'")
  }
  if (length(names(p.RE)) > 0) {
    if (!'sigma.sq.p' %in% names(p.RE)) {
      stop("error: sigma.sq.p must be a tag in p.RE with values for the detection random effect variances")
    }
    if (!'levels' %in% names(p.RE)) {
      stop("error: levels must be a tag in p.RE with the number of random effect levels for each detection random intercept.")
    }
  }

  # range.probs -----------------------
  if (!missing(range.probs)) {
    if (length(range.probs) != N) {
      stop(paste("error: range.probs must be a numeric vector of length ", N, sep = ''))
    }
  } else {
    range.probs <- rep(1, N)
  }

  # Grid for spatial REs that doesn't match the sites ---------------------
  if (!missing(grid) & sp) {
    if (!is.atomic(grid)) {
      stop("grid must be a vector")
    }
    if (length(grid) != J) {
      stop(paste0("grid must be of length ", J))
    }
  } else {
    grid <- 1:J
  }

  # Subroutines -----------------------------------------------------------
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords.full <- as.matrix(expand.grid(s.x, s.y))
  coords <- cbind(tapply(coords.full[, 1], grid, mean),
                  tapply(coords.full[, 2], grid, mean))
  
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
  X.p <- array(NA, dim = c(J, n.rep.max, p.det))
  # Get index of surveyed replicates for each site. 
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- sample(1:n.rep.max, n.rep[j], replace = FALSE)
  }
  X.p[, , 1] <- 1
  if (p.det > 1) {
    for (i in 2:p.det) {
      for (j in 1:J) {
        X.p[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
      } # j
    } # i
  }

  # Simulate latent (spatial) random effect for each species --------------
  # Matrix of spatial locations
  w.star <- vector(mode = "list", length = p.svc)
  w <- vector(mode = "list", length = p.svc)
  lambda <- vector(mode = "list", length = p.svc)
  J.w <- nrow(coords)
  # Form spatial process for each spatially-varying covariate
  for (i in 1:p.svc) {
    w.star[[i]] <- matrix(0, nrow = N, ncol = J.w)
    w[[i]] <- rep(0, J.w)
    if (shared.spatial) {
      if (cov.model == 'matern') {
        theta <- cbind(phi, nu)
      } else {
        theta <- as.matrix(phi)
      }
      Sigma.full <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
      w[[i]] <- rmvn(1, rep(0, J.w), Sigma.full)
      for (ll in 1:N) {
        w.star[[i]][ll, ] <- w[[i]]
      }
    } else {
      if (factor.model) {
        lambda[[i]] <- matrix(rnorm(N * n.factors, 0, 1), N, n.factors) 
        # Set diagonals to 1
        diag(lambda[[i]]) <- 1
        # Set upper tri to 0
        lambda[[i]][upper.tri(lambda[[i]])] <- 0
        w[[i]] <- matrix(NA, n.factors, J.w)
        if (sp) { # sfMsPGOcc
          if (cov.model == 'matern') {
            # Assume all spatial parameters ordered by svc first, then factor
            theta <- cbind(phi[((i - 1) * n.factors + 1):(i * n.factors)], 
          		 nu[((i - 1) * n.factors + 1):(i * n.factors)])
          } else {
            theta <- as.matrix(phi[((i - 1) * n.factors + 1):(i * n.factors)])
          }
          for (ll in 1:n.factors) {
            Sigma <- mkSpCov(coords, as.matrix(1), as.matrix(0), 
                	     theta[ll, ], cov.model)
            w[[i]][ll, ] <- rmvn(1, rep(0, J.w), Sigma)
          }

        } else { # lsMsPGOcc
          for (ll in 1:n.factors) {
            w[[i]][ll, ] <- rnorm(J.w)
          } # ll  
        }
        for (j in 1:J.w) {
          w.star[[i]][, j] <- lambda[[i]] %*% w[[i]][, j]
        }
      } else {
        if (sp) { # spMsPGOcc
          lambda <- NA
          if (cov.model == 'matern') {
            theta <- cbind(phi[((i - 1) * N + 1):(i * N)], 
          		 nu[((i - 1) * N + 1):(i * N)])
          } else {
            theta <- as.matrix(phi[((i - 1) * N + 1):(i * N)])
          }
          # Spatial random effects for each species
          for (ll in 1:N) {
            Sigma <- mkSpCov(coords, as.matrix(sigma.sq[(i - 1) * N + ll]), as.matrix(0), 
                	     theta[ll, ], cov.model)
            w.star[[i]][ll, ] <- rmvn(1, rep(0, J.w), Sigma)
          }
        }
        # For naming consistency
        w <- w.star
        lambda <- NA
      }
    }
  } # i (spatially-varying coefficient)
  # Design matrix for spatially-varying coefficients
  X.w <- X[, svc.cols, drop = FALSE]

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
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
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
    X.p.re <- array(NA, dim = c(J, n.rep.max, p.det.re))
    for (l in 1:p.det.re) {
      X.p.re[, , l] <- matrix(sample(1:p.RE$levels[l], J * n.rep.max, replace = TRUE), 
		              J, n.rep.max)	      
      for (i in 1:N) {
        alpha.star[i, which(alpha.star.indx == l)] <- rnorm(p.RE$levels[l], 0, sqrt(p.RE$sigma.sq.p[l]))
      }
    }
    for (j in 1:J) {
      X.p.re[j, -rep.indx[[j]], ] <- NA
    }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        X.p.re[, , j] <- X.p.re[, , j] + max(X.p.re[, , j - 1], na.rm = TRUE) 
      }
    }
    alpha.star.sites <- array(NA, c(N, J, n.rep.max))
    for (i in 1:N) {
      alpha.star.sites[i, , ] <- apply(X.p.re, c(1, 2), function(a) sum(alpha.star[i, a]))
    }
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }

  # Latent Occupancy Process ----------------------------------------------
  psi <- matrix(0, nrow = N, ncol = J)
  z <- matrix(0, nrow = N, ncol = J)
  range.ind <- matrix(NA, N, J)
  for (i in 1:N) {
    range.ind[i, ] <- rbinom(J, 1, range.probs[i])
    w.star.curr.mat <- sapply(w.star, function(a) a[i, ])
    for (j in 1:J) {  
      if (range.ind[i, j]) {
        if (sp | factor.model) {
          if (length(psi.RE) > 0) {
            psi[i, j] <- logit.inv(X[j, ] %*% as.matrix(beta[i, ]) + 
            		      X.w[j, ] %*% w.star.curr.mat[grid[j], ] + 
            		      beta.star.sites[i, j])
          } else {
            psi[i, j] <- logit.inv(X[j, ] %*% as.matrix(beta[i, ]) +
                                  X.w[j, ] %*% w.star.curr.mat[grid[j], ])
          }
        } else {
          if (length(psi.RE) > 0) {
            psi[i, j] <- logit.inv(X[j, ] %*% as.matrix(beta[i, ]) + beta.star.sites[i, j])
          } else {
            psi[i, j] <- logit.inv(X[j, ] %*% as.matrix(beta[i, ]))
          }
        }
        z[i, j] <- rbinom(1, 1, psi[i, j])
      } else {
        psi[i, j] <- 0
        z[i, j] <- 0
      }
    }
    z[i, ] <- rbinom(J, 1, psi[i, ])
  }

  # Data Formation --------------------------------------------------------
  p <- array(NA, dim = c(N, J, n.rep.max))
  y <- array(NA, dim = c(N, J, n.rep.max))
  for (i in 1:N) {
    for (j in 1:J) {
      if (length(p.RE) > 0) {
        p[i, j, rep.indx[[j]]] <- logit.inv(X.p[j, rep.indx[[j]], ] %*% as.matrix(alpha[i, ]) + alpha.star.sites[i, j, rep.indx[[j]]])
      } else {
        p[i, j, rep.indx[[j]]] <- logit.inv(X.p[j, rep.indx[[j]], ] %*% as.matrix(alpha[i, ]))
      }
 
        y[i, j, rep.indx[[j]]] <- rbinom(n.rep[j], 1, p[i, j, rep.indx[[j]]] * z[i, j]) 
    } # j
  } # i
  return(
    list(X = X, X.p = X.p, coords = coords, coords.full = coords.full, 
	 w = w, psi = psi, z = z, p = p, y = y, X.p.re = X.p.re, 
	 X.re = X.re, alpha.star = alpha.star, beta.star = beta.star, 
	 lambda = lambda, X.w = X.w, range.ind = range.ind)
  )
}
