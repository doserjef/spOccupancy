simIntMsOcc <- function(n.data, J.x, J.y, J.obs, n.rep, N, beta, alpha, psi.RE = list(), 
		        p.RE = list(), sp = FALSE, cov.model, 
		        sigma.sq, phi, nu, factor.model = FALSE, n.factors, ...) {

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
  # n.data -------------------------------
  if (missing(n.data)) {
    stop("error: n.data must be specified")
  }
  if (length(n.data) != 1) {
    stop("error: n.data must be a single numeric value.")
  }
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
  # J.obs -----------------------------
  if (missing(J.obs)) {
    stop("error: J.obs must be specified")
  }
  if (length(J.obs) != n.data) {
    stop(paste("error: J.obs must be a vector of length ", n.data, sep = ''))
  }
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (!is.list(n.rep)) {
    stop(paste("error: n.rep must be a list of ", n.data, " vectors", sep = ''))
  }
  if (length(n.rep) != n.data) {
    stop(paste("error: n.rep must be a list of ", n.data, " vectors", sep = ''))
  }
  for (i in 1:n.data) {
    if (length(n.rep[[i]]) != J.obs[i]) {
      stop(paste("error: n.rep[[", i, "]] must be of length ", J.obs[i], sep = ''))
    }
  }
  # N ---------------------------------
  if (missing(N)) {
    stop("error: N must be specified")
  }
  if (length(N) != n.data) {
    stop(paste("error: N must be a vector of ", n.data, 
	       " values indicating the number of species in each data set", sep = ''))
  }
  N.max <- max(N)
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified")
  }
  if (!is.matrix(beta)) {
    stop(paste("error: beta must be a numeric matrix with ", N.max, " rows", sep = ''))
  }
  if (nrow(beta) != N.max) {
    stop(paste("error: beta must be a numeric matrix with ", N.max, " rows", sep = ''))
  }
  # alpha -----------------------------
  if (missing(alpha)) {
    stop("error: alpha must be specified.")
  }
  if (!is.list(alpha)) {
    stop(paste("error: alpha must be a list with ", n.data, " vectors", sep = ''))
  }
  if (sp & !factor.model) {
    # sigma.sq --------------------------
    if (missing(sigma.sq)) {
      stop("error: sigma.sq must be specified when sp = TRUE")
    }
    if (length(sigma.sq) != N.max) {
      stop(paste("error: sigma.sq must be a vector of length ", N.max, sep = ''))
    }
    # phi -------------------------------
    if(missing(phi)) {
      stop("error: phi must be specified when sp = TRUE")
    }
    if (length(phi) != N.max) {
      stop(paste("error: phi must be a vector of length ", N.max, sep = ''))
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
  if (factor.model) {
    # n.factors -----------------------
    if (missing(n.factors)) {
      stop("error: n.factors must be specified when factor.model = TRUE")
    }
    if (sp) {
      if (!missing(sigma.sq)) {
        message("sigma.sq is specified but will be set to 1 for spatial latent factor model")
      }
      if(missing(phi)) {
        stop("error: phi must be specified when sp = TRUE")
      }
      if (length(phi) != n.factors) {
        stop(paste("error: phi must be a vector of length ", n.factors, sep = ''))
      }
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
  # TODO: detection REs not currently supported.
  # p.RE ----------------------------
  # names(p.RE) <- tolower(names(p.RE))
  # if (!is.list(p.RE)) {
  #   stop("error: if specified, p.RE must be a list with tags 'levels' and 'sigma.sq.p'")
  # }
  # if (length(names(p.RE)) > 0) {
  #   if (!'sigma.sq.p' %in% names(p.RE)) {
  #     stop("error: sigma.sq.p must be a tag in p.RE with values for the detection random effect variances")
  #   }
  #   if (!'levels' %in% names(p.RE)) {
  #     stop("error: levels must be a tag in p.RE with the number of random effect levels for each detection random intercept.")
  #   }
  # }

  # Subroutines -----------------------------------------------------------
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

  # Get site ids for each of the data sets --------------------------------
  # Data sources can be obtained at multiple different sites. 
  sites <- list()
  for (i in 1:n.data) {
    sites[[i]] <- sort(sample(1:J, J.obs[i], replace = FALSE))
  }

  # Form detection covariates (if any) ------------------------------------
  X.p <- list()
  for (i in 1:n.data) {
    n.alpha.curr <- ncol(alpha[[i]])
    K.curr <- n.rep[[i]]
    J.curr <- J.obs[[i]]
    X.p[[i]] <- array(NA, dim = c(J.curr, max(K.curr), n.alpha.curr))
    X.p[[i]][, , 1] <- 1
    if (n.alpha.curr > 1) {
      for (q in 2:n.alpha.curr) {
        for (j in 1:J.curr) {
          X.p[[i]][j, 1:K.curr[j], q] <- rnorm(K.curr[j])
        } # j
      } # q
    }
  } # i

  # Simulate latent (spatial) random effect for each species --------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  w.star <- matrix(0, nrow = N.max, ncol = J)
  if (factor.model) {
    lambda <- matrix(rnorm(N.max * n.factors, 0, 1), N.max, n.factors) 
    # Set diagonals to 1
    diag(lambda) <- 1
    # Set upper tri to 0
    lambda[upper.tri(lambda)] <- 0
    w <- matrix(NA, n.factors, J)
    if (sp) { # sfIntMsPGOcc
      if (cov.model == 'matern') {
        theta <- cbind(phi, nu)
      } else {
        theta <- as.matrix(phi)
      }
      for (ll in 1:n.factors) {
        Sigma <- mkSpCov(coords, as.matrix(1), as.matrix(0), 
            	     theta[ll, ], cov.model)
        w[ll, ] <- rmvn(1, rep(0, J), Sigma)
      }

    } else { # lfIntMsPGOcc
      for (ll in 1:n.factors) {
        w[ll, ] <- rnorm(J)
      } # ll  
    }
    for (j in 1:J) {
      w.star[, j] <- lambda %*% w[, j]
    }
  } else {
    if (sp) { # not supported
      stop("integrated multi-species spatial occupancy models are only supported for a factor model formulation. Set factor.model = TRUE")
    }
    # For naming consistency
    w <- w.star
    lambda <- NA
  }

  # Random effects --------------------------------------------------------
  if (length(psi.RE) > 0) {
    p.occ.re <- length(psi.RE$levels)
    sigma.sq.psi <- rep(NA, p.occ.re)
    n.occ.re.long <- psi.RE$levels
    n.occ.re <- sum(n.occ.re.long)
    beta.star.indx <- rep(1:p.occ.re, n.occ.re.long)
    beta.star <- matrix(0, N.max, n.occ.re)
    X.re <- matrix(NA, J, p.occ.re)
    for (l in 1:p.occ.re) {
      X.re[, l] <- sample(1:psi.RE$levels[l], J, replace = TRUE)         
      for (i in 1:N.max) {
        beta.star[i, which(beta.star.indx == l)] <- rnorm(psi.RE$levels[l], 0, 
							  sqrt(psi.RE$sigma.sq.psi[l]))
      }
    }
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
      }
    } 
    beta.star.sites <- matrix(NA, N.max, J)
    for (i in 1:N.max) {
      beta.star.sites[i, ] <- apply(X.re, 1, function(a) sum(beta.star[i, a]))
    }
  } else {
    X.re <- NA
    beta.star <- NA
  }
  # TODO: detection REs not currently supported. 
  # Detection (not yet implemented)
  # if (length(p.RE) > 0) {
  #   p.det.re <- length(p.RE$levels)
  #   sigma.sq.p <- rep(NA, p.det.re)
  #   n.det.re.long <- p.RE$levels
  #   n.det.re <- sum(n.det.re.long)
  #   alpha.star.indx <- rep(1:p.det.re, n.det.re.long)
  #   alpha.star <- matrix(0, N.max, n.det.re)
  #   X.p.re <- array(NA, dim = c(J, max(n.rep), p.det.re))
  #   for (l in 1:p.det.re) {
  #     X.p.re[, , l] <- matrix(sample(1:p.RE$levels[l], J * max(n.rep), replace = TRUE), 
  #       	              J, max(n.rep))	      
  #     for (i in 1:N.max) {
  #       alpha.star[i, which(alpha.star.indx == l)] <- rnorm(p.RE$levels[l], 0, sqrt(p.RE$sigma.sq.p[l]))
  #     }
  #   }
  #   for (j in 1:J) {
  #     X.p.re[j, -(1:n.rep[j]), ] <- NA
  #   }
  #   if (p.det.re > 1) {
  #     for (j in 2:p.det.re) {
  #       X.p.re[, , j] <- X.p.re[, , j] + max(X.p.re[, , j - 1], na.rm = TRUE) 
  #     }
  #   }
  #   alpha.star.sites <- array(NA, c(N.max, J, max(n.rep)))
  #   for (i in 1:N.max) {
  #     alpha.star.sites[i, , ] <- apply(X.p.re, c(1, 2), function(a) sum(alpha.star[i, a]))
  #   }
  # } else {
    X.p.re <- NA
    alpha.star <- NA
  # }

  # Get species ids for each of the data sets -----------------------------
  # Data sources can sample different amounts of all the species in the community. 
  species <- list()
  for (i in 1:n.data) {
    species[[i]] <- sort(sample(1:N.max, N[i], replace = FALSE))    
  }

  # Latent Occupancy Process ----------------------------------------------
  psi <- matrix(NA, nrow = N.max, ncol = J)
  z <- matrix(NA, nrow = N.max, ncol = J)
  for (i in 1:N.max) {
    if (sp | factor.model) {
      if (length(psi.RE) > 0) {
        psi[i, ] <- logit.inv(X %*% as.matrix(beta[i, ]) + w.star[i, ] + beta.star.sites[i, ])
      } else {
        psi[i, ] <- logit.inv(X %*% as.matrix(beta[i, ]) + w.star[i, ])
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
  p <- list()
  y <- list()
  for (t in 1:n.data) {
    p[[t]] <- array(NA, dim = c(N[t], J.obs[t], max(n.rep[[t]])))
    y[[t]] <- array(NA, dim = c(N[t], J.obs[t], max(n.rep[[t]])))
    for (i in 1:N[t]) {
      for (j in 1:J.obs[t]) {
        # if (length(p.RE) > 0) {
        #   # TODO: detection REs not currently supported 
        #   p[[t]][i, j, 1:n.rep[[t]][j]] <- logit.inv(X.p[[t]][j, 1:n.rep[[t]][j], ] %*% as.matrix(alpha[[t]][i, ]) + alpha.star.sites[i, j, 1:n.rep[[t]][j]])
        # } else {
          p[[t]][i, j, 1:n.rep[[t]][j]] <- logit.inv(X.p[[t]][j, 1:n.rep[[t]][j], ] %*% as.matrix(alpha[[t]][i, ]))
        # }
          y[[t]][i, j, 1:n.rep[[t]][j]] <- rbinom(n.rep[[t]][j], 1, p[[t]][i, j, 1:n.rep[[t]][j]] * z[species[[t]][i], sites[[t]][j]]) 
      } # j
    } # i
  }

  # Split up into observed and predicted ----------------------------------
  sites.obs <- sort(unique(unlist(sites))) 
  sites.pred <- (1:J)[!(1:J %in% sites.obs)]
  X.obs <- X[sites.obs, , drop = FALSE]
  X.pred <- X[sites.pred, , drop = FALSE]
  if (length(psi.RE) > 0) {
    X.re.obs <- X.re[sites.obs, , drop = FALSE]
    X.re.pred <- X.re[sites.pred, , drop = FALSE]
  } else {
    X.re.obs <- NA
    X.re.pred <- NA
  }
  z.obs <- z[, sites.obs]
  z.pred <- z[, sites.pred]
  coords.obs <- coords[sites.obs,, drop = FALSE]
  coords.pred <- coords[sites.pred,, drop = FALSE]
  if (sp) {
    w.obs <- w[, sites.obs, drop = FALSE]
    w.pred <- w[, sites.pred, drop = FALSE]
  } else {
    w.obs <- NA
    w.pred <- NA
  }
  psi.obs <- psi[, sites.obs, drop = FALSE]
  psi.pred <- psi[, sites.pred, drop = FALSE]
  sites.vec <- unlist(sites)
  # Need to get index of each site value in only the sites.obs
  sites.new <- rep(0, length(sites.vec))
  for (i in 1:length(sites.vec)) {
    sites.new[i] <- which(sites.obs == sites.vec[i])
  }
  sites.return <- list()
  indx <- 1
  for (i in 1:n.data) {
    sites.return[[i]] <- sites.new[indx:(indx + J.obs[i] - 1)]
    indx <- indx + J.obs[i]
  }

  return(
    list(X.obs = X.obs, X.pred = X.pred, X.p = X.p, 
	 coords.obs = coords.obs, coords.pred = coords.pred, 
	 w.obs = w.obs, w.pred = w.pred, 
	 psi.obs = psi.obs, psi.pred = psi.pred, z.obs = z.obs, 
	 z.pred = z.pred, p = p, y = y, sites = sites.return,
	 X.p.re = X.p.re, X.re.obs = X.re.obs, 
	 X.re.pred = X.re.pred, alpha.star = alpha.star, 
	 beta.star = beta.star, lambda = lambda, species = species
	)
  )
}
