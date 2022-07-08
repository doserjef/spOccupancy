simTOcc <- function(J.x, J.y, n.time, n.rep, beta, alpha, sp.only = 0, 
		    trend = TRUE, psi.RE = list(), p.RE = list(), sp = FALSE, 
		    cov.model, sigma.sq, phi, nu, ar1 = FALSE, rho, sigma.sq.t, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
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
  # n.time ---------------------------
  if (missing(n.time)) {
    stop("error: n.time must be specified.")
  }
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (!is.matrix(n.rep)) {
    stop(paste("error: n.rep must be a matrix with ", J, " rows and ", max(n.time), " columns", sep = ''))
  }
  if (nrow(n.rep) != J | ncol(n.rep) != max(n.time)) {
    stop(paste("error: n.rep must be a matrix with ", J, " rows and ", max(n.time), " columns", sep = ''))
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
    if (length(beta) <= 1) {
      stop("error: beta must have at least two elements (intercept and trend)")
    }
  }
  # alpha -----------------------------
  if (missing(alpha)) {
    stop("error: alpha must be specified.")
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
  # AR1 -------------------------------
  if (ar1) {
    if (missing(rho)) {
      stop("error: rho must be specified when ar1 = TRUE")
    }
    if (missing(sigma.sq.t)) {
      stop("error: sigma.sq.t must be specified when ar1 = TRUE")
    }
  }
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

  # Occurrence ------------------------------------------------------------
  p.occ <- length(beta)
  n.time.max <- max(n.time, na.rm = TRUE)
  X <- array(NA, dim = c(J, n.time.max, p.occ))
  X[, , 1] <- 1
  if (p.occ > 1) {
    if (trend) { # If simulating data with a trend
      # By default the second simulated covariate is a standardized trend
      X[, , 2] <- scale(c(matrix(rep(1:n.time.max, each = J), nrow = J, ncol = n.time.max)))
      if (p.occ > 2) {
        for (i in 3:p.occ) {
          if (i %in% sp.only) {
            X[, , i] <- rep(rnorm(J), n.time.max)
          } else {
            X[, , i] <- rnorm(J * n.time.max)
          }
        }
      }
    } else { # If not simulating data with a trend
      if (p.occ > 1) {
        for (i in 2:p.occ) {
          if (i %in% sp.only) {
            X[, , i] <- rep(rnorm(J), n.time.max)
          } else {
            X[, , i] <- rnorm(J * n.time.max)
          }
        }
      }
    }
  }

  # Detection -------------------------------------------------------------
  # Time dependent --------------------
  p.det <- length(alpha)
  n.rep.max <- max(n.rep, na.rm = TRUE)
  X.p <- array(NA, dim = c(J, n.time.max, n.rep.max, p.det))
  X.p[, , , 1] <- 1
  if (p.det > 1) {
    for (j in 1:J) {
      for (t in 1:n.time[j]) {
        for (k in 1:n.rep[j, t]) {
          X.p[j, t, k, 2:p.det] <- rnorm(p.det - 1)
        } # k
      } # t
    } # j  
  }

  # Random effects --------------------------------------------------------
  if (length(psi.RE) > 0) {
    p.occ.re <- length(psi.RE$levels)
    sigma.sq.psi <- rep(NA, p.occ.re)
    n.occ.re.long <- psi.RE$levels
    n.occ.re <- sum(n.occ.re.long)
    beta.star.indx <- rep(1:p.occ.re, n.occ.re.long)
    beta.star <- rep(0, n.occ.re)
    X.re <- array(NA, dim = c(J, n.time.max, p.occ.re))
    for (i in 1:p.occ.re) {
      if (length(psi.RE$site.re) == 0) psi.RE$site.re <- FALSE
      if (psi.RE$site.re == TRUE) {
        if (i == 1) {
          site.vals <- 1:J 
          X.re[, , i] <- site.vals
	} else {
          X.re[, , i] <- sample(1:psi.RE$levels[i], J * n.time.max, replace = TRUE)         
	}
      } else {
        X.re[, , i] <- sample(1:psi.RE$levels[i], J * n.time.max, replace = TRUE)         
      }
      beta.star[which(beta.star.indx == i)] <- rnorm(psi.RE$levels[i], 0, sqrt(psi.RE$sigma.sq.psi[i]))
    }
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, , j] <- X.re[, , j] + max(X.re[, , j - 1])
      }
    } 
    beta.star.sites <- apply(X.re, c(1, 2), function(a) sum(beta.star[a]))
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
    X.p.re <- array(NA, dim = c(J, n.time.max, n.rep.max, p.det.re))
    for (i in 1:p.det.re) {
      X.p.re[, , , i] <- sample(1:p.RE$levels[i], J * n.time.max * n.rep.max, replace = TRUE) 
      alpha.star[which(alpha.star.indx == i)] <- rnorm(p.RE$levels[i], 0, sqrt(p.RE$sigma.sq.p[i]))
    }
    # for (j in 1:J) {
    #   X.p.re[j, , -(1:n.rep[j]), ] <- NA
    # }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        X.p.re[, , , j] <- X.p.re[, , , j] + max(X.p.re[, , , j - 1]) 
      }
    }
    alpha.star.sites <- apply(X.p.re, c(1, 2, 3), function(a) sum(alpha.star[a]))
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }

  # Simulate spatial random effect ----------------------------------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  if (sp) {
    if (cov.model == 'matern') {
      theta <- c(phi, nu)
    } else {
      theta <- phi
    }
    Sigma <- mkSpCov(coords, as.matrix(sigma.sq), as.matrix(0), theta, cov.model)
    # Random spatial process
    w <- rmvn(1, rep(0, J), Sigma)
  } else {
    w <- matrix(rep(0, J))
  }

  # Simulate temporal (AR1) random effect ---------------------------------
  if (ar1) {
    exponent <- abs(matrix(1:n.time.max - 1, nrow = n.time.max, 
			   ncol = n.time.max, byrow = TRUE) - (1:n.time.max - 1))
    Sigma.eta <- sigma.sq.t * rho^exponent
    eta <- rmvn(1, rep(0, n.time.max), Sigma.eta)
  } else {
    eta <- matrix(rep(0, n.time.max))
  }

  # Latent Occupancy Process ----------------------------------------------
  psi <- matrix(NA, J, max(n.time))
  z <- matrix(NA, J, max(n.time))
  for (j in 1:J) {
    for (t in 1:n.time.max) {
      if (length(psi.RE) > 0) {
        psi[j, t] <- logit.inv(X[j, t, ] %*% as.matrix(beta) + w[j] + beta.star.sites[j, t] + 
	                       eta[t])
      } else {
        psi[j, t] <- logit.inv(X[j, t, ] %*% as.matrix(beta) + w[j] + eta[t])
      }
      z[j, t] <- rbinom(1, 1, psi[j, t])
    } # t
  } # j

  # Detection Model -------------------------------------------------------
  p <- array(NA, dim = c(J, max(n.time), n.rep.max))
  y <- array(NA, dim = c(J, max(n.time), n.rep.max))
  for (j in 1:J) {
    for (t in 1:n.time[j]) {
      if (length(p.RE) > 0) {
        p[j, t, 1:n.rep[j, t]] <- logit.inv(X.p[j, t, 1:n.rep[j, t], ] %*% as.matrix(alpha) + 
	                                    alpha.star.sites[j, t, 1:n.rep[j, t]])
      } else {
        p[j, t, 1:n.rep[j, t]] <- logit.inv(X.p[j, t, 1:n.rep[j, t], ] %*% as.matrix(alpha))
      }
      y[j, t, 1:n.rep[j, t]] <- rbinom(n.rep[j, t], 1, p[j, t, 1:n.rep[j, t]] * z[j, t]) 
    } # t
  } # j

  return(
    list(X = X, X.p = X.p, coords = coords,
	 psi = psi, z = z, p = p, y = y, w = w,
	 X.p.re = X.p.re, X.re = X.re,
	 alpha.star = alpha.star, beta.star = beta.star, eta = eta)
  )
}
