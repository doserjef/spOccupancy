simTBinom <- function(J.x, J.y, n.time, weights, beta, sp.only = 0, 
		      trend = TRUE, psi.RE = list(), sp = FALSE, 
		      cov.model, sigma.sq, phi, nu, svc.cols = 1, 
                      ar1 = FALSE, rho, sigma.sq.t, x.positive = FALSE, ...) {

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
  # weights -----------------------------
  if (missing(weights)) {
    stop("error: weights must be specified.")
  }
  if (!is.matrix(weights)) {
    stop(paste("error: weights must be a matrix with ", J, " rows and ", max(n.time), " columns", sep = ''))
  }
  if (nrow(weights) != J | ncol(weights) != max(n.time)) {
    stop(paste("error: weights must be a matrix with ", J, " rows and ", max(n.time), " columns", sep = ''))
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
    if (length(beta) <= 1) {
      stop("error: beta must have at least two elements (intercept and trend)")
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
  # Spatial parameters ----------------
  if (length(svc.cols) > 1 & !sp) {
    stop("error: if simulating data with spatially-varying coefficients, set sp = TRUE")
  }
  if (sp) {
    if(missing(sigma.sq)) {
      stop("error: sigma.sq must be specified when sp = TRUE")
    }
    if(missing(phi)) {
      stop("error: phi must be specified when sp = TRUE")
    }
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
    p.svc <- length(svc.cols)
    if (length(phi) != p.svc) {
      stop("error: phi must have the same number of elements as svc.cols")
    }
    if (length(sigma.sq) != p.svc) {
      stop("error: sigma.sq must have the same number of elements as svc.cols")
    }
    if (cov.model == 'matern') {
      if (length(nu) != p.svc) {
        stop("error: nu must have the same number of elements as svc.cols")
      }
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

  if (x.positive) {
    if (p.occ > 1) {
      for (i in 2:p.occ) {
        X[, , i] <- abs(min(X[, , i])) + X[, , i]
      }
    }
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

  # Simulate spatial random effect ----------------------------------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  if (sp) {
    w.mat <- matrix(NA, J, p.svc)
    if (cov.model == 'matern') {
      theta <- cbind(phi, nu)
    } else {
      theta <- as.matrix(phi)
    }
    for (i in 1:p.svc) {
      Sigma <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
      # Random spatial process
      w.mat[, i] <- rmvn(1, rep(0, J), Sigma)
    }
    X.w <- X[, , svc.cols, drop = FALSE]
    w.sites <- matrix(0, J, n.time.max)
    for (j in 1:J) {
      for (t in 1:n.time.max) {
        w.sites[j, t] <- w.mat[j, ] %*% X.w[j, t, ] 
      }
    }
  } else {
    w.mat <- NA
    w.sites <- matrix(0, J, n.time.max) 
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

  # Simulate data ---------------------------------------------------------
  psi <- matrix(NA, J, max(n.time))
  y <- matrix(NA, J, max(n.time))
  for (j in 1:J) {
    for (t in 1:n.time.max) {
      if (length(psi.RE) > 0) {
        psi[j, t] <- logit.inv(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + beta.star.sites[j, t] + 
	                       eta[t])
      } else {
        psi[j, t] <- logit.inv(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + eta[t])
      }
    } # t
    for (t in 1:n.time[j]) {
      y[j, t] <- rbinom(1, weights[j, t], psi[j, t])
    }
  } # j

  return(
    list(X = X, coords = coords, psi = psi, y = y, w = w.mat,
	 X.re = X.re, beta.star = beta.star, eta = eta)
  )
}
