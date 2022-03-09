simOcc <- function(J.x, J.y, n.rep, beta, alpha, psi.RE = list(), p.RE = list(), 
		   sp = FALSE, cov.model, sigma.sq, phi, nu, ...) {

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
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (length(n.rep) != J) {
    stop(paste("error: n.rep must be a vector of length ", J, sep = ''))
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
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
  # Spatial parameters ----------------
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

  # Form occupancy covariates (if any) ------------------------------------
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
    w <- NA
  }

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
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
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
        X.p.re[, , j] <- X.p.re[, , j] + max(X.p.re[, , j - 1], na.rm = TRUE) 
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
    list(X = X, X.p = X.p, coords = coords,
         w = w, psi = psi, z = z, p = p, y = y, 
	 X.p.re = X.p.re, X.re = X.re, 
	 alpha.star = alpha.star, beta.star = beta.star)
  )
}
