simBinom <- function(J.x, J.y, weights, beta, psi.RE = list(), sp = FALSE, 
		     svc.cols = 1, cov.model, sigma.sq, phi, nu, x.positive = FALSE, 
		     ...) {

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
  # weights -----------------------------
  if (missing(weights)) {
    stop("error: weights must be specified.")
  }
  if (length(weights) != J) {
    stop(paste("error: weights must be a vector of length ", J, sep = ''))
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
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
    if (!'beta.indx' %in% names(psi.RE)) {
      psi.RE$beta.indx <- list()
      for (i in 1:length(psi.RE$sigma.sq.psi)) {
        psi.RE$beta.indx[[i]] <- 1
      }
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
    if (x.positive) {
      for (i in 2:n.beta) {
        X[, i] <- runif(J, 0, 1)
        # X[, i] <- abs(min(X[, i])) + X[, i]
      }
    }
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
    X.w <- X[, svc.cols, drop = FALSE]
    # Convert w to a J*ptilde x 1 vector, sorted so that the p.svc values for 
    # each site are given, then the next site, then the next, etc.
    w <- c(t(w.mat))
    # Create X.tilde, which is a J x J*p.tilde matrix. 
    X.tilde <- matrix(0, J, J * p.svc)
    # Fill in the matrix
    for (j in 1:J) {
      X.tilde[j, ((j - 1) * p.svc + 1):(j * p.svc)] <- X.w[j, ]
    }
    w.sites <- X.tilde %*% w
  } else {
    w.mat <- NA
    X.w <- NA
    w.sites <- NA
  }

  # Random effects --------------------------------------------------------
  if (length(psi.RE) > 0) {
    p.occ.re <- length(unlist(psi.RE$beta.indx))
    tmp <- sapply(psi.RE$beta.indx, length)
    re.col.indx <- unlist(lapply(1:length(psi.RE$beta.indx), function(a) rep(a, tmp[a])))
    sigma.sq.psi <- psi.RE$sigma.sq.psi[re.col.indx]
    n.occ.re.long <- psi.RE$levels[re.col.indx]
    n.occ.re <- sum(n.occ.re.long)
    beta.star.indx <- rep(1:p.occ.re, n.occ.re.long)
    beta.star <- rep(0, n.occ.re)
    X.random <- X[, unlist(psi.RE$beta.indx), drop = FALSE]
    n.random <- ncol(X.random)
    X.re <- matrix(NA, J, length(psi.RE$levels))
    for (i in 1:length(psi.RE$levels)) {
      X.re[, i] <- sample(1:psi.RE$levels[i], J, replace = TRUE)  
    }
    indx.mat <- X.re[, re.col.indx, drop = FALSE]
    for (i in 1:p.occ.re) {
      beta.star[which(beta.star.indx == i)] <- rnorm(n.occ.re.long[i], 0, 
						     sqrt(sigma.sq.psi[i]))
    }
    if (length(psi.RE$levels) > 1) {
      for (j in 2:length(psi.RE$levels)) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
      }
    }
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        indx.mat[, j] <- indx.mat[, j] + max(indx.mat[, j - 1], na.rm = TRUE)
      }
    }
    beta.star.sites <- rep(NA, J)
    for (j in 1:J) {
      beta.star.sites[j] <- beta.star[indx.mat[j, , drop = FALSE]] %*% t(X.random[j, , drop = FALSE])
    }
  } else {
    X.re <- NA
    beta.star <- NA
  }
  # Latent Occupancy Process ----------------------------------------------
  if (sp) {
    if (length(psi.RE) > 0) {
      psi <- logit.inv(X %*% as.matrix(beta) + X.tilde %*% w + beta.star.sites)
    } else {
      psi <- logit.inv(X %*% as.matrix(beta) + X.tilde %*% w)
    }
  } else {
    if (length(psi.RE) > 0) {
      psi <- logit.inv(X %*% as.matrix(beta) + beta.star.sites)
    } else {
      psi <- logit.inv(X %*% as.matrix(beta))
    }
  }
  y <- rbinom(J, weights, psi)

  return(
    list(X = X, coords = coords, w = w.mat, psi = psi, y = y, 
	 X.re = X.re, X.w = X.w, beta.star = beta.star)
  )
}
