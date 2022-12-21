simIntOcc <- function(n.data, J.x, J.y, J.obs, n.rep, n.rep.max, beta, alpha, 
		      sp = FALSE, cov.model, sigma.sq, phi, nu, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
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
  if (missing(n.rep.max)) {
    n.rep.max <- sapply(n.rep, max, na.rm = TRUE)
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
  }
  # alpha -----------------------------
  if (missing(alpha)) {
    stop("error: alpha must be specified.")
  }
  if (!is.list(alpha)) {
    stop(paste("error: alpha must be a list with ", n.data, " vectors", sep = ''))
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

  # Get site ids for each of the data sets --------------------------------
  # Data sources can be obtained at multiple different sites. 
  sites <- list()
  for (i in 1:n.data) {
    sites[[i]] <- sort(sample(1:J, J.obs[i], replace = FALSE))
  }

  # Form detection covariates (if any) ------------------------------------
  X.p <- list()
  rep.indx <- list()
  for (i in 1:n.data) {
    rep.indx[[i]] <- list()
    n.alpha.curr <- length(alpha[[i]])
    K.curr <- n.rep[[i]]
    J.curr <- J.obs[[i]]
    for (j in 1:J.curr) {
      rep.indx[[i]][[j]] <- sample(1:n.rep.max[i], K.curr[j], replace = FALSE)
    }
    X.p[[i]] <- array(NA, dim = c(J.curr, n.rep.max[i], n.alpha.curr))
    X.p[[i]][, , 1] <- 1
    if (n.alpha.curr > 1) {
      for (q in 2:n.alpha.curr) {
        for (j in 1:J.curr) {
          X.p[[i]][j, rep.indx[[i]][[j]], q] <- rnorm(K.curr[j])
        } # j
      } # q
    }
  } # i

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

  # Latent Occupancy Process ----------------------------------------------
  if (sp) {
    psi <- logit.inv(X %*% as.matrix(beta) + w)
  } else {
      psi <- logit.inv(X %*% as.matrix(beta))
    }
  z <- rbinom(J, 1, psi)

  # Data Formation --------------------------------------------------------
  p <- list()
  y <- list()
  for (i in 1:n.data) {
    K.curr <- n.rep[[i]]
    J.curr <- J.obs[[i]]
    p[[i]] <- matrix(NA, nrow = J.curr, ncol = n.rep.max[i])
    y[[i]] <- matrix(NA, nrow = J.curr, ncol = n.rep.max[i])
    sites.curr <- sites[[i]]
    X.p.curr <- X.p[[i]]
    alpha.curr <- as.matrix(alpha[[i]])
    for (j in 1:J.curr) {
      p[[i]][j, rep.indx[[i]][[j]]] <- logit.inv(X.p.curr[j, rep.indx[[i]][[j]], ] %*% alpha.curr)
      y[[i]][j, rep.indx[[i]][[j]]] <- rbinom(K.curr[j], 1, p[[i]][j, rep.indx[[i]][[j]]] * z[sites.curr[j]])
    } # j
  } # i

  # Split up into observed and predicted ----------------------------------
  sites.obs <- sort(unique(unlist(sites))) 
  sites.pred <- (1:J)[!(1:J %in% sites.obs)]
  X.obs <- X[sites.obs, , drop = FALSE]
  X.pred <- X[sites.pred, , drop = FALSE]
  z.obs <- z[sites.obs]
  z.pred <- z[sites.pred]
  coords.obs <- coords[sites.obs,, drop = FALSE]
  coords.pred <- coords[sites.pred,, drop = FALSE]
  if (sp) {
    w.obs <- w[sites.obs, , drop = FALSE]
    w.pred <- w[sites.pred, , drop = FALSE]
  } else {
    w.obs <- NA
    w.pred <- NA
  }
  psi.obs <- psi[sites.obs, , drop = FALSE]
  psi.pred <- psi[sites.pred, , drop = FALSE]
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
	 z.pred = z.pred, p = p, y = y, sites = sites.return
	)
  )
}
