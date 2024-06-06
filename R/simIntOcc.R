simIntOcc <- function(n.data, J.x, J.y, J.obs, n.rep, n.rep.max, beta, alpha, 
                      psi.RE = list(), p.RE = list(),
                      sp = FALSE, cov.model, sigma.sq, phi, nu, ...) {

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
  if (!is.list(p.RE)) {
    stop(paste("error: if species, p.RE must be a list with ", n.data, " lists", sep = ''))
  }
  if (length(p.RE) > 0) {
    for (q in 1:n.data) {
      names(p.RE[[q]]) <- tolower(names(p.RE[[q]]))
      if (!is.list(p.RE[[q]])) {
        stop("error: if specified, p.RE[[", q, "]] must be a list with tags 'levels' and 'sigma.sq.p'")
      }
      if (length(names(p.RE[[q]])) > 0) {
        if (!'sigma.sq.p' %in% names(p.RE[[q]])) {
          stop("error: sigma.sq.p must be a tag in p.RE[[", q, "]] with values for the detection random effect variances")
        }
        if (!'levels' %in% names(p.RE[[q]])) {
          stop("error: levels must be a tag in p.RE[[", q, "]] with the number of random effect levels for each detection random intercept.")
        }
      }
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
    w <- as.matrix(rmvn(1, rep(0, J), Sigma))
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
    p.det.re <- list()
    n.det.re.long <- list()
    n.det.re <- list()
    alpha.star.indx <- list()
    alpha.star <- list()
    alpha.star.sites <- list()
    X.p.re <- list()
    for (q in 1:n.data) {
      if (length(p.RE[[q]]) > 0) {
      p.det.re[[q]] <- length(p.RE[[q]]$levels)
      n.det.re.long[[q]] <- p.RE[[q]]$levels
      n.det.re[[q]] <- sum(n.det.re.long[[q]])
      alpha.star.indx[[q]] <- rep(1:p.det.re[[q]], n.det.re.long[[q]])
      alpha.star[[q]] <- rep(0, n.det.re[[q]])
      X.p.re[[q]] <- array(NA, dim = c(J.obs[[q]], max(n.rep[[q]]), p.det.re[[q]]))
      for (i in 1:p.det.re[[q]]) {
        X.p.re[[q]][, , i] <- matrix(sample(1:p.RE[[q]]$levels[i], 
					    J.obs[[q]] * max(n.rep[[q]]), replace = TRUE), 
          	              J.obs[[q]], max(n.rep[[q]]))	      
        alpha.star[[q]][which(alpha.star.indx[[q]] == i)] <- rnorm(p.RE[[q]]$levels[i], 
								   0, 
								   sqrt(p.RE[[q]]$sigma.sq.p[i]))
      }
      for (j in 1:J.obs[[q]]) {
        X.p.re[[q]][j, -rep.indx[[q]][[j]], ] <- NA
      }
      if (p.det.re[[q]] > 1) {
        for (j in 2:p.det.re[[q]]) {
          X.p.re[[q]][, , j] <- X.p.re[[q]][, , j] + max(X.p.re[[q]][, , j - 1], na.rm = TRUE) 
        }
      }
        alpha.star.sites[[q]] <- apply(X.p.re[[q]], c(1, 2), function(a) sum(alpha.star[[q]][a]))
      } else {
        X.p.re[[q]] <- NA
        alpha.star[[q]] <- NA
	alpha.star.sites[[q]] <- NA
      }
    }
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
    if (length(p.RE) > 0) {
      alpha.star.sites.curr <- alpha.star.sites[[i]]
    }
    for (j in 1:J.curr) {
      if (length(p.RE) > 0) { # If any detection random effects
        if (length(p.RE[[i]]) > 0) { # If any detection random effects in this data set
          p[[i]][j, rep.indx[[i]][[j]]] <- logit.inv(X.p.curr[j, rep.indx[[i]][[j]], ] %*% alpha.curr + 
                                              alpha.star.sites.curr[j, rep.indx[[i]][[j]]])
        } else { # Detection random effects, but none in this data set
          p[[i]][j, rep.indx[[i]][[j]]] <- logit.inv(X.p.curr[j, rep.indx[[i]][[j]], ] %*% alpha.curr)
	}
      }	else { # No detection random effects
        p[[i]][j, rep.indx[[i]][[j]]] <- logit.inv(X.p.curr[j, rep.indx[[i]][[j]], ] %*% alpha.curr)
      }
      y[[i]][j, rep.indx[[i]][[j]]] <- rbinom(K.curr[j], 1, p[[i]][j, rep.indx[[i]][[j]]] * z[sites.curr[j]])
    } # j
  } # i

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
	 z.pred = z.pred, p = p, y = y, sites = sites.return, 
	 X.re.obs = X.re.obs, X.re.pred = X.re.pred, beta.star = beta.star, 
	 X.p.re = X.p.re, alpha.star = alpha.star
	)
  )
}

