simTOcc <- function(J.x, J.y, n.time, n.rep, n.rep.max, beta, alpha, sp.only = 0, 
		    trend = TRUE, psi.RE = list(), p.RE = list(), sp = FALSE, 
		    svc.cols = 1, cov.model, sigma.sq, phi, nu, ar1 = FALSE, 
                    rho, sigma.sq.t, x.positive = FALSE, mis.spec.type = 'none', 
		    scale.param = 1, avail, grid, ...) {

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
  if (missing(n.rep.max)) {
    n.rep.max <- max(n.rep, na.rm = TRUE)
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
  if (length(svc.cols) > 1 & !sp) {
    stop("error: if simulating data with spatially-varying coefficients, set sp = TRUE")
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


  # Mis-specification -----------------
  if (!(mis.spec.type %in% c("none", "scale", "line", "probit"))) {
    stop("mis-specification type not allowed.")
  }
  
  if (scale.param <= 0) {
    stop("scale parameter must be greater than zero.")
  }

  # Availability process --------------------------------------------------
  if (missing(avail)) {
    avail <- array(1, dim = c(J, max(n.time, na.rm = TRUE), n.rep.max))  
  } else {
    if (length(dim(avail)) != 3) {
      stop(paste0("avail must be an array with dimensions of ", J, " x ", max(n.time), " x ", max(n.rep), "."))
    }
    if (dim(avail)[1] != J | dim(avail)[2] != max(n.time) | dim(avail)[3] != max(n.rep)) {
      stop(paste0("avail must be an array with dimensions of ", J, " x ", max(n.time), " x ", max(n.rep), "."))
    }
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

  # Occurrence ------------------------------------------------------------
  p.occ <- length(beta)
  n.time.max <- max(n.time, na.rm = TRUE)
  time.indx <- list()
  for (j in 1:J) {
    time.indx[[j]] <- sample(which(!is.na(n.rep[j, ])), n.time[j], replace = FALSE) 
  }
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
        X[, , i] <- runif(J * n.time.max, 0, 1)
      }
    }
  }

  # Detection -------------------------------------------------------------
  # Time dependent --------------------
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- list()
    for (t in time.indx[[j]]) {
      rep.indx[[j]][[t]] <- sample(1:n.rep.max, n.rep[j, t], replace = FALSE)
    }
  }
  p.det <- length(alpha)
  X.p <- array(NA, dim = c(J, n.time.max, n.rep.max, p.det))
  X.p[, , , 1] <- 1
  if (p.det > 1) {
    for (j in 1:J) {
      for (t in time.indx[[j]]) {
        for (k in rep.indx[[j]][[t]]) {
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
  if (sp) {
    w.mat.full <- matrix(NA, J, p.svc)
    if (cov.model == 'matern') {
      theta <- cbind(phi, nu)
    } else {
      theta <- as.matrix(phi)
    }
    w.mat <- matrix(NA, nrow(coords), p.svc)
    for (i in 1:p.svc) {
      Sigma.full <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
      w.mat[, i] <- rmvn(1, rep(0, nrow(Sigma.full)), Sigma.full)
      # Random spatial process
      w.mat.full[, i] <- w.mat[grid, i]
    }
    X.w <- X[, , svc.cols, drop = FALSE]
    w.sites <- matrix(0, J, n.time.max)
    for (j in 1:J) {
      for (t in 1:n.time.max) {
        w.sites[j, t] <- w.mat.full[j, ] %*% X.w[j, t, ] 
      }
    }
  } else {
    w.mat <- NA
    w.mat.full <- NA
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

  # Latent Occupancy Process ----------------------------------------------
  psi <- matrix(NA, J, max(n.time))
  z <- matrix(NA, J, max(n.time))
  for (j in 1:J) {
    for (t in 1:n.time.max) {
      if (mis.spec.type == "none") {
      
      if (length(psi.RE) > 0) {
        psi[j, t] <- logit.inv(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + beta.star.sites[j, t] + 
                                 eta[t])
      } else {
        psi[j, t] <- logit.inv(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + eta[t])
      }
      z[j, t] <- rbinom(1, 1, psi[j, t])
      
      } else if (mis.spec.type == "scale") {
        
        if (length(psi.RE) > 0) {
          psi[j, t] <- scale.param * logit.inv(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + beta.star.sites[j, t] + 
                                   eta[t]) 
          psi[j, t] <- pmax(psi[j, t], 0)
          psi[j, t] <- pmin(psi[j, t], 1)
        } else {
          psi[j, t] <- scale.param * logit.inv(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + eta[t])
          psi[j, t] <- pmax(psi[j, t], 0)
          psi[j, t] <- pmin(psi[j, t], 1)
        }
        z[j, t] <- rbinom(1, 1, psi[j, t])
        
        
      } else if (mis.spec.type == "line") {
        if (length(psi.RE) > 0) {
          psi[j, t] <- X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + beta.star.sites[j, t] + 
                                   eta[t] 
          psi[j, t] <- pmax(psi[j, t], 0)
          psi[j, t] <- pmin(psi[j, t], 1)
        } else {
          psi[j, t] <- X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + eta[t]
          psi[j, t] <- pmax(psi[j, t], 0)
          psi[j, t] <- pmin(psi[j, t], 1)
        }
        z[j, t] <- rbinom(1, 1, psi[j, t])
        
      } else if (mis.spec.type == 'probit') {
      
        if (length(psi.RE) > 0) {
          psi[j, t] <- pnorm(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + beta.star.sites[j, t] + 
                                 eta[t])
        } else {
          psi[j, t] <- pnorm(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + eta[t])
        }
        z[j, t] <- rbinom(1, 1, psi[j, t])
      }
    } # t
  } # j

  # Detection Model -------------------------------------------------------
  p <- array(NA, dim = c(J, max(n.time), n.rep.max))
  y <- array(NA, dim = c(J, max(n.time), n.rep.max))
  for (j in 1:J) {
    for (t in time.indx[[j]]) {
      
      if (mis.spec.type == "none") {
        if (length(p.RE) > 0) {
          p[j, t, rep.indx[[j]][[t]]] <- logit.inv(X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha) + 
                                                     alpha.star.sites[j, t, rep.indx[[j]][[t]]])
        } else {
          p[j, t, rep.indx[[j]][[t]]] <- logit.inv(X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha))
        }
        
      } else if (mis.spec.type == "scale") {
        
        if (length(p.RE) > 0) {
          p[j, t, rep.indx[[j]][[t]]] <- (1/scale.param) * logit.inv(X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha) + 
                                                     alpha.star.sites[j, t, rep.indx[[j]][[t]]]) 
          p[j, t, rep.indx[[j]][[t]]] <- pmax(p[j, t, rep.indx[[j]][[t]]], 0)
          p[j, t, rep.indx[[j]][[t]]] <- pmin(p[j, t, rep.indx[[j]][[t]]], 1)
        } else {
          p[j, t, rep.indx[[j]][[t]]] <- (1/scale.param) * logit.inv(X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha))
          
          p[j, t, rep.indx[[j]][[t]]] <- pmax(p[j, t, rep.indx[[j]][[t]]], 0)
          p[j, t, rep.indx[[j]][[t]]] <- pmin(p[j, t, rep.indx[[j]][[t]]], 1)
        }
        
      } else if (mis.spec.type == "line"){
        if (length(p.RE) > 0) {
          p[j, t, rep.indx[[j]][[t]]] <- X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha) + 
                                                     alpha.star.sites[j, t, rep.indx[[j]][[t]]]
          p[j, t, rep.indx[[j]][[t]]] <- pmax(p[j, t, rep.indx[[j]][[t]]], 0)
          p[j, t, rep.indx[[j]][[t]]] <- pmin(p[j, t, rep.indx[[j]][[t]]], 1)
        } else {
          p[j, t, rep.indx[[j]][[t]]] <- X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha)
          p[j, t, rep.indx[[j]][[t]]] <- pmax(p[j, t, rep.indx[[j]][[t]]], 0)
          p[j, t, rep.indx[[j]][[t]]] <- pmin(p[j, t, rep.indx[[j]][[t]]], 1)
        }
      } else if (mis.spec.type == 'probit') {
        if (length(p.RE) > 0) {
          p[j, t, rep.indx[[j]][[t]]] <- pnorm(X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha) + 
                                                     alpha.star.sites[j, t, rep.indx[[j]][[t]]]) 
        } else {
          p[j, t, rep.indx[[j]][[t]]] <- pnorm(X.p[j, t, rep.indx[[j]][[t]], ] %*% as.matrix(alpha))
        }
      }
      # Allow for closure to be violated if specified
      # tmp <- rbinom(n.rep[j, t], 1, avail[j, t, rep.indx[[j]][[t]]]) 
      # y[j, t, rep.indx[[j]][[t]]] <- rbinom(n.rep[j, t], 1, p[j, t, rep.indx[[j]][[t]]] * tmp * z[j, t]) 
      y[j, t, rep.indx[[j]][[t]]] <- rbinom(n.rep[j, t], 1, p[j, t, rep.indx[[j]][[t]]] * avail[j, t, rep.indx[[j]][[t]]] * z[j, t]) 
    } # t
  } # j

  return(
    list(X = X, X.p = X.p, coords = coords, coords.full = coords.full,
	 psi = psi, z = z, p = p, y = y, w = w.mat.full, w.grid = w.mat, 
	 X.p.re = X.p.re, X.re = X.re,
	 alpha.star = alpha.star, beta.star = beta.star, eta = eta)
  )
}
