spIntPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
		       tuning, cov.model = "exponential", NNGP = TRUE, 
		       n.neighbors = 15, search.type = "cb",
		       n.batch, batch.length, accept.rate = 0.43, 
		       n.omp.threads = 1, verbose = TRUE,  
		       n.report = 100, 
		       n.burn = round(.10 * n.batch * batch.length),
		       n.thin = 1, n.chains = 1, k.fold, k.fold.threads = 1, 
		       k.fold.seed = 100, k.fold.data, k.fold.only = FALSE, ...){

  ptm <- proc.time()

  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
  rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
  }

  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing to run the model\n");
    cat("----------------------------------------\n");
  }
  # Check for unused arguments ------------------------------------------	
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  # Returns a call in which all of the specified arguments are 
  # specified by their full names. 
  cl <- match.call()

  # Some initial checks -------------------------------------------------
  if (missing(data)) {
    stop("error: data must be specified")
  }
  if (!is.list(data)) {
    stop("error: data must be a list")
  }
  names(data) <- tolower(names(data))
  if (missing(occ.formula)) {
    stop("error: occ.formula must be specified")
  }
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  if (!is.list(data$y)) {
    stop("error: y must be a list of detection-nondetection data sets")
  }
  y <- data$y
  n.data <- length(y)
  # Check if individual data sets are provided as vectors for nonreplicated data, 
  # then convert all to matrices to allow for data frames.
  for (q in 1:n.data) {
    if (is.null(dim(y[[q]]))) {
      message(paste("Data source ", q, " is provided as a one-dimensional vector.\nAssuming this is a nonreplicated detection-nondetection data source.\n", sep = ''))
    }
    y[[q]] <- as.matrix(y[[q]])
  }
  if (!'sites' %in% names(data)) {
    stop("error: site ids must be specified in data")
  }
  sites <- data$sites
  # Number of sites with at least one data source
  J <- length(unique(unlist(sites)))
  # Number of sites for each data set
  J.long <- sapply(y, function(a) dim(a)[[1]])

  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- matrix(1, J, 1)
    } else {
      stop("error: occ.covs must be specified in data for an occupancy model with covariates")
    }
  }

  if (!'det.covs' %in% names(data)) {
    data$det.covs <- list()
    for (i in 1:n.data) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model for each data source.\n")
      }
      det.formula.curr <- det.formula[[i]]
      if (det.formula.curr == ~ 1) {
        for (i in 1:n.data) {
          data$det.covs[[i]] <- list(int = matrix(1, dim(y[[i]])[1], dim(y[[i]])[2]))
        }
      } else {
        stop("error: det.covs must be specified in data for a detection model with covariates")
      }
    }
  }

  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  coords <- as.matrix(data$coords)
  if (!missing(k.fold)) {
    if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
      stop("error: k.fold must be a single integer value >= 2")  
    }
  }

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    # Order by x column. Could potentailly allow this to be user defined.
    ord <- order(coords[,1])
    # Reorder everything to align with NN ordering. 
    coords <- coords[ord, ]
    data$occ.covs <- data$occ.covs[ord, , drop = FALSE]
    rownames(data$occ.covs) <- 1:nrow(data$occ.covs)
    sites.orig <- sites
    # Don't need to actually reorder the data, can just reorder the site 
    # indices. 
    for (i in 1:n.data) {
      for (j in 1:length(sites[[i]])) {
      sites[[i]][j] <- which(ord == sites[[i]][j])
      }
    }
  }

  # Check if all detection covariates are at site level for a given
  # data set.
  binom <- rep(FALSE, n.data)
  # Make all covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  for (i in 1:n.data) {
    # Get in required R model format. 
    data$det.covs[[i]] <- data.frame(lapply(data$det.covs[[i]], function(a) unlist(c(a))))
    # Replicate det.covs if only covariates are at the site level. 
    if (nrow(data$det.covs[[i]]) == nrow(y[[i]])) {
      binom[i] <- TRUE
      # Check if there are missing site-level covariates 
      if (sum(is.na(data$det.covs[[i]])) != 0) {
        stop("error: missing values in site-level det.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
      }
      data$det.covs[[i]] <- data.frame(sapply(data$det.covs[[i]], rep,
      					times = dim(y[[i]])[2]))
    }
  }
  data$occ.covs <- as.data.frame(data$occ.covs)

  # Checking missing values ---------------------------------------------
  # y -------------------------------
  for (q in 1:n.data) {
    y.na.test <- apply(y[[q]], 1, function(a) sum(!is.na(a)))
    if (sum(y.na.test == 0) > 0) {
      stop(paste("error: some sites in data source ", q, " in y have all missing detection histories.\n Remove these sites from y and all objects in the 'data' argument if the site is not surveyed by another data source\n, then use 'predict' to obtain predictions at these locations if desired.", sep = ''))
    }
  }
  # occ.covs ------------------------
  if (sum(is.na(data$occ.covs)) != 0) {
    stop("error: missing values in occ.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
  }
  # det.covs ------------------------
  for (q in 1:n.data) {
    if (!binom[q]) {
      for (i in 1:ncol(data$det.covs[[q]])) {
        if (sum(is.na(data$det.covs[[q]][, i])) > sum(is.na(y[[q]]))) {
          stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
        }
      }
      # Misalignment between y and det.covs
      y.missing <- which(is.na(y[[q]]))
      det.covs.missing <- lapply(data$det.covs[[q]], function(a) which(is.na(a)))
      for (i in 1:length(det.covs.missing)) {
        tmp.indx <- !(y.missing %in% det.covs.missing[[i]])
        if (sum(tmp.indx) > 0) {
          if (i == 1 & verbose) {
            message("There are missing values in data$y with corresponding non-missing values in data$det.covs.\nRemoving these site/replicate combinations for fitting the model.")
          }
          data$det.covs[[q]][y.missing, i] <- NA
        }
      }
    }
  }

  # Formula -------------------------------------------------------------
  # Occupancy -----------------------

  if (is(occ.formula, 'formula')) {
    tmp <- parseFormula(occ.formula, data$occ.covs)
    X <- as.matrix(tmp[[1]])
    x.names <- tmp[[2]]
  } else {
    stop("error: occ.formula is misspecified")
  }

  # Detection -----------------------
  if (!is.list(det.formula)) {
    stop(paste("error: det.formula must be a list of ", n.data, " formulas", sep = ''))
  }
  X.p <- list()
  x.p.names <- list()
  for (i in 1:n.data) {
    if (is(det.formula[[i]], 'formula')) {
      tmp <- parseFormula(det.formula[[i]], data$det.covs[[i]])
      X.p[[i]] <- as.matrix(tmp[[1]])
      x.p.names[[i]] <- tmp[[2]]
    } else {
      stop(paste("error: det.formula for data source ", i, " is misspecified", sep = ''))
    }
  }
  x.p.names <- unlist(x.p.names)

    # Get basic info from inputs ------------------------------------------
    # Total number of sites
    J.all <- nrow(X)
    if (length(X.p) != n.data | length(y) != n.data) {
      stop(paste("error: y and X.p must be lists of length ", n.data, ".", sep = ''))
    }
    # Number of occupancy parameters 
    p.occ <- ncol(X)
    # Number of detection parameters for each data set
    p.det.long <- sapply(X.p, function(a) dim(a)[[2]])
    # Total number of detection parameters
    p.det <- sum(p.det.long)
    n.rep <- lapply(y, function(a1) apply(a1, 1, function(a2) sum(!is.na(a2))))
    # Max number of repeat visits for each data set
    K.long.max <- sapply(y, function(a) dim(a)[2])
    # Number of repeat visits for each data set site. 
    K <- unlist(n.rep)
    n.samples <- n.batch * batch.length
    if (missing(n.samples)) {
      stop("error: must specify number of MCMC samples")
    }
    if (n.burn > n.samples) {
      stop("error: n.burn must be less than n.samples")
    }
    if (n.thin > n.samples) {
      stop("error: n.thin must be less than n.samples")
    }
    if (!missing(k.fold)) {
      if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
        stop("error: k.fold must be a single integer value >= 2")  
      }
    }

    # Get indics to map z to y --------------------------------------------
    y.big <- y
    names.long <- list()
    # Remove missing observations when the covariate data are available but
    # there are missing detection-nondetection data
    for (i in 1:n.data) {
      if (nrow(X.p[[i]]) == length(y[[i]])) {
        X.p[[i]] <- X.p[[i]][!is.na(y[[i]]), , drop = FALSE]
      }
      # Need these for later on
      names.long[[i]] <- which(!is.na(y[[i]]))
    }
    n.obs.long <- sapply(X.p, nrow)
    n.obs <- sum(n.obs.long)
    z.long.indx.r <- list()
    for (i in 1:n.data) {
      z.long.indx.r[[i]] <- rep(sites[[i]], K.long.max[i])
      z.long.indx.r[[i]] <- z.long.indx.r[[i]][!is.na(c(y[[i]]))]
    }
    z.long.indx.r <- unlist(z.long.indx.r)
    # Subtract 1 for c indices
    z.long.indx.c <- z.long.indx.r - 1
    y <- unlist(y)
    y <- y[!is.na(y)]
    # Index indicating the data set associated with each data point in y
    data.indx.r <- rep(NA, n.obs)
    indx <- 1
    for (i in 1:n.data) {
      data.indx.r[indx:(indx + n.obs.long[i] - 1)] <- rep(i, n.obs.long[i])
      indx <- indx + n.obs.long[i]
    }
    data.indx.c <- data.indx.r - 1

    X.p.all <- matrix(NA, n.obs, max(p.det.long))
    indx <- 1
    for (i in 1:n.data) {
      X.p.all[indx:(indx + nrow(X.p[[i]]) - 1), 1:p.det.long[i]] <- X.p[[i]] 
      indx <- indx + nrow(X.p[[i]])
    }

    # Priors --------------------------------------------------------------
    if (missing(priors)) {
      priors <- list()
    }
    names(priors) <- tolower(names(priors))
    # beta -----------------------
    if ("beta.normal" %in% names(priors)) {
      if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
        stop("error: beta.normal must be a list of length 2")
      }
      mu.beta <- priors$beta.normal[[1]]
      sigma.beta <- priors$beta.normal[[2]]
      if (length(mu.beta) != p.occ & length(mu.beta) != 1) {
        if (p.occ == 1) {
          stop(paste("error: beta.normal[[1]] must be a vector of length ",
          	     p.occ, " with elements corresponding to betas' mean", sep = ""))
        } else {
          stop(paste("error: beta.normal[[1]] must be a vector of length ",
          	     p.occ, " or 1 with elements corresponding to betas' mean", sep = ""))
        }
      }
      if (length(sigma.beta) != p.occ & length(sigma.beta) != 1) {
        if (p.occ == 1) {
          stop(paste("error: beta.normal[[2]] must be a vector of length ",
        	   p.occ, " with elements corresponding to betas' variance", sep = ""))
        } else {
          stop(paste("error: beta.normal[[2]] must be a vector of length ",
        	   p.occ, " or 1 with elements corresponding to betas' variance", sep = ""))
        }
      }
      if (length(sigma.beta) != p.occ) {
        sigma.beta <- rep(sigma.beta, p.occ)
      }
      if (length(mu.beta) != p.occ) {
        mu.beta <- rep(mu.beta, p.occ)
      }
      Sigma.beta <- sigma.beta * diag(p.occ)
    } else {
      if (verbose) {
        message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
      }
      mu.beta <- rep(0, p.occ)
      sigma.beta <- rep(2.72, p.occ)
      Sigma.beta <- diag(p.occ) * 2.72
    }

    # alpha -----------------------
    if ("alpha.normal" %in% names(priors)) {
      if (!is.list(priors$alpha.normal) | length(priors$alpha.normal) != 2) {
        stop("error: alpha.normal must be a list of length 2")
      }
      mu.alpha <- priors$alpha.normal[[1]]
      sigma.alpha <- priors$alpha.normal[[2]]
      if (length(mu.alpha) != n.data | !is.list(mu.alpha)) {
        stop(paste("error: alpha.normal[[1]] must be a list of length ", 
		   n.data, " with elements corresponding to alphas' mean for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(mu.alpha[[q]]) != p.det.long[q] & length(mu.alpha[[q]]) != 1) {
          if (p.det.long[q] == 1) {
            stop(paste("error: prior means for alpha.normal[[1]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
            stop(paste("error: prior means for alpha.normal[[1]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], "or 1", sep = ""))
          }
        }
        if (length(mu.alpha[[q]]) != p.det.long[q]) {
          mu.alpha[[q]] <- rep(mu.alpha[[q]], p.det.long[q])
        }
      }
      mu.alpha <- unlist(mu.alpha)
      if (length(sigma.alpha) != n.data | !is.list(sigma.alpha)) {
        stop(paste("error: alpha.normal[[2]] must be a list of length ", 
		   n.data, " with elements corresponding to alphas' variance for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(sigma.alpha[[q]]) != p.det.long[q] & length(sigma.alpha[[q]]) != 1) {
          if (p.det.long[q] == 1) {
          stop(paste("error: prior variances for alpha.normal[[2]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
          stop(paste("error: prior variances for alpha.normal[[2]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], " or 1", sep = ""))

          }
        }
      if (length(sigma.alpha[[q]]) != p.det.long[q]) {
        sigma.alpha[[q]] <- rep(sigma.alpha[[q]], p.det.long[q])
      }
      }
      sigma.alpha <- unlist(sigma.alpha)
    } else {
      if (verbose) {
        message("No prior specified for alpha.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
      }
      mu.alpha <- rep(0, p.det)
      sigma.alpha <- rep(2.72, p.det) 
    }

    # phi -----------------------------
    # Get distance matrix which is used if priors are not specified
    if (!NNGP) {
      coords.D <- iDist(coords)
    }
    if ("phi.unif" %in% names(priors)) {
      if (!is.vector(priors$phi.unif) | !is.atomic(priors$phi.unif) | length(priors$phi.unif) != 2) {
        stop("error: phi.unif must be a vector of length 2 with elements corresponding to phi's lower and upper bounds")
      }
      phi.a <- priors$phi.unif[1]
      phi.b <- priors$phi.unif[2]
    } else {
      if (verbose) {
        message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
      }
      if (NNGP) {
        coords.D <- iDist(coords)
      }
      phi.a <- 3 / max(coords.D)
      phi.b <- 3 / sort(unique(c(coords.D)))[2]
    }

    # sigma.sq -----------------------------
    fixed.sigma.sq <- FALSE
    # Check if both an ig and uniform prior are specified
    if (("sigma.sq.ig" %in% names(priors)) & ("sigma.sq.unif" %in% names(priors))) {
      stop("error: cannot specify both an IG and a uniform prior for sigma.sq")
    }
    if ("sigma.sq.ig" %in% names(priors)) { # inverse-gamma prior.
      sigma.sq.ig <- TRUE
      if (priors$sigma.sq.ig[1] == 'fixed') {
        fixed.sigma.sq <- TRUE
        sigma.sq.a <- 1
        sigma.sq.b <- 1
      } else {
        if (!is.vector(priors$sigma.sq.ig) | !is.atomic(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
          stop("error: sigma.sq.ig must be a vector of length 2 with elements corresponding to sigma.sq's shape and scale parameters")
        }
        sigma.sq.a <- priors$sigma.sq.ig[1]
        sigma.sq.b <- priors$sigma.sq.ig[2]
      }
    } else if ("sigma.sq.unif" %in% names(priors)) { # uniform prior
      if (priors$sigma.sq.unif[1] == 'fixed') {
        sigma.sq.ig <- TRUE # This just makes the C++ side a bit easier. 
        fixed.sigma.sq <- TRUE 
        sigma.sq.a <- 1
        sigma.sq.b <- 1
      } else {
        sigma.sq.ig <- FALSE
        if (!is.vector(priors$sigma.sq.unif) | !is.atomic(priors$sigma.sq.unif) | length(priors$sigma.sq.unif) != 2) {
          stop("error: sigma.sq.unif must be a vector of length 2 with elements corresponding to sigma.sq's lower and upper bounds")
        }
        sigma.sq.a <- priors$sigma.sq.unif[1]
        sigma.sq.b <- priors$sigma.sq.unif[2]
      }
       
    } else {
      if (verbose) {
        message("No prior specified for sigma.sq.\nUsing an inverse-Gamma prior with the shape parameter set to 2 and scale parameter to 1.\n")
      }
      sigma.sq.ig <- TRUE
      sigma.sq.a <- 2
      sigma.sq.b <- 1
    }
    # nu -----------------------------
    if (cov.model == 'matern') {
      if (!"nu.unif" %in% names(priors)) {
        stop("error: nu.unif must be specified in priors value list")
      }
      if (!is.vector(priors$nu.unif) | !is.atomic(priors$nu.unif) | length(priors$nu.unif) != 2) {
        stop("error: nu.unif must be a vector of length 2 with elements corresponding to nu's lower and upper bounds")
      }
      nu.a <- priors$nu.unif[1]
      nu.b <- priors$nu.unif[2]
    } else {
      nu.a <- 0
      nu.b <- 0
    }

    # Starting values -----------------------------------------------------
    if (missing(inits)) {
      inits <- list()
    }
    names(inits) <- tolower(names(inits))
    # z -------------------------------
    if ("z" %in% names(inits)) {
      z.inits <- inits$z
      if (!is.vector(z.inits)) {
        stop(paste("error: initial values for z must be a vector of length ",
		   J, sep = ""))
      }
      if (length(z.inits) != J) {
        stop(paste("error: initial values for z must be a vector of length ",
		   J, sep = ""))
      }
      z.test <- tapply(y, z.long.indx.r, max, na.rm = TRUE)	
      init.test <- sum(z.inits < z.test)
      if (init.test > 0) {
        stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site.")
      }
    } else {
      z.inits <- tapply(y, z.long.indx.r, max, na.rm = TRUE)
      if (verbose) {
        message("z is not specified in initial values.\nSetting initial values based on observed data\n")
      }
    }
    # beta -----------------------
    if ("beta" %in% names(inits)) {
      beta.inits <- inits[["beta"]]
      if (length(beta.inits) != p.occ & length(beta.inits) != 1) {
        if (p.occ == 1) {
          stop(paste("error: initial values for beta must be of length ", p.occ,
		     sep = ""))

        } else {
          stop(paste("error: initial values for beta must be of length ", p.occ, " or 1",
	  	     sep = ""))
        }
      }
      if (length(beta.inits) != p.occ) {
        beta.inits <- rep(beta.inits, p.occ)
      }
    } else {
      beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
      if (verbose) {
        message('beta is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
      }
    }
    # alpha -----------------------
    if ("alpha" %in% names(inits)) {
      alpha.inits <- inits[["alpha"]]
      if (length(alpha.inits) != n.data | !is.list(alpha.inits)) {
        stop(paste("error: initial values for alpha must be a list of length ", n.data,
		   sep = ""))
      }
      for (q in 1:n.data) {
        if (length(alpha.inits[[q]]) != p.det.long[q] & length(alpha.inits[[q]]) != 1) {
          if (p.det.long[q] == 1) {
            stop(paste("error: initial values for alpha[[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
            stop(paste("error: initial values for alpha[[", q, "]] must be a vector of length ", 
		       p.det.long[q], " or 1", sep = ""))
          }
        }
        if (length(alpha.inits[[q]]) != p.det.long[q]) {
          alpha.inits[[q]] <- rep(alpha.inits, p.det.long[q])
        }
      }
      alpha.inits <- unlist(alpha.inits)
    } else {
      if (verbose) {
        message("alpha is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
    }

    alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
    alpha.indx.c <- alpha.indx.r - 1

    # phi -----------------------------
    if ("phi" %in% names(inits)) {
      phi.inits <- inits[["phi"]]
      if (length(phi.inits) != 1) {
        stop("error: initial values for phi must be of length 1")
      }
    } else {
      phi.inits <- runif(1, phi.a, phi.b)
      if (verbose) {
        message("phi is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
      }
    }

    # sigma.sq ------------------------
    if ("sigma.sq" %in% names(inits)) {
      sigma.sq.inits <- inits[["sigma.sq"]]
      if (length(sigma.sq.inits) != 1) {
        stop("error: initial values for sigma.sq must be of length 1")
      }
    } else {
      if (sigma.sq.ig) {
        sigma.sq.inits <- rigamma(1, sigma.sq.a, sigma.sq.b)
      } else {
        sigma.sq.inits <- runif(1, sigma.sq.a, sigma.sq.b)
      }
      if (verbose) {
        message("sigma.sq is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
      }
    }
    # w -----------------------------
    if ("w" %in% names(inits)) {
      w.inits <- inits[["w"]]
      if (!is.vector(w.inits)) {
        stop(paste("error: initial values for w must be a vector of length ",
        	   J, sep = ""))
      }
      if (length(w.inits) != J) {
        stop(paste("error: initial values for w must be a vector of length ",
        	   J, sep = ""))
      }
    } else {
      w.inits <- rep(0, J)
      if (verbose) {
        message("w is not specified in initial values.\nSetting initial value to 0\n")
      }
    }
    # nu ------------------------
    if ("nu" %in% names(inits)) {
      nu.inits <- inits[["nu"]]
      if (length(nu.inits) != 1) {
        stop("error: initial values for nu must be of length 1")
      }
    } else {
      if (cov.model == 'matern') {
        if (verbose) {
          message("nu is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
        }
        nu.inits <- runif(1, nu.a, nu.b)
      } else {
        nu.inits <- 0
      }
    }

    # Should initial values be fixed --
    if ("fix" %in% names(inits)) {
      fix.inits <- inits[["fix"]]
      if ((fix.inits != TRUE) & (fix.inits != FALSE)) {
        stop(paste("error: inits$fix must take value TRUE or FALSE"))
      }
    } else {
      fix.inits <- FALSE
    }
    if (verbose & fix.inits & (n.chains > 1)) {
      message("Fixing initial values across all chains\n")
    }

    # Covariance Model ----------------------------------------------------
    # Order must match util.cpp spCor.
    cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
    if(! cov.model %in% cov.model.names){
      stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ", 
	   paste(cov.model.names, collapse=", ", sep="") ,".")}
    # Obo for cov model lookup on c side
    cov.model.indx <- which(cov.model == cov.model.names) - 1
   
    # Get tuning values ---------------------------------------------------
    # Not accessed, but necessary to keep things in line. 
    sigma.sq.tuning <- 0
    phi.tuning <- 0
    nu.tuning <- 0
    if (missing(tuning)) {
      phi.tuning <- 1
      if (cov.model == 'matern') {
        nu.tuning <- 1
      }
      if (!sigma.sq.ig) {
        sigma.sq.tuning <- 1
      }
    } else {
      names(tuning) <- tolower(names(tuning))
      # phi ---------------------------
      if(!"phi" %in% names(tuning)) {
        stop("error: phi must be specified in tuning value list")
      }
      phi.tuning <- tuning$phi
      if (length(phi.tuning) != 1) {
        stop("error: phi tuning must be a single value")
      } 
      if (cov.model == 'matern') {
        # nu --------------------------
        if(!"nu" %in% names(tuning)) {
          stop("error: nu must be specified in tuning value list")
        }
        nu.tuning <- tuning$nu
        if (length(nu.tuning) != 1) {
          stop("error: nu tuning must be a single value")
        } 
      }
      if (!sigma.sq.ig) {
        # sigma.sq --------------------------
        if(!"sigma.sq" %in% names(tuning)) {
          stop("error: sigma.sq must be specified in tuning value list")
        }
        sigma.sq.tuning <- tuning$sigma.sq
        if (length(sigma.sq.tuning) != 1) {
          stop("error: sigma.sq tuning must be a single value")
        } 
      }
    }
    tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))

    # Set model.deviance to NA for returning when no cross-validation
    model.deviance <- NA
    curr.chain <- 1

  if (!NNGP) { 
    # Specify storage modes -----------------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
    storage.mode(coords.D) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.det.long) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(n.obs) <- "integer"
    storage.mode(n.obs.long) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(J.long) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(n.data) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(sigma.alpha) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(sigma.sq.a) <- "double"
    storage.mode(sigma.sq.b) <- "double"
    storage.mode(sigma.sq.ig) <- "integer"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(tuning.c) <- "double"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(n.burn) <- "integer"
    storage.mode(n.thin) <- "integer"
    storage.mode(curr.chain) <- "integer"
    storage.mode(n.chains) <- "integer"
    storage.mode(fixed.sigma.sq) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"

    out.tmp <- list()
    out <- list()
    if (!k.fold.only) {
      for (i in 1:n.chains) {
        # Change initial values if i > 1
        if ((i > 1) & (!fix.inits)) {
          beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
          alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
          if (!fixed.sigma.sq) {
            if (sigma.sq.ig) {
              sigma.sq.inits <- rigamma(1, sigma.sq.a, sigma.sq.b)
            } else {
              sigma.sq.inits <- runif(1, sigma.sq.a, sigma.sq.b)
            }
          }
          phi.inits <- runif(1, phi.a, phi.b)
          if (cov.model == 'matern') {
            nu.inits <- runif(1, nu.a, nu.b)
          }
        }
        storage.mode(curr.chain) <- "integer" 
        # Run the model in C
        out.tmp[[i]] <- .Call("spIntPGOcc", y, X, X.p.all, coords.D, p.occ, p.det, p.det.long, 
          	            J, J.long, K, n.obs, n.obs.long, n.data, 
          	            beta.inits, alpha.inits, z.inits, w.inits, 
          	            phi.inits, sigma.sq.inits, nu.inits, 
          	            z.long.indx.c, data.indx.c, alpha.indx.c, mu.beta, mu.alpha, 
          	            Sigma.beta, sigma.alpha, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	            nu.a, nu.b, tuning.c, cov.model.indx, 
          	            n.batch, batch.length, accept.rate,  
          	            n.omp.threads, verbose, n.report, n.burn, n.thin, n.post.samples, 
          		    curr.chain, n.chains, fixed.sigma.sq, sigma.sq.ig)
        curr.chain <- curr.chain + 1
      }
      # Calculate R-Hat ---------------
      out <- list()
      out$rhat <- list()
      if (n.chains > 1) {
        out$rhat$beta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2]
        out$rhat$alpha <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.samples)))), 
        			      autoburnin = FALSE)$psrf[, 2]
        if (fixed.sigma.sq) {
          out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$theta.samples)))), 
          			      autoburnin = FALSE)$psrf[, 2]
        } else {
          out$rhat$theta <- c(NA, gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$theta.samples[-1, , drop = FALSE])))), 
          			      autoburnin = FALSE)$psrf[, 2])
        }

      } else {
        out$rhat$beta <- rep(NA, p.occ)
        out$rhat$alpha <- rep(NA, p.det)
        out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3, 2))
      }

      out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(do.call(rbind, 
        				lapply(out.tmp, function(a) t(a$alpha.samples))))
      colnames(out$alpha.samples) <- x.p.names
      out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
      if (cov.model != 'matern') {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi')
      } else {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
      }
      out$z.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$z.samples))))
      out$psi.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$psi.samples))))
      out$w.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples))))
      # Calculate effective sample sizes
      out$ESS <- list()
      out$ESS$beta <- effectiveSize(out$beta.samples)
      out$ESS$alpha <- effectiveSize(out$alpha.samples)
      out$ESS$theta <- effectiveSize(out$theta.samples)
      out$X <- X
      out$X.p <- X.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$sites <- sites
      out$cov.model.indx <- cov.model.indx
      out$type <- "GP"
      out$coords <- coords
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$n.chains <- n.chains
    }

    # K-fold cross-validation -------
    if (!missing(k.fold)) {
      if (verbose) {      
        cat("----------------------------------------\n");
        cat("\tCross-validation\n");
        cat("----------------------------------------\n");
        message(paste("Performing ", k.fold, "-fold cross-validation using ", k.fold.threads,
      	      " thread(s).", sep = ''))
      }
      set.seed(k.fold.seed)
      if (missing(k.fold.data)) {
        k.fold.data <- NULL
      }
      # Check to see if only one data source should be used for hold-out evaluations
      if (!is.null(k.fold.data)) {
        if (!is.numeric(k.fold.data) | length(k.fold.data) != 1) {
          stop("error: if specified, k.fold.data must be a single numeric value")
        }
        if (verbose) {
          message(paste("Only holding out data from data source ", k.fold.data, ".", sep = ''))
	}
	sites.random <- sample(sites[[k.fold.data]])
      } else {
        # Number of sites in each hold out data set. 
        sites.random <- sample(1:J)    
      }
      sites.k.fold <- split(sites.random, rep(1:k.fold, length.out = length(sites.random)))
      registerDoParallel(k.fold.threads)
      model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
        curr.set <- sort(sites.k.fold[[i]])
        curr.set.pred <- curr.set
	curr.set.fit <- (1:J)[-curr.set]
	if (!is.null(k.fold.data)) {
          curr.set.fit <- sort(unique(c(curr.set.fit, unlist(sites[-k.fold.data]))))
        }
        y.indx <- !((z.long.indx.r) %in% curr.set)
        if (!is.null(k.fold.data)) {
          y.indx <- ifelse(data.indx.r == k.fold.data, y.indx, TRUE)
        } 
        y.fit <- y[y.indx]
        y.0 <- y[!y.indx]
        z.inits.fit <- z.inits[curr.set.fit]
	w.inits.fit <- w.inits[curr.set.fit]
	coords.fit <- coords[curr.set.fit, , drop = FALSE]
	coords.0 <- coords[curr.set.pred, , drop = FALSE]
	coords.D.fit <- coords.D[curr.set.fit, curr.set.fit, drop = FALSE]
	coords.D.0 <- coords.D[curr.set.pred, curr.set.pred, drop = FALSE]
        X.p.fit <- X.p.all[y.indx, , drop = FALSE]
        X.p.0 <- X.p.all[!y.indx, , drop = FALSE]
        X.fit <- X[curr.set.fit, , drop = FALSE]
        X.0 <- X[curr.set.pred, , drop = FALSE]
        J.fit <- nrow(X.fit)
	# Site indices for fitted data
	sites.fit <- sapply(sites, 
			    function(a) which(as.numeric(row.names(X.fit)) %in% a[a %in% curr.set.fit]))
	tmp <- sapply(sites, function(a) a %in% curr.set.fit)
	if (!is.null(k.fold.data)) {
          sites.fit[[k.fold.data]] <- which(as.numeric(row.names(X.fit)) %in% sites[[k.fold.data]][sites[[k.fold.data]] %in% (1:J)[-curr.set]])
          tmp[[k.fold.data]] <- sites[[k.fold.data]] %in% (1:J)[-curr.set]
        }
        K.fit <- K[unlist(tmp)]

        K.fit <- K[unlist(tmp)]
	y.big.fit <- y.big
	for (q in 1:n.data) {
          for (j in 1:J.long[q]) {
            if (is.null(k.fold.data)) {
              if (sites[[q]][j] %in% curr.set) {
                y.big.fit[[q]][j, ] <- NA
	      }
	    } else {
              if (q == k.fold.data) {
                if (sites[[q]][j] %in% curr.set) {
                  y.big.fit[[q]][j, ] <- NA
	        }
	      }
	    }
	  }
	  y.big.fit[[q]] <- y.big.fit[[q]][which(apply(y.big.fit[[q]], 
						       1, function(a) sum(!is.na(a)) > 0)), ]
	}
	y.big.0 <- y.big
	for (q in 1:n.data) {
          for (j in 1:J.long[q]) {
            if (is.null(k.fold.data)) {
              if (! sites[[q]][j] %in% curr.set) {
                y.big.0[[q]][j, ] <- NA
	      }
	    } else {
              if (q == k.fold.data) {
                if (! sites[[q]][j] %in% curr.set) {
                  y.big.0[[q]][j, ] <- NA
	        }
	      }
	    }
	  }
	  y.big.0[[q]] <- y.big.0[[q]][which(apply(y.big.0[[q]], 
						       1, function(a) sum(!is.na(a)) > 0)), ]
	}

	z.long.indx.fit <- list()
	for (q in 1:n.data) {
          z.long.indx.fit[[q]] <- rep(sites.fit[[q]], K.long.max[q])
	  z.long.indx.fit[[q]] <- z.long.indx.fit[[q]][!is.na(c(y.big.fit[[q]]))]
	}
	z.long.indx.fit <- unlist(z.long.indx.fit)
	z.long.indx.fit <- z.long.indx.fit - 1
	
	# Site indices for hold out data
	sites.0 <- sapply(sites, 
			    function(a) which(as.numeric(row.names(X.0)) %in% a[a %in% curr.set.pred]))
	tmp <- sapply(sites, function(a) a %in% curr.set.pred)
	if (!is.null(k.fold.data)) {
          sites.0[-k.fold.data] <- NA
	  for (q in 1:n.data) {
            if (q != k.fold.data) {
              tmp[[q]] <- rep(FALSE, J.long[q])
            }
          }
        }
        K.0 <- K[unlist(tmp)]
	if (is.null(k.fold.data)) {
	  z.long.indx.0 <- list()
	  for (q in 1:n.data) {
            z.long.indx.0[[q]] <- rep(sites.0[[q]], K.long.max[q])
	    z.long.indx.0[[q]] <- z.long.indx.0[[q]][!is.na(c(y.big.0[[q]]))]
	  }
	  z.long.indx.0 <- unlist(z.long.indx.0) - 1
	} else {
          z.long.indx.0 <- rep(sites.0[[k.fold.data]], K.long.max[k.fold.data])
	  z.long.indx.0 <- z.long.indx.0[!is.na(c(y.big.0[[k.fold.data]]))] - 1
	}
        verbose.fit <- FALSE
        n.omp.threads.fit <- 1
	n.obs.fit <- length(y.fit)
	n.obs.0 <- length(y.0)
	data.indx.c.fit <- data.indx.c[y.indx]
	data.indx.0 <- data.indx.c[!y.indx] + 1
	n.obs.long.fit <- as.vector(table(data.indx.c.fit))
	n.obs.long.0 <- n.obs.long - n.obs.long.fit
	J.long.fit <- as.vector(tapply(z.long.indx.fit, factor(data.indx.c.fit), 
			     FUN = function(a) length(unique(a))))


        storage.mode(y.fit) <- "double"
        storage.mode(z.inits.fit) <- "double"
        storage.mode(X.p.fit) <- "double"
        storage.mode(X.fit) <- "double"
        storage.mode(coords.D.fit) <- "double"
        storage.mode(p.det) <- "integer"
        storage.mode(p.det.long) <- "integer"
        storage.mode(p.occ) <- "integer"
        storage.mode(n.obs.fit) <- "integer"
        storage.mode(n.obs.long.fit) <- "integer"
        storage.mode(J.fit) <- "integer"
        storage.mode(J.long.fit) <- "integer"
        storage.mode(K.fit) <- "integer"
        storage.mode(n.data) <- "integer"
        storage.mode(beta.inits) <- "double"
        storage.mode(alpha.inits) <- "double"
        storage.mode(phi.inits) <- "double"
        storage.mode(nu.inits) <- "double"
        storage.mode(w.inits.fit) <- "double"
        storage.mode(sigma.sq.inits) <- "double"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(data.indx.c.fit) <- "integer"
        storage.mode(alpha.indx.c) <- "integer"
        storage.mode(mu.beta) <- "double"
        storage.mode(Sigma.beta) <- "double"
        storage.mode(mu.alpha) <- "double"
        storage.mode(sigma.alpha) <- "double"
        storage.mode(phi.a) <- "double"
        storage.mode(phi.b) <- "double"
        storage.mode(sigma.sq.a) <- "double"
        storage.mode(sigma.sq.b) <- "double"
        storage.mode(nu.a) <- "double"
        storage.mode(nu.b) <- "double"
        storage.mode(tuning.c) <- "double"
        storage.mode(cov.model.indx) <- "integer"
        storage.mode(n.batch) <- "integer"
        storage.mode(batch.length) <- "integer"
        storage.mode(accept.rate) <- "double"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(n.burn) <- "integer"
        storage.mode(n.thin) <- "integer"
	curr.chain <- 1
	storage.mode(curr.chain) <- "integer"

        out.fit <- .Call("spIntPGOcc", y.fit, X.fit, X.p.fit, coords.D.fit, p.occ, p.det, p.det.long, 
		         J.fit, J.long.fit, K.fit, n.obs.fit, n.obs.long.fit, n.data, 
		         beta.inits, alpha.inits, z.inits.fit, w.inits.fit, 
		         phi.inits, sigma.sq.inits, nu.inits, 
		         z.long.indx.fit, data.indx.c.fit, alpha.indx.c, mu.beta, mu.alpha, 
		         Sigma.beta, sigma.alpha, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		         nu.a, nu.b, tuning.c, cov.model.indx, 
		         n.batch, batch.length, accept.rate,  
		         n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, n.post.samples, 
			 curr.chain, n.chains, fixed.sigma.sq, sigma.sq.ig)

        out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
        colnames(out.fit$beta.samples) <- x.names
        out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
        if (cov.model != 'matern') {
          colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi')
        } else {
          colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi', 'nu')
        }
        out.fit$w.samples <- mcmc(t(out.fit$w.samples))
	out.fit$z.samples <- mcmc(t(out.fit$z.samples))
	out.fit$psi.samples <- mcmc(t(out.fit$psi.samples))
        out.fit$X <- X.fit
        out.fit$X.p <- X.p.fit
        out.fit$call <- cl
        out.fit$n.samples <- batch.length * n.batch
        out.fit$cov.model.indx <- cov.model.indx
        out.fit$type <- "GP"
        out.fit$coords <- coords.fit
        out.fit$n.post <- n.post.samples
        out.fit$n.thin <- n.thin
        out.fit$n.burn <- n.burn
	out.fit$n.chains <- 1
	out.fit$pRE <- FALSE
	out.fit$psiRE <- FALSE
	class(out.fit) <- "spPGOcc"

        out.pred <- predict.spPGOcc(out.fit, X.0, coords.0, verbose = FALSE)
        # Detection 
	p.0.samples <- matrix(NA, n.post.samples, nrow(X.p.0))
        like.samples <- rep(NA, nrow(X.p.0))
        for (j in 1:nrow(X.p.0)) {
          p.0.samples[, j] <- logit.inv(X.p.0[j, 1:sum(alpha.indx.r == data.indx.0[j])] %*% out.fit$alpha.samples[which(alpha.indx.r == data.indx.0[j]), ])
          like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * out.pred$z.0.samples[, z.long.indx.0[j] + 1]))
        }
	as.vector(tapply(like.samples, data.indx.0, function(a) sum(log(a))))
      }
      model.deviance <- -2 * model.deviance
      # Return objects from cross-validation
      out$k.fold.deviance <- model.deviance
      stopImplicitCluster()
    }
  } else {

    # Nearest Neighbor Search ---------------------------------------------
    if(verbose){
      cat("----------------------------------------\n");
      cat("\tBuilding the neighbor list\n");
      cat("----------------------------------------\n");
    }

    search.type.names <- c("brute", "cb")
    
    if(!search.type %in% search.type.names){
      stop("error: specified search.type '",search.type,
	   "' is not a valid option; choose from ", 
	   paste(search.type.names, collapse=", ", sep="") ,".")
    }
    
    ## Indexes
    if(search.type == "brute"){
      indx <- mkNNIndx(coords, n.neighbors, n.omp.threads)
    } else{
      indx <- mkNNIndxCB(coords, n.neighbors, n.omp.threads)
    }
    
    nn.indx <- indx$nnIndx
    nn.indx.lu <- indx$nnIndxLU
    nn.indx.run.time <- indx$run.time
    

    if(verbose){
      cat("----------------------------------------\n");
      cat("Building the neighbors of neighbors list\n");
      cat("----------------------------------------\n");
    }
    
    indx <- mkUIndx(J, n.neighbors, nn.indx, nn.indx.lu, u.search.type)
    
    u.indx <- indx$u.indx
    u.indx.lu <- indx$u.indx.lu
    ui.indx <- indx$ui.indx
    u.indx.run.time <- indx$run.time

    # Specify storage modes -----------------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.det.long) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(n.obs) <- "integer"
    storage.mode(n.obs.long) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(J.long) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(n.data) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(sigma.alpha) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(sigma.sq.a) <- "double"
    storage.mode(sigma.sq.b) <- "double"
    storage.mode(sigma.sq.ig) <- "integer"
    storage.mode(tuning.c) <- "double"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.indx) <- "integer"
    storage.mode(u.indx.lu) <- "integer"
    storage.mode(ui.indx) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(n.burn) <- "integer"
    storage.mode(n.thin) <- "integer"
    storage.mode(curr.chain) <- "integer"
    storage.mode(n.chains) <- "integer"
    storage.mode(fixed.sigma.sq) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"

    out.tmp <- list()
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
        alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
	if (!fixed.sigma.sq) {
          if (sigma.sq.ig) {
            sigma.sq.inits <- rigamma(1, sigma.sq.a, sigma.sq.b)
	  } else {
            sigma.sq.inits <- runif(1, sigma.sq.a, sigma.sq.b)
	  }
	}
        phi.inits <- runif(1, phi.a, phi.b)
	if (cov.model == 'matern') {
          nu.inits <- runif(1, nu.a, nu.b)
	}
      }
      storage.mode(curr.chain) <- "integer" 
      # Run the model in C
      out.tmp[[i]] <- .Call("spIntPGOccNNGP", y, X, X.p.all, coords, p.occ, p.det, p.det.long, 
		            J, J.long, K, n.obs, n.obs.long, n.data, 
		            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
		            beta.inits, alpha.inits, z.inits, w.inits, 
		            phi.inits, sigma.sq.inits, nu.inits, 
		            z.long.indx.c, data.indx.c, alpha.indx.c, mu.beta, mu.alpha, 
		            Sigma.beta, sigma.alpha, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		            nu.a, nu.b, tuning.c, cov.model.indx,
		            n.batch, batch.length, accept.rate,  
		            n.omp.threads, verbose, n.report, n.burn, n.thin, n.post.samples, 
			    curr.chain, n.chains, fixed.sigma.sq, sigma.sq.ig)
      curr.chain <- curr.chain + 1
    }
    # Calculate R-Hat ---------------
    out <- list()
    out$rhat <- list()
    if (n.chains > 1) {
      out$rhat$beta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$beta.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2]
      out$rhat$alpha <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$alpha.samples)))), 
      			      autoburnin = FALSE)$psrf[, 2]
      if (fixed.sigma.sq) {
        out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$theta.samples)))), 
        			      autoburnin = FALSE)$psrf[, 2]
      } else {
        out$rhat$theta <- c(NA, gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$theta.samples[-1, , drop = FALSE])))), 
        			      autoburnin = FALSE)$psrf[, 2])
      }

    } else {
      out$rhat$beta <- rep(NA, p.occ)
      out$rhat$alpha <- rep(NA, p.det)
      out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3, 2))
    }

    # Get everything back in the original order
    out$coords <- coords[order(ord), ]
    out$z.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$z.samples))))
    out$z.samples <- mcmc(out$z.samples[, order(ord), drop = FALSE])
    out$X <- X[order(ord), , drop = FALSE]
    out$w.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples))))
    out$w.samples <- mcmc(out$w.samples[, order(ord), drop = FALSE])
    out$psi.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$psi.samples))))
    out$psi.samples <- mcmc(out$psi.samples[, order(ord), drop = FALSE])
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- x.names
    out$alpha.samples <- mcmc(do.call(rbind, 
      				lapply(out.tmp, function(a) t(a$alpha.samples))))
    colnames(out$alpha.samples) <- x.p.names
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    if (cov.model != 'matern') {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi')
    } else {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
    }
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$alpha <- effectiveSize(out$alpha.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    out$X.p <- X.p
    out$y <- y.big
    out$call <- cl
    out$n.samples <- batch.length * n.batch
    out$sites <- sites.orig
    out$n.neighbors <- n.neighbors
    out$cov.model.indx <- cov.model.indx
    out$type <- "NNGP"
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains

    # K-fold cross-validation -------
    if (!missing(k.fold)) {
      if (verbose) {      
        cat("----------------------------------------\n");
        cat("\tCross-validation\n");
        cat("----------------------------------------\n");
        message(paste("Performing ", k.fold, "-fold cross-validation using ", k.fold.threads,
      	      " thread(s).", sep = ''))
      }
      # Currently implemented without parellization. 
      set.seed(k.fold.seed)
      if (missing(k.fold.data)) {
        k.fold.data <- NULL
      }
      # Check to see if only one data source should be used for hold-out evaluations
      if (!is.null(k.fold.data)) {
        if (!is.numeric(k.fold.data) | length(k.fold.data) != 1) {
          stop("error: if specified, k.fold.data must be a single numeric value")
        }
        if (verbose) {
          message(paste("Only holding out data from data source ", k.fold.data, ".", sep = ''))
	}
	sites.random <- sample(sites[[k.fold.data]])
      } else {
        # Number of sites in each hold out data set. 
        sites.random <- sample(1:J)    
      }
      sites.k.fold <- split(sites.random, rep(1:k.fold, length.out = length(sites.random)))
      registerDoParallel(k.fold.threads)
      model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
        curr.set <- sort(sites.k.fold[[i]])
        curr.set.pred <- curr.set
	# curr.set.fit are the rows of X.fit that are included in the fitted model
	curr.set.fit <- (1:J)[-curr.set]
	if (!is.null(k.fold.data)) {
          curr.set.fit <- sort(unique(c(curr.set.fit, unlist(sites[-k.fold.data]))))
        }
        y.indx <- !((z.long.indx.r) %in% curr.set)
        if (!is.null(k.fold.data)) {
          y.indx <- ifelse(data.indx.r == k.fold.data, y.indx, TRUE)
        } 
        y.fit <- y[y.indx]
        y.0 <- y[!y.indx]
        z.inits.fit <- z.inits[curr.set.fit]
	w.inits.fit <- w.inits[curr.set.fit]
	coords.fit <- coords[curr.set.fit, , drop = FALSE]
	coords.0 <- coords[curr.set.pred, , drop = FALSE]
        X.p.fit <- X.p.all[y.indx, , drop = FALSE]
        X.p.0 <- X.p.all[!y.indx, , drop = FALSE]
        X.fit <- X[curr.set.fit, , drop = FALSE]
        X.0 <- X[curr.set.pred, , drop = FALSE]
        J.fit <- nrow(X.fit)
	# Site indices for fitted data
	sites.fit <- sapply(sites, 
			    function(a) which(as.numeric(row.names(X.fit)) %in% a[a %in% curr.set.fit]))
	tmp <- sapply(sites, function(a) a %in% curr.set.fit)
	if (!is.null(k.fold.data)) {
          sites.fit[[k.fold.data]] <- which(as.numeric(row.names(X.fit)) %in% sites[[k.fold.data]][sites[[k.fold.data]] %in% (1:J)[-curr.set]])
          tmp[[k.fold.data]] <- sites[[k.fold.data]] %in% (1:J)[-curr.set]
        }
        K.fit <- K[unlist(tmp)]
	y.big.fit <- y.big
	for (q in 1:n.data) {
          for (j in 1:J.long[q]) {
            if (is.null(k.fold.data)) {
              if (sites[[q]][j] %in% curr.set) {
                y.big.fit[[q]][j, ] <- NA
	      }
	    } else {
              if (q == k.fold.data) {
                if (sites[[q]][j] %in% curr.set) {
                  y.big.fit[[q]][j, ] <- NA
	        }
	      }
	    }
	  }
	  y.big.fit[[q]] <- y.big.fit[[q]][which(apply(y.big.fit[[q]], 
						       1, function(a) sum(!is.na(a)) > 0)), ]
	}
	y.big.0 <- y.big
	for (q in 1:n.data) {
          for (j in 1:J.long[q]) {
            if (is.null(k.fold.data)) {
              if (! sites[[q]][j] %in% curr.set) {
                y.big.0[[q]][j, ] <- NA
	      }
	    } else {
              if (q == k.fold.data) {
                if (! sites[[q]][j] %in% curr.set) {
                  y.big.0[[q]][j, ] <- NA
	        }
	      }
	    }
	  }
	  y.big.0[[q]] <- y.big.0[[q]][which(apply(y.big.0[[q]], 
						       1, function(a) sum(!is.na(a)) > 0)), ]
	}

	z.long.indx.fit <- list()
	for (q in 1:n.data) {
          z.long.indx.fit[[q]] <- rep(sites.fit[[q]], K.long.max[q])
	  z.long.indx.fit[[q]] <- z.long.indx.fit[[q]][!is.na(c(y.big.fit[[q]]))]
	}
	z.long.indx.fit <- unlist(z.long.indx.fit)
	z.long.indx.fit <- z.long.indx.fit - 1

	# Site indices for hold out data
	sites.0 <- sapply(sites, 
			    function(a) which(as.numeric(row.names(X.0)) %in% a[a %in% curr.set.pred]))
	tmp <- sapply(sites, function(a) a %in% curr.set.pred)
	if (!is.null(k.fold.data)) {
          sites.0[-k.fold.data] <- NA
	  for (q in 1:n.data) {
            if (q != k.fold.data) {
              tmp[[q]] <- rep(FALSE, J.long[q])
            }
          }
        }
        K.0 <- K[unlist(tmp)]
	if (is.null(k.fold.data)) {
	  z.long.indx.0 <- list()
	  for (q in 1:n.data) {
            z.long.indx.0[[q]] <- rep(sites.0[[q]], K.long.max[q])
	    z.long.indx.0[[q]] <- z.long.indx.0[[q]][!is.na(c(y.big.0[[q]]))]
	  }
	  z.long.indx.0 <- unlist(z.long.indx.0) - 1
	} else {
          z.long.indx.0 <- rep(sites.0[[k.fold.data]], K.long.max[k.fold.data])
	  z.long.indx.0 <- z.long.indx.0[!is.na(c(y.big.0[[k.fold.data]]))] - 1
	}
	z.long.indx.0 <- unlist(z.long.indx.0)
        verbose.fit <- FALSE
        n.omp.threads.fit <- 1
	n.obs.fit <- length(y.fit)
	n.obs.0 <- length(y.0)
	data.indx.c.fit <- data.indx.c[y.indx]
	data.indx.0 <- data.indx.c[!y.indx] + 1
	n.obs.long.fit <- as.vector(table(data.indx.c.fit))
	n.obs.long.0 <- n.obs.long - n.obs.long.fit
	J.long.fit <- as.vector(tapply(z.long.indx.fit, factor(data.indx.c.fit), 
			     FUN = function(a) length(unique(a))))

        # Nearest Neighbor Search ---
        ## Indexes
        if(search.type == "brute"){
          indx <- mkNNIndx(coords.fit, n.neighbors, n.omp.threads.fit)
        } else{
          indx <- mkNNIndxCB(coords.fit, n.neighbors, n.omp.threads.fit)
        }
        
        nn.indx.fit <- indx$nnIndx
        nn.indx.lu.fit <- indx$nnIndxLU
        
        indx <- mkUIndx(J.fit, n.neighbors, nn.indx.fit, 
	      	  nn.indx.lu.fit, u.search.type)
        
        u.indx.fit <- indx$u.indx
        u.indx.lu.fit <- indx$u.indx.lu
        ui.indx.fit <- indx$ui.indx

        storage.mode(y.fit) <- "double"
        storage.mode(z.inits.fit) <- "double"
        storage.mode(X.p.fit) <- "double"
        storage.mode(X.fit) <- "double"
        storage.mode(coords.fit) <- "double"
        storage.mode(p.det) <- "integer"
        storage.mode(p.det.long) <- "integer"
        storage.mode(p.occ) <- "integer"
        storage.mode(n.obs.fit) <- "integer"
        storage.mode(n.obs.long.fit) <- "integer"
        storage.mode(J.fit) <- "integer"
        storage.mode(J.long.fit) <- "integer"
        storage.mode(K.fit) <- "integer"
        storage.mode(n.data) <- "integer"
        storage.mode(beta.inits) <- "double"
        storage.mode(alpha.inits) <- "double"
        storage.mode(phi.inits) <- "double"
        storage.mode(sigma.sq.inits) <- "double"
        storage.mode(nu.inits) <- "double"
        storage.mode(w.inits.fit) <- "double"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(data.indx.c.fit) <- "integer"
        storage.mode(alpha.indx.c) <- "integer"
        storage.mode(mu.beta) <- "double"
        storage.mode(Sigma.beta) <- "double"
        storage.mode(mu.alpha) <- "double"
        storage.mode(sigma.alpha) <- "double"
        storage.mode(phi.a) <- "double"
        storage.mode(phi.b) <- "double"
        storage.mode(nu.a) <- "double"
        storage.mode(nu.b) <- "double"
        storage.mode(sigma.sq.a) <- "double"
        storage.mode(sigma.sq.b) <- "double"
        storage.mode(tuning.c) <- "double"
        storage.mode(n.batch) <- "integer"
        storage.mode(batch.length) <- "integer"
        storage.mode(accept.rate) <- "double"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(cov.model.indx) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(nn.indx.fit) <- "integer"
        storage.mode(nn.indx.lu.fit) <- "integer"
        storage.mode(u.indx.fit) <- "integer"
        storage.mode(u.indx.lu.fit) <- "integer"
        storage.mode(ui.indx.fit) <- "integer"
        storage.mode(n.neighbors) <- "integer"
        storage.mode(n.burn) <- "integer"
        storage.mode(n.thin) <- "integer"
	curr.chain <- 1
	storage.mode(curr.chain) <- "integer"

        out.fit <- .Call("spIntPGOccNNGP", y.fit, X.fit, X.p.fit, coords.fit, p.occ, p.det, p.det.long, 
		         J.fit, J.long.fit, K.fit, n.obs.fit, n.obs.long.fit, n.data, 
		         n.neighbors, nn.indx.fit, nn.indx.lu.fit, 
			 u.indx.fit, u.indx.lu.fit, ui.indx.fit, 
		         beta.inits, alpha.inits, z.inits.fit, w.inits.fit, 
		         phi.inits, sigma.sq.inits, nu.inits, 
		         z.long.indx.fit, data.indx.c.fit, alpha.indx.c, mu.beta, mu.alpha, 
		         Sigma.beta, sigma.alpha, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		         nu.a, nu.b, tuning.c, cov.model.indx,
		         n.batch, batch.length, accept.rate,  
		         n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, n.post.samples, 
			 curr.chain, n.chains, fixed.sigma.sq, sigma.sq.ig)

        out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
        colnames(out.fit$beta.samples) <- x.names
        out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
        if (cov.model != 'matern') {
          colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi')
        } else {
          colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi', 'nu')
        }
	# Don't need to reorder these, because nothing is re-ordered and it doesn't matter
	# since you are just returning the model deviance value.  
	out.fit$z.samples <- mcmc(t(out.fit$z.samples))
	out.fit$psi.samples <- mcmc(t(out.fit$psi.samples))
        out.fit$w.samples <- mcmc(t(out.fit$w.samples))
        out.fit$X <- X.fit
        out.fit$X.p <- X.p.fit
        out.fit$call <- cl
        out.fit$n.samples <- batch.length * n.batch
        out.fit$cov.model.indx <- cov.model.indx
        out.fit$type <- "NNGP"
	out.fit$n.neighbors <- n.neighbors
        out.fit$coords <- coords.fit
        out.fit$n.post <- n.post.samples
        out.fit$n.thin <- n.thin
        out.fit$n.burn <- n.burn
	out.fit$n.chains <- 1
	out.fit$pRE <- FALSE
	out.fit$psiRE <- FALSE
	class(out.fit) <- "spPGOcc"

        out.pred <- predict.spPGOcc(out.fit, X.0, coords.0, verbose = FALSE)
        # Detection 
	p.0.samples <- matrix(NA, n.post.samples, nrow(X.p.0))
        like.samples <- rep(NA, nrow(X.p.0))
        for (j in 1:nrow(X.p.0)) {
          p.0.samples[, j] <- logit.inv(X.p.0[j, 1:sum(alpha.indx.r == data.indx.0[j])] %*% out.fit$alpha.samples[which(alpha.indx.r == data.indx.0[j]), ])
          like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * out.pred$z.0.samples[, z.long.indx.0[j] + 1]))
        }
	as.vector(tapply(like.samples, data.indx.0, function(a) sum(log(a))))
      }
      model.deviance <- -2 * model.deviance
      # Return objects from cross-validation
      out$k.fold.deviance <- model.deviance
      stopImplicitCluster()
    }
  }
  class(out) <- "spIntPGOcc"
  out$run.time <- proc.time() - ptm
  # TODO: placeholder until you add in random effects for the integrated 
  #       models. 
  out$pRE <- FALSE
  out$psiRE <- FALSE
  out
}
