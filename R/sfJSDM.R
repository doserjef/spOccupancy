sfJSDM <- function(formula, data, inits, priors, 
		   tuning, cov.model = 'exponential', NNGP = TRUE, 
		   n.neighbors = 15, search.type = "cb", n.factors, 
		   n.batch, batch.length, accept.rate = 0.43,
		   n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		   n.burn = round(.10 * n.batch * batch.length), 
		   n.thin = 1, n.chains = 1, k.fold, k.fold.threads = 1, 
		   k.fold.seed = 100, k.fold.only = FALSE, monitors, 
		   keep.only.mean.95, ...){

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
  # Only implemented for NNGP
  if (!NNGP) {
    stop("error: sfJSDM is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 
  }
  if (missing(data)) {
    stop("error: data must be specified")
  }
  if (!is.list(data)) {
    stop("error: data must be a list")
  }
  names(data) <- tolower(names(data))
  data.orig <- data
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  if (length(dim(data$y)) != 2) {
    stop("error: detection-nondetection data y must be a two-dimensional array with dimensions corresponding to sites and replicates.")
  }
  y <- data$y
  sp.names <- attr(y, 'dimnames')[[1]]
  if (!'covs' %in% names(data)) {
    if (formula == ~ 1) {
      if (verbose) {
        message("covariates (covs) not specified in data.\nAssuming intercept only model.\n")
      }
      data$covs <- matrix(1, dim(y)[2], 1)
    } else {
      stop("error: covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  coords <- as.matrix(data$coords)
  # Check if all spatial coordinates are unique. 
  unique.coords <- unique(data$coords)
  if (nrow(unique.coords) < nrow(data$coords)) {
    stop("coordinates provided in coords are not all unique. spOccupancy requires each site to have its own unique pair of spatial coordinates. This may be the result of an error in preparing the data list, or you will need to change what you consider a 'site' in order to meet this requirement.") 
  }
  if (!missing(k.fold)) {
    if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
      stop("error: k.fold must be a single integer value >= 2")  
    }
  }
  if (missing(n.factors)) {
    stop("error: n.factors must be specified for a spatial factor occupancy model")
  }

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[, ord, drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Occupancy covariates
    data$covs <- data$covs[ord, , drop = FALSE]
  }

  data$covs <- as.data.frame(data$covs)

  # Checking missing values ---------------------------------------------
  if (sum(is.na(y) != 0)) {
    stop("error: some sites in y have missing data. Remove these from the data, and subsequently predict at those sites if interested.")
  }
  # occ.covs ------------------------
  if (sum(is.na(data$covs)) != 0) {
    stop("error: missing values in covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).")
  }

  # Check whether random effects are sent in as numeric, and
  # return error if they are. 
  # Occurrence ----------------------
  if (!is.null(findbars(formula))) {
    re.names <- sapply(findbars(formula), all.vars)
    for (i in 1:length(re.names)) {
      if (is(data$covs[, re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$covs[, re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }

  # Formula -------------------------------------------------------------
  # Occupancy -----------------------
  if (missing(formula)) {
    stop("error: formula must be specified")
  }

  if (is(formula, 'formula')) {
    tmp <- parseFormula(formula, data$covs)
    X <- as.matrix(tmp[[1]])
    X.re <- as.matrix(tmp[[4]])
    x.re.names <- colnames(X.re)
    x.names <- tmp[[2]]
  } else {
    stop("error: formula is misspecified")
  }
  # Get RE level names
  re.level.names <- lapply(data$covs[, x.re.names, drop = FALSE],
			   function (a) sort(unique(a)))

  # Extract data from inputs --------------------------------------------
  # Number of species 
  N <- dim(y)[1]
  # Number of latent factors
  q <- n.factors
  # Number of occupancy parameters 
  p.occ <- ncol(X)
  # Number of occupancy random effect parameters
  p.occ.re <- ncol(X.re)
  # Number of latent occupancy random effect values
  n.occ.re <- length(unlist(apply(X.re, 2, unique)))
  n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of sites
  J <- nrow(X)
  if (missing(n.batch)) {
    stop("error: must specify number of MCMC batches")
  }
  if (missing(batch.length)) {
    stop("error: must specify length of each MCMC batch")
  }
  n.samples <- n.batch * batch.length
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

  y.big <- y
  y <- c(y)
  y <- y[!is.na(y)]

  # Get random effect matrices all set ----------------------------------
  if (p.occ.re > 1) {
    for (j in 2:p.occ.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }
  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))

  # beta.comm -----------------------
  if ("beta.comm.normal" %in% names(priors)) {
    if (!is.list(priors$beta.comm.normal) | length(priors$beta.comm.normal) != 2) {
      stop("error: beta.comm.normal must be a list of length 2")
    }
    mu.beta.comm <- priors$beta.comm.normal[[1]]
    sigma.beta.comm <- priors$beta.comm.normal[[2]]
    if (length(mu.beta.comm) != p.occ & length(mu.beta.comm) != 1) {
      if (p.occ == 1) {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ",
        	     p.occ, " with elements corresponding to beta.comms' mean", sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ",
        	     p.occ, " or 1 with elements corresponding to beta.comms' mean", sep = ""))
      }
    }
    if (length(sigma.beta.comm) != p.occ & length(sigma.beta.comm) != 1) {
      if (p.occ == 1) {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ",
      	   p.occ, " with elements corresponding to beta.comms' variance", sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ",
      	   p.occ, " or 1 with elements corresponding to beta.comms' variance", sep = ""))
      }
    }
    if (length(sigma.beta.comm) != p.occ) {
      sigma.beta.comm <- rep(sigma.beta.comm, p.occ)
    }
    if (length(mu.beta.comm) != p.occ) {
      mu.beta.comm <- rep(mu.beta.comm, p.occ)
    }
    Sigma.beta.comm <- sigma.beta.comm * diag(p.occ)
  } else {
    if (verbose) {
      message("No prior specified for beta.comm.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.beta.comm <- rep(0, p.occ)
    sigma.beta.comm <- rep(2.72, p.occ)
    Sigma.beta.comm <- diag(p.occ) * 2.72
  }
  # tau.sq.beta -----------------------
  if ("tau.sq.beta.ig" %in% names(priors)) {
    if (!is.list(priors$tau.sq.beta.ig) | length(priors$tau.sq.beta.ig) != 2) {
      stop("error: tau.sq.beta.ig must be a list of length 2")
    }
    tau.sq.beta.a <- priors$tau.sq.beta.ig[[1]]
    tau.sq.beta.b <- priors$tau.sq.beta.ig[[2]]
    if (length(tau.sq.beta.a) != p.occ & length(tau.sq.beta.a) != 1) {
      if (p.occ == 1) {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ", 
      	   p.occ, " with elements corresponding to tau.sq.betas' shape", sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ", 
      	   p.occ, " or 1 with elements corresponding to tau.sq.betas' shape", sep = ""))
      }
    }
    if (length(tau.sq.beta.b) != p.occ & length(tau.sq.beta.b) != 1) {
      if (p.occ == 1) {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ", 
      	   p.occ, " with elements corresponding to tau.sq.betas' scale", sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ", 
      	   p.occ, " or 1 with elements corresponding to tau.sq.betas' scale", sep = ""))
      }
    }
    if (length(tau.sq.beta.a) != p.occ) {
      tau.sq.beta.a <- rep(tau.sq.beta.a, p.occ)
    }
    if (length(tau.sq.beta.b) != p.occ) {
      tau.sq.beta.b <- rep(tau.sq.beta.b, p.occ)
    }
  } else {
    if (verbose) {	    
      message("No prior specified for tau.sq.beta.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.beta.a <- rep(0.1, p.occ)
    tau.sq.beta.b <- rep(0.1, p.occ)
  }
  # phi -----------------------------
  if (!NNGP) {
    coords.D <- iDist(coords)
  }
  # Get distance matrix which is used if priors are not specified
  if ("phi.unif" %in% names(priors)) {
    if (!is.list(priors$phi.unif) | length(priors$phi.unif) != 2) {
      stop("error: phi.unif must be a list of length 2")
    }
    phi.a <- priors$phi.unif[[1]]
    phi.b <- priors$phi.unif[[2]]
    if (length(phi.a) != q & length(phi.a) != 1) {
      stop(paste("error: phi.unif[[1]] must be a vector of length ", 
      	   q, " or 1 with elements corresponding to phis' lower bound for each latent factor", sep = ""))
    }
    if (length(phi.b) != q & length(phi.b) != 1) {
      stop(paste("error: phi.unif[[2]] must be a vector of length ", 
      	   q, " or 1 with elements corresponding to phis' upper bound for each latent factor", sep = ""))
    }
    if (length(phi.a) != q) {
      phi.a <- rep(phi.a, q)
    }
    if (length(phi.b) != q) {
      phi.b <- rep(phi.b, q)
    }
  } else {
    if (verbose) {
    message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
    }
    if (NNGP) {
      coords.D <- iDist(coords)
    }
    phi.a <- rep(3 / max(coords.D), q)
    phi.b <- rep(3 / sort(unique(c(coords.D)))[2], q)
  }

  # nu -----------------------------
  if (cov.model == "matern") {
    if (!"nu.unif" %in% names(priors)) {
      stop("error: nu.unif must be specified in priors value list")
    }
    nu.a <- priors$nu.unif[[1]]
    nu.b <- priors$nu.unif[[2]]
    if (!is.list(priors$nu.unif) | length(priors$nu.unif) != 2) {
      stop("error: nu.unif must be a list of length 2")
    }
    if (length(nu.a) != q & length(nu.a) != 1) {
      stop(paste("error: nu.unif[[1]] must be a vector of length ", 
      	   q, " or 1 with elements corresponding to nus' lower bound for each latent factor", sep = ""))
    }
    if (length(nu.b) != q & length(nu.b) != 1) {
      stop(paste("error: nu.unif[[2]] must be a vector of length ", 
      	   q, " or 1 with elements corresponding to nus' upper bound for each latent factor", sep = ""))
    }
    if (length(nu.a) != q) {
      nu.a <- rep(nu.a, q)
    }
    if (length(nu.b) != q) {
      nu.b <- rep(nu.b, q)
    }
  } else {
    nu.a <- rep(0, q)
    nu.b <- rep(0, q)
  }

  # sigma.sq.psi --------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.psi.ig) | length(priors$sigma.sq.psi.ig) != 2) {
        stop("error: sigma.sq.psi.ig must be a list of length 2")
      }
      sigma.sq.psi.a <- priors$sigma.sq.psi.ig[[1]]
      sigma.sq.psi.b <- priors$sigma.sq.psi.ig[[2]]
      if (length(sigma.sq.psi.a) != p.occ.re & length(sigma.sq.psi.a) != 1) {
        if (p.occ.re == 1) {
        stop(paste("error: sigma.sq.psi.ig[[1]] must be a vector of length ", 
        	   p.occ.re, " with elements corresponding to sigma.sq.psis' shape", sep = ""))
        } else {
        stop(paste("error: sigma.sq.psi.ig[[1]] must be a vector of length ", 
        	   p.occ.re, " or 1 with elements corresponding to sigma.sq.psis' shape", sep = ""))
        }
      }
      if (length(sigma.sq.psi.b) != p.occ.re & length(sigma.sq.psi.b) != 1) {
        if (p.occ.re == 1) {
          stop(paste("error: sigma.sq.psi.ig[[2]] must be a vector of length ", 
        	   p.occ.re, " with elements corresponding to sigma.sq.psis' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.psi.ig[[2]] must be a vector of length ", 
        	   p.occ.re, " or 1with elements corresponding to sigma.sq.psis' scale", sep = ""))
        }
      }
      if (length(sigma.sq.psi.a) != p.occ.re) {
        sigma.sq.psi.a <- rep(sigma.sq.psi.a, p.occ.re)
      }
      if (length(sigma.sq.psi.b) != p.occ.re) {
        sigma.sq.psi.b <- rep(sigma.sq.psi.b, p.occ.re)
      }
  }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.psi.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.psi.a <- rep(0.1, p.occ.re)
      sigma.sq.psi.b <- rep(0.1, p.occ.re)
    }
  } else {
    sigma.sq.psi.a <- 0
    sigma.sq.psi.b <- 0
  }
  # Initial values --------------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # beta.comm -----------------------
  # ORDER: a p.occ vector ordered by the effects in the formula.
  if ("beta.comm" %in% names(inits)) {
    beta.comm.inits <- inits[["beta.comm"]]
    if (length(beta.comm.inits) != p.occ & length(beta.comm.inits) != 1) {
      if (p.occ == 1) {
        stop(paste("error: initial values for beta.comm must be of length ", p.occ, 
      	   sep = ""))
      } else {
        stop(paste("error: initial values for beta.comm must be of length ", p.occ, 
      	   , " or 1", sep = ""))
      }
    }
    if (length(beta.comm.inits) != p.occ) {
      beta.comm.inits <- rep(beta.comm.inits, p.occ)
    }
  } else {
    beta.comm.inits <- rnorm(p.occ, mu.beta.comm, sqrt(sigma.beta.comm))
    if (verbose) {
      message('beta.comm is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
    }
  }
  # tau.sq.beta ------------------------
  # ORDER: a p.occ vector ordered by the effects in the occurrence formula
  if ("tau.sq.beta" %in% names(inits)) {
    tau.sq.beta.inits <- inits[["tau.sq.beta"]]
    if (length(tau.sq.beta.inits) != p.occ & length(tau.sq.beta.inits) != 1) {
      if (p.occ == 1) {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p.occ, 
      	   sep = ""))
      } else {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p.occ, 
      	   " or 1", sep = ""))
      }
    }
    if (length(tau.sq.beta.inits) != p.occ) {
      tau.sq.beta.inits <- rep(tau.sq.beta.inits, p.occ)
    }
  } else {
    tau.sq.beta.inits <- runif(p.occ, 0.5, 10)
    if (verbose) {
      message('tau.sq.beta is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n')
    }
  }
  # beta ----------------------------
  # ORDER: N x p.occ matrix sent in as a column-major vector ordered by 
  #        parameter then species within parameter. 
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (is.matrix(beta.inits)) {
      if (ncol(beta.inits) != p.occ | nrow(beta.inits) != N) {
        stop(paste("error: initial values for beta must be a matrix with dimensions ", 
        	   N, "x", p.occ, " or a single numeric value", sep = ""))
      }
    }
    if (!is.matrix(beta.inits) & length(beta.inits) != 1) {
      stop(paste("error: initial values for beta must be a matrix with dimensions ", 
      	   N, " x ", p.occ, " or a single numeric value", sep = ""))
    }
    if (length(beta.inits) == 1) {
      beta.inits <- matrix(beta.inits, N, p.occ)
    }
  } else {
      beta.inits <- matrix(rnorm(N * p.occ, beta.comm.inits, sqrt(tau.sq.beta.inits)), N, p.occ)
      if (verbose) {
        message('beta is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n')
      }
  }
  # Create a N * p.occ x 1 matrix of the species-level regression coefficients. 
  # This is ordered by parameter, then species within a parameter. 
  beta.inits <- c(beta.inits)
  # phi -----------------------------
  # ORDER: a length N vector ordered by species in the detection-nondetection data.
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != q & length(phi.inits) != 1) {
      stop(paste("error: initial values for phi must be of length ", q, " or 1", 
      	   sep = ""))
    }
    if (length(phi.inits) != q) {
      phi.inits <- rep(phi.inits, q)
    }
  } else {
    phi.inits <- runif(q, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random values from the prior distribution\n")
    }
  }
  # lambda ----------------------------
  # ORDER: an N x q matrix sent in as a column-major vector, which is ordered by 
  #        factor, then species within factor. 
  if ("lambda" %in% names(inits)) {
    lambda.inits <- inits[["lambda"]]
    if (!is.matrix(lambda.inits)) {
      stop(paste("error: initial values for lambda must be a matrix with dimensions ",
		 N, " x ", q, sep = ""))
    }
    if (nrow(lambda.inits) != N | ncol(lambda.inits) != q) {
      stop(paste("error: initial values for lambda must be a matrix with dimensions ",
		 N, " x ", q, sep = ""))
    }
    if (!all.equal(diag(lambda.inits), rep(1, q))) {
      stop("error: diagonal of inits$lambda matrix must be all 1s")
    }
    if (sum(lambda.inits[upper.tri(lambda.inits)]) != 0) {
      stop("error: upper triangle of inits$lambda must be all 0s")
    }
  } else {
    lambda.inits <- matrix(0, N, q)
    diag(lambda.inits) <- 1
    lambda.inits[lower.tri(lambda.inits)] <- 0
    if (verbose) {
      message("lambda is not specified in initial values.\nSetting initial values of the lower triangle to 0\n")
    }
    # lambda.inits are organized by factor, then by species. This is necessary for working
    # with dgemv.  
    lambda.inits <- c(lambda.inits)
  }
  # nu ------------------------
  if ("nu" %in% names(inits)) {
    nu.inits <- inits[["nu"]]
    if (length(nu.inits) != q & length(nu.inits) != 1) {
      stop(paste("error: initial values for nu must be of length ", q,  " or 1",
      	   sep = ""))
    }
    if (length(nu.inits) != q) {
      nu.inits <- rep(nu.inits, q)
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      nu.inits <- runif(q, nu.a, nu.b)
    } else {
      nu.inits <- rep(0, q)
    }
  }

  if (p.occ.re > 0) {
  # sigma.sq.psi ------------------
  # ORDER: a length p.occ.re vector ordered by the random effects in the formula.
    if ("sigma.sq.psi" %in% names(inits)) {
      sigma.sq.psi.inits <- inits[["sigma.sq.psi"]]
      if (length(sigma.sq.psi.inits) != p.occ.re & length(sigma.sq.psi.inits) != 1) {
        if (p.occ.re == 1) {
          stop(paste("error: initial values for sigma.sq.psi must be of length ", p.occ.re, 
      	     sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.psi must be of length ", p.occ.re, 
      	     " or 1", sep = ""))
        }
      }
      if (length(sigma.sq.psi.inits) != p.occ.re) {
        sigma.sq.psi.inits <- rep(sigma.sq.psi.inits, p.occ.re)
      }
    } else {
      sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.psi is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    beta.star.indx <- rep(0:(p.occ.re - 1), n.occ.re.long)
  # beta.star -------------------------
  # ORDER: an N x n.occ.re matrix of random effects values for different levels.
  if ("beta.star" %in% names(inits)) {
    beta.star.inits <- inits[["beta.star"]]
    if (is.matrix(beta.star.inits)) {
      if (ncol(beta.star.inits) != n.occ.re | nrow(beta.star.inits) != N) {
        stop(paste("error: initial values for beta.star must be a matrix with dimensions ", 
        	   N, "x", n.occ.re, " or a single numeric value", sep = ""))
      }
    }
    if (!is.matrix(beta.star.inits) & length(beta.star.inits) != 1) {
      stop(paste("error: initial values for beta.star must be a matrix with dimensions ", 
      	   N, " x ", n.occ.re, " or a single numeric value", sep = ""))
    }
    if (length(beta.star.inits) == 1) {
      beta.star.inits <- matrix(beta.star.inits, N, n.occ.re)
    }
    beta.star.inits <- t(beta.star.inits)
  } else {
      beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
      beta.star.inits <- rep(beta.star.inits, N)
      if (verbose) {
        message('beta.star is not specified in initial values.\nSetting initial values to random values from the random effects variance\n')
      }
  }
  beta.star.inits <- c(beta.star.inits)
  } else {
    sigma.sq.psi.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
  }
  # w -----------------------------
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (!is.matrix(w.inits)) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   q, " x ", J, sep = ""))
    }
    if (nrow(w.inits) != q | ncol(w.inits) != J) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   q, " x ", J, sep = ""))
    }
    if (NNGP) {
      w.inits <- w.inits[, ord]
    }
  } else {
    w.inits <- matrix(0, q, J)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
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
  # Not accessed, but necessary to keep things in line with the underlying functions. 
  sigma.sq.tuning <- rep(0, q)
  phi.tuning <- rep(0, q)
  nu.tuning <- rep(0, q)
  if (missing(tuning)) {
    phi.tuning <- rep(1, q)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, q)
    }
  } else {
    names(tuning) <- tolower(names(tuning))
    # phi ---------------------------
    if(!"phi" %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (length(phi.tuning) == 1) {
      phi.tuning <- rep(tuning$phi, q)
    } else if (length(phi.tuning) != q) {
      stop(paste("error: phi tuning must be either a single value or a vector of length ",
      	   q, sep = ""))
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) == 1) {
        nu.tuning <- rep(tuning$nu, q)
      } else if (length(nu.tuning) != q) {
        stop(paste("error: nu tuning must be either a single value or a vector of length ",
        	   q, sep = ""))
      }
    }
  }
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA
  curr.chain <- 1

  if (!NNGP) {

    stop("error: sfJSDM is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 

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

    # Determine parameters that are monitored -----------------------------
    # The monitored parameters will be directly manipulated on the C side 
    n.track <- 11
    beta.comm.monitor <- 1
    tau.sq.beta.monitor <- 2
    beta.monitor <- 3
    z.monitor <- 4
    psi.monitor <- 5
    lambda.monitor <- 6
    theta.monitor <- 7
    w.monitor <- 8
    like.monitor <- 9
    beta.star.monitor <- 10 
    sigma.sq.psi.monitor <- 11
    if (missing(monitors)) {
      monitors <- rep(1, n.track) 
    } else {
      monitors.input <- monitors
      monitors <- rep(0, n.track)
      if ('beta.comm' %in% monitors.input) {
        monitors[beta.comm.monitor] <- 1
      }
      if ('tau.sq.beta' %in% monitors.input) {
        monitors[tau.sq.beta.monitor] <- 1
      }
      if ('beta' %in% monitors.input) {
        monitors[beta.monitor] <- 1
      }
      if ('z' %in% monitors.input) {
        monitors[z.monitor] <- 1
      }
      if ('psi' %in% monitors.input) {
        monitors[psi.monitor] <- 1
      }
      if ('lambda' %in% monitors.input) {
        monitors[lambda.monitor] <- 1
      }
      if ('theta' %in% monitors.input) {
        monitors[theta.monitor] <- 1
      }
      if ('w' %in% monitors.input) {
        monitors[w.monitor] <- 1
      }
      if ('like' %in% monitors.input) {
        monitors[like.monitor] <- 1
      }
      if ('beta.star' %in% monitors.input) {
        monitors[beta.star.monitor] <- 1
      }
      if ('sigma.sq.psi' %in% monitors.input) {
        monitors[sigma.sq.psi.monitor] <- 1
      }
    }
    
    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(coords) <- "double"
    # consts order: N, J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, n.det.re, q)
    consts <- c(N, J, p.occ, p.occ.re, n.occ.re, q)
    storage.mode(consts) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(beta.comm.inits) <- "double"
    storage.mode(tau.sq.beta.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(lambda.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(tau.sq.beta.a) <- "double"
    storage.mode(tau.sq.beta.b) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(tuning.c) <- "double"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.indx) <- "integer"
    storage.mode(u.indx.lu) <- "integer"
    storage.mode(ui.indx) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    # chain.info order: current chain, total number of chains
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    # samples.info order: burn-in, thinning rate, number of posterior samples
    samples.info <- c(n.burn, n.thin, n.post.samples)
    storage.mode(samples.info) <- "integer"
    # For occurrence random effects
    storage.mode(X.re) <- "integer"
    beta.level.indx <- sort(unique(c(X.re)))
    storage.mode(beta.level.indx) <- "integer"
    storage.mode(sigma.sq.psi.inits) <- "double"
    storage.mode(sigma.sq.psi.a) <- "double"
    storage.mode(sigma.sq.psi.b) <- "double"
    storage.mode(beta.star.inits) <- "double"
    storage.mode(beta.star.indx) <- "integer"
    # Monitors
    storage.mode(monitors) <- "integer"
    # Initial seed
    if (! exists(".Random.seed")) runif(1)
    init.seed <- .Random.seed

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    # Random seed information for each chain of the model. 
    seeds.list <- list()
    out <- list()
    if (!k.fold.only) {
      for (i in 1:n.chains) {
        # Change initial values if i > 1
        if ((i > 1) & (!fix.inits)) {
          beta.comm.inits <- rnorm(p.occ, mu.beta.comm, sqrt(sigma.beta.comm))
          tau.sq.beta.inits <- runif(p.occ, 0.5, 10)
          beta.inits <- matrix(rnorm(N * p.occ, beta.comm.inits, 
                		     sqrt(tau.sq.beta.inits)), N, p.occ)
          beta.inits <- c(beta.inits)
          lambda.inits <- matrix(0, N, q)
          diag(lambda.inits) <- 1
          lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
          lambda.inits <- c(lambda.inits)
          phi.inits <- runif(q, phi.a, phi.b)
          if (cov.model == 'matern') {
            nu.inits <- runif(q, nu.a, nu.b)
          }
          if (p.occ.re > 0) {
            sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
            beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
            beta.star.inits <- rep(beta.star.inits, N)
          }
        }

        storage.mode(chain.info) <- "integer"
        # Run the model in C
        out.tmp[[i]] <- .Call("sfJSDMNNGP", y, X, coords, X.re, consts, n.occ.re.long, 
          	            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
          	            beta.inits, beta.comm.inits, tau.sq.beta.inits, phi.inits, 
          	            lambda.inits, nu.inits, sigma.sq.psi.inits, beta.star.inits, w.inits,
          		    beta.star.indx, beta.level.indx, mu.beta.comm, Sigma.beta.comm, 
          	            tau.sq.beta.a, tau.sq.beta.b, phi.a, phi.b,
          	            nu.a, nu.b, sigma.sq.psi.a, sigma.sq.psi.b, 
          		    tuning.c, cov.model.indx, n.batch, 
          	            batch.length, accept.rate, n.omp.threads, verbose, n.report, 
          	            samples.info, chain.info, monitors)
        chain.info[1] <- chain.info[1] + 1
	seeds.list[[i]] <- .Random.seed
      }
      # Calculate R-Hat ---------------
      out <- list()
      out$rhat <- list()
      if (n.chains > 1) {
        # as.vector removes the "Upper CI" when there is only 1 variable. 
        if (monitors[beta.comm.monitor]) {
          out$rhat$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$beta.comm.samples)))), 
          			     autoburnin = FALSE)$psrf[, 2])
        }
        if (monitors[tau.sq.beta.monitor]) {
          out$rhat$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$tau.sq.beta.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        }
        if (monitors[beta.monitor]) {
          out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					         mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        }
        if (monitors[theta.monitor]) {
          out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$theta.samples)))), 
        			      autoburnin = FALSE)$psrf[, 2]
        }
        if (monitors[lambda.monitor]) {
          lambda.mat <- matrix(lambda.inits, N, q)
          out$rhat$lambda.lower.tri <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					       mcmc(t(a$lambda.samples[c(lower.tri(lambda.mat)), ])))), 
          					       autoburnin = FALSE)$psrf[, 2])
        }
        if (p.occ.re > 0) {
          if (monitors[sigma.sq.psi.monitor]) {
            out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$sigma.sq.psi.samples)))), 
          			     autoburnin = FALSE)$psrf[, 2])
          }
        }
      } else {
        out$rhat$beta.comm <- rep(NA, p.occ)
        out$rhat$tau.sq.beta <- rep(NA, p.occ)
        out$rhat$beta <- rep(NA, p.occ * N)
        out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 2 * q, q))
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
        }
      }

      # Put everything into MCMC objects
      if (monitors[beta.comm.monitor]) {
        out$beta.comm.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.comm.samples))))
        colnames(out$beta.comm.samples) <- x.names
      }
      if (monitors[tau.sq.beta.monitor]) {
        out$tau.sq.beta.samples <- mcmc(do.call(rbind, 
          				lapply(out.tmp, function(a) t(a$tau.sq.beta.samples))))
        colnames(out$tau.sq.beta.samples) <- x.names
      }

      if (is.null(sp.names)) {
        sp.names <- paste('sp', 1:N, sep = '')
      }
      coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
      if (monitors[beta.monitor]) {
        out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
        colnames(out$beta.samples) <- coef.names
      }
      if (p.occ.re > 0) {
        if (monitors[sigma.sq.psi.monitor]) {
          out$sigma.sq.psi.samples <- mcmc(
            do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.psi.samples))))
          colnames(out$sigma.sq.psi.samples) <- x.re.names
        }
        if (monitors[beta.star.monitor]) {
          out$beta.star.samples <- mcmc(
            do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
          tmp.names <- unlist(re.level.names)
          beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
          beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
          colnames(out$beta.star.samples) <- beta.star.names
        }
        out$re.level.names <- re.level.names
      }
      loadings.names <- paste(rep(sp.names, times = n.factors), rep(1:n.factors, each = N), sep = '-')
      if (monitors[lambda.monitor]) {
        out$lambda.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$lambda.samples))))
        colnames(out$lambda.samples) <- loadings.names
      }
      if (cov.model != 'matern') {
        theta.names <- paste(rep(c('phi'), each = q), 1:q, sep = '-')
      } else {
        theta.names <- paste(rep(c('phi', 'nu'), each = q), 1:q, sep = '-')
      } 
      if (monitors[theta.monitor]) {
        out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
        colnames(out$theta.samples) <- theta.names
      }

      # Return things back in the original order. 
      if (monitors[w.monitor]) {
        out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
          								dim = c(q, J, n.post.samples))))
        out$w.samples <- out$w.samples[, order(ord), , drop = FALSE]
        out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
      }
      if (monitors[psi.monitor]) {
        out$psi.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$psi.samples, 
          								dim = c(N, J, n.post.samples))))
        out$psi.samples <- out$psi.samples[, order(ord), ]
        out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      }
      if (monitors[like.monitor]) {
        out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples, 
          								dim = c(N, J, n.post.samples))))
        out$like.samples <- out$like.samples[, order(ord), ]
        out$like.samples <- aperm(out$like.samples, c(3, 1, 2))
      }
      if (monitors[z.monitor]) {
        out$z.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$z.samples, 
          								dim = c(N, J, n.post.samples))))
        out$z.samples <- out$z.samples[, order(ord), ]
        out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      }
      out$X.re <- X.re[order(ord), , drop = FALSE]
      # Calculate effective sample sizes
      out$ESS <- list()
      if (monitors[beta.comm.monitor]) {
        out$ESS$beta.comm <- effectiveSize(out$beta.comm.samples)
      }
      if (monitors[tau.sq.beta.monitor]) {
        out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
      }
      if (monitors[beta.monitor]) {
        out$ESS$beta <- effectiveSize(out$beta.samples)
      }
      if (monitors[theta.monitor]) {
        out$ESS$theta <- effectiveSize(out$theta.samples)
      }
      if (monitors[lambda.monitor]) {
        out$ESS$lambda <- effectiveSize(out$lambda.samples)
      }
      if (p.occ.re > 0) {
        if (monitors[sigma.sq.psi.monitor]) {
          out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
        }
      }
      # Only save 95% CIs and means for certain variables
      get.vals <- function(a) {
        tmp <- rbind(apply(a, 2, mean), apply(a, 2, quantile, c(0.025, 0.975)))
        row.names(tmp) <- c('mean', '0.025', '0.975')
        tmp
      }
      get.vals.big <- function(a) {
        tmp <- abind(array(apply(a, c(2, 3), mean), dim = c(1, N, J)), 
          	   apply(a, c(2, 3), quantile, c(0.025, 0.975)), along = 1)
        row.names(tmp) <- c('mean', '0.025', '0.975')
        tmp
      }
      if (!missing(keep.only.mean.95)) {
        if('beta.comm' %in% keep.only.mean.95) {
          out$beta.comm.samples <- get.vals(out$beta.comm.samples)
        }  
        if('tau.sq.beta' %in% keep.only.mean.95) {
          out$tau.sq.beta.samples <- get.vals(out$tau.sq.beta.samples)
        }  
        if('beta' %in% keep.only.mean.95) {
          out$beta.samples <- get.vals(out$beta.samples)
        }  
        if('z' %in% keep.only.mean.95) {
          out$z.samples <- get.vals.big(out$z.samples)
        }  
        if('psi' %in% keep.only.mean.95) {
          out$psi.samples <- get.vals.big(out$psi.samples)
        }  
        if('lambda' %in% keep.only.mean.95) {
          out$lambda.samples <- get.vals(out$lambda.samples)
        }  
        if('theta' %in% keep.only.mean.95) {
          out$theta.samples <- get.vals(out$theta.samples)
        }  
        if('w' %in% keep.only.mean.95) {
          out$w.samples <- get.vals.big(out$w.samples)
        }  
        if('like' %in% keep.only.mean.95) {
          out$like.samples <- get.vals.big(out$like.samples)
        }  
        if('sigma.sq.psi' %in% keep.only.mean.95) {
          out$sigma.sq.psi.samples <- get.vals(out$sigma.sq.psi.samples)
        }  
        if('beta.star' %in% keep.only.mean.95) {
          out$beta.star.samples <- get.vals(out$beta.star.samples)
        }  
      }
      out$X <- X[order(ord), , drop = FALSE]
      out$y <- y.big[, order(ord), drop = FALSE]
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$theta.names <- theta.names
      out$type <- "NNGP"
      out$coords <- coords[order(ord), ]
      out$cov.model.indx <- cov.model.indx
      out$n.neighbors <- n.neighbors
      out$q <- q
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$n.chains <- n.chains
      out$monitors <- monitors
      # Send out objects needed for updateMCMC 
      update.list <- list()
      tmp.val <- ifelse(cov.model == 'matern', q * 3, q * 2)
      update.list$tuning <- matrix(NA, tmp.val, n.chains)
      for (i in 1:n.chains) {
        update.list$tuning[, i] <- exp(out.tmp[[i]]$tuning)
      }
      update.list$accept.rate <- accept.rate
      update.list$n.batch <- n.batch
      update.list$batch.length <- batch.length
      update.list$n.omp.threads <- n.omp.threads
      update.list$data <- data.orig
      update.list$cov.model <- cov.model
      update.list$priors <- priors
      update.list$search.type <- search.type
      update.list$formula <- formula
      # Random seed to have for updating. 
      update.list$final.seed <- seeds.list
      out$update <- update.list
      if (p.occ.re > 0) {
        out$psiRE <- TRUE
      } else {
        out$psiRE <- FALSE
      }
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
      # Currently implemented without parellization. 
      # TODO: don't think this works when updating. 
      set.seed(k.fold.seed)
      # Number of sites in each hold out data set. 
      sites.random <- sample(1:J)    
      sites.k.fold <- split(sites.random, sites.random %% k.fold)
      registerDoParallel(k.fold.threads)
      model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
        curr.set <- sort(sites.random[sites.k.fold[[i]]])
        y.fit <- c(y.big[, -curr.set, drop = FALSE])
        y.fit <- y.fit[!is.na(y.fit)]
        y.0 <- c(y.big[, curr.set, drop = FALSE])
        y.0 <- y.0[!is.na(y.0)]
        y.big.fit <- y.big[, -curr.set, drop = FALSE]
        y.big.0 <- y.big[, curr.set, drop = FALSE]
        X.fit <- X[-curr.set, , drop = FALSE]
        X.0 <- X[curr.set, , drop = FALSE]
	w.inits.fit <- w.inits[, curr.set, drop = FALSE]
        coords.fit <- coords[-curr.set, , drop = FALSE]
        coords.0 <- coords[curr.set, , drop = FALSE]
        J.fit <- nrow(X.fit)
        J.0 <- nrow(X.0)
	X.re.fit <- X.re[-curr.set, , drop = FALSE]
	X.re.0 <- X.re[curr.set, , drop = FALSE]
	n.occ.re.fit <- length(unique(c(X.re.fit)))
	if (n.occ.re.fit == 0) {
          n.occ.re.long.fit <- 0
	} else {
          n.occ.re.long.fit <- apply(X.re.fit, 2, function(a) length(unique(a)))
	}
        if (p.occ.re > 0) {	
          beta.star.indx.fit <- rep(0:(p.occ.re - 1), n.occ.re.long.fit)
	  beta.level.indx.fit <- sort(unique(c(X.re.fit)))
          beta.star.inits.fit <- rnorm(n.occ.re.fit, 
	  			      sqrt(sigma.sq.psi.inits[beta.star.indx.fit + 1]))
          beta.star.inits.fit <- rep(beta.star.inits.fit, N)
          re.level.names.fit <- list()
	  for (t in 1:p.occ.re) {
            tmp.indx <- beta.level.indx.fit[beta.star.indx.fit == t - 1]
            re.level.names.fit[[t]] <- unlist(re.level.names)[tmp.indx + 1]    
          }
	} else {
          beta.star.indx.fit <- beta.star.indx
	  beta.level.indx.fit <- beta.level.indx
	  beta.star.inits.fit <- beta.star.inits
	  re.level.names.fit <- re.level.names
	}
        verbose.fit <- FALSE
        n.omp.threads.fit <- 1
        # Don't need to reorder things, since they are already sorted by 
        # the horizontal location in the coordinates. 

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
        storage.mode(X.fit) <- "double"
        storage.mode(coords.fit) <- "double"
        consts.fit <- c(N, J.fit, p.occ, p.occ.re, n.occ.re.fit, q)
        storage.mode(consts.fit) <- "integer"	
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(nn.indx.fit) <- "integer"
        storage.mode(nn.indx.lu.fit) <- "integer"
        storage.mode(u.indx.fit) <- "integer"
        storage.mode(u.indx.lu.fit) <- "integer"
        storage.mode(ui.indx.fit) <- "integer"
        storage.mode(X.re.fit) <- "integer"
        storage.mode(n.occ.re.long.fit) <- "integer"
        storage.mode(beta.star.inits.fit) <- "double"
        storage.mode(beta.star.indx.fit) <- "integer"
	storage.mode(beta.level.indx.fit) <- "integer"
        chain.info[1] <- 1
        storage.mode(chain.info) <- "integer"
	monitors.fit <- rep(1, n.track)
	storage.mode(monitors.fit) <- "integer"

      out.fit <- .Call("sfJSDMNNGP", y.fit, X.fit, coords.fit, 
		       X.re.fit, consts.fit, n.occ.re.long.fit, 
        	       n.neighbors, nn.indx.fit, nn.indx.lu.fit, u.indx.fit, 
		       u.indx.lu.fit, ui.indx.fit, beta.inits, 
        	       beta.comm.inits, tau.sq.beta.inits, phi.inits, 
        	       lambda.inits, nu.inits, sigma.sq.psi.inits,
		       beta.star.inits.fit, w.inits.fit, beta.star.indx.fit, beta.level.indx.fit, 
		       mu.beta.comm, Sigma.beta.comm,
        	       tau.sq.beta.a, tau.sq.beta.b, phi.a, phi.b,
        	       nu.a, nu.b, sigma.sq.psi.a, sigma.sq.psi.b, 
		       tuning.c, cov.model.indx, n.batch, 
        	       batch.length, accept.rate, n.omp.threads.fit, 
		       verbose.fit, n.report, 
        	       samples.info, chain.info, monitors.fit)

        if (is.null(sp.names)) {
          sp.names <- paste('sp', 1:N, sep = '')
        }
        coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
        out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
        colnames(out.fit$beta.samples) <- coef.names
        out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
        if (cov.model != 'matern') {
          theta.names <- paste(rep(c('phi'), each = q), 1:q, sep = '-')
        } else {
          theta.names <- paste(rep(c('phi', 'nu'), each = q), 1:q, sep = '-')
        } 
        colnames(out.fit$theta.samples) <- theta.names
        loadings.names <- paste(rep(sp.names, times = n.factors), 
				rep(1:n.factors, each = N), sep = '-')
        out.fit$lambda.samples <- mcmc(t(out.fit$lambda.samples))
        colnames(out.fit$lambda.samples) <- loadings.names
        out.fit$w.samples <- array(out.fit$w.samples, dim = c(q, J, n.post.samples))
        out.fit$w.samples <- aperm(out.fit$w.samples, c(3, 1, 2))
        out.fit$X <- X.fit
        out.fit$y <- y.big
        out.fit$call <- cl
        out.fit$n.samples <- n.samples
        out.fit$type <- "NNGP"
        out.fit$n.neighbors <- n.neighbors
        out.fit$q <- q
        out.fit$coords <- coords.fit
        out.fit$cov.model.indx <- cov.model.indx
        out.fit$n.post <- n.post.samples
        out.fit$n.thin <- n.thin
        out.fit$n.burn <- n.burn
        out.fit$n.chains <- 1
        if (p.occ.re > 0) {
          out.fit$sigma.sq.psi.samples <- mcmc(t(out.fit$sigma.sq.psi.samples))
          colnames(out.fit$sigma.sq.psi.samples) <- x.re.names
          out.fit$beta.star.samples <- mcmc(t(out.fit$beta.star.samples))
          tmp.names <- unlist(re.level.names.fit)
          beta.star.names <- paste(rep(x.re.names, n.occ.re.long.fit), tmp.names, sep = '-')
          beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re.fit), 
				   sep = '-')
          colnames(out.fit$beta.star.samples) <- beta.star.names
          out.fit$re.level.names <- re.level.names.fit
	  out.fit$X.re <- X.re.fit
        }
	if (p.occ.re > 0) {
          out.fit$psiRE <- TRUE
	} else {
          out.fit$psiRE <- FALSE	
	}
        class(out.fit) <- "sfJSDM"

        # Predict occurrence at new sites. 
	if (p.occ.re > 0) {
	  tmp <- unlist(re.level.names)
	  X.re.0 <- matrix(tmp[c(X.re.0 + 1)], nrow(X.re.0), ncol(X.re.0))
	  colnames(X.re.0) <- x.re.names
	}
	if (p.occ.re > 0) {X.0 <- cbind(X.0, X.re.0)}
        out.pred <- predict.sfJSDM(out.fit, X.0, 
				      coords.0, verbose = FALSE, ignore.RE = FALSE)
	like.samples <- matrix(NA, N, J.0)
	for (q in 1:N) {
          for (j in 1:J.0) {
            like.samples[q, j] <- mean(dbinom(y.big.0[q, j], 1, out.pred$psi.0.samples[, q, j]))
	  } # j
	} # q

        apply(like.samples, 1, function(a) sum(log(a), na.rm = TRUE))
      }
      model.deviance <- -2 * model.deviance
      # Return objects from cross-validation
      out$k.fold.deviance <- model.deviance
      stopImplicitCluster()
    }
   
    class(out) <- "sfJSDM"
  }

  out$run.time <- proc.time() - ptm
  return(out)
}
