svcPGBinom <- function(formula, data, inits, priors, tuning, 
                       svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
                       n.neighbors = 15, search.type = 'cb', n.batch, 
                       batch.length, accept.rate = 0.43,
                       n.omp.threads = 1, verbose = TRUE, n.report = 100, 
                       n.burn = round(.10 * n.batch * batch.length), 
                       n.thin = 1, n.chains = 1, 
                       k.fold, k.fold.threads = 1, k.fold.seed = 100, 
                       k.fold.only = FALSE, ...){

  ptm <- proc.time()

  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing to run the model\n");
    cat("----------------------------------------\n");
  }

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
  rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
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
  if (missing(formula)) {
    stop("error: formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  y <- c(data$y)
  if (!'weights' %in% names(data)) {
    stop("error: weights (binomial denominator) must be specified in data")
  }
  weights <- data$weights
  if (!'covs' %in% names(data)) {
    if (formula == ~ 1) {
      if (verbose) {
        message("covariates (covs) not specified in data.\nAssuming intercept only model.\n")
      }
      data$covs <- matrix(1, length(y), 1)
    } else {
      stop("error: covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!is.matrix(data$covs) & !is.data.frame(data$covs)) {
    stop("error: covs must be a matrix or data frame")
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("error: coords must be a matrix or data frame")
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

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[ord, drop = FALSE]
    weights <- weights[ord, drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Occupancy covariates
    data$covs <- data$covs[ord, , drop = FALSE]
  }

  data$covs <- as.data.frame(data$covs)
  
  # Checking missing values ---------------------------------------------
  # y -------------------------------
  if (sum(is.na(y) > 0)) {
    stop("error: some sites in y have missing values. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
  }
  # covs ------------------------
  if (sum(is.na(data$covs)) != 0) {
    stop("error: missing values in covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
  }

  # Check whether random effects are sent in as numeric, and
  # return error if they are. 
  # Occurrence ----------------------
  if (!is.null(findbars(formula))) {
    occ.re.names <- sapply(findbars(formula), all.vars)
    for (i in 1:length(occ.re.names)) {
      if (is(data$covs[, occ.re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", occ.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$covs[, occ.re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", occ.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }

  # Formula -------------------------------------------------------------
  # Occupancy -----------------------
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

  # Get basic info from inputs ------------------------------------------
  # Number of sites
  J <- nrow(coords)
  # Number of occupancy parameters 
  p <- ncol(X)
  # Number of occurrence random effect parameters
  p.re <- ncol(X.re)
  # Number of latent occupancy random effect values
  n.re <- length(unlist(apply(X.re, 2, unique)))
  n.re.long <- apply(X.re, 2, function(a) length(unique(a)))
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
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
  }
  # Check SVC columns -----------------------------------------------------
  if (is.character(svc.cols)) {
    # Check if all column names in svc are in covs
    if (!all(svc.cols %in% x.names)) {
        missing.cols <- svc.cols[!(svc.cols %in% x.names)]
        stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not inurrence covariates", sep=""))
    }
    # Convert desired column names into the numeric column index
    svc.cols <- (1:p)[x.names %in% svc.cols]
    
  } else if (is.numeric(svc.cols)) {
    # Check if all column indices are in 1:p
    if (!all(svc.cols %in% 1:p)) {
        missing.cols <- svc.cols[!(svc.cols %in% (1:p))]
        stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
    }
  }
  p.svc <- length(svc.cols)

  # Get random effect matrices all set ----------------------------------
  if (p.re > 1) {
    for (j in 2:p.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }

  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))
  # Logical vector indicating what parameters are estimated and what 
  # parameters are fixed. 6 is the total number of parameter types that 
  # can be estimated here. Note that phi and nu are both fixed if phi.unif = 'fixed' 
  all.params <- c('beta', 'alpha', 'phi', 'sigma.sq', 
		  'sigma.sq.psi', 'sigma.sq.p')
  n.params <- length(all.params)
  fixed.params <- rep(FALSE, n.params)
  # beta -----------------------
  if ("beta.normal" %in% names(priors)) {
    if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
      stop("error: beta.normal must be a list of length 2")
    }
    mu.beta <- priors$beta.normal[[1]]
    sigma.beta <- priors$beta.normal[[2]]
    if (length(mu.beta) != p & length(mu.beta) != 1) {
      if (p == 1) {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p, " with elements corresponding to betas' mean", sep = ""))
      } else {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p, " or 1 with elements corresponding to betas' mean", sep = ""))
      }
    }
    if (length(sigma.beta) != p & length(sigma.beta) != 1) {
      if (p == 1) {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
      	   p, " with elements corresponding to betas' variance", sep = ""))
      } else {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
      	   p, " or 1 with elements corresponding to betas' variance", sep = ""))
      }
    }
    if (length(sigma.beta) != p) {
      sigma.beta <- rep(sigma.beta, p)
    }
    if (length(mu.beta) != p) {
      mu.beta <- rep(mu.beta, p)
    }
    Sigma.beta <- sigma.beta * diag(p)
  } else {
    if (verbose) {
      message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.beta <- rep(0, p)
    sigma.beta <- rep(2.72, p)
    Sigma.beta <- diag(p) * 2.72
  }
  # phi -----------------------------
  # Get distance matrix which is used if priors are not specified
  if ("phi.unif" %in% names(priors)) {
    if (priors$phi.unif[1] == 'fixed') {
      fixed.params[which(all.params == 'phi')] <- TRUE 
      phi.a <- rep(1, p.svc)
      phi.b <- rep(1, p.svc)
      if (cov.model == 'matern') {
        message("phi is specified as fixed in priors$phi.unif. This will also fix nu at its initial value. Fixing phi without nu (or vice versa) is not supported.")
      }
    } else {
      if (!is.list(priors$phi.unif) | length(priors$phi.unif) != 2) {
        stop("error: phi.unif must be a list of length 2")
      }
      phi.a <- priors$phi.unif[[1]]
      phi.b <- priors$phi.unif[[2]]
      if (length(phi.a) != p.svc & length(phi.a) != 1) {
        stop(paste("error: phi.unif[[1]] must be a vector of length ", 
        	   p.svc, 
             " or 1 with elements corresponding to phis' lower bound for each covariate with spatially-varying coefficients",
             sep = ""))
      }
      if (length(phi.b) != p.svc & length(phi.b) != 1) {
        stop(paste("error: phi.unif[[2]] must be a vector of length ", 
        	   p.svc, 
             " or 1 with elements corresponding to phis' upper bound for each covariate with spatially-varying coefficients", sep = ""))
      }
      if (length(phi.a) != p.svc) {
        phi.a <- rep(phi.a, p.svc)
      }
      if (length(phi.b) != p.svc) {
        phi.b <- rep(phi.b, p.svc)
      }
    }
  } else {
    if (verbose) {
    message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
    }
    coords.D <- iDist(coords)
    phi.a <- rep(3 / max(coords.D), p.svc)
    phi.b <- rep(3 / sort(unique(c(coords.D)))[2], p.svc)
  }
  # sigma.sq -----------------------------
  if (("sigma.sq.ig" %in% names(priors)) & ("sigma.sq.unif" %in% names(priors))) {
    stop("error: cannot specify both an IG and a uniform prior for sigma.sq")
  }
  if ("sigma.sq.ig" %in% names(priors)) {
    sigma.sq.ig <- TRUE
    if (priors$sigma.sq.ig[1] == 'fixed') { # inverse-Gamma
      fixed.params[which(all.params == 'sigma.sq')] <- TRUE 
      sigma.sq.a <- rep(1, p.svc)
      sigma.sq.b <- rep(1, p.svc)
    } else {
      if (!is.list(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
        stop("error: sigma.sq.ig must be a list of length 2")
      }
      sigma.sq.a <- priors$sigma.sq.ig[[1]]
      sigma.sq.b <- priors$sigma.sq.ig[[2]]
      if (length(sigma.sq.a) != p.svc & length(sigma.sq.a) != 1) {
        stop(paste("error: sigma.sq.ig[[1]] must be a vector of length ", 
        	   p.svc, " or 1 with elements corresponding to sigma.sqs' shape for each covariate with spatially-varying coefficients", sep = ""))
      }
      if (length(sigma.sq.b) != p.svc & length(sigma.sq.b) != 1) {
        stop(paste("error: sigma.sq.ig[[2]] must be a vector of length ", 
        	   p.svc, " or 1 with elements corresponding to sigma.sqs' scale for each covariate with spatially-varying coefficients", sep = ""))
      }
      if (length(sigma.sq.a) != p.svc) {
        sigma.sq.a <- rep(sigma.sq.a, p.svc)
      }
      if (length(sigma.sq.b) != p.svc) {
        sigma.sq.b <- rep(sigma.sq.b, p.svc)
      }
    }
  } else if ("sigma.sq.unif" %in% names(priors)) { # uniform prior
    if (priors$sigma.sq.unif[1] == 'fixed') {
      sigma.sq.ig <- TRUE # This just makes the C++ side a bit easier
      fixed.params[which(all.params == 'sigma.sq')] <- TRUE
      sigma.sq.a <- 1
      sigma.sq.b <- 1
    } else {
      sigma.sq.ig <- FALSE
      if (!is.list(priors$sigma.sq.unif) | length(priors$sigma.sq.unif) != 2) {
        stop("error: sigma.sq.unif must be a list of length 2")
      }
      sigma.sq.a <- priors$sigma.sq.unif[[1]]
      sigma.sq.b <- priors$sigma.sq.unif[[2]]
      if (length(sigma.sq.a) != p.svc & length(sigma.sq.a) != 1) {
        stop(paste("error: sigma.sq.unif[[1]] must be a vector of length ", 
        	   p.svc, " or 1 with elements corresponding to sigma.sqs' shape for each covariate with spatially-varying coefficients", sep = ""))
      }
      if (length(sigma.sq.b) != p.svc & length(sigma.sq.b) != 1) {
        stop(paste("error: sigma.sq.unif[[2]] must be a vector of length ", 
        	   p.svc, " or 1 with elements corresponding to sigma.sqs' scale for each covariate with spatially-varying coefficients", sep = ""))
      }
      if (length(sigma.sq.a) != p.svc) {
        sigma.sq.a <- rep(sigma.sq.a, p.svc)
      }
      if (length(sigma.sq.b) != p.svc) {
        sigma.sq.b <- rep(sigma.sq.b, p.svc)
      }
    }
  } else {
    if (verbose) {
      message("No prior specified for sigma.sq.\nUsing an inverse-Gamma prior with the shape parameter set to 2 and scale parameter to 1.\n")
    }
    sigma.sq.ig <- TRUE
    sigma.sq.a <- rep(2, p.svc)
    sigma.sq.b <- rep(1, p.svc)
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
    if (length(nu.a) != p.svc & length(nu.a) != 1) {
      stop(paste("error: nu.unif[[1]] must be a vector of length ", 
      	   p.svc, " or 1 with elements corresponding to nus' lower bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(nu.b) != p.svc & length(nu.b) != 1) {
      stop(paste("error: nu.unif[[2]] must be a vector of length ", 
      	   p.svc, " or 1 with elements corresponding to nus' upper bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(nu.a) != p.svc) {
      nu.a <- rep(nu.a, p.svc)
    }
    if (length(nu.b) != p.svc) {
      nu.b <- rep(nu.b, p.svc)
    }
  } else {
    nu.a <- rep(0, p.svc)
    nu.b <- rep(0, p.svc)
  }

  # sigma.sq.psi --------------------
  if (p.re > 0) {
    if ("sigma.sq.psi.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.psi.ig) | length(priors$sigma.sq.psi.ig) != 2) {
        stop("error: sigma.sq.psi.ig must be a list of length 2")
      }
      sigma.sq.psi.a <- priors$sigma.sq.psi.ig[[1]]
      sigma.sq.psi.b <- priors$sigma.sq.psi.ig[[2]]
      if (length(sigma.sq.psi.a) != p.re & length(sigma.sq.psi.a) != 1) {
        if (p.re == 1) {
        stop(paste("error: sigma.sq.psi.ig[[1]] must be a vector of length ", 
        	   p.re, " with elements corresponding to sigma.sq.psis' shape", sep = ""))
        } else {
        stop(paste("error: sigma.sq.psi.ig[[1]] must be a vector of length ", 
        	   p.re, " or 1 with elements corresponding to sigma.sq.psis' shape", sep = ""))
        }
      }
      if (length(sigma.sq.psi.b) != p.re & length(sigma.sq.psi.b) != 1) {
        if (p.re == 1) {
          stop(paste("error: sigma.sq.psi.ig[[2]] must be a vector of length ", 
        	   p.re, " with elements corresponding to sigma.sq.psis' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.psi.ig[[2]] must be a vector of length ", 
        	   p.re, " or 1with elements corresponding to sigma.sq.psis' scale", sep = ""))
        }
      }
      if (length(sigma.sq.psi.a) != p.re) {
        sigma.sq.psi.a <- rep(sigma.sq.psi.a, p.re)
      }
      if (length(sigma.sq.psi.b) != p.re) {
        sigma.sq.psi.b <- rep(sigma.sq.psi.b, p.re)
      }
  }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.psi.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.psi.a <- rep(0.1, p.re)
      sigma.sq.psi.b <- rep(0.1, p.re)
    }
  } else {
    sigma.sq.psi.a <- 0
    sigma.sq.psi.b <- 0
  }

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # beta -----------------------
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (length(beta.inits) != p & length(beta.inits) != 1) {
      if (p == 1) {
        stop(paste("error: initial values for beta must be of length ", p,
      	     sep = ""))

      } else {
        stop(paste("error: initial values for beta must be of length ", p, " or 1",
        	     sep = ""))
      }
    }
    if (length(beta.inits) != p) {
      beta.inits <- rep(beta.inits, p)
    }
  } else {
    beta.inits <- rnorm(p, mu.beta, sqrt(sigma.beta))
    if (verbose) {
      message('beta is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
    }
  }
  # phi -----------------------------
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != p.svc & length(phi.inits) != 1) {
      stop(paste("error: initial values for phi must be of length ", p.svc, " or 1", 
      	   sep = ""))
    }
    if (length(phi.inits) != p.svc) {
      phi.inits <- rep(phi.inits, p.svc)
    }
  } else {
    phi.inits <- runif(p.svc, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random values from the prior distribution\n")
    }
  }
  # sigma.sq ------------------------
  if ("sigma.sq" %in% names(inits)) {
    sigma.sq.inits <- inits[["sigma.sq"]]
    if (length(sigma.sq.inits) != p.svc & length(sigma.sq.inits) != 1) {
      stop(paste("error: initial values for sigma.sq must be of length ", p.svc,  " or 1",
      	   sep = ""))
    }
    if (length(sigma.sq.inits) != p.svc) {
      sigma.sq.inits <- rep(sigma.sq.inits, p.svc)
    }
  } else {
    if (sigma.sq.ig) {
      sigma.sq.inits <- rigamma(p.svc, sigma.sq.a, sigma.sq.b)
    } else {
      sigma.sq.inits <- runif(p.svc, sigma.sq.a, sigma.sq.b)
    }
    if (verbose) {
      message("sigma.sq is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  # w -----------------------------
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (!is.matrix(w.inits)) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   p.svc, " x ", J, sep = ""))
    }
    if (nrow(w.inits) != p.svc | ncol(w.inits) != J) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   p.svc, " x ", J, sep = ""))
    }
    if (NNGP) {
      w.inits <- w.inits[, ord, drop = FALSE]
    }
  } else {
    w.inits <- matrix(0, p.svc, J)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
    }
  }
  # nu ------------------------
  if ("nu" %in% names(inits)) {
    nu.inits <- inits[["nu"]]
    if (length(nu.inits) != p.svc & length(nu.inits) != 1) {
      stop(paste("error: initial values for nu must be of length ", p.svc,  " or 1",
      	   sep = ""))
    }
    if (length(nu.inits) != p.svc) {
      nu.inits <- rep(nu.inits, p.svc)
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      nu.inits <- runif(p.svc, nu.a, nu.b)
    } else {
      nu.inits <- rep(0, p.svc)
    }
  }
  # sigma.sq.psi -------------------
  if (p.re > 0) {
    if ("sigma.sq.psi" %in% names(inits)) {
      sigma.sq.psi.inits <- inits[["sigma.sq.psi"]]
      if (length(sigma.sq.psi.inits) != p.re & length(sigma.sq.psi.inits) != 1) {
        if (p.re == 1) {
          stop(paste("error: initial values for sigma.sq.psi must be of length ", p.re, 
      	       sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.psi must be of length ", p.re, 
      	       " or 1", sep = ""))
        }
      }
      if (length(sigma.sq.psi.inits) != p.re) {
        sigma.sq.psi.inits <- rep(sigma.sq.psi.inits, p.re)
      }
    } else {
      sigma.sq.psi.inits <- runif(p.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.psi is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    beta.star.indx <- rep(0:(p.re - 1), n.re.long)
    beta.star.inits <- rnorm(n.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
  } else {
    sigma.sq.psi.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
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
  storage.mode(cov.model.indx) <- "integer"

  # Prep for SVCs ---------------------------------------------------------
  X.w <- X[, svc.cols, drop = FALSE]

  # Get tuning values ---------------------------------------------------
  # Not accessed, but necessary to keep things in line. 
  sigma.sq.tuning <- rep(0, p.svc)
  phi.tuning <- rep(0, p.svc)
  nu.tuning <- rep(0, p.svc)
  if (missing(tuning)) {
    phi.tuning <- rep(1, p.svc)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, p.svc)
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
    if (length(phi.tuning) == 1) {
      phi.tuning <- rep(tuning$phi, p.svc)
    } else if (length(phi.tuning) != p.svc) {
      stop(paste("error: phi tuning must be either a single value or a vector of length ",
      	   p.svc, sep = ""))
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) == 1) {
        nu.tuning <- rep(tuning$nu, p.svc)
      } else if (length(nu.tuning) != p.svc) {
        stop(paste("error: nu tuning must be either a single value or a vector of length ",
        	   p.svc, sep = ""))
      }
    }
    if (!sigma.sq.ig) {
      # sigma.sq --------------------------
      if(!"sigma.sq" %in% names(tuning)) {
        stop("error: sigma.sq must be specified in tuning value list")
      }
      sigma.sq.tuning <- tuning$sigma.sq
      if (length(sigma.sq.tuning) == 1) {
        sigma.sq.tuning <- rep(tuning$sigma.sq, p.svc)
      } else if (length(sigma.sq.tuning) != p.svc) {
        stop(paste("error: sigma.sq tuning must be either a single value or a vector of length ",
		   p.svc, sep = ""))
      }
    }
  }
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA
  curr.chain <- 1

  if (!NNGP) {
    stop("error: svcPGBinom is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 

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
    
    storage.mode(n.neighbors) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    ## Indexes
    if(search.type == "brute"){
      indx <- mkNNIndx(coords, n.neighbors, n.omp.threads)
    } else{
      indx <- mkNNIndxCB(coords, n.neighbors, n.omp.threads)
    }
    
    nn.indx <- indx$nnIndx
    nn.indx.lu <- indx$nnIndxLU
    nn.indx.run.time <- indx$run.time
    
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.search.type) <- "integer"
    storage.mode(J) <- "integer"

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

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(X.w) <- "double"
    consts <- c(J, p, p.re, n.re, p.svc)
    storage.mode(consts) <- "integer"
    storage.mode(weights) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
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
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.indx) <- "integer"
    storage.mode(u.indx.lu) <- "integer"
    storage.mode(ui.indx) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    storage.mode(fixed.params) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"
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

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    out <- list()
    if (!k.fold.only) {
      for (i in 1:n.chains) {
        # Change initial values if i > 1
        if ((i > 1) & (!fix.inits)) {
          beta.inits <- rnorm(p, mu.beta, sqrt(sigma.beta))
          if (!fixed.params[which(all.params == 'sigma.sq')]) {
            if (sigma.sq.ig) {
              sigma.sq.inits <- rigamma(p.svc, sigma.sq.a, sigma.sq.b)
            } else {
              sigma.sq.inits <- runif(p.svc, sigma.sq.a, sigma.sq.b)
            }
          }
          phi.inits <- runif(p.svc, phi.a, phi.b)
          if (cov.model == 'matern') {
            nu.inits <- runif(p.svc, nu.a, nu.b)
          }
          if (p.re > 0) {
            sigma.sq.psi.inits <- runif(p.re, 0.5, 10)
            beta.star.inits <- rnorm(n.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
          }
        }
        storage.mode(chain.info) <- "integer"
        # Run the model in C    
        out.tmp[[i]] <- .Call("svcPGBinomNNGP", y, X, X.w, coords, X.re, consts, 
                              weights, n.re.long, 
                              n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
                              beta.inits, sigma.sq.psi.inits, beta.star.inits, 
                              w.inits, phi.inits, sigma.sq.inits, nu.inits, 
                              beta.star.indx, beta.level.indx, mu.beta, 
                              Sigma.beta, phi.a, phi.b, 
                              sigma.sq.a, sigma.sq.b, nu.a, nu.b, 
                              sigma.sq.psi.a, sigma.sq.psi.b, 
                              tuning.c, cov.model.indx,
                              n.batch, batch.length, 
                              accept.rate, n.omp.threads, verbose, n.report, 
                              samples.info, chain.info, fixed.params, sigma.sq.ig)
        chain.info[1] <- chain.info[1] + 1
      }
      # Calculate R-Hat ---------------
      out$rhat <- list()
      if (n.chains > 1) {
        # as.vector removes the "Upper CI" when there is only 1 variable. 
        out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					         mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        if (!fixed.params[which(all.params == 'sigma.sq')] & 
            !fixed.params[which(all.params == 'phi')]) { # none are fixed
          out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					        mcmc(t(a$theta.samples)))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
        } else if (fixed.params[which(all.params == 'sigma.sq')] & 
          	 !fixed.params[which(all.params == 'phi')]) { # sigma.sq is fixed
          out$rhat$theta <- c(NA, gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					        mcmc(t(a$theta.samples[-1, , drop = FALSE])))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        } else if (!fixed.params[which(all.params == 'sigma.sq')] & 
          	 fixed.params[which(all.params == 'phi')]) { # phi/nu is fixed
          tmp <- ifelse(cov.model == 'matern', NA, c(NA, NA))
          out$rhat$theta <- c(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					        mcmc(t(a$theta.samples[1, , drop = FALSE])))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2], tmp)
        } else { # both are fixed
          out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3, 2))
        } 
        if (p.re > 0) {
          out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$sigma.sq.psi.samples)))), 
          			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        }
      } else {
        out$rhat$beta <- rep(NA, p)
        out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3 * p.svc, 2 * p.svc))
        if (p.re > 0) {
          out$rhat$sigma.sq.psi <- rep(NA, p.re)
        }
      }
      out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
      colnames(out$beta.samples) <- x.names
      out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
      if (cov.model != 'matern') {
        theta.names <- paste(rep(c('sigma.sq', 'phi'), each = p.svc), x.names[svc.cols], sep = '-')
      } else {
        theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = p.svc), x.names[svc.cols], sep = '-')
      } 
      colnames(out$theta.samples) <- theta.names
      # Get everything back in the original order
      out$coords <- coords[order(ord), ]
      out$y.rep.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.samples))))
      out$y.rep.samples <- mcmc(out$y.rep.samples[, order(ord), drop = FALSE])
      out$X <- X[order(ord), , drop = FALSE]
      out$X.re <- X.re[order(ord), , drop = FALSE]
      out$X.w <- X.w[order(ord), , drop = FALSE]
      # Account for case when intercept only spatial model. 
      if (p.svc == 1) {
        tmp <- do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples)))
        tmp <- tmp[, order(ord), drop = FALSE]
        out$w.samples <- array(NA, dim = c(p.svc, J, n.post.samples * n.chains))
        out$w.samples[1, , ] <- t(tmp)
      } else {
        out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
          								dim = c(p.svc, J, n.post.samples))))
        out$w.samples <- out$w.samples[, order(ord), ]
      }
      out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
      out$psi.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$psi.samples))))
      out$psi.samples <- mcmc(out$psi.samples[, order(ord), drop = FALSE])
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
      out$like.samples <- mcmc(out$like.samples[, order(ord), drop = FALSE])
      out$y <- y[order(ord)]
      if (p.re > 0) {
        out$sigma.sq.psi.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.psi.samples))))
        colnames(out$sigma.sq.psi.samples) <- x.re.names
        out$beta.star.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
        tmp.names <- unlist(re.level.names)
        beta.star.names <- paste(rep(x.re.names, n.re.long), tmp.names, sep = '-')
        colnames(out$beta.star.samples) <- beta.star.names
        out$re.level.names <- re.level.names
      }
      # Calculate effective sample sizes
      out$ESS <- list()
      out$ESS$beta <- effectiveSize(out$beta.samples)
      out$ESS$theta <- effectiveSize(out$theta.samples)
      if (p.re > 0) {
        out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
      }
      out$call <- cl
      out$n.samples <- batch.length * n.batch
      out$n.neighbors <- n.neighbors
      out$cov.model.indx <- cov.model.indx
      out$svc.cols <- svc.cols
      out$type <- "NNGP"
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$n.chains <- n.chains
      if (p.re > 0) {
        out$psiRE <- TRUE
      } else {
        out$psiRE <- FALSE
      }
    }
    # K-fold cross-validation ---------
    if (!missing(k.fold)) {
      if (verbose) {
        cat("----------------------------------------\n");
        cat("\tCross-validation\n");
        cat("----------------------------------------\n");
        message(paste("Performing ", k.fold, "-fold cross-validation using ", k.fold.threads,
        	        " thread(s).", sep = ''))
      }
      set.seed(k.fold.seed)
      # Number of sites in each hold out data set. 
      sites.random <- sample(1:J)    
      sites.k.fold <- split(sites.random, sites.random %% k.fold)
      registerDoParallel(k.fold.threads)
      model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
        curr.set <- sort(sites.random[sites.k.fold[[i]]])
        y.fit <- y[-curr.set]
	y.0 <- y[curr.set]
	X.fit <- X[-curr.set, , drop = FALSE]
	X.0 <- X[curr.set, , drop = FALSE]
	X.w.fit <- X.w[-curr.set, , drop = FALSE]
	X.w.0 <- X.w[curr.set, , drop = FALSE]
	J.fit <- nrow(X.fit)
	J.0 <- nrow(X.0)
        coords.fit <- coords[-curr.set, , drop = FALSE]
        coords.0 <- coords[curr.set, , drop = FALSE]
        weights.fit <- weights[-curr.set]
        weights.0 <- weights[curr.set]
	# Random Occurrence Effects
        X.re.fit <- X.re[-curr.set, , drop = FALSE]
        X.re.0 <- X.re[curr.set, , drop = FALSE]
        n.re.fit <- length(unique(c(X.re.fit)))
        n.re.long.fit <- apply(X.re.fit, 2, function(a) length(unique(a)))
        if (p.re > 0) {	
          beta.star.indx.fit <- rep(0:(p.re - 1), n.re.long.fit)
          beta.level.indx.fit <- sort(unique(c(X.re.fit)))
          beta.star.inits.fit <- rnorm(n.re.fit, 0,
          			      sqrt(sigma.sq.psi.inits[beta.star.indx.fit + 1]))
          re.level.names.fit <- list()
          for (t in 1:p.re) {
            tmp.indx <- beta.level.indx.fit[beta.star.indx.fit == t - 1]
            re.level.names.fit[[t]] <- unlist(re.level.names)[tmp.indx + 1]    
          }
        } else {
          beta.star.indx.fit <- beta.star.indx
          beta.level.indx.fit <- beta.level.indx
          beta.star.inits.fit <- beta.star.inits
          re.level.names.fit <- re.level.names
        }
        # Nearest Neighbor Search ---
	verbose.fit <- FALSE
	n.omp.threads.fit <- 1
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
	storage.mode(X.w.fit) <- "double"
        storage.mode(coords.fit) <- "double"
	storage.mode(weights.fit) <- "double"
        consts.fit <- c(J.fit, p, p.re, n.re.fit, p.svc)
        storage.mode(consts.fit) <- "integer"
        storage.mode(nn.indx.fit) <- "integer"
        storage.mode(nn.indx.lu.fit) <- "integer"
        storage.mode(u.indx.fit) <- "integer"
        storage.mode(u.indx.lu.fit) <- "integer"
        storage.mode(ui.indx.fit) <- "integer"
        storage.mode(n.samples) <- "integer"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(X.re.fit) <- "integer"
        storage.mode(n.re.long.fit) <- "integer"
        storage.mode(beta.star.inits.fit) <- "double"
        storage.mode(beta.star.indx.fit) <- "integer"
        storage.mode(beta.level.indx.fit) <- "integer"
        chain.info[1] <- 1
        storage.mode(chain.info) <- "integer"
        # Run the model in C
        out.fit <- .Call("svcPGBinomNNGP", y.fit, X.fit, X.w.fit, coords.fit, X.re.fit, 
			 consts.fit, weights.fit, n.re.long.fit, 
        	         n.neighbors, nn.indx.fit, nn.indx.lu.fit, u.indx.fit, 
			 u.indx.lu.fit, ui.indx.fit, beta.inits, 
			 sigma.sq.psi.inits, beta.star.inits.fit, 
			 w.inits, phi.inits, sigma.sq.inits, 
			 nu.inits, beta.star.indx.fit, 
			 beta.level.indx.fit, mu.beta, 
			 Sigma.beta, phi.a, phi.b, 
			 sigma.sq.a, sigma.sq.b, nu.a, nu.b, sigma.sq.psi.a, sigma.sq.psi.b, 
			 tuning.c, cov.model.indx, 
			 n.batch, batch.length, accept.rate, n.omp.threads.fit, verbose.fit, 
			 n.report, samples.info, chain.info, fixed.params, sigma.sq.ig)
        out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
        colnames(out.fit$beta.samples) <- x.names
        out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
        # Account for case when intercept only spatial model. 
        if (p.svc == 1) {
          tmp <- t(out.fit$w.samples)
          out.fit$w.samples <- array(NA, dim = c(p.svc, J.fit, n.post.samples * n.chains))
          out.fit$w.samples[1, , ] <- tmp
        } else {
          out.fit$w.samples <- array(out.fit$w.samples, dim = c(p.svc, J.fit, n.post.samples))
        }
        out.fit$w.samples <- aperm(out.fit$w.samples, c(3, 1, 2))
        out.fit$X <- X.fit
        out.fit$y <- y.fit
        out.fit$call <- cl
        out.fit$type <- "NNGP"
        out.fit$n.neighbors <- n.neighbors
        out.fit$n.samples <- n.samples
        out.fit$coords <- coords.fit
        out.fit$cov.model.indx <- cov.model.indx
	out.fit$svc.cols <- svc.cols
        out.fit$n.post <- n.post.samples
        out.fit$n.thin <- n.thin
        out.fit$n.burn <- n.burn
        out.fit$n.chains <- 1
        if (p.re > 0) {
          out.fit$sigma.sq.psi.samples <- mcmc(t(out.fit$sigma.sq.psi.samples))
          colnames(out.fit$sigma.sq.psi.samples) <- x.re.names
          out.fit$beta.star.samples <- mcmc(t(out.fit$beta.star.samples))
          tmp.names <- unlist(re.level.names.fit)
          beta.star.names <- paste(rep(x.re.names, n.re.long.fit), tmp.names, sep = '-')
          colnames(out.fit$beta.star.samples) <- beta.star.names
          out.fit$re.level.names <- re.level.names.fit
          out.fit$X.re <- X.re.fit
        }
        if (p.re > 0) {
          out.fit$psiRE <- TRUE
        } else {
          out.fit$psiRE <- FALSE	
        }
        class(out.fit) <- "svcPGBinom"

	# Predict occurrence at new sites
	if (p.re > 0) {
	  tmp <- unlist(re.level.names)
	  X.re.0 <- matrix(tmp[c(X.re.0 + 1)], nrow(X.re.0), ncol(X.re.0))
	  colnames(X.re.0) <- x.re.names
	}
        if (p.re > 0) {
	  X.0 <- cbind(X.0, X.re.0)
	}
        out.pred <- predict.svcPGBinom(out.fit, X.0, coords.0, weights.0, verbose = FALSE)

        like.samples <- rep(NA, J.0)
	for (j in 1:J.0) {
          like.samples[j] <- mean(dbinom(y.0[j], weights.0[j], out.pred$psi.0.samples[, j]))
	}
	sum(log(like.samples), na.rm = TRUE)
      }
      model.deviance <- -2 * model.deviance
      # Return objects from cross-validation
      out$k.fold.deviance <- model.deviance
      stopImplicitCluster()
    } # cross-validation
  } # NNGP or GP
  class(out) <- "svcPGBinom"
  out$run.time <- proc.time() - ptm
  return(out)
}

