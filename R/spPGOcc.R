spPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
		    tuning, cov.model = 'exponential', NNGP = TRUE, 
		    n.neighbors = 15, search.type = 'cb', n.batch, 
		    batch.length, accept.rate = 0.43,
		    n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		    n.burn = round(.10 * n.batch * batch.length), 
		    n.thin = 1, k.fold, k.fold.threads = 1, 
		    k.fold.seed = 100, ...){

  ptm <- proc.time()

  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing the data\n");
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
  if (missing(occ.formula)) {
    stop("error: occ.formula must be specified")
  }
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  y <- as.matrix(data$y)
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- matrix(1, dim(y)[1], 1)
    } else {
      stop("error: occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!'det.covs' %in% names(data)) {
    if (det.formula == ~ 1) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
      }
      data$det.covs <- list(int = matrix(1, dim(y)[1], dim(y)[2]))
    } else {
      stop("error: det.covs must be specified in data for a detection model with covariates")
    }
  }
  if (!is.matrix(data$occ.covs) & !is.data.frame(data$occ.covs)) {
    stop("error: occ.covs must be a matrix or data frame")
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("error: coords must be a matrix or data frame")
  }
  coords <- as.matrix(data$coords)
  if (!is.list(data$det.covs)) {
    stop("error: det.covs must be a list of matrices, data frames, and/or vectors")
  }
  if (!missing(k.fold)) {
    if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
      stop("error: k.fold must be a single integer value >= 2")  
    }
  }

  # Checking missing values ---------------------------------------------
  y.na.test <- apply(y, 1, function(a) sum(!is.na(a)))
  if (sum(y.na.test == 0) > 0) {
    stop("error: some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
  }


  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[ord, , drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Occupancy covariates
    data$occ.covs <- data$occ.covs[ord, , drop = FALSE]
    for (i in 1:length(data$det.covs)) {
      if (!is.null(dim(data$det.covs[[i]]))) {
        data$det.covs[[i]] <- data$det.covs[[i]][ord, , drop = FALSE]
      } else {
        data$det.covs[[i]] <- data$det.covs[[i]][ord]
      }
    }
  }

  # First subset detection covariates to only use those that are included in the analysis. 
  data$det.covs <- data$det.covs[names(data$det.covs) %in% all.vars(det.formula)]
  # Null model support
  if (length(data$det.covs) == 0) {
    data$det.covs <- list(int = rep(1, dim(y)[1]))
  }
  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  data$det.covs <- data.frame(lapply(data$det.covs, function(a) unlist(c(a))))
  binom <- FALSE
  # Check if all detection covariates are at site level, and simplify the data
  # if necessary
  y.big <- y
  if (nrow(data$det.covs) == nrow(y)) {
   # Convert data to binomial form
   y <- apply(y, 1, sum, na.rm = TRUE) 
   binom <- TRUE
  }
  data$occ.covs <- as.data.frame(data$occ.covs)

  # Formula -------------------------------------------------------------
  # Occupancy -----------------------
  if (class(occ.formula) == 'formula') {
    tmp <- parseFormula(occ.formula, data$occ.covs)
    X <- as.matrix(tmp[[1]])
    x.names <- tmp[[2]]
  } else {
    stop("error: occ.formula is misspecified")
  }

  # Detection -----------------------
  if (class(det.formula) == 'formula') {
    tmp <- parseFormula(det.formula, data$det.covs)
    X.p <- as.matrix(tmp[[1]])
    x.p.names <- tmp[[2]]
    X.p.re <- as.matrix(tmp[[4]])
    x.p.re.names <- colnames(X.p.re)
  } else {
    stop("error: det.formula is misspecified")
  }

  # Get basic info from inputs ------------------------------------------
  # Number of sites
  J <- nrow(y.big)
  # Number of occupancy parameters 
  p.occ <- ncol(X)
  # Number of detection parameters
  p.det <- dim(X.p)[2]
  # Number of detection random effect parameters
  p.det.re <- ncol(X.p.re)
  # Number of latent detection random effect values
  n.det.re <- length(unlist(apply(X.p.re, 2, unique)))
  n.det.re.long <- apply(X.p.re, 2, function(a) length(unique(a)))
  # Number of replicates at each site
  n.rep <- apply(y.big, 1, function(a) sum(!is.na(a)))
  # Max number of repeat visits
  K.max <- max(n.rep)
  # Because I like K better than n.rep
  K <- n.rep
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

  # Get indices to map z to y -------------------------------------------
  if (!binom) {
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y.big))]
    # Subtract 1 for indices in C
    z.long.indx <- z.long.indx - 1
  } else {
    z.long.indx <- 0:(J - 1)
  }
  # y is stored in the following order: species, site, visit
  y <- c(y)
  names.long <- which(!is.na(y))
  # Remove missing observations when the covariate data are available but
  # there are missing detection-nondetection data. 
  if (nrow(X.p) == length(y)) {
    X.p <- X.p[!is.na(y), , drop = FALSE]  
  }
  if (nrow(X.p.re) == length(y) & p.det.re > 0) {
    X.p.re <- X.p.re[!is.na(y), , drop = FALSE]
  }
  y <- y[!is.na(y)]
  # Number of pseudoreplicates
  n.obs <- nrow(X.p)

  # Get random effect matrices all set ----------------------------------
  if (p.det.re > 1) {
    for (j in 2:p.det.re) {
      X.p.re[, j] <- X.p.re[, j] + max(X.p.re[, j - 1]) + 1
    }
  }
  lambda.p <- matrix(0, n.obs, n.det.re)
  if (p.det.re > 0) {
    for (i in 1:n.det.re) {
      lambda.p[which(X.p.re == (i - 1), arr.ind = TRUE)[, 1], i] <- 1
    }
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
    if (length(mu.alpha) != p.det & length(mu.alpha) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ",
        	     p.det, " with elements corresponding to alphas' mean", sep = ""))
      } else {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ",
        	     p.det, " or 1 with elements corresponding to alphas' mean", sep = ""))
      }
    }
    if (length(sigma.alpha) != p.det & length(sigma.alpha) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ",
      	   p.det, " with elements corresponding to alphas' variance", sep = ""))
      } else {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ",
      	   p.det, " or 1 with elements corresponding to alphas' variance", sep = ""))
      }
    }
    if (length(sigma.alpha) != p.det) {
      sigma.alpha <- rep(sigma.alpha, p.det)
    }
    if (length(mu.alpha) != p.det) {
      mu.alpha <- rep(mu.alpha, p.det)
    }
    Sigma.alpha <- sigma.alpha * diag(p.det)
  } else {
    if (verbose) {
      message("No prior specified for alpha.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.alpha <- rep(0, p.det)
    sigma.alpha <- rep(2.72, p.det)
    Sigma.alpha <- diag(p.det) * 2.72
  }

  # phi -----------------------------
  # Get distance matrix which is used if priors are not specified
  coords.D <- iDist(coords)
  lower.unif <- 3 / max(coords.D)
  upper.unif <- 3 / sort(unique(c(coords.D)))[2]
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
    phi.a <- 3 / max(coords.D)
    phi.b <- 3 / sort(unique(c(coords.D)))[2]
  }
  # sigma.sq -----------------------------
  if ("sigma.sq.ig" %in% names(priors)) {
    if (!is.vector(priors$sigma.sq.ig) | !is.atomic(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
      stop("error: sigma.sq.ig must be a vector of length 2 with elements corresponding to sigma.sq's shape and scale parameters")
    }
    sigma.sq.a <- priors$sigma.sq.ig[1]
    sigma.sq.b <- priors$sigma.sq.ig[2]
  } else {
    if (verbose) {
      message("No prior specified for sigma.sq.ig.\nSetting the shape and scale parameters to 2.\n")
    }
    sigma.sq.a <- 2
    sigma.sq.b <- 2
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

  # sigma.sq.p --------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.p.ig) | length(priors$sigma.sq.p.ig) != 2) {
        stop("error: sigma.sq.p.ig must be a list of length 2")
      }
      sigma.sq.p.a <- priors$sigma.sq.p.ig[[1]]
      sigma.sq.p.b <- priors$sigma.sq.p.ig[[2]]
      if (length(sigma.sq.p.a) != p.det.re & length(sigma.sq.p.a) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: sigma.sq.p.ig[[1]] must be a vector of length ", 
        	   p.det.re, " with elements corresponding to sigma.sq.ps' shape", sep = ""))
        } else {
          stop(paste("error: sigma.sq.p.ig[[1]] must be a vector of length ", 
        	   p.det.re, " or 1 with elements corresponding to sigma.sq.ps' shape", sep = ""))
        }
      }
      if (length(sigma.sq.p.b) != p.det.re & length(sigma.sq.p.b) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: sigma.sq.p.ig[[2]] must be a vector of length ", 
        	     p.det.re, " with elements corresponding to sigma.sq.ps' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.p.ig[[2]] must be a vector of length ", 
        	     p.det.re, " or 1 with elements corresponding to sigma.sq.ps' scale", sep = ""))
        }
      }
      if (length(sigma.sq.p.a) != p.det.re) {
        sigma.sq.p.a <- rep(sigma.sq.p.a, p.det.re)
      }
      if (length(sigma.sq.p.b) != p.det.re) {
        sigma.sq.p.b <- rep(sigma.sq.p.b, p.det.re)
      }
  }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.p.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.p.a <- rep(0.1, p.det.re)
      sigma.sq.p.b <- rep(0.1, p.det.re)
    }
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
    # Reorder the user supplied inits values
    if (NNGP) {
      z.inits <- z.inits[ord]
    }
    z.test <- apply(y.big, 1, max, na.rm = TRUE)
    init.test <- sum(z.inits < z.test)
    if (init.test > 0) {
      stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site.")
    }
  } else {
    # In correct order since you reordered y.
    z.inits <- apply(y.big, 1, max, na.rm = TRUE)
    if (verbose) {
      message("z.inits is not specified in initial values.\nSetting initial values based on observed data\n")
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
    if (length(alpha.inits) != p.det & length(alpha.inits) != 1) {
      if (p.det == 1) {
      stop(paste("error: initial values for alpha must be of length ", p.det,
      	   sep = ""))
      } else {
        stop(paste("error: initial values for alpha must be of length ", p.det, " or 1",
      	     sep = ""))
      }
    }
    if (length(alpha.inits) != p.det) {
      alpha.inits <- rep(alpha.inits, p.det)
    }
  } else {
    alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
    if (verbose) {
      message("alpha is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
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
    sigma.sq.inits <- rigamma(1, sigma.sq.a, sigma.sq.b)
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
  # sigma.sq.p ------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p" %in% names(inits)) {
      sigma.sq.p.inits <- inits[["sigma.sq.p"]]
      if (length(sigma.sq.p.inits) != p.det.re & length(sigma.sq.p.inits) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: initial values for sigma.sq.p must be of length ", p.det.re, 
      	     sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.p must be of length ", p.det.re, 
      	       " or 1", sep = ""))
          
        }
      }
      if (length(sigma.sq.p.inits) != p.det.re) {
        sigma.sq.p.inits <- rep(sigma.sq.p.inits, p.det.re)
      }
    } else {
      sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.p is not specified in initial values.\nSetting initial value to random values between 0.5 and 10\n")
      }
    }
    alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
    alpha.star.inits <- rnorm(n.det.re, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
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
  }
  # Log the tuning values since they are used in the AMCMC. 
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA

  if (!NNGP) {

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(n.obs) <- "integer"
    storage.mode(K) <- "double"
    storage.mode(coords.D) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(Sigma.alpha) <- "double"
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
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(n.burn) <- "integer"
    storage.mode(n.thin) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"

    if (p.det.re > 0) {

      storage.mode(p.det.re) <- "integer"
      storage.mode(X.p.re) <- "integer"
      storage.mode(n.det.re) <- "integer"
      storage.mode(n.det.re.long) <- "integer"
      storage.mode(sigma.sq.p.inits) <- "double"
      storage.mode(sigma.sq.p.a) <- "double"
      storage.mode(sigma.sq.p.b) <- "double"
      storage.mode(alpha.star.inits) <- "double"
      storage.mode(alpha.star.indx) <- "integer"
      storage.mode(lambda.p) <- "double"

      # Run the model in C    
      out <- .Call("spPGOccRE", y, X, X.p, coords.D, X.p.re, lambda.p, p.occ, p.det, 
		   p.det.re, J, n.obs, K, n.det.re, n.det.re.long, 
                   beta.inits, alpha.inits, sigma.sq.p.inits, 
		   alpha.star.inits, z.inits,
                   w.inits, phi.inits, sigma.sq.inits, nu.inits, z.long.indx, 
                   alpha.star.indx, mu.beta, mu.alpha, 
                   Sigma.beta, Sigma.alpha, phi.a, phi.b, 
                   sigma.sq.a, sigma.sq.b, nu.a, nu.b, sigma.sq.p.a, sigma.sq.p.b, 
		   tuning.c, cov.model.indx,
                   n.batch, batch.length, 
                   accept.rate, n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	   n.post.samples)

      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi')
      } else {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
      }
      out$z.samples <- mcmc(t(out$z.samples))
      out$psi.samples <- mcmc(t(out$psi.samples))
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$w.samples <- mcmc(t(out$w.samples))
      out$X <- X
      out$X.p <- X.p
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- batch.length * n.batch
      out$cov.model.indx <- cov.model.indx
      out$type <- "GP"
      out$coords <- coords
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE

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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
          }
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  y.big.fit <- y.big[-curr.set, , drop = FALSE]
	  y.big.0 <- y.big[curr.set, , drop = FALSE]
	  z.inits.fit <- z.inits[-curr.set]
	  w.inits.fit <- w.inits[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  J.fit <- nrow(X.fit)
	  J.0 <- nrow(X.0)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  n.obs.fit <- nrow(X.p.fit)
	  n.obs.0 <- nrow(X.p.0)
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  coords.D.fit <- coords.D[-curr.set, -curr.set, drop = FALSE]
	  coords.D.0 <- coords.D[curr.set, curr.set, drop = FALSE]
	  lambda.p.fit <- lambda.p[y.indx, , drop = FALSE]
	  lambda.p.0 <- lambda.p[!y.indx, , drop = FALSE]
	  X.p.re.fit <- X.p.re[y.indx, , drop = FALSE]
	  X.p.re.0 <- X.p.re[!y.indx, , drop = FALSE]
	  # Gotta be a better way, but will do for now. 
	  if (binom) {
            z.long.indx.fit <- 0:(J.fit - 1)
            z.0.long.indx <- 1:J.0
          } else {
	    z.long.indx.fit <- matrix(NA, J.fit, max(K.fit))
	    for (j in 1:J.fit) {
              z.long.indx.fit[j, 1:K.fit[j]] <- j  
            }
	    z.long.indx.fit <- c(z.long.indx.fit)
	    z.long.indx.fit <- z.long.indx.fit[!is.na(z.long.indx.fit)] - 1
	    z.0.long.indx <- matrix(NA, nrow(X.0), max(K.0))
	    for (j in 1:nrow(X.0)) {
              z.0.long.indx[j, 1:K.0[j]] <- j  
            }
	    z.0.long.indx <- c(z.0.long.indx)
	    z.0.long.indx <- z.0.long.indx[!is.na(z.0.long.indx)] 
	  }
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1

          storage.mode(y.fit) <- "double"
          storage.mode(z.inits.fit) <- "double"
          storage.mode(X.p.fit) <- "double"
          storage.mode(X.fit) <- "double"
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(n.obs.fit) <- "integer"
          storage.mode(K.fit) <- "double"
          storage.mode(coords.D.fit) <- "double"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(phi.inits) <- "double"
          storage.mode(sigma.sq.inits) <- "double"
          storage.mode(nu.inits) <- "double"
          storage.mode(w.inits.fit) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
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
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"
          storage.mode(p.det.re) <- "integer"
          storage.mode(X.p.re.fit) <- "integer"
          storage.mode(n.det.re) <- "integer"
          storage.mode(n.det.re.long) <- "integer"
          storage.mode(sigma.sq.p.inits) <- "double"
          storage.mode(sigma.sq.p.a) <- "double"
          storage.mode(sigma.sq.p.b) <- "double"
          storage.mode(alpha.star.inits) <- "double"
          storage.mode(alpha.star.indx) <- "integer"
          storage.mode(lambda.p.fit) <- "double"
      
	  out.fit <- .Call("spPGOccRE", y.fit, X.fit, X.p.fit, coords.D.fit, 
			   X.p.re.fit, lambda.p.fit, p.occ, p.det, 
		           p.det.re, J.fit, n.obs.fit, K.fit, n.det.re, n.det.re.long, 
                           beta.inits, alpha.inits, sigma.sq.p.inits, 
		           alpha.star.inits, z.inits.fit,
                           w.inits.fit, phi.inits, sigma.sq.inits, nu.inits, z.long.indx.fit, 
                           alpha.star.indx, mu.beta, mu.alpha, 
                           Sigma.beta, Sigma.alpha, phi.a, phi.b, 
                           sigma.sq.a, sigma.sq.b, nu.a, nu.b, sigma.sq.p.a, sigma.sq.p.b, 
		           tuning.c, cov.model.indx,
                           n.batch, batch.length, 
                           accept.rate, n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
          	           n.post.samples)
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- x.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi')
          } else {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi', 'nu')
          }
          out.fit$w.samples <- mcmc(t(out.fit$w.samples))
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
          out.fit$pRE <- TRUE
	  class(out.fit) <- "spPGOcc"

	  # Predict occurrence at new sites
	  out.pred <- predict.spPGOcc(out.fit, X.0, coords.0, verbose = FALSE, 
				      n.omp.threads = 1)
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples + 
				     lambda.p.0 %*% out.fit$alpha.star.samples))
	  if (binom) {
            like.samples <- matrix(NA, nrow(y.big.0), ncol(y.big.0))
            for (j in 1:nrow(X.p.0)) {
              for (k in 1:K.0[j]) {
                like.samples[j, k] <- mean(dbinom(y.big.0[j, k], 1,
					         p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
	      }
            }
          } else {
	    like.samples <- rep(NA, nrow(X.p.0))
	    for (j in 1:nrow(X.p.0)) {
              like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
            }
	  }
	  sum(log(like.samples), na.rm = TRUE)
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }
      class(out) <- "spPGOcc"
    }

    if (p.det.re == 0) {   

      # Run the model in C    
      out <- .Call("spPGOcc", y, X, X.p, coords.D, p.occ, p.det, J, n.obs, K, 
                   beta.inits, alpha.inits, z.inits,
                   w.inits, phi.inits, sigma.sq.inits, nu.inits, z.long.indx, 
                   mu.beta, mu.alpha, 
                   Sigma.beta, Sigma.alpha, phi.a, phi.b, 
                   sigma.sq.a, sigma.sq.b, nu.a, nu.b, tuning.c, cov.model.indx,
                   n.batch, batch.length, 
                   accept.rate, n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	 n.post.samples)

      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi')
      } else {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
      }
      out$z.samples <- mcmc(t(out$z.samples))
      out$psi.samples <- mcmc(t(out$psi.samples))
      out$w.samples <- mcmc(t(out$w.samples))
      out$X <- X
      out$X.p <- X.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- batch.length * n.batch
      out$cov.model.indx <- cov.model.indx
      out$type <- "GP"
      out$coords <- coords
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE

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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
          }
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  y.big.fit <- y.big[-curr.set, , drop = FALSE]
	  y.big.0 <- y.big[curr.set, , drop = FALSE]
	  z.inits.fit <- z.inits[-curr.set]
	  w.inits.fit <- w.inits[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  n.obs.fit <- nrow(X.p.fit)
	  n.obs.0 <- nrow(X.p.0)
	  J.fit <- nrow(X.fit)
	  J.0 <- nrow(X.0)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  coords.D.fit <- coords.D[-curr.set, -curr.set, drop = FALSE]
	  coords.D.0 <- coords.D[curr.set, curr.set, drop = FALSE]
	  # Gotta be a better way, but will do for now. 
	  if (binom) {
            z.long.indx.fit <- 0:(J.fit - 1)
            z.0.long.indx <- 1:J.0
          } else {
	    z.long.indx.fit <- matrix(NA, J.fit, max(K.fit))
	    for (j in 1:J.fit) {
              z.long.indx.fit[j, 1:K.fit[j]] <- j  
            }
	    z.long.indx.fit <- c(z.long.indx.fit)
	    z.long.indx.fit <- z.long.indx.fit[!is.na(z.long.indx.fit)] - 1
	    z.0.long.indx <- matrix(NA, nrow(X.0), max(K.0))
	    for (j in 1:nrow(X.0)) {
              z.0.long.indx[j, 1:K.0[j]] <- j  
            }
	    z.0.long.indx <- c(z.0.long.indx)
	    z.0.long.indx <- z.0.long.indx[!is.na(z.0.long.indx)] 
	  }
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1

          storage.mode(y.fit) <- "double"
          storage.mode(z.inits.fit) <- "double"
          storage.mode(X.p.fit) <- "double"
          storage.mode(X.fit) <- "double"
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(n.obs.fit) <- "integer"
          storage.mode(K.fit) <- "double"
          storage.mode(coords.D.fit) <- "double"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(phi.inits) <- "double"
          storage.mode(sigma.sq.inits) <- "double"
          storage.mode(nu.inits) <- "double"
          storage.mode(w.inits.fit) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
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
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"

          # Run the model in C
          out.fit <- .Call("spPGOcc", y.fit, X.fit, X.p.fit, coords.D.fit, p.occ, p.det, J.fit, 
			   n.obs.fit, K.fit, beta.inits, alpha.inits, z.inits.fit,
                           w.inits.fit, phi.inits, sigma.sq.inits, nu.inits, z.long.indx.fit, 
                           mu.beta, mu.alpha, 
                           Sigma.beta, Sigma.alpha, phi.a, phi.b, 
                           sigma.sq.a, sigma.sq.b, nu.a, nu.b, tuning.c, cov.model.indx,
                           n.batch, batch.length, 
                           accept.rate, n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
              	           n.post.samples)
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- x.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi')
          } else {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi', 'nu')
          }
          out.fit$w.samples <- mcmc(t(out.fit$w.samples))
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
          out.fit$pRE <- FALSE
	  class(out.fit) <- "spPGOcc"

	  # Predict occurrence at new sites
	  out.pred <- predict.spPGOcc(out.fit, X.0, coords.0, verbose = FALSE, 
				      n.omp.threads = 1)
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples))
	  if (binom) {
            like.samples <- matrix(NA, nrow(y.big.0), ncol(y.big.0))
            for (j in 1:nrow(X.p.0)) {
              for (k in 1:K.0[j]) {
                like.samples[j, k] <- mean(dbinom(y.big.0[j, k], 1,
					         p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
	      }
            }
          } else {
	    like.samples <- rep(NA, nrow(X.p.0))
	    for (j in 1:nrow(X.p.0)) {
              like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
            }
	  }
	  sum(log(like.samples), na.rm = TRUE)
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }
      class(out) <- "spPGOcc"
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
    
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"

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
    storage.mode(z.inits) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(n.obs) <- "integer"
    storage.mode(K) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(Sigma.alpha) <- "double"
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
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.indx) <- "integer"
    storage.mode(u.indx.lu) <- "integer"
    storage.mode(ui.indx) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(n.burn) <- "integer"
    storage.mode(n.thin) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"

    if (p.det.re > 0) {

      # Run the model in C --------------------------------------------------

      storage.mode(p.det.re) <- "integer"
      storage.mode(X.p.re) <- "integer"
      storage.mode(n.det.re) <- "integer"
      storage.mode(n.det.re.long) <- "integer"
      storage.mode(sigma.sq.p.inits) <- "double"
      storage.mode(sigma.sq.p.a) <- "double"
      storage.mode(sigma.sq.p.b) <- "double"
      storage.mode(alpha.star.inits) <- "double"
      storage.mode(alpha.star.indx) <- "integer"
      storage.mode(lambda.p) <- "double"

      out <- .Call("spPGOccNNGPRE", y, X, X.p, coords, X.p.re, lambda.p, 
		   p.occ, p.det, p.det.re, J, n.obs, K, n.det.re, n.det.re.long, 
          	   n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
      	           beta.inits, alpha.inits, sigma.sq.p.inits, 
		   alpha.star.inits, z.inits,
      	           phi.inits, sigma.sq.inits, nu.inits, 
          	   z.long.indx, alpha.star.indx, mu.beta, mu.alpha, 
      	           Sigma.beta, Sigma.alpha, phi.a, phi.b, 
      	           sigma.sq.a, sigma.sq.b, nu.a, nu.b, sigma.sq.p.a, 
		   sigma.sq.p.b, tuning.c, 
          	   cov.model.indx, n.batch, batch.length, 
      	           accept.rate, n.omp.threads, verbose, n.report, 
          	   n.burn, n.thin, n.post.samples)

      # Get everything back in the original order
      out$coords <- coords[order(ord), ]
      out$z.samples <- mcmc(t(out$z.samples[order(ord), , drop = FALSE]))
      out$X <- X[order(ord), , drop = FALSE]
      out$w.samples <- mcmc(t(out$w.samples[order(ord), , drop = FALSE]))
      out$psi.samples <- mcmc(t(out$psi.samples[order(ord), , drop = FALSE]))
      # Get detection covariate stuff in right order. Method of doing this
      # depends on if there are observation level covariates or not. 
      if (!binom) {
        tmp <- matrix(NA, J * K.max, p.det)
        tmp[names.long, ] <- X.p
        tmp <- array(tmp, dim = c(J, K.max, p.det))
        tmp <- tmp[order(ord), , ]
        out$X.p <- matrix(tmp, J * K.max, p.det)
        out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
        tmp <- matrix(NA, J * K.max, p.det.re)
        tmp[names.long, ] <- X.p.re
        tmp <- array(tmp, dim = c(J, K.max, p.det.re))
        tmp <- tmp[order(ord), , ]
        out$X.p.re <- matrix(tmp, J * K.max, p.det.re)
        out$X.p.re <- out$X.p.re[apply(out$X.p.re, 1, function(a) sum(is.na(a))) == 0, ]
        tmp <- matrix(NA, J * K.max, n.det.re)
        tmp[names.long, ] <- lambda.p
        tmp <- array(tmp, dim = c(J, K.max, n.det.re))
        tmp <- tmp[order(ord), , ]
        out$lambda.p <- matrix(tmp, J * K.max, n.det.re)
        out$lambda.p <- out$lambda.p[apply(out$lambda.p, 1, function(a) sum(is.na(a))) == 0, ]
      } else {
        out$X.p <- X.p[order(ord), , drop = FALSE]
        out$lambda.p <- lambda.p[order(ord), , drop = FALSE]
	out$X.p.re <- X.p.re[order(ord), , drop = FALSE]
      }
      out$y <- y.big[order(ord), , drop = FALSE]
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names

      # Return all the other good stuff. 
      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi')
      } else {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
      }
      out$call <- cl
      out$n.samples <- batch.length * n.batch
      out$n.neighbors <- n.neighbors
      out$cov.model.indx <- cov.model.indx
      out$type <- "NNGP"
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE

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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
          }
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  y.big.fit <- y.big[-curr.set, , drop = FALSE]
	  y.big.0 <- y.big[curr.set, , drop = FALSE]
	  z.inits.fit <- z.inits[-curr.set]
	  w.inits.fit <- w.inits[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  J.fit <- nrow(X.fit)
	  J.0 <- nrow(X.0)
	  n.obs.fit <- nrow(X.p.fit)
	  n.obs.0 <- nrow(X.p.0)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  lambda.p.fit <- lambda.p[y.indx, , drop = FALSE]
	  lambda.p.0 <- lambda.p[!y.indx, , drop = FALSE]
	  X.p.re.fit <- X.p.re[y.indx, , drop = FALSE]
	  X.p.re.0 <- X.p.re[!y.indx, , drop = FALSE]
	  # Don't need to reorder things, since they are already sorted by 
	  # the horizontal location in the coordinates. 
	  # Gotta be a better way, but will do for now. 
	  if (binom) {
            z.long.indx.fit <- 0:(J.fit - 1)
            z.0.long.indx <- 1:J.0
          } else {
	    z.long.indx.fit <- matrix(NA, J.fit, max(K.fit))
	    for (j in 1:J.fit) {
              z.long.indx.fit[j, 1:K.fit[j]] <- j  
            }
	    z.long.indx.fit <- c(z.long.indx.fit)
	    z.long.indx.fit <- z.long.indx.fit[!is.na(z.long.indx.fit)] - 1
	    z.0.long.indx <- matrix(NA, nrow(X.0), max(K.0))
	    for (j in 1:nrow(X.0)) {
              z.0.long.indx[j, 1:K.0[j]] <- j  
            }
	    z.0.long.indx <- c(z.0.long.indx)
	    z.0.long.indx <- z.0.long.indx[!is.na(z.0.long.indx)] 
	  }
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1
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
          storage.mode(n.obs.fit) <- "integer"
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "double"
          storage.mode(coords.fit) <- "double"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(phi.inits) <- "double"
          storage.mode(sigma.sq.inits) <- "double"
          storage.mode(nu.inits) <- "double"
          storage.mode(w.inits) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
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
          storage.mode(nn.indx.fit) <- "integer"
          storage.mode(nn.indx.lu.fit) <- "integer"
          storage.mode(u.indx.fit) <- "integer"
          storage.mode(u.indx.lu.fit) <- "integer"
          storage.mode(ui.indx.fit) <- "integer"
          storage.mode(n.neighbors) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"
          storage.mode(p.det.re) <- "integer"
          storage.mode(X.p.re.fit) <- "integer"
          storage.mode(n.det.re) <- "integer"
          storage.mode(n.det.re.long) <- "integer"
          storage.mode(sigma.sq.p.inits) <- "double"
          storage.mode(sigma.sq.p.a) <- "double"
          storage.mode(sigma.sq.p.b) <- "double"
          storage.mode(alpha.star.inits) <- "double"
          storage.mode(alpha.star.indx) <- "integer"
          storage.mode(lambda.p.fit) <- "double"
     
	  # Run the model in C 
          out.fit <- .Call("spPGOccNNGPRE", y.fit, X.fit, X.p.fit, 
			   coords.fit, X.p.re.fit, lambda.p.fit, 
		           p.occ, p.det, p.det.re, J.fit, n.obs.fit, 
			   K.fit, n.det.re, n.det.re.long, 
          	           n.neighbors, nn.indx.fit, nn.indx.lu.fit, 
			   u.indx.fit, u.indx.lu.fit, ui.indx.fit, 
      	                   beta.inits, alpha.inits, sigma.sq.p.inits, 
		           alpha.star.inits, z.inits.fit,
      	                   phi.inits, sigma.sq.inits, nu.inits, 
          	           z.long.indx.fit, alpha.star.indx, mu.beta, mu.alpha, 
      	                   Sigma.beta, Sigma.alpha, phi.a, phi.b, 
      	                   sigma.sq.a, sigma.sq.b, nu.a, nu.b, sigma.sq.p.a, 
		           sigma.sq.p.b, tuning.c, 
          	           cov.model.indx, n.batch, batch.length, 
      	                   accept.rate, n.omp.threads.fit, verbose.fit, n.report, 
          	           n.burn, n.thin, n.post.samples)
	  # Don't need to reorder anything since you're not returning 
	  # any of this stuff. 
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- x.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi')
          } else {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi', 'nu')
          }
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
          out.fit$pRE <- TRUE
	  class(out.fit) <- "spPGOcc"

	  # Predict occurrence at new sites
	  out.pred <- predict.spPGOcc(out.fit, X.0, coords.0, verbose = FALSE, 
				      n.omp.threads = 1)
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples + 
				     lambda.p.0 %*% out.fit$alpha.star.samples))
	  if (binom) {
            like.samples <- matrix(NA, nrow(y.big.0), ncol(y.big.0))
            for (j in 1:nrow(X.p.0)) {
              for (k in 1:K.0[j]) {
                like.samples[j, k] <- mean(dbinom(y.big.0[j, k], 1,
					         p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
	      }
            }
          } else {
	    like.samples <- rep(NA, nrow(X.p.0))
	    for (j in 1:nrow(X.p.0)) {
              like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
            }
	  }
	  sum(log(like.samples), na.rm = TRUE)
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "spPGOcc"
      
    }

    if (p.det.re == 0) {

      # Run the model in C --------------------------------------------------
      out <- .Call("spPGOccNNGP", y, X, X.p, coords, p.occ, p.det, J, n.obs, K, 
          	 n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
      	         beta.inits, alpha.inits, z.inits,
      	         phi.inits, sigma.sq.inits, nu.inits, 
          	 z.long.indx, mu.beta, mu.alpha, 
      	         Sigma.beta, Sigma.alpha, phi.a, phi.b, 
      	         sigma.sq.a, sigma.sq.b, nu.a, nu.b, tuning.c, 
          	 cov.model.indx, n.batch, batch.length, 
      	         accept.rate, n.omp.threads, verbose, n.report, 
          	 n.burn, n.thin, n.post.samples)

      # Get everything back in the original order
      out$coords <- coords[order(ord), ]
      out$z.samples <- mcmc(t(out$z.samples[order(ord), , drop = FALSE]))
      out$X <- X[order(ord), , drop = FALSE]
      out$w.samples <- mcmc(t(out$w.samples[order(ord), , drop = FALSE]))
      out$psi.samples <- mcmc(t(out$psi.samples[order(ord), , drop = FALSE]))
      if (!binom) {
        tmp <- matrix(NA, J * K.max, p.det)
        tmp[names.long, ] <- X.p
        tmp <- array(tmp, dim = c(J, K.max, p.det))
        tmp <- tmp[order(ord), , ]
        out$X.p <- matrix(tmp, J * K.max, p.det)
        out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
      } else {
        out$X.p <- X.p[order(ord), , drop = FALSE]
      }
      out$y <- y.big[order(ord), , drop = FALSE]

      # Return all the other good stuff. 
      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi')
      } else {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
      }
      out$call <- cl
      out$n.samples <- batch.length * n.batch
      out$n.neighbors <- n.neighbors
      out$cov.model.indx <- cov.model.indx
      out$type <- "NNGP"
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE

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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
          }
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  y.big.fit <- y.big[-curr.set, , drop = FALSE]
	  y.big.0 <- y.big[curr.set, , drop = FALSE]
	  z.inits.fit <- z.inits[-curr.set]
	  w.inits.fit <- w.inits[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  n.obs.fit <- nrow(X.p.fit)
	  n.obs.0 <- nrow(X.p.0)
	  J.fit <- nrow(X.fit)
	  J.0 <- nrow(X.0)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  # Don't need to reorder things, since they are already sorted by 
	  # the horizontal location in the coordinates. 

	  # Gotta be a better way, but will do for now. 
	  if (binom) {
            z.long.indx.fit <- 0:(J.fit - 1)
            z.0.long.indx <- 1:J.0
          } else {
	    z.long.indx.fit <- matrix(NA, J.fit, max(K.fit))
	    for (j in 1:J.fit) {
              z.long.indx.fit[j, 1:K.fit[j]] <- j  
            }
	    z.long.indx.fit <- c(z.long.indx.fit)
	    z.long.indx.fit <- z.long.indx.fit[!is.na(z.long.indx.fit)] - 1
	    z.0.long.indx <- matrix(NA, nrow(X.0), max(K.0))
	    for (j in 1:nrow(X.0)) {
              z.0.long.indx[j, 1:K.0[j]] <- j  
            }
	    z.0.long.indx <- c(z.0.long.indx)
	    z.0.long.indx <- z.0.long.indx[!is.na(z.0.long.indx)] 
	  }
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1
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
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(n.obs.fit) <- "integer"
          storage.mode(K.fit) <- "double"
          storage.mode(coords.fit) <- "double"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(phi.inits) <- "double"
          storage.mode(sigma.sq.inits) <- "double"
          storage.mode(nu.inits) <- "double"
          storage.mode(w.inits) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
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
          storage.mode(nn.indx.fit) <- "integer"
          storage.mode(nn.indx.lu.fit) <- "integer"
          storage.mode(u.indx.fit) <- "integer"
          storage.mode(u.indx.lu.fit) <- "integer"
          storage.mode(ui.indx.fit) <- "integer"
          storage.mode(n.neighbors) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"
     
	  # Run the model in C 
	  out.fit <- .Call("spPGOccNNGP", y.fit, X.fit, X.p.fit, coords.fit, 
			   p.occ, p.det, J.fit, n.obs.fit, K.fit, 
          	           n.neighbors, nn.indx.fit, nn.indx.lu.fit, u.indx.fit, 
			   u.indx.lu.fit, ui.indx.fit, 
      	                   beta.inits, alpha.inits, z.inits.fit,
      	                   phi.inits, sigma.sq.inits, nu.inits, 
          	           z.long.indx.fit, mu.beta, mu.alpha, 
      	                   Sigma.beta, Sigma.alpha, phi.a, phi.b, 
      	                   sigma.sq.a, sigma.sq.b, nu.a, nu.b, tuning.c, 
          	           cov.model.indx, n.batch, batch.length, 
      	                   accept.rate, n.omp.threads.fit, verbose.fit, n.report, 
          	           n.burn, n.thin, n.post.samples)

	  # Don't need to reorder anything since you're not returning 
	  # any of this stuff. 
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- x.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi')
          } else {
            colnames(out.fit$theta.samples) <- c('sigma.sq', 'phi', 'nu')
          }
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
          out.fit$pRE <- FALSE
	  class(out.fit) <- "spPGOcc"

	  # Predict occurrence at new sites
	  out.pred <- predict.spPGOcc(out.fit, X.0, coords.0, verbose = FALSE, 
				      n.omp.threads = 1)
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples))
	  if (binom) {
            like.samples <- matrix(NA, nrow(y.big.0), ncol(y.big.0))
            for (j in 1:nrow(X.p.0)) {
              for (k in 1:K.0[j]) {
                like.samples[j, k] <- mean(dbinom(y.big.0[j, k], 1,
					         p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
	      }
            }
          } else {
	    like.samples <- rep(NA, nrow(X.p.0))
	    for (j in 1:nrow(X.p.0)) {
              like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
            }
	  }
	  sum(log(like.samples), na.rm = TRUE)
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "spPGOcc"
      
    }
  }
  out$run.time <- proc.time() - ptm
  return(out)
}
