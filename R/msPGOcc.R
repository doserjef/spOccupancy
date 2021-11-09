msPGOcc <- function(occ.formula, det.formula, data, inits, priors,  
		    n.samples, n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		    n.burn = round(.10 * n.samples), n.thin = 1, 
		    k.fold, k.fold.threads = 1, k.fold.seed = 100, ...){

    ptm <- proc.time()

    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
    rigamma <- function(n, a, b){
      1/rgamma(n = n, shape = a, rate = b)
    }
 
    # Make it look nice
    if (verbose) {
      cat("----------------------------------------\n");
      cat("\tPreparing the data\n");
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
    if (!'y' %in% names(data)) {
      stop("error: detection-nondetection data y must be specified in data")
    }
    if (length(dim(data$y)) != 3) {
      stop("error: detection-nondetection data y must be a three-dimensional array with dimensions corresponding to species, sites, and replicates.")
    }
    y <- data$y
    sp.names <- attr(y, 'dimnames')[[1]]
    if (!'occ.covs' %in% names(data)) {
      if (occ.formula == ~ 1) {
        if (verbose) {
          message("occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
        }
        data$occ.covs <- matrix(1, dim(y)[2], 1)
      } else {
        stop("error: occ.covs must be specified in data for an occupancy model with covariates")
      }
    }
    if (!'det.covs' %in% names(data)) {
      if (det.formula == ~ 1) {
        if (verbose) {
          message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
	}
        data$det.covs <- list(int = matrix(1, dim(y)[2], dim(y)[3]))
      } else {
        stop("error: det.covs must be specified in data for a detection model with covariates")
      }
    }
    if (!missing(k.fold)) {
      if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
        stop("error: k.fold must be a single integer value >= 2")  
      }
    }

    # First subset detection covariates to only use those that are included in the analysis. 
    data$det.covs <- data$det.covs[names(data$det.covs) %in% all.vars(det.formula)]
    # Null model support
    if (length(data$det.covs) == 0) {
      data$det.covs <- list(int = rep(1, dim(y)[2]))
    }
    # Make both covariates a data frame. Unlist is necessary for when factors
    # are supplied. 
    data$det.covs <- data.frame(lapply(data$det.covs, function(a) unlist(c(a))))
    binom <- FALSE
    # Check if all detection covariates are at site level, and simplify the data
    # if necessary
    y.big <- y
    if (nrow(data$det.covs) == dim(y)[2]) {
     # Convert data to binomial form
     y <- apply(y, c(1, 2), sum, na.rm = TRUE) 
     binom <- TRUE
    }
    data$occ.covs <- as.data.frame(data$occ.covs)

    # Checking missing values ---------------------------------------------
    y.na.test <- apply(y.big, c(1, 2), function(a) sum(!is.na(a)))
    if (sum(y.na.test == 0) > 0) {
      stop("error: some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
    }

    # Formula -------------------------------------------------------------
    # Occupancy -----------------------
    if (missing(occ.formula)) {
      stop("error: occ.formula must be specified")
    }

    if (class(occ.formula) == 'formula') {
      tmp <- parseFormula(occ.formula, data$occ.covs)
      X <- as.matrix(tmp[[1]])
      X.re <- as.matrix(tmp[[4]])
      x.re.names <- colnames(X.re)
      x.names <- tmp[[2]]
    } else {
      stop("error: occ.formula is misspecified")
    }

    # Detection -----------------------
    if (missing(det.formula)) {
      stop("error: det.formula must be specified")
    }

    if (class(det.formula) == 'formula') {
      tmp <- parseFormula(det.formula, data$det.covs)
      X.p <- as.matrix(tmp[[1]])
      X.p.re <- as.matrix(tmp[[4]])
      x.p.re.names <- colnames(X.p.re)
      x.p.names <- tmp[[2]]
    } else {
      stop("error: det.formula is misspecified")
    }

    # Extract data from inputs --------------------------------------------
    # Number of species 
    N <- dim(y)[1]
    # Number of occupancy parameters 
    p.occ <- ncol(X)
    # Number of occupancy random effect parameters
    p.occ.re <- ncol(X.re)
    # Number of detection parameters
    p.det <- ncol(X.p)
    # Number of detection random effect parameters
    p.det.re <- ncol(X.p.re)
    # Number of latent occupancy random effect values
    n.occ.re <- length(unlist(apply(X.re, 2, unique)))
    n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
    # Number of latent detection random effect values
    n.det.re <- length(unlist(apply(X.p.re, 2, unique)))
    n.det.re.long <- apply(X.p.re, 2, function(a) length(unique(a)))
    # Number of sites
    J <- nrow(X)
    # Number of repeat visits
    # Note this assumes equivalent detection histories for all species. 
    # May want to change this at some point. 
    n.rep <- apply(y.big[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
    K.max <- max(n.rep)
    # Because I like K better than n.rep
    K <- n.rep
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
    # Get indices to map z to y -------------------------------------------
    if (!binom) {
      z.long.indx <- rep(1:J, K.max)
      z.long.indx <- z.long.indx[!is.na(c(y.big[1, , ]))]
      # Subtract 1 for indices in C
      z.long.indx <- z.long.indx - 1
    } else {
      z.long.indx <- 0:(J - 1)
    }
    # y is stored in the following order: species, site, visit
    y <- c(y)
    # Assumes the missing data are constant across species, which seems likely, 
    # but may eventually need some updating. 
    names.long <- which(!is.na(c(y.big[1, , ])))
    # Only need to check this when there are observation level covariates. 
    if (nrow(X.p) == length(y) / N) {
      if (!binom) {
        X.p <- X.p[!is.na(c(y.big[1, , ])), , drop = FALSE]
      }
    }
    if (nrow(X.p.re) == length(y) / N & p.det.re > 0) {
      if (!binom) {
        X.p.re <- X.p.re[!is.na(c(y.big[1, , ])), , drop = FALSE]
      }
    }
    y <- y[!is.na(y)]
    # Number of pseudoreplicates
    n.obs <- nrow(X.p)

    # Get random effect matrices all set ----------------------------------
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
      }
    }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        X.p.re[, j] <- X.p.re[, j] + max(X.p.re[, j - 1]) + 1
      }
    }
    lambda.psi <- matrix(0, J, n.occ.re)
    lambda.p <- matrix(0, n.obs, n.det.re)
    if (p.occ.re > 0) {
      for (i in 1:n.occ.re) {
        lambda.psi[which(X.re == (i - 1), arr.ind = TRUE)[, 1], i] <- 1
      }
    }
    if (p.det.re > 0) {
      for (i in 1:n.det.re) {
        lambda.p[which(X.p.re == (i - 1), arr.ind = TRUE)[, 1], i] <- 1
      }
    }

    # Separate out priors -------------------------------------------------
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

    # alpha.comm -----------------------
    if ("alpha.comm.normal" %in% names(priors)) {
      if (!is.list(priors$alpha.comm.normal) | length(priors$alpha.comm.normal) != 2) {
        stop("error: alpha.comm.normal must be a list of length 2")
      }
      mu.alpha.comm <- priors$alpha.comm.normal[[1]]
      sigma.alpha.comm <- priors$alpha.comm.normal[[2]]
      if (length(mu.alpha.comm) != p.det & length(mu.alpha.comm) != 1) {
        if (p.det == 1) {
          stop(paste("error: alpha.comm.normal[[1]] must be a vector of length ",
	  	     p.det, " with elements corresponding to alpha.comms' mean", sep = ""))
	} else {
          stop(paste("error: alpha.comm.normal[[1]] must be a vector of length ",
	  	     p.det, " or 1 with elements corresponding to alpha.comms' mean", sep = ""))
        }
      }
      if (length(sigma.alpha.comm) != p.det & length(sigma.alpha.comm) != 1) {
        if (p.det == 1) {
          stop(paste("error: alpha.comm.normal[[2]] must be a vector of length ",
		   p.det, " with elements corresponding to alpha.comms' variance", sep = ""))
        } else {
          stop(paste("error: alpha.comm.normal[[2]] must be a vector of length ",
		   p.det, " or 1 with elements corresponding to alpha.comms' variance", sep = ""))
        }
      }
      if (length(sigma.alpha.comm) != p.det) {
        sigma.alpha.comm <- rep(sigma.alpha.comm, p.det)
      }
      if (length(mu.alpha.comm) != p.det) {
        mu.alpha.comm <- rep(mu.alpha.comm, p.det)
      }
      Sigma.alpha.comm <- sigma.alpha.comm * diag(p.det)
    } else {
      if (verbose) {
        message("No prior specified for alpha.comm.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
      }
      mu.alpha.comm <- rep(0, p.det)
      sigma.alpha.comm <- rep(2.72, p.det)
      Sigma.alpha.comm <- diag(p.det) * 2.72
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

    # tau.sq.alpha -----------------------
    if ("tau.sq.alpha.ig" %in% names(priors)) {
      if (!is.list(priors$tau.sq.alpha.ig) | length(priors$tau.sq.alpha.ig) != 2) {
        stop("error: tau.sq.alpha.ig must be a list of length 2")
      }
      tau.sq.alpha.a <- priors$tau.sq.alpha.ig[[1]]
      tau.sq.alpha.b <- priors$tau.sq.alpha.ig[[2]]
      if (length(tau.sq.alpha.a) != p.det & length(tau.sq.alpha.a) != 1) {
        if (p.det == 1) {
          stop(paste("error: tau.sq.alpha.ig[[1]] must be a vector of length ", 
		   p.det, " with elements corresponding to tau.sq.alphas' shape", sep = ""))
	} else {
          stop(paste("error: tau.sq.alpha.ig[[1]] must be a vector of length ", 
		   p.det, " or 1 with elements corresponding to tau.sq.alphas' shape", sep = ""))
        }
      }
      if (length(tau.sq.alpha.b) != p.det & length(tau.sq.alpha.b) != 1) {
        if (p.det == 1) {
          stop(paste("error: tau.sq.alpha.ig[[2]] must be a vector of length ", 
		   p.det, " with elements corresponding to tau.sq.alphas' scale", sep = ""))
	} else {
          stop(paste("error: tau.sq.alpha.ig[[2]] must be a vector of length ", 
		   p.det, " or 1 with elements corresponding to tau.sq.alphas' scale", sep = ""))
        }
      }
      if (length(tau.sq.alpha.a) != p.det) {
        tau.sq.alpha.a <- rep(tau.sq.alpha.a, p.det)
      }
      if (length(tau.sq.alpha.b) != p.det) {
        tau.sq.alpha.b <- rep(tau.sq.alpha.b, p.det)
      }
    } else {
      if (verbose) {	    
        message("No prior specified for tau.sq.alpha.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      tau.sq.alpha.a <- rep(0.1, p.det)
      tau.sq.alpha.b <- rep(0.1, p.det)
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
      if (!is.matrix(z.inits)) {
        stop(paste("error: initial values for z must be a matrix with dimensions ", 
		   N, " x ", J, sep = ""))
      }
      if (nrow(z.inits) != N | ncol(z.inits) != J) {
        stop(paste("error: initial values for z must be a matrix with dimensions ", 
		   N, " x ", J, sep = ""))
      }
      z.test <- apply(y.big, c(1, 2), max, na.rm = TRUE)
      init.test <- sum(z.inits < z.test)
      if (init.test > 0) {
        stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site.")
      }
    } else {
      z.inits <- apply(y.big, c(1, 2), max, na.rm = TRUE)
      if (verbose) {
        message("z is not specified in initial values.\nSetting initial values based on observed data\n")
      }
    }
    # beta.comm -----------------------
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
    # alpha.comm -----------------------
    if ("alpha.comm" %in% names(inits)) {
      alpha.comm.inits <- inits[["alpha.comm"]]
      if (length(alpha.comm.inits) != p.det & length(alpha.comm.inits) != 1) {
        if (p.det == 1) {
          stop(paste("error: initial values for alpha.comm must be of length ", p.det, 
		   sep = ""))
	} else {
          stop(paste("error: initial values for alpha.comm must be of length ", p.det, 
		   , " or 1", sep = ""))
        }
      }
      if (length(alpha.comm.inits) != p.det) {
        alpha.comm.inits <- rep(alpha.comm.inits, p.det)
      }
    } else {
      alpha.comm.inits <- rnorm(p.det, mu.alpha.comm, sqrt(sigma.alpha.comm))
      if (verbose) {
        message('alpha.comm is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
      }
    }
    # tau.sq.beta ------------------------
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
    # tau.sq.alpha -----------------------
    if ("tau.sq.alpha" %in% names(inits)) {
      tau.sq.alpha.inits <- inits[["tau.sq.alpha"]]
      if (length(tau.sq.alpha.inits) != p.det & length(tau.sq.alpha.inits) != 1) {
        if (p.det == 1) {
          stop(paste("error: initial values for tau.sq.alpha must be of length ", p.det, 
		   sep = ""))
	} else {
          stop(paste("error: initial values for tau.sq.alpha must be of length ", p.det, 
		   " or 1", sep = ""))
        }
      }
      if (length(tau.sq.alpha.inits) != p.det) {
        tau.sq.alpha.inits <- rep(tau.sq.alpha.inits, p.det)
      }
    } else {
      tau.sq.alpha.inits <- runif(p.det, 0.5, 10)
      if (verbose) {
        message('tau.sq.alpha is not specified in initial values.\nSetting to initial values to random values between 0.5 and 10\n')
      }
    }
    # beta ----------------------------
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
    # alpha ----------------------------
    if ("alpha" %in% names(inits)) {
      alpha.inits <- inits[["alpha"]]
      if (is.matrix(alpha.inits)) {
        if (ncol(alpha.inits) != p.det | nrow(alpha.inits) != N) {
          stop(paste("error: initial values for alpha must be a matrix with dimensions ", 
          	   N, "x", p.det, " or a single numeric value", sep = ""))
        }
      }
      if (!is.matrix(alpha.inits) & length(alpha.inits) != 1) {
        stop(paste("error: initial values for alpha must be a matrix with dimensions ", 
		   N, " x ", p.det, " or a single numeric value", sep = ""))
      }
      if (length(alpha.inits) == 1) {
        alpha.inits <- matrix(alpha.inits, N, p.det)
      }
    } else {
      alpha.inits <- matrix(rnorm(N * p.det, alpha.comm.inits, sqrt(tau.sq.alpha.inits)), N, p.det)
      if (verbose) {
        message('alpha is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n')
      }
    }

    # sigma.sq.psi -------------------
    if (p.occ.re > 0) {
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
      beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
      # Starting values for all species 
      beta.star.inits <- rep(beta.star.inits, N)
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
          message("sigma.sq.p is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
        }
      }
      alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
      alpha.star.inits <- rnorm(n.det.re, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
      alpha.star.inits <- rep(alpha.star.inits, N)
    }

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
    storage.mode(N) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(beta.comm.inits) <- "double"
    storage.mode(alpha.comm.inits) <- "double"
    storage.mode(tau.sq.beta.inits) <- "double"
    storage.mode(tau.sq.alpha.inits) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(mu.alpha.comm) <- "double"
    storage.mode(Sigma.alpha.comm) <- "double"
    storage.mode(tau.sq.beta.a) <- "double"
    storage.mode(tau.sq.beta.b) <- "double"
    storage.mode(tau.sq.alpha.a) <- "double"
    storage.mode(tau.sq.alpha.b) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(n.burn) <- "integer"
    storage.mode(n.thin) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"

    # Set model.deviance to NA for returning when no cross-validation
    model.deviance <- NA

    if (p.occ.re > 0 & p.det.re == 0) {
      storage.mode(p.occ.re) <- "integer"
      storage.mode(X.re) <- "integer"
      storage.mode(n.occ.re) <- "integer"
      storage.mode(n.occ.re.long) <- "integer"
      storage.mode(sigma.sq.psi.inits) <- "double"
      storage.mode(sigma.sq.psi.a) <- "double"
      storage.mode(sigma.sq.psi.b) <- "double"
      storage.mode(beta.star.inits) <- "double"
      storage.mode(beta.star.indx) <- "integer"
      storage.mode(lambda.psi) <- "double"

      out <- .Call("msPGOccREOcc", y, X, X.p, X.re, 
		   lambda.psi, p.occ, p.det, p.occ.re, 
		   J, n.obs, K, N, n.occ.re, n.occ.re.long, 
          	   beta.inits, alpha.inits, z.inits,
          	   beta.comm.inits, 
          	   alpha.comm.inits, tau.sq.beta.inits, 
          	   tau.sq.alpha.inits, sigma.sq.psi.inits, 
		   beta.star.inits, 
		   z.long.indx, beta.star.indx, mu.beta.comm, 
          	   mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	   tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	   tau.sq.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
		   n.samples, n.omp.threads, 
		   verbose, n.report, n.burn, n.thin, n.post.samples)

      out$beta.comm.samples <- mcmc(t(out$beta.comm.samples))
      colnames(out$beta.comm.samples) <- x.names
      out$alpha.comm.samples <- mcmc(t(out$alpha.comm.samples))
      colnames(out$alpha.comm.samples) <- x.p.names
      out$tau.sq.beta.samples <- mcmc(t(out$tau.sq.beta.samples))
      colnames(out$tau.sq.beta.samples) <- x.names
      out$tau.sq.alpha.samples <- mcmc(t(out$tau.sq.alpha.samples))
      colnames(out$tau.sq.alpha.samples) <- x.p.names
      if (is.null(sp.names)) {
        sp.names <- paste('sp', 1:N, sep = '')
      }
      coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- coef.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
      colnames(out$alpha.samples) <- coef.names.det
      out$sigma.sq.psi.samples <- mcmc(t(out$sigma.sq.psi.samples))
      colnames(out$sigma.sq.psi.samples) <- x.re.names
      out$beta.star.samples <- mcmc(t(out$beta.star.samples))
      tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
      beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
      beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      out$X <- X
      out$X.p <- X.p
      out$X.re <- X.re
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE
      out$psiRE <- TRUE

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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
	    y.fit <- y[rep(y.indx, N), drop = FALSE]
	    y.0 <- y[rep(y.indx, N), drop = FALSE]
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
	    y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	    y.fit <- y.fit[!is.na(y.fit)]
	    y.0 <- c(y.big[, curr.set, , drop = FALSE])
	    y.0 <- y.0[!is.na(y.0)]
          }
	  z.inits.fit <- z.inits[, -curr.set]
	  y.big.fit <- y.big[, -curr.set, , drop = FALSE]
	  y.big.0 <- y.big[, curr.set, , drop = FALSE]
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
	  lambda.psi.fit <- lambda.psi[-curr.set, , drop = FALSE]
	  X.re.fit <- X.re[-curr.set, , drop = FALSE]
	  X.re.0 <- X.re[curr.set, , drop = FALSE]
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
          storage.mode(K.fit) <- "double"
	  storage.mode(n.obs.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(beta.comm.inits) <- "double"
          storage.mode(alpha.comm.inits) <- "double"
          storage.mode(tau.sq.beta.inits) <- "double"
          storage.mode(tau.sq.alpha.inits) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta.comm) <- "double"
          storage.mode(Sigma.beta.comm) <- "double"
          storage.mode(mu.alpha.comm) <- "double"
          storage.mode(Sigma.alpha.comm) <- "double"
          storage.mode(tau.sq.beta.a) <- "double"
          storage.mode(tau.sq.beta.b) <- "double"
          storage.mode(tau.sq.alpha.a) <- "double"
          storage.mode(tau.sq.alpha.b) <- "double"
          storage.mode(n.samples) <- "integer"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"
          storage.mode(p.occ.re) <- "integer"
          storage.mode(X.re.fit) <- "integer"
          storage.mode(n.occ.re) <- "integer"
          storage.mode(n.occ.re.long) <- "integer"
          storage.mode(sigma.sq.psi.inits) <- "double"
          storage.mode(sigma.sq.psi.a) <- "double"
          storage.mode(sigma.sq.psi.b) <- "double"
          storage.mode(beta.star.inits) <- "double"
          storage.mode(beta.star.indx) <- "integer"
          storage.mode(lambda.psi.fit) <- "double"

          # Run the model in C
          out.fit <- .Call("msPGOccREOcc", y.fit, X.fit, X.p.fit, X.re.fit, 
		           lambda.psi.fit, p.occ, p.det, p.occ.re, 
		           J.fit, n.obs.fit, K.fit, N, n.occ.re, n.occ.re.long, 
          	           beta.inits, alpha.inits, z.inits.fit,
          	           beta.comm.inits, 
          	           alpha.comm.inits, tau.sq.beta.inits, 
          	           tau.sq.alpha.inits, sigma.sq.psi.inits, 
		           beta.star.inits, 
		           z.long.indx.fit, beta.star.indx, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
		           n.samples, n.omp.threads.fit, 
		           verbose.fit, n.report, n.burn, n.thin, n.post.samples)

          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$sigma.sq.psi.samples <- mcmc(t(out.fit$sigma.sq.psi.samples))
          colnames(out.fit$sigma.sq.psi.samples) <- x.re.names
          out.fit$beta.star.samples <- mcmc(t(out.fit$beta.star.samples))
          tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
          beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
          beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
          colnames(out.fit$beta.star.samples) <- beta.star.names
          out.fit$X <- X
	  out.fit$X.re <- X.re
          out.fit$y <- y.big
          out.fit$n.post <- n.post.samples
          out.fit$psiRE <- TRUE
	  class(out.fit) <- "msPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.msPGOcc(out.fit, cbind(X.0, X.re.0))

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
	  if (binom) {
            like.samples <- array(NA, c(N, nrow(X.p.0), dim(y.big.0)[3]))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
              for (j in 1:nrow(X.p.0)) {
                for (k in 1:K.0[j]) {
                  like.samples[q, j, k] <- mean(dbinom(y.big.0[q, j, k], 1,
	          			         p.0.samples[j, q, k] * out.pred$z.0.samples[, q, z.0.long.indx[j]]))
	        }
              }
	    }
          } else {
	    like.samples <- matrix(NA, N, nrow(X.p.0))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
	      for (j in 1:nrow(X.p.0)) {
                like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
	          				p.0.samples[j, q, ] * 
	          			        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
              }
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a), na.rm = TRUE))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "msPGOcc"

    }

    if (p.occ.re == 0 & p.det.re > 0) {
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

      out <- .Call("msPGOccREDet", y, X, X.p, X.p.re, 
		   lambda.p, p.occ, p.det, p.det.re, 
		   J, n.obs, K, N, n.det.re, n.det.re.long,
          	   beta.inits, alpha.inits, z.inits,
          	   beta.comm.inits, 
          	   alpha.comm.inits, tau.sq.beta.inits, 
          	   tau.sq.alpha.inits,  
		   sigma.sq.p.inits, alpha.star.inits, 
		   z.long.indx,
		   alpha.star.indx, mu.beta.comm, 
          	   mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	   tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	   tau.sq.alpha.b, sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads, 
		   verbose, n.report, n.burn, n.thin, n.post.samples)

      out$beta.comm.samples <- mcmc(t(out$beta.comm.samples))
      colnames(out$beta.comm.samples) <- x.names
      out$alpha.comm.samples <- mcmc(t(out$alpha.comm.samples))
      colnames(out$alpha.comm.samples) <- x.p.names
      out$tau.sq.beta.samples <- mcmc(t(out$tau.sq.beta.samples))
      colnames(out$tau.sq.beta.samples) <- x.names
      out$tau.sq.alpha.samples <- mcmc(t(out$tau.sq.alpha.samples))
      colnames(out$tau.sq.alpha.samples) <- x.p.names
      if (is.null(sp.names)) {
        sp.names <- paste('sp', 1:N, sep = '')
      }
      coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- coef.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
      colnames(out$alpha.samples) <- coef.names.det
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      out$X <- X
      out$X.p <- X.p
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE
      out$psiRE <- FALSE

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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
	    y.fit <- y[rep(y.indx, N), drop = FALSE]
	    y.0 <- y[rep(y.indx, N), drop = FALSE]
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
	    y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	    y.fit <- y.fit[!is.na(y.fit)]
	    y.0 <- c(y.big[, curr.set, , drop = FALSE])
	    y.0 <- y.0[!is.na(y.0)]
          }
	  z.inits.fit <- z.inits[, -curr.set]
	  y.big.fit <- y.big[, -curr.set, , drop = FALSE]
	  y.big.0 <- y.big[, curr.set, , drop = FALSE]
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
          storage.mode(K.fit) <- "double"
	  storage.mode(n.obs.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(beta.comm.inits) <- "double"
          storage.mode(alpha.comm.inits) <- "double"
          storage.mode(tau.sq.beta.inits) <- "double"
          storage.mode(tau.sq.alpha.inits) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta.comm) <- "double"
          storage.mode(Sigma.beta.comm) <- "double"
          storage.mode(mu.alpha.comm) <- "double"
          storage.mode(Sigma.alpha.comm) <- "double"
          storage.mode(tau.sq.beta.a) <- "double"
          storage.mode(tau.sq.beta.b) <- "double"
          storage.mode(tau.sq.alpha.a) <- "double"
          storage.mode(tau.sq.alpha.b) <- "double"
          storage.mode(n.samples) <- "integer"
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

          out.fit <- .Call("msPGOccREDet", y.fit, X.fit, X.p.fit, X.p.re.fit, 
		           lambda.p.fit, p.occ, p.det, p.det.re, 
		           J.fit, n.obs.fit, K.fit, N, n.det.re, n.det.re.long,
          	           beta.inits, alpha.inits, z.inits.fit,
          	           beta.comm.inits, 
          	           alpha.comm.inits, tau.sq.beta.inits, 
          	           tau.sq.alpha.inits,  
		           sigma.sq.p.inits, alpha.star.inits, 
		           z.long.indx.fit,
		           alpha.star.indx, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads.fit, 
		           verbose.fit, n.report, n.burn, n.thin, n.post.samples)

          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$X <- X
          out.fit$y <- y.big
          out.fit$n.post <- n.post.samples
          out.fit$psiRE <- FALSE
	  class(out.fit) <- "msPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.msPGOcc(out.fit, X.0)

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
          sp.re.indx <- rep(1:N, each = nrow(out.fit$alpha.star.samples) / N)
	  if (binom) {
            like.samples <- array(NA, c(N, nrow(X.p.0), dim(y.big.0)[3]))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] + 
					      lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
              for (j in 1:nrow(X.p.0)) {
                for (k in 1:K.0[j]) {
                  like.samples[q, j, k] <- mean(dbinom(y.big.0[q, j, k], 1,
	          			         p.0.samples[j, q, k] * out.pred$z.0.samples[, q, z.0.long.indx[j]]))
	        }
              }
	    }
          } else {
	    like.samples <- matrix(NA, N, nrow(X.p.0))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] + 
					      lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
	      for (j in 1:nrow(X.p.0)) {
                like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
	          				p.0.samples[j, q, ] * 
	          			        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
              }
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a), na.rm = TRUE))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "msPGOcc"

    }

    if (p.occ.re > 0 & p.det.re > 0) {

      storage.mode(p.occ.re) <- "integer"
      storage.mode(p.det.re) <- "integer"
      storage.mode(X.re) <- "integer"
      storage.mode(X.p.re) <- "integer"
      storage.mode(n.occ.re) <- "integer"
      storage.mode(n.det.re) <- "integer"
      storage.mode(n.occ.re.long) <- "integer"
      storage.mode(n.det.re.long) <- "integer"
      storage.mode(sigma.sq.psi.inits) <- "double"
      storage.mode(sigma.sq.p.inits) <- "double"
      storage.mode(sigma.sq.psi.a) <- "double"
      storage.mode(sigma.sq.psi.b) <- "double"
      storage.mode(sigma.sq.p.a) <- "double"
      storage.mode(sigma.sq.p.b) <- "double"
      storage.mode(beta.star.inits) <- "double"
      storage.mode(beta.star.indx) <- "integer"
      storage.mode(alpha.star.inits) <- "double"
      storage.mode(alpha.star.indx) <- "integer"
      storage.mode(lambda.psi) <- "double"
      storage.mode(lambda.p) <- "double"

      out <- .Call("msPGOccREBoth", y, X, X.p, X.re, X.p.re, 
		   lambda.psi, lambda.p, p.occ, p.det, p.occ.re, p.det.re, 
		   J, n.obs, K, N, n.occ.re, n.det.re, n.occ.re.long, n.det.re.long,
          	   beta.inits, alpha.inits, z.inits,
          	   beta.comm.inits, 
          	   alpha.comm.inits, tau.sq.beta.inits, 
          	   tau.sq.alpha.inits, sigma.sq.psi.inits, 
		   sigma.sq.p.inits, beta.star.inits, 
		   alpha.star.inits, z.long.indx, beta.star.indx, 
		   alpha.star.indx, mu.beta.comm, 
          	   mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	   tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	   tau.sq.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
		   sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads, 
		   verbose, n.report, n.burn, n.thin, n.post.samples)

      out$beta.comm.samples <- mcmc(t(out$beta.comm.samples))
      colnames(out$beta.comm.samples) <- x.names
      out$alpha.comm.samples <- mcmc(t(out$alpha.comm.samples))
      colnames(out$alpha.comm.samples) <- x.p.names
      out$tau.sq.beta.samples <- mcmc(t(out$tau.sq.beta.samples))
      colnames(out$tau.sq.beta.samples) <- x.names
      out$tau.sq.alpha.samples <- mcmc(t(out$tau.sq.alpha.samples))
      colnames(out$tau.sq.alpha.samples) <- x.p.names
      if (is.null(sp.names)) {
        sp.names <- paste('sp', 1:N, sep = '')
      }
      coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- coef.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
      colnames(out$alpha.samples) <- coef.names.det
      out$sigma.sq.psi.samples <- mcmc(t(out$sigma.sq.psi.samples))
      colnames(out$sigma.sq.psi.samples) <- x.re.names
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$beta.star.samples <- mcmc(t(out$beta.star.samples))
      tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
      beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
      beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      out$X <- X
      out$X.p <- X.p
      out$X.re <- X.re
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE
      out$psiRE <- TRUE


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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
	    y.fit <- y[rep(y.indx, N), drop = FALSE]
	    y.0 <- y[rep(y.indx, N), drop = FALSE]
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
	    y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	    y.fit <- y.fit[!is.na(y.fit)]
	    y.0 <- c(y.big[, curr.set, , drop = FALSE])
	    y.0 <- y.0[!is.na(y.0)]
          }
	  z.inits.fit <- z.inits[, -curr.set]
	  y.big.fit <- y.big[, -curr.set, , drop = FALSE]
	  y.big.0 <- y.big[, curr.set, , drop = FALSE]
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
	  lambda.psi.fit <- lambda.psi[-curr.set, , drop = FALSE]
	  lambda.psi.0 <- lambda.psi[curr.set, , drop = FALSE]
	  X.re.fit <- X.re[-curr.set, , drop = FALSE]
	  X.re.0 <- X.re[curr.set, , drop = FALSE]
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
          storage.mode(K.fit) <- "double"
	  storage.mode(n.obs.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(beta.comm.inits) <- "double"
          storage.mode(alpha.comm.inits) <- "double"
          storage.mode(tau.sq.beta.inits) <- "double"
          storage.mode(tau.sq.alpha.inits) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta.comm) <- "double"
          storage.mode(Sigma.beta.comm) <- "double"
          storage.mode(mu.alpha.comm) <- "double"
          storage.mode(Sigma.alpha.comm) <- "double"
          storage.mode(tau.sq.beta.a) <- "double"
          storage.mode(tau.sq.beta.b) <- "double"
          storage.mode(tau.sq.alpha.a) <- "double"
          storage.mode(tau.sq.alpha.b) <- "double"
          storage.mode(n.samples) <- "integer"
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
          storage.mode(sigma.sq.psi.inits) <- "double"
          storage.mode(sigma.sq.psi.a) <- "double"
          storage.mode(sigma.sq.psi.b) <- "double"
          storage.mode(beta.star.inits) <- "double"
          storage.mode(beta.star.indx) <- "integer"
          storage.mode(lambda.psi.fit) <- "double"

          # Run the model in C
          out.fit <- .Call("msPGOccREBoth", y.fit, X.fit, X.p.fit, X.re.fit, X.p.re.fit, 
		           lambda.psi.fit, lambda.p.fit, p.occ, p.det, p.occ.re, p.det.re, 
		           J.fit, n.obs.fit, K.fit, N, n.occ.re, n.det.re, n.occ.re.long, n.det.re.long,
          	           beta.inits, alpha.inits, z.inits.fit,
          	           beta.comm.inits, 
          	           alpha.comm.inits, tau.sq.beta.inits, 
          	           tau.sq.alpha.inits, sigma.sq.psi.inits, 
		           sigma.sq.p.inits, beta.star.inits, 
		           alpha.star.inits, z.long.indx.fit, beta.star.indx, 
		           alpha.star.indx, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
		           sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads.fit, 
		           verbose.fit, n.report, n.burn, n.thin, n.post.samples)

          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$sigma.sq.psi.samples <- mcmc(t(out.fit$sigma.sq.psi.samples))
          colnames(out.fit$sigma.sq.psi.samples) <- x.re.names
          out.fit$beta.star.samples <- mcmc(t(out.fit$beta.star.samples))
          tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
          beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
          beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
          colnames(out.fit$beta.star.samples) <- beta.star.names
          out.fit$X <- X
	  out.fit$X.re <- X.re
          out.fit$y <- y.big
          out.fit$n.post <- n.post.samples
          out.fit$psiRE <- TRUE
	  class(out.fit) <- "msPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.msPGOcc(out.fit, cbind(X.0, X.re.0))

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
          sp.re.indx <- rep(1:N, each = nrow(out.fit$alpha.star.samples) / N)
	  if (binom) {
            like.samples <- array(NA, c(N, nrow(X.p.0), dim(y.big.0)[3]))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] + 
					      lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
              for (j in 1:nrow(X.p.0)) {
                for (k in 1:K.0[j]) {
                  like.samples[q, j, k] <- mean(dbinom(y.big.0[q, j, k], 1,
	          			         p.0.samples[j, q, k] * out.pred$z.0.samples[, q, z.0.long.indx[j]]))
	        }
              }
	    }
          } else {
	    like.samples <- matrix(NA, N, nrow(X.p.0))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] + 
					      lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
	      for (j in 1:nrow(X.p.0)) {
                like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
	          				p.0.samples[j, q, ] * 
	          			        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
              }
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a), na.rm = TRUE))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "msPGOcc"

    }

    if (p.occ.re == 0 & p.det.re == 0) {
      # Run the model in C    
      out <- .Call("msPGOcc", y, X, X.p, p.occ, p.det, J, n.obs, K, N, 
          	 beta.inits, alpha.inits, z.inits,
          	 beta.comm.inits, 
          	 alpha.comm.inits, tau.sq.beta.inits, 
          	 tau.sq.alpha.inits, z.long.indx, mu.beta.comm, 
          	 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	 tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	 tau.sq.alpha.b, n.samples, n.omp.threads, verbose, n.report, 
          	 n.burn, n.thin, n.post.samples)

      out$beta.comm.samples <- mcmc(t(out$beta.comm.samples))
      colnames(out$beta.comm.samples) <- x.names
      out$alpha.comm.samples <- mcmc(t(out$alpha.comm.samples))
      colnames(out$alpha.comm.samples) <- x.p.names
      out$tau.sq.beta.samples <- mcmc(t(out$tau.sq.beta.samples))
      colnames(out$tau.sq.beta.samples) <- x.names
      out$tau.sq.alpha.samples <- mcmc(t(out$tau.sq.alpha.samples))
      colnames(out$tau.sq.alpha.samples) <- x.p.names
      if (is.null(sp.names)) {
        sp.names <- paste('sp', 1:N, sep = '')
      }
      coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- coef.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
      colnames(out$alpha.samples) <- coef.names.det
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      out$X <- X
      out$X.p <- X.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE
      out$psiRE <- FALSE

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
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.random[sites.k.fold[[i]]])
          if (binom) {
            y.indx <- !(1:J %in% curr.set)
	    y.fit <- y[rep(y.indx, N), drop = FALSE]
	    y.0 <- y[rep(y.indx, N), drop = FALSE]
          } else {
	    y.indx <- !((z.long.indx + 1) %in% curr.set)
	    y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	    y.fit <- y.fit[!is.na(y.fit)]
	    y.0 <- c(y.big[, curr.set, , drop = FALSE])
	    y.0 <- y.0[!is.na(y.0)]
          }
	  z.inits.fit <- z.inits[, -curr.set]
	  y.big.fit <- y.big[, -curr.set, , drop = FALSE]
	  y.big.0 <- y.big[, curr.set, , drop = FALSE]
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
          storage.mode(K.fit) <- "double"
	  storage.mode(n.obs.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.inits) <- "double"
          storage.mode(alpha.inits) <- "double"
          storage.mode(beta.comm.inits) <- "double"
          storage.mode(alpha.comm.inits) <- "double"
          storage.mode(tau.sq.beta.inits) <- "double"
          storage.mode(tau.sq.alpha.inits) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta.comm) <- "double"
          storage.mode(Sigma.beta.comm) <- "double"
          storage.mode(mu.alpha.comm) <- "double"
          storage.mode(Sigma.alpha.comm) <- "double"
          storage.mode(tau.sq.beta.a) <- "double"
          storage.mode(tau.sq.beta.b) <- "double"
          storage.mode(tau.sq.alpha.a) <- "double"
          storage.mode(tau.sq.alpha.b) <- "double"
          storage.mode(n.samples) <- "integer"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"

          # Run the model in C
          out.fit <- .Call("msPGOcc", y.fit, X.fit, X.p.fit, p.occ, p.det, J.fit, n.obs.fit, 
			   K.fit, N, beta.inits, alpha.inits, z.inits.fit,
          	           beta.comm.inits, 
          	           alpha.comm.inits, tau.sq.beta.inits, 
          	           tau.sq.alpha.inits, z.long.indx.fit, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, n.samples, n.omp.threads.fit, verbose.fit, n.report, 
          	           n.burn, n.thin, n.post.samples)
          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$X <- X
          out.fit$y <- y.big
          out.fit$n.post <- n.post.samples
          out.fit$psiRE <- FALSE
	  class(out.fit) <- "msPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.msPGOcc(out.fit, X.0)

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
	  if (binom) {
            like.samples <- array(NA, c(N, nrow(X.p.0), dim(y.big.0)[3]))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
              for (j in 1:nrow(X.p.0)) {
                for (k in 1:K.0[j]) {
                  like.samples[q, j, k] <- mean(dbinom(y.big.0[q, j, k], 1,
	          			         p.0.samples[j, q, k] * out.pred$z.0.samples[, q, z.0.long.indx[j]]))
	        }
              }
	    }
          } else {
	    like.samples <- matrix(NA, N, nrow(X.p.0))
	    for (q in 1:N) {
              p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
	      for (j in 1:nrow(X.p.0)) {
                like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
	          				p.0.samples[j, q, ] * 
	          			        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
              }
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a), na.rm = TRUE))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "msPGOcc"

    }
    out$run.time <- proc.time() - ptm
    out
}
