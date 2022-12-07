PGOcc <- function(occ.formula, det.formula, data, inits, priors, 
		  n.samples, n.omp.threads = 1, verbose = TRUE,
		  n.report = 100, n.burn = round(.10 * n.samples), n.thin = 1, 
		  n.chains = 1, k.fold, k.fold.threads = 1, k.fold.seed = 100, 
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
    if (!is.matrix(data$occ.covs) & !is.data.frame(data$occ.covs)) {
      stop("error: occ.covs must be a matrix or data frame")
    }
    if (sum(is.na(data$occ.covs)) > 0) {
      stop("error: missing covariate values in data$occ.covs. Remove these sites from all data or impute non-missing values.")
    }
    if (!'det.covs' %in% names(data)) {
      if (det.formula == ~ 1) {
        if (verbose) {
          message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
	}
        data$det.covs <- list(int = rep(1, dim(y)[1]))
      } else {
        stop("error: det.covs must be specified in data for a detection model with covariates")
      }
    }
    if (!is.list(data$det.covs)) {
      stop("error: det.covs must be a list of matrices, data frames, and/or vectors")
    }
    if (missing(n.samples)) {
      stop("error: n.samples must be specified")
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

    # Check whether random effects are sent in as numeric, and
    # return error if they are. 
    # Occurrence ----------------------
    if (!is.null(findbars(occ.formula))) {
      occ.re.names <- sapply(findbars(occ.formula), all.vars)
      for (i in 1:length(occ.re.names)) {
        if (is(data$occ.covs[, occ.re.names[i]], 'factor')) {
          stop(paste("error: random effect variable ", occ.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
        } 
        if (is(data$occ.covs[, occ.re.names[i]], 'character')) {
          stop(paste("error: random effect variable ", occ.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
        }
      }
    }
    # Detection -----------------------
    if (!is.null(findbars(det.formula))) {
      det.re.names <- sapply(findbars(det.formula), all.vars)
      for (i in 1:length(det.re.names)) {
        if (is(data$det.covs[, det.re.names[i]], 'factor')) {
          stop(paste("error: random effect variable ", det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
        } 
        if (is(data$det.covs[, det.re.names[i]], 'character')) {
          stop(paste("error: random effect variable ", det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
        }
      }
    }

    # Checking missing values ---------------------------------------------
    # y -------------------------------
    y.na.test <- apply(y.big, 1, function(a) sum(!is.na(a)))
    if (sum(y.na.test == 0) > 0) {
      stop("error: some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
    }
    # occ.covs ------------------------
    if (sum(is.na(data$occ.covs)) != 0) {
      stop("error: missing values in occ.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
    }
    # det.covs ------------------------
    if (!binom) {
      for (i in 1:ncol(data$det.covs)) {
        if (sum(is.na(data$det.covs[, i])) > sum(is.na(y.big))) {
          stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
        }
      }
      # Misalignment between y and det.covs
      y.missing <- which(is.na(data$y))
      det.covs.missing <- lapply(data$det.covs, function(a) which(is.na(a)))
      for (i in 1:length(det.covs.missing)) {
        tmp.indx <- !(y.missing %in% det.covs.missing[[i]])
        if (sum(tmp.indx) > 0) {
          if (i == 1 & verbose) {
            message("There are missing values in data$y with corresponding non-missing values in data$det.covs.\nRemoving these site/replicate combinations for fitting the model.")
          }
          data$det.covs[y.missing, i] <- NA
        }
      }
    }
    # det.covs when binom == TRUE -----
    if (binom) {
      if (sum(is.na(data$det.covs)) != 0) {
        stop("error: missing values in site-level det.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
      }
    }

    # Formula -------------------------------------------------------------
    # Occupancy -----------------------
    if (is(occ.formula, 'formula')) {
      tmp <- parseFormula(occ.formula, data$occ.covs)
      X <- as.matrix(tmp[[1]])
      X.re <- as.matrix(tmp[[4]])
      x.re.names <- colnames(X.re)
      x.names <- tmp[[2]]
    } else {
      stop("error: occ.formula is misspecified")
    }
    # Get RE level names
    re.level.names <- lapply(data$occ.covs[, x.re.names, drop = FALSE],
			     function (a) sort(unique(a)))

    # Detection -----------------------
    if (is(det.formula, 'formula')) {
      tmp <- parseFormula(det.formula, data$det.covs)
      X.p <- as.matrix(tmp[[1]])
      X.p.re <- as.matrix(tmp[[4]])
      x.p.re.names <- colnames(X.p.re)
      x.p.names <- tmp[[2]]
    } else {
      stop("error: det.formula is misspecified")
    }
    p.re.level.names <- lapply(data$det.covs[, x.p.re.names, drop = FALSE],
			       function (a) sort(unique(a)))

    # Get basic info from inputs ------------------------------------------
    # Number of sites
    J <- nrow(y.big)
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
    if (p.det.re == 0) n.det.re.long <- 0
    # Number of replicates at each site
    n.rep <- apply(y.big, 1, function(a) sum(!is.na(a)))
    # Max number of repeat visits
    K.max <- max(n.rep)
    # Because I like K better than n.rep
    K <- n.rep

    # Get indices to map z to y -------------------------------------------
    if (!binom) {
      z.long.indx <- rep(1:J, dim(y.big)[2])
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
    # Number of data points for the y vector
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
    } else {
      sigma.sq.p.a <- 0
      sigma.sq.p.b <- 0
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
      z.test <- apply(y.big, 1, max, na.rm = TRUE)
      init.test <- sum(z.inits < z.test)
      if (init.test > 0) {
        stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site.")
      }
    } else {
      # In correct order since you reordered y.
      z.inits <- apply(y.big, 1, max, na.rm = TRUE)
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
    } else {
      sigma.sq.psi.inits <- 0
      beta.star.indx <- 0
      beta.star.inits <- 0
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
    } else {
      sigma.sq.p.inits <- 0
      alpha.star.indx <- 0
      alpha.star.inits <- 0
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

    curr.chain <- 1
    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, n.det.re)
    storage.mode(consts) <- "integer"
    storage.mode(K) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(Sigma.alpha) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"
    samples.info <- c(n.burn, n.thin, n.post.samples)
    storage.mode(samples.info) <- "integer"
    # For detection random effects
    storage.mode(X.p.re) <- "integer"
    alpha.level.indx <- sort(unique(c(X.p.re)))
    storage.mode(alpha.level.indx) <- "integer"
    storage.mode(n.det.re.long) <- "integer"
    storage.mode(sigma.sq.p.inits) <- "double"
    storage.mode(sigma.sq.p.a) <- "double"
    storage.mode(sigma.sq.p.b) <- "double"
    storage.mode(alpha.star.inits) <- "double"
    storage.mode(alpha.star.indx) <- "integer"
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
          beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
          alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
          if (p.occ.re > 0) {
            sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
            beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
          }
          if (p.det.re > 0) {
            sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
            alpha.star.inits <- rnorm(n.det.re, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
          }
        }
        storage.mode(chain.info) <- "integer"
        # Run the model in C
        out.tmp[[i]] <- .Call("PGOcc", y, X, X.p, X.re, X.p.re, consts, 
          		    K, n.occ.re.long, n.det.re.long, beta.inits, alpha.inits, 
          		    sigma.sq.psi.inits, sigma.sq.p.inits, beta.star.inits, 
          		    alpha.star.inits, z.inits, z.long.indx, beta.star.indx, 
          		    beta.level.indx, alpha.star.indx, alpha.level.indx, mu.beta, 
          		    mu.alpha, Sigma.beta, Sigma.alpha, sigma.sq.psi.a, sigma.sq.psi.b, 
          		    sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads, verbose, 
          		    n.report, samples.info, chain.info)
        chain.info[1] <- chain.info[1] + 1
      } # i   
      # Calculate R-Hat ---------------
      out <- list()
      out$rhat <- list()
      if (n.chains > 1) {
        # as.vector removes the "Upper CI" when there is only 1 variable. 
        out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        out$rhat$alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.samples)))), 
        			      autoburnin = FALSE)$psrf[, 2])
        if (p.det.re > 0) {
        out$rhat$sigma.sq.p <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$sigma.sq.p.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        }
        if (p.occ.re > 0) {
        out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$sigma.sq.psi.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        }
      } else {
        out$rhat$beta <- rep(NA, p.occ)
        out$rhat$alpha <- rep(NA, p.det)
        if (p.det.re > 0) {
          out$rhat$sigma.sq.p <- rep(NA, p.det.re)
        }
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
        }
      }
      # Put everything into MCMC objects
      out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(do.call(rbind, 
        				lapply(out.tmp, function(a) t(a$alpha.samples))))
      colnames(out$alpha.samples) <- x.p.names
      out$z.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$z.samples))))
      out$psi.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$psi.samples))))
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
      if (p.occ.re > 0) {
        out$sigma.sq.psi.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.psi.samples))))
        colnames(out$sigma.sq.psi.samples) <- x.re.names
        out$beta.star.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
        tmp.names <- unlist(re.level.names)
        beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
        colnames(out$beta.star.samples) <- beta.star.names
        out$re.level.names <- re.level.names
      }
      if (p.det.re > 0) {
        out$sigma.sq.p.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.p.samples))))
        colnames(out$sigma.sq.p.samples) <- x.p.re.names
        out$alpha.star.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.star.samples))))
        tmp.names <- unlist(p.re.level.names)
        alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
        colnames(out$alpha.star.samples) <- alpha.star.names
        out$p.re.level.names <- p.re.level.names
      }
      # Calculate effective sample sizes
      out$ESS <- list()
      out$ESS$beta <- effectiveSize(out$beta.samples)
      out$ESS$alpha <- effectiveSize(out$alpha.samples)
      if (p.det.re > 0) {
        out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
      }
      if (p.occ.re > 0) {
        out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
      }
      out$X <- X
      out$X.p <- X.p
      out$X.re <- X.re
      out$X.p.re <- X.p.re
      out$y <- y.big
      out$n.samples <- n.samples
      out$call <- cl
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$n.chains <- n.chains
      if (p.det.re > 0) {
        out$pRE <- TRUE
      } else {
        out$pRE <- FALSE
      }
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
        X.p.re.fit <- X.p.re[y.indx, , drop = FALSE]
        X.p.re.0 <- X.p.re[!y.indx, , drop = FALSE]
        n.det.re.fit <- length(unique(c(X.p.re.fit)))
        n.det.re.long.fit <- apply(X.p.re.fit, 2, function(a) length(unique(a)))
        if (p.det.re > 0) {	
          alpha.star.indx.fit <- rep(0:(p.det.re - 1), n.det.re.long.fit)
          alpha.level.indx.fit <- sort(unique(c(X.p.re.fit)))
          alpha.star.inits.fit <- rnorm(n.det.re.fit, 
          			      sqrt(sigma.sq.p.inits[alpha.star.indx.fit + 1]))
          p.re.level.names.fit <- list()
          for (t in 1:p.det.re) {
            tmp.indx <- alpha.level.indx.fit[alpha.star.indx.fit == t - 1]
            p.re.level.names.fit[[t]] <- unlist(p.re.level.names)[tmp.indx + 1]    
          }
        } else {
          alpha.star.indx.fit <- alpha.star.indx
          alpha.level.indx.fit <- alpha.level.indx
          alpha.star.inits.fit <- alpha.star.inits
        }
	# Random Occurrence Effects
        X.re.fit <- X.re[-curr.set, , drop = FALSE]
        X.re.0 <- X.re[curr.set, , drop = FALSE]
        n.occ.re.fit <- length(unique(c(X.re.fit)))
        n.occ.re.long.fit <- apply(X.re.fit, 2, function(a) length(unique(a)))
        if (p.occ.re > 0) {	
          beta.star.indx.fit <- rep(0:(p.occ.re - 1), n.occ.re.long.fit)
          beta.level.indx.fit <- sort(unique(c(X.re.fit)))
          beta.star.inits.fit <- rnorm(n.occ.re.fit, 
          			      sqrt(sigma.sq.psi.inits[beta.star.indx.fit + 1]))
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
        storage.mode(K.fit) <- "double"
        consts.fit <- c(J.fit, n.obs.fit, p.occ, p.occ.re, n.occ.re.fit, 
                        p.det, p.det.re, n.det.re.fit)
        storage.mode(consts.fit) <- "integer"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(n.samples) <- "integer"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(X.p.re.fit) <- "integer"
        storage.mode(n.det.re.long.fit) <- "integer"
        storage.mode(alpha.star.inits.fit) <- "double"
        storage.mode(alpha.star.indx.fit) <- "integer"
        storage.mode(alpha.level.indx.fit) <- "integer"
        storage.mode(X.re.fit) <- "integer"
        storage.mode(n.occ.re.long.fit) <- "integer"
        storage.mode(beta.star.inits.fit) <- "double"
        storage.mode(beta.star.indx.fit) <- "integer"
        storage.mode(beta.level.indx.fit) <- "integer"
        chain.info[1] <- 1
        storage.mode(chain.info) <- "integer"
        # Run the model in C
        out.fit <- .Call("PGOcc", y.fit, X.fit, X.p.fit, X.re.fit, X.p.re.fit, consts.fit, 
			 K.fit, n.occ.re.long.fit, n.det.re.long.fit, beta.inits, alpha.inits, 
			 sigma.sq.psi.inits, sigma.sq.p.inits, beta.star.inits.fit, 
			 alpha.star.inits.fit, z.inits.fit, z.long.indx.fit, beta.star.indx.fit, 
			 beta.level.indx.fit, alpha.star.indx.fit, alpha.level.indx.fit, mu.beta, 
			 mu.alpha, Sigma.beta, Sigma.alpha, sigma.sq.psi.a, sigma.sq.psi.b, 
			 sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads.fit, verbose.fit, 
			 n.report, samples.info, chain.info)
        out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
        colnames(out.fit$beta.samples) <- x.names
        out.fit$alpha.samples <- mcmc(t(out.fit$alpha.samples))
        colnames(out.fit$alpha.samples) <- x.p.names
        out.fit$X <- X.fit
        out.fit$y <- y.big.fit
        out.fit$X.p <- X.p.fit
        out.fit$call <- cl
        out.fit$n.samples <- n.samples
        out.fit$n.post <- n.post.samples
        out.fit$n.thin <- n.thin
        out.fit$n.burn <- n.burn
        out.fit$n.chains <- 1
        if (p.det.re > 0) {
          out.fit$pRE <- TRUE
        } else {
          out.fit$pRE <- FALSE
        }
        if (p.occ.re > 0) {
          out.fit$sigma.sq.psi.samples <- mcmc(t(out.fit$sigma.sq.psi.samples))
          colnames(out.fit$sigma.sq.psi.samples) <- x.re.names
          out.fit$beta.star.samples <- mcmc(t(out.fit$beta.star.samples))
          tmp.names <- unlist(re.level.names.fit)
          beta.star.names <- paste(rep(x.re.names, n.occ.re.long.fit), tmp.names, sep = '-')
          colnames(out.fit$beta.star.samples) <- beta.star.names
          out.fit$re.level.names <- re.level.names.fit
          out.fit$X.re <- X.re.fit
        }
        if (p.det.re > 0) {
          out.fit$sigma.sq.p.samples <- mcmc(t(out.fit$sigma.sq.p.samples))
          colnames(out.fit$sigma.sq.p.samples) <- x.p.re.names
          out.fit$alpha.star.samples <- mcmc(t(out.fit$alpha.star.samples))
          tmp.names <- unlist(p.re.level.names.fit)
          alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long.fit), tmp.names, sep = '-')
          colnames(out.fit$alpha.star.samples) <- alpha.star.names
          out.fit$p.re.level.names <- p.re.level.names.fit
          out.fit$X.p.re <- X.p.re.fit
        }
        if (p.occ.re > 0) {
          out.fit$psiRE <- TRUE
        } else {
          out.fit$psiRE <- FALSE	
        }
        class(out.fit) <- "PGOcc"

	# Get RE levels correct for when they aren't supplied at values starting at 1.
	if (p.occ.re > 0) {
	  tmp <- unlist(re.level.names)
	  X.re.0 <- matrix(tmp[c(X.re.0 + 1)], nrow(X.re.0), ncol(X.re.0))
	  colnames(X.re.0) <- x.re.names
	}

	# Predict occurrence at new sites
        if (p.occ.re > 0) {X.0 <- cbind(X.0, X.re.0)}
        out.pred <- predict.PGOcc(out.fit, X.0)

	# Detection 
        # Generate detection values
	# Get RE levels correct for when they aren't supplied at values starting at 1.
	if (p.det.re > 0) {
	  tmp <- unlist(p.re.level.names)
	  X.p.re.0 <- matrix(tmp[c(X.p.re.0 + 1)], nrow(X.p.re.0), ncol(X.p.re.0))
	  colnames(X.p.re.0) <- x.p.re.names
	}
        if (p.det.re > 0) {X.p.0 <- cbind(X.p.0, X.p.re.0)}
        out.p.pred <- predict.PGOcc(out.fit, X.p.0, type = 'detection')

	if (binom) {
          like.samples <- matrix(NA, nrow(y.big.0), ncol(y.big.0))
          for (j in 1:nrow(X.p.0)) {
            for (k in 1:K.0[j]) {
              like.samples[j, k] <- mean(dbinom(y.big.0[j, k], 1,
	      			         out.p.pred$p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
	    }
          }
        } else {
	  like.samples <- rep(NA, nrow(X.p.0))
	  for (j in 1:nrow(X.p.0)) {
            like.samples[j] <- mean(dbinom(y.0[j], 1, 
					   out.p.pred$p.0.samples[, j] * out.pred$z.0.samples[, z.0.long.indx[j]]))
          }
	}
	sum(log(like.samples), na.rm = TRUE)
      }
      model.deviance <- -2 * model.deviance
      # Return objects from cross-validation
      out$k.fold.deviance <- model.deviance
      stopImplicitCluster()
    }
    class(out) <- "PGOcc"
    out$run.time <- proc.time() - ptm
    out
  }
