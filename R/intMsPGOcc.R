intMsPGOcc <- function(occ.formula, det.formula, data, inits, priors,  
		    n.samples, n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		    n.burn = round(.10 * n.samples), n.thin = 1, n.chains = 1,
		    k.fold, k.fold.threads = 1, k.fold.seed = 100, 
		    k.fold.only = FALSE, ...){

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
    # y -------------------------------
    if (!'y' %in% names(data)) {
      stop("error: detection-nondetection data y must be specified in data")
    }
    if (!is.list(data$y)) {
      stop("error: y must be a list of detection-nondetection data sets")
    }
    y <- data$y
    # Make sure all data sources are provided as 3-D arrays, even if only one
    # species in the data set or only one replicate. 
    n.data <- length(y)
    for (q in 1:n.data) {
      if (length(dim(y[[q]])) != 3) {
        stop(paste("Data source ", q, " is not provided as a 3-D array. Each data source must be formatted as a three-dimensional array with dimensions of species, site, and replicate survey. This is required even if the data source only has one species and/or one replicate.", sep = ''))
      }
    }
    # Species -------------------------
    if (!'species' %in% names(data)) {
      stop("error: species (list of species names/indexes for each data source) must be specified in data")
    }
    if (!is.list(data$species)) {
      stop("error: species must be a list of species indices for each data set")
    }
    sp.names <- unique(unlist(data$species))
    # Note the -1 for C
    sp.indx.c <- lapply(data$species, function(a) which(sp.names %in% a) - 1)
    sp.names.long <- data$species
    # Sites ---------------------------
    if (!'sites' %in% names(data)) {
      stop("error: site ids must be specified in data")
    }
    sites <- data$sites
    # Number of sites with at least one data source
    J <- length(unique(unlist(sites)))
    # Number of sites for each data set
    J.long <- sapply(y, function(a) dim(a)[[2]])
    # Occurrence covariates -----------
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
    # Detection covariates ------------
    if (!'det.covs' %in% names(data)) {
      data$det.covs <- list()
      for (i in 1:n.data) {
        if (verbose) {
          message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model for each data source.\n")
        }
        det.formula.curr <- det.formula[[i]]
        if (det.formula.curr == ~ 1) {
          for (i in 1:n.data) {
            data$det.covs[[i]] <- list(int = matrix(1, dim(y[[i]])[2], dim(y[[i]])[3]))
          }
        } else {
          stop("error: det.covs must be specified in data for a detection model with covariates")
        }
      }
    }
    # TODO: add in cross-validation eventually
    if (!missing(k.fold)) {
      message("k-fold cross-validation is not currently supported for integrated multi-species occupancy models. Please keep an eye out for future updates with this functionality.")
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
      if (nrow(data$det.covs[[i]]) == ncol(y[[i]])) {
        binom[i] <- TRUE
        # Check if there are missing site-level covariates 
        if (sum(is.na(data$det.covs[[i]])) != 0) {
          stop("error: missing values in site-level det.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
        }
        data$det.covs[[i]] <- data.frame(sapply(data$det.covs[[i]], rep,
						times = dim(y[[i]])[3]))
      }
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
    # TODO: 
    # Detection -----------------------
    # if (!is.null(findbars(det.formula))) {
    #   det.re.names <- sapply(findbars(det.formula), all.vars)
    #   for (i in 1:length(det.re.names)) {
    #     if (is(data$det.covs[, det.re.names[i]], 'factor')) {
    #       stop(paste("error: random effect variable ", det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
    #     } 
    #     if (is(data$det.covs[, det.re.names[i]], 'character')) {
    #       stop(paste("error: random effect variable ", det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
    #     }
    #   }
    # }

    # Checking missing values ---------------------------------------------
    # y -------------------------------
    for (q in 1:n.data) {
      y.na.test <- apply(y[[q]], c(1, 2), function(a) sum(!is.na(a)))
      if (sum(y.na.test == 0) > 0) {
        stop(paste("error: some sites in data source ", q, " have all missing detection histories. Remove these sites from all objects in the 'data' argument if the site is not surveyed by another data source, then use 'predict' to obtain predictions at these locations if desired.", sep = ''))
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
          if (sum(is.na(data$det.covs[[q]][, i])) > sum(is.na(y[[q]][1, , ]))) {
            stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
          }
        }
        # Misalignment between y and det.covs
        y.missing <- which(is.na(y[[q]][1, , ]))
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
      } else {
        if (sum(is.na(data$det.covs[[q]])) != 0) {
          stop("error: missing values in site-level det.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
        }

      }
    }

    # Formula -------------------------------------------------------------
    # Occupancy -----------------------
    if (missing(occ.formula)) {
      stop("error: occ.formula must be specified")
    }

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
    if (missing(det.formula)) {
      stop("error: det.formula must be specified")
    }
    if (!is.list(det.formula)) {
      stop(paste("error: det.formula must be a list of ", n.data, " formulas", sep = ''))
    }
    X.p <- list()
    X.p.re <- list()
    x.p.names <- list()
    # x.p.re.names <- list()
    # p.re.level.names <- list()
    for (i in 1:n.data) {
      if (is(det.formula[[i]], 'formula')) {
        tmp <- parseFormula(det.formula[[i]], data$det.covs[[i]])
        X.p[[i]] <- as.matrix(tmp[[1]])
        # X.p.re[[i]] <- as.matrix(tmp[[4]])
        # x.p.re.names[[i]] <- colnames(X.p.re)
        x.p.names[[i]] <- tmp[[2]]
	# p.re.level.names[[i]] <- lapply(data$det.covs[[i]][, x.p.re.names[[i]], drop = FALSE],
	# 		                function (a) sort(unique(a)))
      } else {
        stop(paste("error: det.formula for data source ", i, " is misspecified", sep = ''))
      }
    }
    x.p.names <- unlist(x.p.names)

    # Extract data from inputs --------------------------------------------
    # Total number of sites
    J.all <- nrow(X)
    if (length(X.p) != n.data | length(y) != n.data) {
      stop(paste("error: y and X.p must be lists of length ", n.data, ".", sep = ''))
    }
    # Number of species in each data set
    N.long <- lapply(y, nrow)
    # Total number of species
    N <- length(sp.names)
    # Number of occupancy parameters 
    p.occ <- ncol(X)
    # Number of occupancy random effect parameters
    p.occ.re <- ncol(X.re)
    # Number of detection parameters for each data set
    p.det.long <- sapply(X.p, function(a) dim(a)[[2]])
    # Total number of detection parameters
    p.det <- sum(p.det.long)
    # Number of latent occupancy random effect values
    n.occ.re <- length(unlist(apply(X.re, 2, unique)))
    n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
    n.rep <- lapply(y, function(a1) apply(matrix(a1[1, , ], dim(a1)[2], dim(a1)[3]), 
					  1, function(a2) sum(!is.na(a2))))
    # Max number of repeat visits for each data set
    K.long.max <- sapply(n.rep, max)
    # Number of repeat visits for each data set site. 
    K <- unlist(n.rep)
    if (missing(n.samples)) {
      stop("error: must specify number of MCMC samples")
    }
    if (n.burn > n.samples) {
      stop("error: n.burn must be less than n.samples")
    }
    if (n.thin > n.samples) {
      stop("error: n.thin must be less than n.samples")
    }
    # Get indices to map z to y --------------------------------------------
    y.big <- y
    names.long <- list()
    # Remove missing observations when the covariate data are available but
    # there are missing detection-nondetection data
    for (i in 1:n.data) {
      if (nrow(X.p[[i]]) == length(y[[i]][1, , ])) {
        X.p[[i]] <- X.p[[i]][!is.na(y[[i]][1, , ]), , drop = FALSE]
      }
      # Need these for later on
      names.long[[i]] <- which(!is.na(y[[i]][1, , ]))
    }
    n.obs.long <- sapply(X.p, nrow)
    n.obs <- sum(n.obs.long)
    z.long.indx.r <- list()
    for (i in 1:n.data) {
      z.long.indx.r[[i]] <- rep(sites[[i]], K.long.max[i])
      z.long.indx.r[[i]] <- z.long.indx.r[[i]][!is.na(c(y[[i]][1, , ]))]
      z.long.indx.r[[i]] <- rep(z.long.indx.r[[i]], each = N.long[[i]])
    }
    z.long.indx.r <- unlist(z.long.indx.r)
    # Subtract 1 for c indices
    z.long.indx.c <- z.long.indx.r - 1
    # Get species index
    sp.long.indx.c <- list()
    for (i in 1:n.data) {
      sp.long.indx.c[[i]] <- rep(sp.indx.c[[i]], n.obs.long[[i]])
    }
    sp.long.indx.c <- unlist(sp.long.indx.c)
    obs.long.indx.c <- list()
    start <- 1
    for (i in 1:n.data) {
      obs.long.indx.c[[i]] <- rep(start:(sum(!is.na(y[[i]][1, , ])) + start - 1), each = N.long[[i]])
      start <- start + sum(!is.na(y[[i]][1, , ]))
    }
    sp.site.indx.c <- matrix(0, N, J)
    for (j in 1:J) {
      tmp <- which(z.long.indx.r == j)
      sp.site.indx.c[unique(sp.long.indx.c[tmp] + 1), j] <- 1
    }
    obs.long.indx.c <- unlist(obs.long.indx.c) - 1
    y <- unlist(y)
    y <- y[!is.na(y)]
    # Index indicating the data set associated with each data point in y
    data.indx.r <- list()
    for (i in 1:n.data) {
      data.indx.r[[i]] <- rep(i, n.obs.long[[i]] * N.long[[i]])
    }
    data.indx.r <- unlist(data.indx.r)
    data.indx.c <- data.indx.r - 1
    sp.dat.long.indx.c <- rep(NA, length(data.indx.c))
    for (i in 1:length(sp.dat.long.indx.c)) {
      sp.dat.long.indx.c[i] <- which(sp.indx.c[[data.indx.r[i]]] == sp.long.indx.c[i])
    }
    sp.dat.long.indx.c <- sp.dat.long.indx.c - 1

    X.p.all <- matrix(NA, n.obs, max(p.det.long))
    indx <- 1
    for (i in 1:n.data) {
      X.p.all[indx:(indx + nrow(X.p[[i]]) - 1), 1:p.det.long[i]] <- X.p[[i]] 
      indx <- indx + nrow(X.p[[i]])
    }

    # Get random effect matrices all set ----------------------------------
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
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

    # alpha.comm ----------------------
    if ("alpha.comm.normal" %in% names(priors)) {
      if (!is.list(priors$alpha.comm.normal) | length(priors$alpha.comm.normal) != 2) {
        stop("error: alpha.comm.normal must be a list of length 2")
      }
      mu.alpha.comm <- priors$alpha.comm.normal[[1]]
      sigma.alpha.comm <- priors$alpha.comm.normal[[2]]
      if (length(mu.alpha.comm) != n.data | !is.list(mu.alpha.comm)) {
        stop(paste("error: alpha.comm.normal[[1]] must be a list of length ", 
		   n.data, " with elements corresponding to alphas.comms' mean for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(mu.alpha.comm[[q]]) != p.det.long[q] & length(mu.alpha.comm[[q]]) != 1) {
          if (p.det.long[q] == 1) {
            stop(paste("error: prior means for alpha.comm.normal[[1]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
            stop(paste("error: prior means for alpha.comm.normal[[1]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], "or 1", sep = ""))
          }
        }
        if (length(mu.alpha.comm[[q]]) != p.det.long[q]) {
          mu.alpha.comm[[q]] <- rep(mu.alpha.comm[[q]], p.det.long[q])
        }
      }
      mu.alpha.comm <- unlist(mu.alpha.comm)
      if (length(sigma.alpha.comm) != n.data | !is.list(sigma.alpha.comm)) {
        stop(paste("error: alpha.comm.normal[[2]] must be a list of length ", 
		   n.data, " with elements corresponding to alphas.comms' variance for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(sigma.alpha.comm[[q]]) != p.det.long[q] & length(sigma.alpha.comm[[q]]) != 1) {
          if (p.det.long[q] == 1) {
          stop(paste("error: prior variances for alpha.comm.normal[[2]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
          stop(paste("error: prior variances for alpha.comm.normal[[2]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], " or 1", sep = ""))

          }
        }
        if (length(sigma.alpha.comm[[q]]) != p.det.long[q]) {
          sigma.alpha.comm[[q]] <- rep(sigma.alpha.comm[[q]], p.det.long[q])
        }
      }
      sigma.alpha.comm <- unlist(sigma.alpha.comm)
    } else {
      if (verbose) {
        message("No prior specified for alpha.comm.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
      }
      mu.alpha.comm <- rep(0, p.det)
      sigma.alpha.comm <- rep(2.72, p.det) 
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
    
    # tau.sq.alpha --------------------   
    if ("tau.sq.alpha.ig" %in% names(priors)) {
      if (!is.list(priors$tau.sq.alpha.ig) | length(priors$tau.sq.alpha.ig) != 2) {
        stop("error: tau.sq.alpha.ig must be a list of length 2")
      }
      tau.sq.alpha.a <- priors$tau.sq.alpha.ig[[1]]
      tau.sq.alpha.b <- priors$tau.sq.alpha.ig[[2]]
      if (length(tau.sq.alpha.a) != n.data | !is.list(tau.sq.alpha.a)) {
        stop(paste("error: tau.sq.alpha.ig[[1]] must be a list of length ", 
		   n.data, " with elements corresponding to tau.sq.alpha's shape for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(tau.sq.alpha.a[[q]]) != p.det.long[q] & length(tau.sq.alpha.a[[q]]) != 1) {
          if (p.det.long[q] == 1) {
            stop(paste("error: prior scales for tau.sq.alpha.ig[[1]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
            stop(paste("error: prior scales for tau.sq.alpha.ig[[1]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], "or 1", sep = ""))
          }
        }
        if (length(tau.sq.alpha.a[[q]]) != p.det.long[q]) {
          tau.sq.alpha.a[[q]] <- rep(tau.sq.alpha.a[[q]], p.det.long[q])
        }
      }
      tau.sq.alpha.a <- unlist(tau.sq.alpha.a)
      if (length(tau.sq.alpha.b) != n.data | !is.list(tau.sq.alpha.b)) {
        stop(paste("error: tau.sq.alpha.ig[[2]] must be a list of length ", 
		   n.data, " with elements corresponding to tau.sq.alpha's scale for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(tau.sq.alpha.b[[q]]) != p.det.long[q] & length(tau.sq.alpha.b[[q]]) != 1) {
          if (p.det.long[q] == 1) {
          stop(paste("error: prior variances for tau.sq.alpha.ig[[2]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
          stop(paste("error: prior variances for tau.sq.alpha.ig[[2]][[", q, "]] must be a vector of length ", 
		     p.det.long[q], " or 1", sep = ""))

          }
        }
        if (length(tau.sq.alpha.b[[q]]) != p.det.long[q]) {
          tau.sq.alpha.b[[q]] <- rep(tau.sq.alpha.b[[q]], p.det.long[q])
        }
      }
      tau.sq.alpha.b <- unlist(tau.sq.alpha.b)
    } else {
      if (verbose) {
        message("No prior specified for tau.sq.alpha.ig.\nSetting prior shape and scale to 0.1\n")
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
    } else {
      sigma.sq.psi.a <- 0
      sigma.sq.psi.b <- 0
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
      z.test <- matrix(0, N, J)
      for (i in 1:n.data) {
        z.test[sp.indx.c[[i]] + 1, sites[[i]]] <- apply(y.big[[i]], c(1, 2), max, na.rm = TRUE)
      }
      init.test <- sum(z.inits < z.test)
      if (init.test > 0) {
        stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site.")
      }
    } else {
      z.inits <- matrix(0, N, J)
      for (i in 1:n.data) {
        z.inits[sp.indx.c[[i]] + 1, sites[[i]]] <- apply(y.big[[i]], c(1, 2), max, na.rm = TRUE)
      }
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
    # alpha.comm ----------------------
    if ("alpha" %in% names(inits)) {
      alpha.comm.inits <- inits[["alpha.comm"]]
      if (length(alpha.comm.inits) != n.data | !is.list(alpha.comm.inits)) {
        stop(paste("error: initial values for alpha.comm must be a list of length ", n.data,
		   sep = ""))
      }
      for (q in 1:n.data) {
        if (length(alpha.comm.inits[[q]]) != p.det.long[q] & length(alpha.comm.inits[[q]]) != 1) {
          if (p.det.long[q] == 1) {
            stop(paste("error: initial values for alpha.comm[[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
            stop(paste("error: initial values for alpha.comm[[", q, "]] must be a vector of length ", 
		       p.det.long[q], " or 1", sep = ""))
          }
        }
        if (length(alpha.comm.inits[[q]]) != p.det.long[q]) {
          alpha.comm.inits[[q]] <- rep(alpha.comm.inits[[q]], p.det.long[q])
        }
      }
      alpha.comm.inits <- unlist(alpha.comm.inits)
    } else {
      if (verbose) {
        message("alpha.comm is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      alpha.comm.inits <- rnorm(p.det, mu.alpha.comm, sqrt(sigma.alpha.comm))
    }

    alpha.comm.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
    alpha.comm.indx.c <- alpha.comm.indx.r - 1

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
    # tau.sq.alpha --------------------
    if ("tau.sq.alpha" %in% names(inits)) {
      tau.sq.alpha.inits <- inits[["tau.sq.alpha"]]
      if (length(tau.sq.alpha.inits) != n.data | !is.list(tau.sq.alpha.inits)) {
        stop(paste("error: initial values for tau.sq.alpha must be a list of length ", n.data,
		   sep = ""))
      }
      for (q in 1:n.data) {
        if (length(tau.sq.alpha.inits[[q]]) != p.det.long[q] & length(tau.sq.alpha.inits[[q]]) != 1) {
          if (p.det.long[q] == 1) {
            stop(paste("error: initial values for tau.sq.alpha[[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
	  } else {
            stop(paste("error: initial values for tau.sq.alpha[[", q, "]] must be a vector of length ", 
		       p.det.long[q], " or 1", sep = ""))
          }
        }
        if (length(tau.sq.alpha.inits[[q]]) != p.det.long[q]) {
          tau.sq.alpha.inits[[q]] <- rep(tau.sq.alpha.inits[[q]], p.det.long[q])
        }
      }
      tau.sq.alpha.inits <- unlist(tau.sq.alpha.inits)
    } else {
      if (verbose) {
        message("tau.sq.beta is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
      tau.sq.alpha.inits <- runif(p.det, 0.5, 10)
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
    # Index of species x detection parameter
    alpha.indx <- list()
    for (q in 1:n.data) {
      alpha.indx[[q]] <- matrix(0, N, p.det.long[q])
      alpha.indx[[q]][sp.indx.c[[q]] + 1, ] <- 1
    }
    alpha.indx <- unlist(alpha.indx)

    if ("alpha" %in% names(inits)) {
      alpha.inits <- inits[["alpha"]]
      if (length(alpha.inits) != n.data | !is.list(alpha.inits)) {
        stop(paste("error: initial values for alpha must be a list of length ", n.data,
		   sep = ""))
      }
      for (q in 1:n.data) {
        if (is.matrix(alpha.inits[[q]])) {
          if (ncol(alpha.inits[[q]]) != p.det.long[q] | nrow(alpha.inits[[q]]) != N.long[[q]]) {
            stop(paste("error: initial values for alpha must be a matrix with dimensions ", 
            	   N.long[[q]], "x", p.det.long[q], " or a single numeric value", sep = ""))
          }
        }
        if (!is.matrix(alpha.inits[[q]]) & length(alpha.inits[[q]]) != 1) {
          stop(paste("error: initial values for alpha must be a matrix with dimensions ", 
          	   N.long[[q]], " x ", p.det.long[q], " or a single numeric value", sep = ""))
        }
        if (length(alpha.inits[[q]]) == 1) {
          alpha.inits[[q]] <- matrix(alpha.inits[[q]], N.long[[q]], p.det.long[q])
        }
      }
      alpha.inits <- unlist(alpha.inits)
    } else {
      alpha.inits <- rnorm(N * p.det, alpha.comm.inits, sqrt(tau.sq.alpha.inits))
      alpha.inits <- alpha.inits[alpha.indx == 1]
      if (verbose) {
        message('alpha is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n')
      }
    }
    # Index to determine the species associated with each alpha value
    alpha.sp.indx.r <- rep(1:N, times = p.det)[alpha.indx == 1]
    alpha.sp.indx.c <- alpha.sp.indx.r - 1
    # Index to get the data set associated with each alpha value. 
    tmp <- list()
    for (i in 1:n.data) {
      tmp[[i]] <- rep(i, p.det.long[i] * N.long[[i]])
    }
    alpha.dat.indx.c <- unlist(tmp) - 1
    # Index to get the unique alpha.comm parameter associated with the given alpha
    alpha.p.det.indx.c <- list()
    for (i in 1:p.det) {
      alpha.p.det.indx.c[[i]] <- rep(i - 1, N.long[alpha.comm.indx.r[i]])
    }
    alpha.p.det.indx.c <- unlist(alpha.p.det.indx.c)

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

    curr.chain <- 1
    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
    storage.mode(K) <- "double"
    # Total number of observation x species data points. 
    n.obs.full <- length(y)
    consts <- c(N, J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, n.data, n.obs.full)
    storage.mode(consts) <- "integer"
    storage.mode(p.det.long) <- "integer"
    storage.mode(n.obs.long) <- "integer"
    storage.mode(n.occ.re.long) <- 'integer'
    storage.mode(J.long) <- "integer"
    storage.mode(N.long) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(beta.comm.inits) <- "double"
    storage.mode(alpha.comm.inits) <- "double"
    storage.mode(tau.sq.beta.inits) <- "double"
    storage.mode(tau.sq.alpha.inits) <- "double"
    storage.mode(alpha.comm.indx.c) <- "integer"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.sp.indx.c) <- "integer"
    storage.mode(alpha.dat.indx.c) <- "integer" 
    storage.mode(alpha.p.det.indx.c) <- 'integer'
    storage.mode(sp.dat.long.indx.c) <- "integer"
    storage.mode(sp.long.indx.c) <- 'integer'
    storage.mode(obs.long.indx.c) <- 'integer'
    storage.mode(sp.site.indx.c) <- 'integer'
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(mu.alpha.comm) <- "double"
    storage.mode(sigma.alpha.comm) <- "double"
    storage.mode(tau.sq.beta.a) <- "double"
    storage.mode(tau.sq.beta.b) <- "double"
    storage.mode(tau.sq.alpha.a) <- "double"
    storage.mode(tau.sq.alpha.b) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
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

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    out <- list()
    if (!k.fold.only) {
      for (i in 1:n.chains) {
        # Change initial values if i > 1
        if ((i > 1) & (!fix.inits)) {
          beta.comm.inits <- rnorm(p.occ, mu.beta.comm, sqrt(sigma.beta.comm))
          alpha.comm.inits <- rnorm(p.det, mu.alpha.comm, sqrt(sigma.alpha.comm))
          tau.sq.beta.inits <- runif(p.occ, 0.5, 10)
          tau.sq.alpha.inits <- runif(p.det, 0.5, 10)
          beta.inits <- matrix(rnorm(N * p.occ, beta.comm.inits, 
                		     sqrt(tau.sq.beta.inits)), N, p.occ)
          alpha.inits <- rnorm(N * p.det, alpha.comm.inits, sqrt(tau.sq.alpha.inits))
          alpha.inits <- alpha.inits[alpha.indx == 1]
          if (p.occ.re > 0) {
            sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
            beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
            beta.star.inits <- rep(beta.star.inits, N)
          }
        }

        storage.mode(chain.info) <- "integer"
        out.tmp[[i]] <- .Call("intMsPGOcc", y, X, X.p.all, X.re, consts, 
        	              n.occ.re.long, p.det.long, n.obs.long, N.long,
			      beta.inits, alpha.inits, z.inits, beta.comm.inits, 
          	              alpha.comm.inits, tau.sq.beta.inits, tau.sq.alpha.inits, 
          		      sigma.sq.psi.inits, 
        	              beta.star.inits, z.long.indx.c, 
			      alpha.comm.indx.c, data.indx.c, alpha.sp.indx.c, 
			      alpha.dat.indx.c, alpha.p.det.indx.c, 
			      sp.long.indx.c, obs.long.indx.c, sp.site.indx.c,
			      sp.dat.long.indx.c, beta.star.indx, beta.level.indx, 
			      mu.beta.comm, mu.alpha.comm, 
          		      Sigma.beta.comm, sigma.alpha.comm, 
          	              tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	              tau.sq.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
        	              n.samples, n.omp.threads, 
        	              verbose, n.report, samples.info, chain.info)
        chain.info[1] <- chain.info[1] + 1
      }
      # Calculate R-Hat ---------------
      out$rhat <- list()
      if (n.chains > 1) {
        # as.vector removes the "Upper CI" when there is only 1 variable. 
        out$rhat$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$beta.comm.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        out$rhat$alpha.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.comm.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        out$rhat$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$tau.sq.beta.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        out$rhat$tau.sq.alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$tau.sq.alpha.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					         mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
        out$rhat$alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.samples)))), 
        			      autoburnin = FALSE)$psrf[, 2])
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$sigma.sq.psi.samples)))), 
          			     autoburnin = FALSE)$psrf[, 2])
        }
      } else {
        out$rhat$beta.comm <- rep(NA, p.occ)
        out$rhat$alpha.comm <- rep(NA, p.det)
        out$rhat$tau.sq.beta <- rep(NA, p.occ)
        out$rhat$tau.sq.alpha <- rep(NA, p.det)
        out$rhat$beta <- rep(NA, p.occ * N)
        out$rhat$alpha <- rep(NA, sum(p.det.long * N.long))
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
        }
      }
      # Put everything into MCMC objects
      out$beta.comm.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.comm.samples))))
      colnames(out$beta.comm.samples) <- x.names
      out$alpha.comm.samples <- mcmc(do.call(rbind, 
        				lapply(out.tmp, function(a) t(a$alpha.comm.samples))))
      colnames(out$alpha.comm.samples) <- x.p.names
      out$tau.sq.beta.samples <- mcmc(do.call(rbind, 
        				lapply(out.tmp, function(a) t(a$tau.sq.beta.samples))))
      colnames(out$tau.sq.beta.samples) <- x.names
      out$tau.sq.alpha.samples <- mcmc(do.call(rbind, 
        				lapply(out.tmp, function(a) t(a$tau.sq.alpha.samples))))
      colnames(out$tau.sq.alpha.samples) <- x.p.names

      if (is.null(sp.names)) {
        sp.names <- paste('sp', 1:N, sep = '')
      }
      coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
      out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
      colnames(out$beta.samples) <- coef.names
      out$alpha.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.samples))))
      coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
      coef.names.det <- coef.names.det[alpha.indx == 1]
      colnames(out$alpha.samples) <- coef.names.det
      out$z.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$z.samples, 
        								dim = c(N, J, n.post.samples))))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$psi.samples, 
        								dim = c(N, J, n.post.samples))))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples, 
        								dim = c(N, J, n.post.samples))))
      out$like.samples <- aperm(out$like.samples, c(3, 1, 2))
      if (p.occ.re > 0) {
        out$sigma.sq.psi.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.psi.samples))))
        colnames(out$sigma.sq.psi.samples) <- x.re.names
        out$beta.star.samples <- mcmc(
          do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
        tmp.names <- unlist(re.level.names)
        beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
        beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
        colnames(out$beta.star.samples) <- beta.star.names
        out$re.level.names <- re.level.names
      }
      # Calculate effective sample sizes
      out$ESS <- list()
      out$ESS$beta.comm <- effectiveSize(out$beta.comm.samples)
      out$ESS$alpha.comm <- effectiveSize(out$alpha.comm.samples)
      out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
      out$ESS$tau.sq.alpha <- effectiveSize(out$tau.sq.alpha.samples)
      out$ESS$beta <- effectiveSize(out$beta.samples)
      out$ESS$alpha <- effectiveSize(out$alpha.samples)
      if (p.occ.re > 0) {
        out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
      }
      out$X <- X
      out$X.p <- X.p
      out$X.p.re <- X.p.re
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
      out$n.chains <- n.chains
      out$alpha.indx <- alpha.dat.indx.c + 1
      out$alpha.comm.indx <- alpha.comm.indx.r
      out$N <- N
      out$sites <- sites
      out$species <- data$species
      # TODO: 
      # if (p.det.re > 0) {
      #   out$pRE <- TRUE
      # } else {
        out$pRE <- FALSE
      # }
      if (p.occ.re > 0) {
        out$psiRE <- TRUE
      } else {
        out$psiRE <- FALSE
      }
    }
    class(out) <- "intMsPGOcc"
    out$run.time <- proc.time() - ptm
    out
  }
