intPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
		     n.samples, n.omp.threads = 1, verbose = TRUE, 
		     n.report = 1000, n.burn = round(.10 * n.samples), 
		     n.thin = 1, k.fold, k.fold.threads = 1, 
		     k.fold.seed = 100, k.fold.data, ...){

    ptm <- proc.time()

    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

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
    if (!missing(k.fold)) {
      if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
        stop("error: k.fold must be a single integer value >= 2")  
      }
    }

    # Make all covariates a data frame. Unlist is necessary for when factors
    # are supplied. 
    for (i in 1:n.data) {
      # Get in required R model format. 
      data$det.covs[[i]] <- data.frame(lapply(data$det.covs[[i]], function(a) unlist(c(a))))
      # Replicate det.covs if only covariates are at the site level. 
      if (nrow(data$det.covs[[i]]) == nrow(y[[i]])) {
        data$det.covs[[i]] <- data.frame(sapply(data$det.covs[[i]], rep, times = dim(y[[i]])[2]))
      }
    }
    data$occ.covs <- as.data.frame(data$occ.covs)

    # Checking missing values ---------------------------------------------
    for (q in 1:n.data) {
      y.na.test <- apply(y[[q]], 1, function(a) sum(!is.na(a)))
      if (sum(y.na.test == 0) > 0) {
        stop(paste("error: some sites in data source ", q, " in y have all missing detection histories.\n Remove these sites from y and all objects in the 'data' argument if the site is not surveyed by another data source\n, then use 'predict' to obtain predictions at these locations if desired.", sep = ''))
      }
    }

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

    # Formula -------------------------------------------------------------
    # Occupancy -----------------------
    if (missing(occ.formula)) {
      stop("error: occ.formula must be specified")
    }

    if (class(occ.formula) == 'formula') {
      tmp <- parseFormula(occ.formula, data$occ.covs)
      X <- as.matrix(tmp[[1]])
      x.names <- tmp[[2]]
    } else {
      stop("error: occ.formula is misspecified")
    }

    # Detection -----------------------
    if (missing(det.formula)) {
      stop("error: det.formula must be specified")
    }
    if (!is.list(det.formula)) {
      stop(paste("error: det.formula must be a list of ", n.data, " formulas", sep = ''))
    }
    X.p <- list()
    x.p.names <- list()
    for (i in 1:n.data) {
      if (class(det.formula[[i]]) == 'formula') {
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
    K.long.max <- sapply(n.rep, max)
    # Number of repeat visits for each data set site. 
    K <- unlist(n.rep)

    # Get indics to map z to y --------------------------------------------
    X.p.orig <- X.p
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

    # Set model.deviance to NA for returning when no cross-validation
    model.deviance <- NA

    # Specify storage modes -----------------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
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
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(sigma.alpha) <- "double"
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

    # Run the model in C    
    out <- .Call("intPGOcc", y, X, X.p.all, p.occ, p.det, p.det.long, 
		 J, J.long, K, n.obs, n.obs.long, n.data, 
		 beta.inits, alpha.inits, z.inits,
		 z.long.indx.c, data.indx.c, alpha.indx.c, mu.beta, mu.alpha, 
		 Sigma.beta, sigma.alpha, n.samples, 
		 n.omp.threads, verbose, n.report, n.burn, n.thin, 
		 n.post.samples)

    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- x.names
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    colnames(out$alpha.samples) <- x.p.names
    out$z.samples <- mcmc(t(out$z.samples))
    out$psi.samples <- mcmc(t(out$psi.samples))
    # y.rep.samples is returned as a list, where each element 
    # corresponds to a different data set. 
    tmp <- list()
    indx <- 1
    for (q in 1:n.data) {
      tmp[[q]] <- array(NA, dim = c(J.long[q] * K.long.max[q], n.post.samples))
      tmp[[q]][names.long[[q]], ] <- out$y.rep.samples[indx:(indx + n.obs.long[q] - 1), ] 
      tmp[[q]] <- array(tmp[[q]], dim = c(J.long[q], K.long.max[q], n.post.samples))
      tmp[[q]] <- aperm(tmp[[q]], c(3, 1, 2))
      indx <- indx + n.obs.long[q]
    }
    out$y.rep.samples <- tmp
    out$X <- X
    out$X.p <- X.p.orig
    out$y <- y.big
    out$n.samples <- n.samples
    out$call <- cl
    out$sites <- sites
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn

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
	z.long.indx.fit <- list()
	tmp.indx <- 1
	for (q in 1:n.data) {
          z.long.indx.fit[[q]] <- matrix(NA, length(sites.fit[[q]]), K.long.max[q])
          for (j in 1:length(sites.fit[[q]])) {
            z.long.indx.fit[[q]][j, 1:K.fit[tmp.indx]] <- sites.fit[[q]][j]
	    tmp.indx <- tmp.indx + 1
          }
          z.long.indx.fit[[q]] <- c(z.long.indx.fit[[q]])
          z.long.indx.fit[[q]] <- z.long.indx.fit[[q]][!is.na(z.long.indx.fit[[q]])] - 1
        }
	z.long.indx.fit <- unlist(z.long.indx.fit)
	
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
	z.long.indx.0 <- list()
	tmp.indx <- 1
	for (q in 1:n.data) {
          if (!is.na(sites.0[[q]][1])) {
            z.long.indx.0[[q]] <- matrix(NA, length(sites.0[[q]]), K.long.max[q])
            for (j in 1:length(sites.0[[q]])) {
              z.long.indx.0[[q]][j, 1:K.0[tmp.indx]] <- sites.0[[q]][j]
	      tmp.indx <- tmp.indx + 1
            }
            z.long.indx.0[[q]] <- c(z.long.indx.0[[q]])
            z.long.indx.0[[q]] <- z.long.indx.0[[q]][!is.na(z.long.indx.0[[q]])] - 1
	  }
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

        storage.mode(y.fit) <- "double"
        storage.mode(z.inits.fit) <- "double"
        storage.mode(X.p.fit) <- "double"
        storage.mode(X.fit) <- "double"
        storage.mode(p.det) <- "integer"
        storage.mode(p.det.long) <- "integer"
        storage.mode(p.occ) <- "integer"
        storage.mode(n.obs.fit) <- "integer"
        storage.mode(n.obs.long.fit) <- "integer"
        storage.mode(J.fit) <- "integer"
        storage.mode(J.long.fit) <- "integer"
        storage.mode(K.fit) <- "integer"
	storage.mode(n.obs.fit) <- "integer"
	storage.mode(n.obs.long.fit) <- "integer"
        storage.mode(n.data) <- "integer"
        storage.mode(beta.inits) <- "double"
        storage.mode(alpha.inits) <- "double"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(data.indx.c.fit) <- "integer"
        storage.mode(alpha.indx.c) <- "integer"
        storage.mode(mu.beta) <- "double"
        storage.mode(Sigma.beta) <- "double"
        storage.mode(mu.alpha) <- "double"
        storage.mode(sigma.alpha) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(n.burn) <- "integer"
        storage.mode(n.thin) <- "integer"
    
	out.fit <- .Call("intPGOcc", y.fit, X.fit, X.p.fit, p.occ, p.det, p.det.long, 
		         J.fit, J.long.fit, K.fit, n.obs.fit, n.obs.long.fit, n.data, 
		         beta.inits, alpha.inits, z.inits.fit,
		         z.long.indx.fit, data.indx.c.fit, alpha.indx.c, mu.beta, mu.alpha, 
		         Sigma.beta, sigma.alpha, n.samples, 
		         n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
		         n.post.samples)

        # Predict occurrence at new sites. 
        psi.0.samples <- logit.inv(t(X.0 %*% out.fit$beta.samples))
        # Check this if you encounter a bug. 
        z.0.samples <- matrix(NA, n.post.samples, nrow(X.0))
        for (j in 1:nrow(X.0)) {
          z.0.samples[, j] <- rbinom(n.post.samples, 1, psi.0.samples[, j])
        }
        # Detection 

	p.0.samples <- matrix(NA, n.post.samples, nrow(X.p.0))
        like.samples <- rep(NA, nrow(X.p.0))
        for (j in 1:nrow(X.p.0)) {
          p.0.samples[, j] <- logit.inv(X.p.0[j, 1:sum(alpha.indx.r == data.indx.0[j])] %*% out.fit$alpha.samples[which(alpha.indx.r == data.indx.0[j]), ])
          like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * z.0.samples[, z.long.indx.0[j] + 1]))
        }
	as.vector(tapply(like.samples, data.indx.0, function(a) sum(log(a))))
      }
      model.deviance <- -2 * model.deviance
      # Return objects from cross-validation
      out$k.fold.deviance <- model.deviance
      stopImplicitCluster()
    }

    class(out) <- "intPGOcc"

    out$run.time <- proc.time() - ptm
    out
}
