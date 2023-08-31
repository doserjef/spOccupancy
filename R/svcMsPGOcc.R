svcMsPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
		      tuning, svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
		      n.neighbors = 15, search.type = "cb", std.by.sp = FALSE, n.factors, 
		      n.batch, batch.length, accept.rate = 0.43,
		      n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		      n.burn = round(.10 * n.batch * batch.length), 
		      n.thin = 1, n.chains = 1, ...){

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
    stop("error: svcMsPGOcc is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 
  }
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
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  coords <- as.matrix(data$coords)
  # Check if all spatial coordinates are unique. 
  unique.coords <- unique(data$coords)
  if (nrow(unique.coords) < nrow(data$coords)) {
    stop("coordinates provided in coords are not all unique. spOccupancy requires each site to have its own unique pair of spatial coordinates. This may be the result of an error in preparing the data list, or you will need to change what you consider a 'site' in order to meet this requirement.") 
  }
  if (missing(n.factors)) {
    stop("error: n.factors must be specified for a spatial factor occupancy model")
  }

  if (std.by.sp != FALSE & std.by.sp != TRUE) {
    stop("error: std.by.sp must be either TRUE or FALSE")
  }
  if ('range.ind' %in% names(data)) {
    range.ind <- data$range.ind
    if (!is.matrix(range.ind)) {
      stop(paste("error: if specified, range.ind must be a matrix of 0s and 1s with ", 
                 nrow(y), " rows and ", ncol(y), " columns.", sep = ''))
    }
    if (nrow(range.ind) != nrow(y) | ncol(range.ind) != ncol(y)) {
      stop(paste("error: if specified, range.ind must be a matrix of 0s and 1s with ", 
                 nrow(y), " rows and ", ncol(y), " columns.", sep = ''))
    }
  } else {
    range.ind <- matrix(1, nrow(y), ncol(y))
  }

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[, ord, , drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    range.ind <- range.ind[, ord, drop = FALSE]
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
  # y -------------------------------
  y.na.test <- apply(y.big, c(1, 2), function(a) sum(!is.na(a)))
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
      # Note that this assumes the same detection history for each species.  
      if (sum(is.na(data$det.covs[, i])) > sum(is.na(y.big[1, , ]))) {
        stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
      }
    }
    # Misalignment between y and det.covs
    y.missing <- which(is.na(y[1, , ]))
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

  # Get a separate X for each species for standardization within species
  # if desired.
  X.big <- array(NA, dim = c(nrow(X), ncol(X), dim(y)[1]))
  species.sds <- matrix(NA, nrow = dim(y)[1], ncol = ncol(X))
  species.means <- matrix(NA, nrow = dim(y)[1], ncol = ncol(X))
  for (i in 1:nrow(y)) {
    curr.indx <- which(range.ind[i, ] == 1)
    X.big[curr.indx, , i] <- X[curr.indx, , drop = FALSE]
    if (std.by.sp) {
      for (r in 1:ncol(X)) {
        if (sd(X[, r], na.rm = TRUE) != 0) {
          species.sds[i, r] <- sd(c(X.big[curr.indx, r, i]), na.rm = TRUE)
          species.means[i, r] <- mean(c(X.big[curr.indx, r, i]), na.rm = TRUE)
          X.big[curr.indx, r, i] <- c((X.big[curr.indx, r, i] - species.means[i, r])) / species.sds[i, r]
	}
      }
    }
  }
  X.big <- ifelse(is.na(X.big), 0, X.big)

  # Detection -----------------------
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }

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

  # Extract data from inputs --------------------------------------------
  # Number of species 
  N <- dim(y)[1]
  # Number of latent factors
  q <- n.factors
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
  # Number of sites
  J <- nrow(X)
  # Number of repeat visits
  n.rep <- apply(y.big[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- which(!is.na(y.big[1, j, ]))
  }
  K.max <- dim(y.big)[3]
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

  # Check SVC columns -----------------------------------------------------
  if (is.character(svc.cols)) {
    # Check if all column names in svc are in occ.covs
    if (!all(svc.cols %in% x.names)) {
        missing.cols <- svc.cols[!(svc.cols %in% x.names)]
        stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
    }
    # Convert desired column names into the numeric column index
    svc.cols <- (1:p.occ)[x.names %in% svc.cols]
    
  } else if (is.numeric(svc.cols)) {
    # Check if all column indices are in 1:p.occ
    if (!all(svc.cols %in% 1:p.occ)) {
        missing.cols <- svc.cols[!(svc.cols %in% (1:p.occ))]
        stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
    }
  }
  p.svc <- length(svc.cols)
  q.p.svc <- q * p.svc

  # Get indices to map z to y -------------------------------------------
  if (!binom) {
    z.long.indx <- rep(1:J, dim(y.big)[3])
    z.long.indx <- z.long.indx[!is.na(c(y.big[1, , ]))]
    # Subtract 1 for indices in C
    z.long.indx <- z.long.indx - 1
  } else {
    z.long.indx <- 0:(J - 1)
  }
  # y is order as follows: sorted by visit, site within visit, species within site
  # This will match up with z, which is sorted by site, then species within site. 
  # This matches up with spMsPGOcc.
  y <- c(y)
  # Assumes the missing data are constant across species, which seems likely, 
  # but may eventually need some updating. 
  # Removing missing observations when covariate data are available but 
  # there are missing detection-nondetection data. 
  names.long <- which(!is.na(c(y.big[1, , ])))
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
    if (length(phi.a) != q.p.svc & length(phi.a) != 1) {
      stop(paste("error: phi.unif[[1]] must be a vector of length ", 
      	   q.p.svc, ", a matrix with ", q, " rows and ", p.svc, 
	   " columns, or a vector of length 1 with elements corresponding to phis' lower bound for each latent factor and spatially-varying coefficient", sep = ""))
    }
    if (length(phi.b) != q.p.svc & length(phi.b) != 1) {
      stop(paste("error: phi.unif[[2]] must be a vector of length ", 
      	   q.p.svc, ", a matrix with ", q, " rows and ", p.svc, 
	   " columns, or a vector of length 1 with elements corresponding to phis' upper bound for each latent factor and spatially-varying coefficient", sep = ""))
    }
    if (length(phi.a) != q.p.svc) {
      phi.a <- rep(phi.a, q.p.svc)
    }
    if (length(phi.b) != q.p.svc) {
      phi.b <- rep(phi.b, q.p.svc)
    }
  } else {
    if (verbose) {
    message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
    }
    if (NNGP) {
      coords.D <- iDist(coords)
    }
    phi.a <- rep(3 / max(coords.D), q.p.svc)
    phi.b <- rep(3 / sort(unique(c(coords.D)))[2], q.p.svc)
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
    if (length(nu.a) != q.p.svc & length(nu.a) != 1) {
      stop(paste("error: nu.unif[[1]] must be a vector of length ", 
      	   q.p.svc, ", a matrix with ", q, " rows and ", p.svc, 
	   " columns, or a vector of length 1 with elements corresponding to nus' lower bound for each latent factor and spatially-varying coefficient", sep = ""))
    }
    if (length(nu.b) != q & length(nu.b) != 1) {
      stop(paste("error: nu.unif[[2]] must be a vector of length ", 
      	   q.p.svc, ", a matrix with ", q, " rows and ", p.svc, 
	   " columns, or a vector of length 1 with elements corresponding to nus' upper bound for each latent factor and spatially-varying coefficient", sep = ""))
    }
    if (length(nu.a) != q.p.svc) {
      nu.a <- rep(nu.a, q.p.svc)
    }
    if (length(nu.b) != q.p.svc) {
      nu.b <- rep(nu.b, q.p.svc)
    }
  } else {
    nu.a <- rep(0, q.p.svc)
    nu.b <- rep(0, q.p.svc)
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

  # Initial values --------------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # z -------------------------------
  # ORDER: an N x J matrix sent in as a column-major vector ordered by site, then species
  #        within site. 
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
    # Reorder the user supplied inits values for NNGP models
    if (NNGP) {
      z.inits <- z.inits[, ord]
    }
    z.test <- apply(y.big, c(1, 2), max, na.rm = TRUE)
    init.test <- sum(z.inits < z.test)
    if (init.test > 0) {
      stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site.")
    }
  } else {
    # In correct order since you reordered y for NNGP. 
    z.inits <- apply(y.big, c(1, 2), max, na.rm = TRUE)
    if (verbose) {
      message("z is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
  # z.inits is ordered by site then species within site
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
  # alpha.comm -----------------------
  # ORDER: a p.det vector ordered by the effects in the detection formula.
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
  # tau.sq.alpha -----------------------
  # ORDER: a p.det vector ordered by the effects in the detection formula.
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
      message('tau.sq.alpha is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n')
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
  # alpha ----------------------------
  # ORDER: N x p.det matrix sent in as a column-major vector ordered by 
  #        parameter then species within parameter. 
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
  # Create a N * p.det x 1 matrix of the species-level regression coefficients. 
  # This is ordered by parameter, then species within parameter. 
  alpha.inits <- c(alpha.inits)
  # phi -----------------------------
  # ORDER: a q x p.svc matrix sent in as a column-major vector sorted first by 
  #        the spatially-varying coefficient, then latent factor within svc.
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != q.p.svc & length(phi.inits) != 1) {
      stop(paste("error: initial values for phi must be of length ", q.p.svc, " or 1", 
      	   sep = ""))
    }
    if (length(phi.inits) != q.p.svc) {
      phi.inits <- rep(phi.inits, q.p.svc)
    }
  } else {
    phi.inits <- runif(q.p.svc, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random values from the prior distribution\n")
    }
  }
  # lambda ----------------------------
  # ORDER: an p.svc N x q matrices sent in as a list, each in 
  #        column-major vector, which is ordered by factor, 
  #        then species within factor. Eventually sent into 
  #        C++ as a stacked p.svc x N x q matrix, ordered by 
  #        svc, then factor within svc, then species within factor.
  if ("lambda" %in% names(inits)) {
    lambda.inits <- inits[["lambda"]]
    if (!is.list(lambda.inits)) {
      stop(paste("error: initial values for lambda must be a list comprised of ", 
		 p.svc, " matrices, each with dimensions ", N, " x ", q, sep = ""))
    }
    for (i in 1:p.svc) {
      if (nrow(lambda.inits[[i]]) != N | ncol(lambda.inits[[i]]) != q) {
        stop(paste("error: initial values for lambda[[", i, 
		   "]] must be a matrix with dimensions ", N, " x ", q, sep = ""))
      }
      if (!all.equal(diag(lambda.inits[[i]]), rep(1, q))) {
        stop("error: diagonal of inits$lambda[[", i, "]] matrix must be all 1s")
      }
      if (sum(lambda.inits[[i]][upper.tri(lambda.inits[[i]])]) != 0) {
        stop("error: upper triangle of inits$lambda[[", i, "]] must be all 0s")
      }
    }
    lambda.inits <- unlist(lambda.inits)
  } else {
    lambda.inits <- list()
    for (i in 1:p.svc) {
      lambda.inits[[i]] <- matrix(0, N, q)
      diag(lambda.inits[[i]]) <- 1
      lambda.inits[[i]][lower.tri(lambda.inits[[i]])] <- rnorm(sum(lower.tri(lambda.inits[[i]])))
    }
    if (verbose) {
      message("lambda is not specified in initial values.\nSetting initial values of the lower triangle to random values from a standard normal\n")
    }
    lambda.inits <- unlist(lambda.inits)
  }
  # nu ------------------------
  if ("nu" %in% names(inits)) {
    nu.inits <- inits[["nu"]]
    if (length(nu.inits) != q.p.svc & length(nu.inits) != 1) {
      stop(paste("error: initial values for nu must be of length ", q.p.svc,  " or 1",
      	   sep = ""))
    }
    if (length(nu.inits) != q.p.svc) {
      nu.inits <- rep(nu.inits, q.p.svc)
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      nu.inits <- runif(q.p.svc, nu.a, nu.b)
    } else {
      nu.inits <- rep(0, q.p.svc)
    }
  }
  # sigma.sq.psi ------------------
  # ORDER: a length p.occ.re vector ordered by the random effects in the formula.
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
    beta.star.inits <- rep(beta.star.inits, N)
  } else {
    sigma.sq.psi.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
  }

  # sigma.sq.p ------------------
  # ORDER: a length p.det.re vector ordered by the random effects in the formula.
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
  } else {
    sigma.sq.p.inits <- 0
    alpha.star.indx <- 0
    alpha.star.inits <- 0
  }
  # w -----------------------------
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (!is.list(w.inits)) {
      stop(paste("error: initial values for w must be a list comprised of ", 
		 p.svc, " matrices, each with dimensions ", q, " x ", J, sep = ""))
    }
    for (i in 1:p.svc) {
      if (!is.matrix(w.inits[[i]])) {
        stop(paste("error: initial values for w must be a matrix with dimensions ",
        	   q, " x ", J, sep = ""))
      }
      if (nrow(w.inits[[i]]) != q | ncol(w.inits[[i]]) != J) {
        stop(paste("error: initial values for w must be a matrix with dimensions ",
        	   q, " x ", J, sep = ""))
      }
      if (NNGP) {
        w.inits[[i]] <- w.inits[[i]][, ord]
      }
    }
    w.inits <- unlist(w.inits)
  } else {
    w.inits <- list()
    for (i in 1:p.svc) {
      w.inits[[i]] <- matrix(0, q, J)
    }
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
    }
    w.inits <- unlist(w.inits)
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

  # Prep for SVCs ---------------------------------------------------------
  X.w <- X[, svc.cols, drop = FALSE]
  x.w.names <- colnames(X.w)
  X.w.big <- X.big[, svc.cols, , drop = FALSE]

  # Get tuning values ---------------------------------------------------
  # Not accessed, but necessary to keep things in line with the underlying functions. 
  sigma.sq.tuning <- rep(0, q.p.svc)
  phi.tuning <- rep(0, q.p.svc)
  nu.tuning <- rep(0, q.p.svc)
  if (missing(tuning)) {
    phi.tuning <- rep(1, q.p.svc)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, q.p.svc)
    }
  } else {
    names(tuning) <- tolower(names(tuning))
    # phi ---------------------------
    if(!"phi" %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (length(phi.tuning) == 1) {
      phi.tuning <- rep(tuning$phi, q.p.svc)
    } else if (length(phi.tuning) != q.p.svc) {
      stop(paste("error: phi tuning must be either a single value or a vector of length ",
      	   q.p.svc, sep = ""))
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) == 1) {
        nu.tuning <- rep(tuning$nu, q.p.svc)
      } else if (length(nu.tuning) != q.p.svc) {
        stop(paste("error: nu tuning must be either a single value or a vector of length ",
        	   q.p.svc, sep = ""))
      }
    }
  }
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA
  curr.chain <- 1

  # Names for spatial parameters
  if (cov.model != 'matern') {
    theta.names <- paste(rep(c('phi'), each = q), 1:q, sep = '-')
    theta.names <- paste(rep(theta.names, times = p.svc), 
      		   rep(x.w.names, each = q), sep = '-')
  } else {
    theta.names <- paste(rep(c('phi', 'nu'), each = q), 1:q, sep = '-')
    theta.names <- paste(rep(theta.names, times = p.svc), 
      		   rep(x.w.names, each = 2 * q), sep = '-')
  } 

  if (!NNGP) {

    stop("error: svcMsPGOcc is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 

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
    
    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X.big) <- "double"
    storage.mode(X.w.big) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(K) <- "double"
    storage.mode(range.ind) <- "double"
    # consts order: N, J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, n.det.re,
    #               q, p.svc
    consts <- c(N, J, n.obs, p.occ, p.occ.re, n.occ.re, 
		p.det, p.det.re, n.det.re, q, p.svc)
    storage.mode(consts) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(beta.comm.inits) <- "double"
    storage.mode(alpha.comm.inits) <- "double"
    storage.mode(tau.sq.beta.inits) <- "double"
    storage.mode(tau.sq.alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(lambda.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(mu.alpha.comm) <- "double"
    storage.mode(Sigma.alpha.comm) <- "double"
    storage.mode(tau.sq.beta.a) <- "double"
    storage.mode(tau.sq.beta.b) <- "double"
    storage.mode(tau.sq.alpha.a) <- "double"
    storage.mode(tau.sq.alpha.b) <- "double"
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
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        beta.comm.inits <- rnorm(p.occ, mu.beta.comm, sqrt(sigma.beta.comm))
        alpha.comm.inits <- rnorm(p.det, mu.alpha.comm, sqrt(sigma.alpha.comm))
        tau.sq.beta.inits <- runif(p.occ, 0.5, 10)
        tau.sq.alpha.inits <- runif(p.det, 0.5, 10)
        beta.inits <- matrix(rnorm(N * p.occ, beta.comm.inits, 
              		     sqrt(tau.sq.beta.inits)), N, p.occ)
        beta.inits <- c(beta.inits)
        alpha.inits <- matrix(rnorm(N * p.det, alpha.comm.inits, 
              		      sqrt(tau.sq.alpha.inits)), N, p.det)
        alpha.inits <- c(alpha.inits)
        lambda.inits <- list()
        for (j in 1:p.svc) {
          lambda.inits[[j]] <- matrix(0, N, q)
          diag(lambda.inits[[j]]) <- 1
          lambda.inits[[j]][lower.tri(lambda.inits[[j]])] <- rnorm(sum(lower.tri(lambda.inits[[j]])))
        }
        lambda.inits <- unlist(lambda.inits)
        phi.inits <- runif(q.p.svc, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(q.p.svc, nu.a, nu.b)
        }
        if (p.det.re > 0) {
          sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
          alpha.star.inits <- rnorm(n.det.re, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
          alpha.star.inits <- rep(alpha.star.inits, N)
        }
        if (p.occ.re > 0) {
          sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
          beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
          beta.star.inits <- rep(beta.star.inits, N)
        }
      }

      storage.mode(chain.info) <- "integer"
      # Run the model in C
      # Getting real close to 65 arguments....
      out.tmp[[i]] <- .Call("svcMsPGOccNNGP", y, X.big, X.w.big, X.p, coords, X.re, X.p.re, 
      		      range.ind, consts, K, n.occ.re.long, n.det.re.long, 
        	            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
        	            beta.inits, alpha.inits, z.inits,
        	            beta.comm.inits, 
        	            alpha.comm.inits, tau.sq.beta.inits, 
        	            tau.sq.alpha.inits, phi.inits, 
        	            lambda.inits, w.inits, nu.inits, 
      		    sigma.sq.psi.inits, sigma.sq.p.inits, 
        		    beta.star.inits, alpha.star.inits, z.long.indx, 
        		    beta.star.indx, beta.level.indx, alpha.star.indx, 
        		    alpha.level.indx, mu.beta.comm, 
        	            mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
        	            tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
        	            tau.sq.alpha.b, phi.a, phi.b,
        	            nu.a, nu.b, sigma.sq.psi.a, sigma.sq.psi.b, 
        		    sigma.sq.p.a, sigma.sq.p.b, 
        		    tuning.c, cov.model.indx, n.batch, 
        	            batch.length, accept.rate, n.omp.threads, verbose, n.report, 
        	            samples.info, chain.info)
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
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$theta.samples)))), 
      			      autoburnin = FALSE)$psrf[, 2]
      out$rhat$lambda.lower.tri <- list()
      for (j in 1:p.svc) {
        lambda.mat <- matrix(0, N, q)
        indx <- (((j - 1) * N * q + 1):(j * N * q))[c(lower.tri(lambda.mat))]
        out$rhat$lambda.lower.tri[[j]] <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					       mcmc(t(a$lambda.samples[indx, ])))), 
          					       autoburnin = FALSE)$psrf[, 2])
      }
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
      out$rhat$beta.comm <- rep(NA, p.occ)
      out$rhat$alpha.comm <- rep(NA, p.det)
      out$rhat$tau.sq.beta <- rep(NA, p.occ)
      out$rhat$tau.sq.alpha <- rep(NA, p.det)
      out$rhat$beta <- rep(NA, p.occ * N)
      out$rhat$alpha <- rep(NA, p.det * N)
      out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 2 * q.p.svc, q.p.svc))
      if (p.det.re > 0) {
        out$rhat$sigma.sq.p <- rep(NA, p.det.re)
      }
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
    colnames(out$alpha.samples) <- coef.names.det
    if (p.det.re > 0) {
      out$sigma.sq.p.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.p.samples))))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.star.samples))))
      tmp.names <- unlist(p.re.level.names)
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$p.re.level.names <- p.re.level.names
    }
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
    loadings.names <- paste(rep(sp.names, times = q), rep(1:q, each = N), sep = '-')
    loadings.names <- paste(rep(loadings.names, times = p.svc), 
        		    rep(x.w.names, each = N * q), sep = '-')
    out$lambda.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$lambda.samples))))
    colnames(out$lambda.samples) <- loadings.names
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    colnames(out$theta.samples) <- theta.names

    # Return things back in the original order. 
    out$z.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$z.samples, 
      								dim = c(N, J, n.post.samples))))
    out$z.samples <- out$z.samples[, order(ord), ]
    out$z.samples <- aperm(out$z.samples, c(3, 1, 2))

    # Account for case when there is only 1 svc. 
    if (p.svc == 1) {
      tmp <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
        						      dim = c(q, J, n.post.samples))))
      tmp <- tmp[, order(ord), , drop = FALSE]
      out$w.samples <- array(NA, dim = c(q, J, p.svc, n.post.samples * n.chains))
      out$w.samples[, , 1, ] <- tmp
    } else {
      out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
        								dim = c(q, J, p.svc, n.post.samples))))
      out$w.samples <- out$w.samples[, order(ord), , , drop = FALSE]
    }
    out$w.samples <- aperm(out$w.samples, c(4, 1, 2, 3))
    out$psi.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$psi.samples, 
      								dim = c(N, J, n.post.samples))))
    out$psi.samples <- out$psi.samples[, order(ord), ]
    out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
    out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples, 
      								dim = c(N, J, n.post.samples))))
    out$like.samples <- out$like.samples[, order(ord), ]
    out$like.samples <- aperm(out$like.samples, c(3, 1, 2))
    if (!binom) {
     tmp <- matrix(NA, J * K.max, p.det)
     tmp[names.long, ] <- X.p
     tmp <- array(tmp, dim = c(J, K.max, p.det))
     tmp <- tmp[order(ord), , ]
     out$X.p <- matrix(tmp, J * K.max, p.det)
     out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
     colnames(out$X.p) <- x.p.names
     tmp <- matrix(NA, J * K.max, p.det.re)
     tmp[names.long, ] <- X.p.re
     tmp <- array(tmp, dim = c(J, K.max, p.det.re))
     tmp <- tmp[order(ord), , ]
     out$X.p.re <- matrix(tmp, J * K.max, p.det.re)
     out$X.p.re <- out$X.p.re[apply(out$X.p.re, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
     colnames(out$X.p.re) <- x.p.re.names
     tmp <- matrix(NA, J * K.max, n.det.re)
    } else {
      out$X.p <- X.p[order(ord), , drop = FALSE]
      out$X.p.re <- X.p.re[order(ord), , drop = FALSE]
    }
    out$X.re <- X.re[order(ord), , drop = FALSE]
    out$X.w <- X.w[order(ord), , drop = FALSE]
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta.comm <- effectiveSize(out$beta.comm.samples)
    out$ESS$alpha.comm <- effectiveSize(out$alpha.comm.samples)
    out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
    out$ESS$tau.sq.alpha <- effectiveSize(out$tau.sq.alpha.samples)
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$alpha <- effectiveSize(out$alpha.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    out$ESS$lambda <- effectiveSize(out$lambda.samples)
    if (p.det.re > 0) {
      out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
    }
    if (p.occ.re > 0) {
      out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
    }
    out$X <- X[order(ord), , drop = FALSE]
    out$X.big <- X.big[order(ord), , , drop = FALSE]
    out$y <- y.big[, order(ord), , drop = FALSE]
    out$call <- cl
    out$n.samples <- n.samples
    out$x.names <- x.names
    out$sp.names <- sp.names
    out$x.p.names <- x.p.names
    out$theta.names <- theta.names
    out$type <- "NNGP"
    out$coords <- coords[order(ord), ]
    out$cov.model.indx <- cov.model.indx
    out$svc.cols <- svc.cols
    out$n.neighbors <- n.neighbors
    out$q <- q
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$species.sds <- species.sds
    out$species.means <- species.means
    out$std.by.sp <- std.by.sp
    out$range.ind <- range.ind[, order(ord)]
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
   
    class(out) <- "svcMsPGOcc"
  }

  out$run.time <- proc.time() - ptm
  return(out)
}
