stPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
		     tuning, cov.model = 'exponential', NNGP = TRUE, 
		     n.neighbors = 15, search.type = 'cb', n.batch, 
		     batch.length, accept.rate = 0.43, n.omp.threads = 1, 
		     verbose = TRUE, ar1 = FALSE, n.report = 100, 
		     n.burn = round(.10 * n.batch * batch.length), 
		     n.thin = 1, n.chains = 1, k.fold, k.fold.threads = 1, 
		     k.fold.seed = 100, k.fold.only = FALSE, ...){

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
  y <- as.array(data$y)
  # Check occupancy covariates
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("Occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- matrix(1, dim(y)[1], dim(y)[2])
    } else {
      stop("error: occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!is.list(data$occ.covs)) {
    stop("error: occ.covs must be a list of matrices, data frames, and/or vectors")
  }
  # Check detection covariates
  if (!'det.covs' %in% names(data)) {
    if (det.formula == ~ 1) {
      if (verbose) {
        message("Detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
      }
      data$det.covs <- list(int = array(1, dim = dim(y)))
    } else {
      stop("error: det.covs must be specified in data for a detection model with covariates")
    }
  }
  if (!is.list(data$det.covs)) {
    stop("error: det.covs must be a list of arrays, matrices, data frames, and/or vectors")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("error: coords must be a matrix or data frame")
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
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[ord, , , drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Occupancy covariates
    for (i in 1:length(data$occ.covs)) {
      if (!is.null(dim(data$occ.covs[[i]]))) { # Time/space varying
        data$occ.covs[[i]] <- data$occ.covs[[i]][ord, , drop = FALSE]
      } else { # Space-varying
        data$occ.covs[[i]] <- data$occ.covs[[i]][ord]
      }
    } 
    for (i in 1:length(data$det.covs)) {
      if (!is.null(dim(data$det.covs[[i]]))) { 
        if (length(dim(data$det.covs[[i]])) == 2) { # Time/space varying
          data$det.covs[[i]] <- data$det.covs[[i]][ord, , drop = FALSE]
	}
        if (length(dim(data$det.covs[[i]])) == 3) { # Time/space/rep varying
          data$det.covs[[i]] <- data$det.covs[[i]][ord, , , drop = FALSE]
	}	
      } else { # Space-varying
        data$det.covs[[i]] <- data$det.covs[[i]][ord]
      }
    }
  }

  # Reformat covariates ---------------------------------------------------
  # Get detection covariates in proper format
  # First subset detection covariates to only use those that are included in the analysis.
  data$det.covs <- data$det.covs[names(data$det.covs) %in% all.vars(det.formula)]
  # Null model support
  if (length(data$det.covs) == 0) {
    data$det.covs <- list(int = array(1, dim = dim(y)))
  }
  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  # Ordered by visit, year, then site. 
  data$det.covs <- data.frame(lapply(data$det.covs, function(a) unlist(c(a))))
  # Get detection covariates in site x year x replicate format
  if (nrow(data$det.covs) == dim(y)[1]) { # if only site-level covariates. 
    data$det.covs <- as.data.frame(mapply(rep, data$det.covs, dim(y)[2] * dim(y)[3]))
  } else if (nrow(data$det.covs) == dim(y)[1] * dim(y)[2]) { # if only site/year level covariates
    data$det.covs <- as.data.frame(mapply(rep, data$det.covs, dim(y)[3]))
  }
  y.big <- y
  # Get occurrence covariates in proper format
  # Subset covariates to only use those that are included in the analysis
  data$occ.covs <- data$occ.covs[names(data$occ.covs) %in% all.vars(occ.formula)]
  # Null model support
  if (length(data$occ.covs) == 0) {
    data$occ.covs <- list(int = matrix(1, nrow = dim(y)[1], ncol = dim(y)[2]))
  }
  # Ordered by year, then site within year. 
  data$occ.covs <- data.frame(lapply(data$occ.covs, function(a) unlist(c(a))))
  # Check if only site-level covariates are included
  if (nrow(data$occ.covs) == dim(y)[1]) {
    data$occ.covs <- as.data.frame(mapply(rep, data$occ.covs, dim(y)[2]))
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
  for (i in 1:ncol(data$det.covs)) {
    if (sum(is.na(data$det.covs[, i])) > sum(is.na(y.big))) {
      stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  if (det.formula != ~ 1) {
    # Misalignment between y and det.covs
    y.missing <- which(is.na(y))
    det.covs.missing <- lapply(data$det.covs, function(a) which(is.na(a)))
    for (i in 1:length(det.covs.missing)) {
      tmp.indx <- !(y.missing %in% det.covs.missing[[i]])
      if (sum(tmp.indx) > 0) {
        if (i == 1 & verbose) {
          message("There are missing values in data$y with corresponding non-missing values in data$det.covs.\nRemoving these site/time/replicate combinations for fitting the model.")
        }
        data$det.covs[y.missing, i] <- NA
      }
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

  # Check ar1 parameter ---------------------------------------------------
  if (!(ar1 %in% c(TRUE, FALSE))) {
    stop("error: ar1 must be either TRUE or FALSE")
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
    x.p.names <- tmp[[2]]
    X.p.re <- as.matrix(tmp[[4]])
    x.p.re.names <- colnames(X.p.re)
  } else {
    stop("error: det.formula is misspecified")
  }
  p.re.level.names <- lapply(data$det.covs[, x.p.re.names, drop = FALSE],
			     function (a) sort(unique(a)))

  # Get basic info from inputs ------------------------------------------
  # Number of sites
  J <- dim(y.big)[1]
  # Total number of years
  n.years.max <- dim(y.big)[2]
  # Number of years for each site
  n.years <- rep(NA, J)
  for (j in 1:J) {
    n.years[j] <- sum(apply(y.big[j, , ], 1, function(a) sum(!is.na(a))) != 0)
  }
  # Number of occupancy effects
  p.occ <- ncol(X)
  # Number of detection effects
  p.det <- ncol(X.p)
  # Number of occurrence random effect parameters
  p.occ.re <- ncol(X.re)
  # Number of detection random effect parameters
  p.det.re <- ncol(X.p.re)
  # Number of latent occupancy random effect values
  n.occ.re <- length(unlist(apply(X.re, 2, unique)))
  n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of latent detection random effect values
  n.det.re <- length(unlist(apply(X.p.re, 2, unique)))
  n.det.re.long <- apply(X.p.re, 2, function(a) length(unique(a)))
  if (p.det.re == 0) n.det.re.long <- 0
  # Number of replicates at each site/year combo
  # This also inherently contains info on which sites were sampled
  # in which years. 
  n.rep <- apply(y.big, c(1, 2), function(a) sum(!is.na(a)))
  # Max number of repeat visits
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

  # Get indices for mapping different values in Z. 
  z.site.indx <- rep(1:J, n.years.max) - 1
  z.year.indx <- rep(1:n.years.max, each = J) - 1
  z.dat.indx <- c(ifelse(K > 0, 1, 0))

  # Get indices to map z to y -------------------------------------------
  # This corresponds to the specific value in the n.years.max * J length 
  # vector of latent occurrence values, and corresponds with the z.site.indx
  # and z.year.indx
  z.long.indx <- rep(1:(J * n.years.max), dim(y.big)[3])
  # Order of c(y.big): All sites year 1 v 1, All sites year 2 v 1, ....
  z.long.indx <- z.long.indx[!is.na(c(y.big))]
  # Subtract 1 for indices in C
  z.long.indx <- z.long.indx - 1
  # Index that links observations to sites. 
  z.long.site.indx <- rep(rep(1:J, n.years.max), K.max)
  z.long.site.indx <- z.long.site.indx[!is.na(c(y.big))]
  # Subtract 1 for indices in C
  z.long.site.indx <- z.long.site.indx - 1
  # Order: visit, year within visit, site within year
  y <- c(y)
  # Need these later for naming things
  names.long <- which(!is.na(y))
  # Null model support for missing values
  if (nrow(X.p) != length(names.long)) {
    X.p <- X.p[which(!is.na(y)), , drop = FALSE]
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
    if (priors$beta.normal[1] == 'fixed') {
      fixed.params[which(all.params == 'beta')] <- TRUE 
      mu.beta <- rep(0, p.occ)
      Sigma.beta <- diag(p.occ)
    } else {
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
    }
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
    if (priors$alpha.normal[1] == 'fixed') {
      fixed.params[which(all.params == 'alpha')] <- TRUE 
      mu.alpha <- rep(0, p.det)
      Sigma.alpha <- diag(p.det)
    } else {
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
    }
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
  if (!NNGP) {
    coords.D <- iDist(coords)
  }
  if ("phi.unif" %in% names(priors)) {
    if (priors$phi.unif[1] == 'fixed') {
      fixed.params[which(all.params == 'phi')] <- TRUE 
      phi.a <- 1
      phi.b <- 1
      if ((cov.model == 'matern') & (priors$nu.unif[1] != 'fixed')) {
        message("phi is specified as fixed in priors$phi.unif but nu is not, which is not allowed. Continuing to fit the model with both phi and nu fixed.")
      }
    } else {
      if (!is.vector(priors$phi.unif) | !is.atomic(priors$phi.unif) | length(priors$phi.unif) != 2) {
        stop("error: phi.unif must be a vector of length 2 with elements corresponding to phi's lower and upper bounds")
      }
      phi.a <- priors$phi.unif[1]
      phi.b <- priors$phi.unif[2]
    }
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
  # Check if both an ig and uniform prior are specified
  if (("sigma.sq.ig" %in% names(priors)) & ("sigma.sq.unif" %in% names(priors))) {
    stop("error: cannot specify both an IG and a uniform prior for sigma.sq")
  }
  if ("sigma.sq.ig" %in% names(priors)) { # inverse-gamma prior.
    sigma.sq.ig <- TRUE
    if (priors$sigma.sq.ig[1] == 'fixed') {
      fixed.params[which(all.params == 'sigma.sq')] <- TRUE 
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
      fixed.params[which(all.params == 'sigma.sq')] <- TRUE 
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
    if (priors$nu.unif[1] == 'fixed') {
      fixed.params[which(all.params == 'phi')] <- TRUE 
      nu.a <- 1
      nu.b <- 1
      if ((cov.model == 'matern') & (priors$phi.unif[1] != 'fixed')) {
        message("nu is specified as fixed in priors$nu.unif but phi is not, which is not allowed. Continuing to fit the model with both phi and nu fixed.")
      }
    } else {
      if (!is.vector(priors$nu.unif) | !is.atomic(priors$nu.unif) | length(priors$nu.unif) != 2) {
        stop("error: nu.unif must be a vector of length 2 with elements corresponding to nu's lower and upper bounds")
      }
      nu.a <- priors$nu.unif[1]
      nu.b <- priors$nu.unif[2]
    }
  } else {
    nu.a <- 0
    nu.b <- 0
  }
  # sigma.sq.psi --------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi.ig" %in% names(priors)) {
      if (priors$sigma.sq.psi.ig[1] == 'fixed') {
        fixed.params[which(all.params == 'sigma.sq.psi')] <- TRUE
        sigma.sq.psi.a <- rep(1, p.occ.re)
	sigma.sq.psi.b <- rep(1, p.occ.re)
      } else {
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
      if (priors$sigma.sq.p.ig[1] == 'fixed') {
        fixed.params[which(all.params == 'sigma.sq.p')] <- TRUE
        sigma.sq.p.a <- rep(1, p.det.re)
	sigma.sq.p.b <- rep(1, p.det.re)
      } else {
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
  if (ar1) {
    # rho -----------------------------
    if ("rho.unif" %in% names(priors)) {
      if (!is.vector(priors$rho.unif) | !is.atomic(priors$rho.unif) | length(priors$rho.unif) != 2) {
        stop("error: rho.unif must be a vector of length 2 with elements corresponding to rho's lower and upper bounds")
      }
      rho.a <- priors$rho.unif[1]
      rho.b <- priors$rho.unif[2]
    } else {
      if (verbose) {
        message("No prior specified for rho.unif.\nSetting uniform bounds to -1 and 1.\n")
      }
      rho.a <- -1 
      rho.b <- 1 
    }
    
    # sigma.sq.t.t ----------------------
    if ("sigma.sq.t.ig" %in% names(priors)) { 
      if (!is.vector(priors$sigma.sq.t.ig) | !is.atomic(priors$sigma.sq.t.ig) | length(priors$sigma.sq.t.ig) != 2) {
        stop("error: sigma.sq.t.ig must be a vector of length 2 with elements corresponding to sigma.sq.t's shape and scale parameters")
      }
      sigma.sq.t.a <- priors$sigma.sq.t.ig[1]
      sigma.sq.t.b <- priors$sigma.sq.t.ig[2]
    } else {
      if (verbose) {
        message("No prior specified for sigma.sq.t.\nUsing an inverse-Gamma prior with the shape parameter set to 2 and scale parameter to 0.5.\n")
      }
      sigma.sq.t.a <- 2
      sigma.sq.t.b <- 0.5
    }
  } else {
    rho.a <- 0
    rho.b <- 0
    sigma.sq.t.a <- 0
    sigma.sq.t.b <- 0
  }

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # z -------------------------------
  #ORDER: stored in column-major order as a vector, with values sorted by year, then 
  #       site within year. 
  if ("z" %in% names(inits)) {
    z.inits <- inits$z
    if (!is.matrix(z.inits) & !is.data.frame(z.inits)) {
      stop(paste("error: initial values for z must be a matrix or data frame with ",
      	   J, " rows and ", n.years.max, " columns.", sep = ""))
    }
    if (nrow(z.inits) != J | ncol(z.inits) != n.years.max) {
      stop(paste("error: initial values for z must be a matrix or data frame with ",
      	   J, " rows and ", n.years.max, " columns.", sep = ""))
    }
    # Reorder the user supplied inits values for NNGP models
    if (NNGP) {
      z.inits <- z.inits[ord, ]
    }
    z.test <- apply(y.big, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    init.test <- sum(z.inits < z.test, na.rm = TRUE)
    if (init.test > 0) {
      stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site/year combination.")
    }
  } else {
    z.inits <- apply(y.big, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
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
  # w ---------------------------------
  # Just set initial W values to 0. 
  w.inits <- rep(0, J)
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
  if (ar1) {
    # rho -----------------------------
    if ("rho" %in% names(inits)) {
      rho.inits <- inits[["rho"]]
      if (length(rho.inits) != 1) {
        stop("error: initial values for rho must be of length 1")
      }
    } else {
      rho.inits <- runif(1, rho.a, rho.b)
      if (verbose) {
        message("rho is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
      }
    }
    # sigma.sq.t ------------------------
    if ("sigma.sq.t" %in% names(inits)) {
      sigma.sq.t.inits <- inits[["sigma.sq.t"]]
      if (length(sigma.sq.t.inits) != 1) {
        stop("error: initial values for sigma.sq.t must be of length 1")
      }
    } else {
      sigma.sq.t.inits <- runif(1, 0.5, 10)
      if (verbose) {
        message("sigma.sq.t is not specified in initial values.\nSetting initial value to random value between 0.5 and 10\n")
      }
    }
  } else {
    rho.inits <- 0
    sigma.sq.t.inits <- 0
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
  sigma.sq.tuning <- 0
  phi.tuning <- 0
  nu.tuning <- 0
  rho.tuning <- 0
  sigma.sq.t.tuning <- 0
  if (missing(tuning)) {
    phi.tuning <- 1
    rho.tuning <- 1
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
    if (ar1) {
      # rho ---------------------------
      if(!"rho" %in% names(tuning)) {
        stop("error: rho must be specified in tuning value list")
      }
      rho.tuning <- tuning$rho
      if (length(rho.tuning) != 1) {
        stop("error: rho tuning must be a single value")
      } 
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
  # Log the tuning values since they are used in the AMCMC.
  # Need to shift the order depending on what's in the model. 
  if (ar1) {
    if (cov.model == 'matern') {
      tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning, 
			sigma.sq.t.tuning, rho.tuning))
    } else {
      tuning.c <- log(c(sigma.sq.tuning, phi.tuning,  
			sigma.sq.t.tuning, rho.tuning))
    }
  } else {
    tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  }
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA
  curr.chain <- 1

  if (!NNGP) {
    stop("error: stPGOcc is currently only implemented for NNGP models. Please set NNGP = TRUE") 
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
    consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, 
		n.det.re, n.years.max)
    storage.mode(consts) <- "integer"
    storage.mode(coords) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(z.year.indx) <- "integer"
    storage.mode(z.dat.indx) <- "integer"
    storage.mode(z.long.site.indx) <- "integer"
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
    storage.mode(n.burn) <- "integer"
    storage.mode(n.thin) <- "integer"
    storage.mode(curr.chain) <- "integer"
    storage.mode(n.chains) <- "integer"
    storage.mode(fixed.params) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"
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
    storage.mode(n.occ.re.long) <- "integer"
    storage.mode(beta.star.inits) <- "double"
    storage.mode(beta.star.indx) <- "integer"  
    # AR1 parameters
    storage.mode(ar1) <- "integer"
    ar1.vals <- c(rho.a, rho.b, sigma.sq.t.a, sigma.sq.t.b, 
                  rho.inits, sigma.sq.t.inits)
    storage.mode(ar1.vals) <- "double"
                                               
    # stPGOccNNGP                             
    out.tmp <- list()                          
    out <- list()                              
    if (!k.fold.only) {                        
      for (i in 1:n.chains) {                  
        # Change initial values if i > 1       
        if ((i > 1) & (!fix.inits)) {
          beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
          alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
          if (sigma.sq.ig) {
            sigma.sq.inits <- rigamma(1, sigma.sq.a, sigma.sq.b)
	  } else {
            sigma.sq.inits <- runif(1, sigma.sq.a, sigma.sq.b)
	  }
          phi.inits <- runif(1, phi.a, phi.b)
          if (cov.model == 'matern') {
            nu.inits <- runif(1, nu.a, nu.b)
          }
          if (p.det.re > 0) {
            sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
            alpha.star.inits <- rnorm(n.det.re, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
          }
          if (p.occ.re > 0) {
            sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
            beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
          }
	  if (ar1) {
            ar1.vals[5] <- runif(1, rho.a, rho.b)
            ar1.vals[6] <- runif(1, 0.5, 10)	
	  }
        }
        storage.mode(curr.chain) <- "integer"
        # Note that the upper limit on the number of arguments is 65, which 
        # you're getting close to. 
        out.tmp[[i]] <- .Call("stPGOccNNGP", y, X, X.p, coords, X.re, X.p.re, 
                              consts, n.occ.re.long, n.det.re.long,
                              n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
                              beta.inits, alpha.inits, sigma.sq.psi.inits, 
                              sigma.sq.p.inits, beta.star.inits, alpha.star.inits, 
          		      phi.inits, sigma.sq.inits, nu.inits,
                              w.inits, z.inits, z.long.indx, z.year.indx,
                              z.dat.indx, z.long.site.indx, 
                              beta.star.indx, beta.level.indx,
                              alpha.star.indx, alpha.level.indx,
                              mu.beta, Sigma.beta, mu.alpha, Sigma.alpha, 
                              phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b,
                              sigma.sq.psi.a, sigma.sq.psi.b, sigma.sq.p.a, sigma.sq.p.b,
			      ar1, ar1.vals,
                              tuning.c, cov.model.indx, n.batch, batch.length, accept.rate, 
                              n.omp.threads, verbose, n.report,  
                              n.burn, n.thin, n.post.samples, curr.chain, n.chains, sigma.sq.ig)
        curr.chain <- curr.chain + 1
      }
      out <- list()

      # Calculate R-Hat ---------------
      if (ar1) {
        if (cov.model == 'matern') {
          n.theta <- 5
          theta.names <- c('sigma.sq', 'phi', 'nu', 'sigma.sq.t', 'rho')
	} else {
          n.theta <- 4
          theta.names <- c('sigma.sq', 'phi', 'sigma.sq.t', 'rho')
	}
      } else {
        if (cov.model == 'matern') {
          n.theta <- 3
          theta.names <- c('sigma.sq', 'phi', 'nu')
	} else {
          n.theta <- 2
          theta.names <- c('sigma.sq', 'phi')
	}
      }
      out <- list()
      out$rhat <- list()
      if (n.chains > 1) {
        out$rhat$beta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2]
        out$rhat$alpha <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2]
        out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$theta.samples)))), 
          			      autoburnin = FALSE)$psrf[, 2]
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
        out$rhat$theta <- rep(NA, n.theta)
        if (p.det.re > 0) {
          out$rhat$sigma.sq.p <- rep(NA, p.det.re)
        }
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
        }
      }

      # Put everything into an MCMC objects
      out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(do.call(rbind, 
        				lapply(out.tmp, function(a) t(a$alpha.samples))))
      colnames(out$alpha.samples) <- x.p.names
      out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
      colnames(out$theta.samples) <- theta.names
      if (ar1) {
        out$eta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$eta.samples))))
      } 
      # Return things back in the original order
      out$coords <- coords[order(ord), ]
      out$X <- array(X, dim = c(J, n.years.max, p.occ))
      out$X <- out$X[order(ord), , , drop = FALSE]
      dimnames(out$X)[[3]] <- x.names
      out$X.re <- array(X.re, dim = c(J, n.years.max, p.occ.re))
      out$X.re <- out$X.re[order(ord), , , drop = FALSE]
      dimnames(out$X.re)[[3]] <- x.re.names
      out$w.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples))))
      out$w.samples <- mcmc(out$w.samples[, order(ord), drop = FALSE])
      out$z.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$z.samples, 
        								dim = c(J, n.years.max, n.post.samples))))
      out$z.samples <- out$z.samples[order(ord), , ]
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$psi.samples, 
        								dim = c(J, n.years.max, n.post.samples))))
      out$psi.samples <- out$psi.samples[order(ord), , ]
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples, 
        								dim = c(J, n.years.max, n.post.samples))))
      out$like.samples <- out$like.samples[order(ord), , ]
      out$like.samples <- aperm(out$like.samples, c(3, 1, 2))
      out$y <- y.big[order(ord), , , drop = FALSE]
      # Get all detection covariate stuff in the right order
      tmp <- matrix(NA, J * K.max * n.years.max, p.det)
      tmp[names.long, ] <- X.p
      tmp <- array(tmp, dim = c(J, n.years.max, K.max, p.det))
      tmp <- tmp[order(ord), , , , drop = FALSE]
      out$X.p <- matrix(tmp, J * K.max * n.years.max, p.det)
      out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
      colnames(out$X.p) <- x.p.names
      tmp <- matrix(NA, J * K.max * n.years.max, p.det.re)
      tmp[names.long, ] <- X.p.re
      tmp <- array(tmp, dim = c(J, n.years.max, K.max, p.det.re))
      tmp <- tmp[order(ord), , , , drop = FALSE]
      out$X.p.re <- matrix(tmp, J * K.max * n.years.max, p.det.re)
      out$X.p.re <- out$X.p.re[apply(out$X.p.re, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
      colnames(out$X.p.re) <- x.p.re.names
      tmp <- matrix(NA, J * K.max * n.years.max, n.det.re)

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
      out$ESS$theta <- effectiveSize(out$theta.samples)
      if (p.det.re > 0) {
        out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
      }
      if (p.occ.re > 0) {
        out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
      }
      out$call <- cl
      out$n.samples <- batch.length * n.batch
      out$n.neighbors <- n.neighbors
      out$cov.model.indx <- cov.model.indx
      out$type <- "NNGP"
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$n.chains <- n.chains
      out$ar1 <- as.logical(ar1)
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
        y.indx <- !((z.long.site.indx + 1) %in% curr.set)
	year.indx <- !((z.site.indx + 1) %in% curr.set)
        y.fit <- y[y.indx]
	y.0 <- y[!y.indx]
	y.big.fit <- y.big[-curr.set, , , drop = FALSE]
	y.big.0 <- y.big[curr.set, , , drop = FALSE]
	z.inits.fit <- z.inits[-curr.set, ]
	X.p.fit <- X.p[y.indx, , drop = FALSE]
	X.p.0 <- X.p[!y.indx, , drop = FALSE]
	X.fit <- X[year.indx, , drop = FALSE]
	X.0 <- X[!year.indx, , drop = FALSE]
        coords.fit <- coords[-curr.set, , drop = FALSE]
        coords.0 <- coords[curr.set, , drop = FALSE]
	J.fit <- nrow(coords.fit)
	J.0 <- nrow(coords.0)
	K.fit <- apply(K, 2, function(a) a[-curr.set])
	K.0 <- apply(K, 2, function(a) a[curr.set])
	n.obs.fit <- nrow(X.p.fit)
	n.obs.0 <- nrow(X.p.0)
	# Random Detection Effects
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
        X.re.fit <- X.re[year.indx, , drop = FALSE]
        X.re.0 <- X.re[!year.indx, , drop = FALSE]
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
	# Get all the indices for linking to the new fitted data. 
        z.long.indx.fit <- rep(1:(J.fit * n.years.max), K.max)
        z.0.long.indx <- rep(1:(J.0 * n.years.max), K.max)
        # Order of c(y.big): All sites year 1 v 1, All sites year 2 v 1, ....
        z.long.indx.fit <- z.long.indx.fit[!is.na(c(y.big.fit))]
	z.0.long.indx <- z.0.long.indx[!is.na(c(y.big.0))]
        # Subtract 1 for indices in C (only do this for the fitted ones)
        z.long.indx.fit <- z.long.indx.fit - 1
        # Index that links observations to sites. 
        z.long.site.indx.fit <- rep(rep(1:J.fit, n.years.max), K.max)
        z.long.site.indx.fit <- z.long.site.indx.fit[!is.na(c(y.big.fit))]
        # Subtract 1 for indices in C
        z.long.site.indx.fit <- z.long.site.indx.fit - 1
        z.site.indx.fit <- rep(1:J.fit, n.years.max) - 1
        z.year.indx.fit <- rep(1:n.years.max, each = J.fit) - 1
        z.dat.indx.fit <- c(ifelse(K.fit > 0, 1, 0))

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
        storage.mode(z.inits.fit) <- "double"
        storage.mode(X.p.fit) <- "double"
        storage.mode(X.fit) <- "double"
        consts.fit <- c(J.fit, n.obs.fit, p.occ, p.occ.re, n.occ.re.fit, 
			p.det, p.det.re, n.det.re.fit, n.years.max)
        storage.mode(consts.fit) <- "integer"
        storage.mode(coords.fit) <- "double"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(z.year.indx.fit) <- "integer"
        storage.mode(z.dat.indx.fit) <- "integer"
        storage.mode(z.long.site.indx.fit) <- "integer"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(nn.indx.fit) <- "integer"
        storage.mode(nn.indx.lu.fit) <- "integer"
        storage.mode(u.indx.fit) <- "integer"
        storage.mode(u.indx.lu.fit) <- "integer"
        storage.mode(ui.indx.fit) <- "integer"
	curr.chain <- 1
        storage.mode(curr.chain) <- "integer"
        # For detection random effects
        storage.mode(X.p.re.fit) <- "integer"
        storage.mode(alpha.level.indx.fit) <- "integer"
        storage.mode(n.det.re.long.fit) <- "integer"
        storage.mode(alpha.star.inits.fit) <- "double"
        storage.mode(alpha.star.indx.fit) <- "integer"
        # For occurrence random effects
        storage.mode(X.re.fit) <- "integer"
	storage.mode(n.occ.re.long.fit) <- "integer"
        storage.mode(beta.level.indx.fit) <- "integer"
        storage.mode(beta.star.inits.fit) <- "double"
        storage.mode(beta.star.indx.fit) <- "integer"

        out.fit <- .Call("stPGOccNNGP", y.fit, X.fit, X.p.fit, coords.fit, 
			 X.re.fit, X.p.re.fit, 
                         consts.fit, n.occ.re.long.fit, n.det.re.long.fit,
                         n.neighbors, nn.indx.fit, nn.indx.lu.fit, u.indx.fit, 
			 u.indx.lu.fit, ui.indx.fit,
                         beta.inits, alpha.inits, sigma.sq.psi.inits, 
                         sigma.sq.p.inits, beta.star.inits.fit, alpha.star.inits.fit, 
			 phi.inits, sigma.sq.inits, nu.inits,
                         w.inits, z.inits.fit, z.long.indx.fit, z.year.indx.fit,
                         z.dat.indx.fit, z.long.site.indx.fit, 
                         beta.star.indx.fit, beta.level.indx.fit,
                         alpha.star.indx.fit, alpha.level.indx.fit,
                         mu.beta, Sigma.beta, mu.alpha, Sigma.alpha, 
                         phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b,
                         sigma.sq.psi.a, sigma.sq.psi.b, sigma.sq.p.a, sigma.sq.p.b,
			 ar1, ar1.vals,
                         tuning.c, cov.model.indx, n.batch, batch.length, accept.rate, 
                         n.omp.threads.fit, verbose.fit, n.report,  
                         n.burn, n.thin, n.post.samples, curr.chain, n.chains, sigma.sq.ig)

        out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
        colnames(out.fit$beta.samples) <- x.names
        out.fit$alpha.samples <- mcmc(t(out.fit$alpha.samples))
        colnames(out.fit$alpha.samples) <- x.p.names
        out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
	# colnames(out.fit$theta.samples) <- theta.names
        out.fit$w.samples <- mcmc(t(out.fit$w.samples))
        if (ar1) {
          out.fit$eta.samples <- mcmc(t(out.fit$eta.samples))
	  out.fit$ar1 <- TRUE
        } else {
          out.fit$ar1 <- FALSE
	}
        out.fit$X <- array(X.fit, dim = c(J.fit, n.years.max, p.occ))
	dimnames(out.fit$X)[[3]] <- x.names
        out.fit$y <- y.big.fit
        out.fit$X.p <- X.p.fit
        out.fit$call <- cl
        out.fit$type <- "NNGP"
        out.fit$n.neighbors <- n.neighbors
        out.fit$n.samples <- n.samples
        out.fit$coords <- coords.fit
        out.fit$cov.model.indx <- cov.model.indx
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
          out.fit$X.re <- array(X.re.fit, dim = c(J.fit, n.years.max, p.occ.re))
	  dimnames(out.fit$X.re)[[3]] <- x.re.names
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
        class(out.fit) <- "stPGOcc"

	# Predict occurrence at new sites
        if (p.occ.re > 0) {
          tmp <- unlist(re.level.names)
          X.re.0 <- matrix(tmp[c(X.re.0 + 1)], nrow(X.re.0), ncol(X.re.0))
          colnames(X.re.0) <- x.re.names
          X.0 <- cbind(X.0, X.re.0)
        }
	tmp.names <- colnames(X.0)
	X.0 <- array(X.0, dim = c(J.0, n.years.max, ncol(X.0)))
	dimnames(X.0)[[3]] <- tmp.names
        out.pred <- predict.stPGOcc(out.fit, X.0, coords.0, 
				     t.cols = 1:n.years.max, 
				     verbose = FALSE)

        # Predict detection values 
        if (p.det.re > 0) {
          tmp <- unlist(p.re.level.names)
          X.p.re.0 <- matrix(tmp[c(X.p.re.0 + 1)], nrow(X.p.re.0), ncol(X.p.re.0))
          colnames(X.p.re.0) <- x.p.re.names
        }
        if (p.det.re > 0) {X.p.0 <- cbind(X.p.0, X.p.re.0)}
        tmp.names <- colnames(X.p.0)
        X.p.0 <- array(X.p.0, dim = c(nrow(X.p.0), 1, ncol(X.p.0)))
        dimnames(X.p.0)[[3]] <- tmp.names
        out.p.pred <- predict.stPGOcc(out.fit, X.p.0, t.cols = 1, type = 'detection')

        like.samples <- rep(NA, nrow(X.p.0))
        out.pred$z.0.samples <- aperm(out.pred$z.0.samples, c(2, 3, 1))
        out.pred$z.0.samples <- matrix(out.pred$z.0.samples, J.0 * n.years.max, n.post.samples)
        for (j in 1:nrow(X.p.0)) {
          like.samples[j] <- mean(dbinom(y.0[j], 1, 
              			   out.p.pred$p.0.samples[, j, 1] * out.pred$z.0.samples[z.0.long.indx[j], ]))
        }
        sum(log(like.samples), na.rm = TRUE)
      } # parallel loop
      model.deviance <- -2 * model.deviance
      # Return objects from cross-validation
      out$k.fold.deviance <- model.deviance
      stopImplicitCluster()
    } # cross-validation
  } # NNGP
  class(out) <- "stPGOcc"
  out$run.time <- proc.time() - ptm
  out
}
