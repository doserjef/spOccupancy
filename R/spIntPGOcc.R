spIntPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
                       tuning, cov.model = "exponential", NNGP = TRUE, 
                       n.neighbors = 15, search.type = "cb",
                       n.batch, batch.length, accept.rate = 0.43, 
                       n.omp.threads = 1, verbose = TRUE,  
                       n.report = 100, 
                       n.burn = round(.10 * n.batch * batch.length),
                       n.thin = 1, n.chains = 1, 
                       k.fold, k.fold.threads = 1, k.fold.seed = 100, 
                       k.fold.data, k.fold.only = FALSE, ...){

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
  for (q in 1:n.data) {
    if (!is.null(findbars(det.formula[[q]]))) {
      det.re.names <- sapply(findbars(det.formula[[q]]), all.vars)
      for (i in 1:length(det.re.names)) {
        if (is(data$det.covs[[q]][, det.re.names[i]], 'factor')) {
          stop(paste("error: random effect variable ", det.re.names[i], " in data source ", q, " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
        } 
        if (is(data$det.covs[[q]][, det.re.names[i]], 'character')) {
          stop(paste("error: random effect variable ", det.re.names[i], "in data source ", q, " specified as character. Random effect variables must be specified as numeric.", sep = ''))
        }
      }
    }
  }

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
  if (!is.list(det.formula)) {
    stop(paste("error: det.formula must be a list of ", n.data, " formulas", sep = ''))
  }
  X.p <- list()
  X.p.re <- list()
  x.p.names <- list()
  x.p.re.names <- list()
  p.re.level.names <- list()
  for (i in 1:n.data) {
    if (is(det.formula[[i]], 'formula')) {
      tmp <- parseFormula(det.formula[[i]], data$det.covs[[i]])
      X.p[[i]] <- as.matrix(tmp[[1]])
      x.p.names[[i]] <- tmp[[2]]
      if (ncol(tmp[[4]]) > 0) {
        X.p.re[[i]] <- as.matrix(tmp[[4]])
        x.p.re.names[[i]] <- colnames(X.p.re[[i]])
        p.re.level.names[[i]] <- lapply(data$det.covs[[i]][, x.p.re.names[[i]], drop = FALSE],
       		                function (a) sort(unique(a)))
      } else {
        X.p.re[[i]] <- matrix(NA, 0, 0)
        x.p.re.names[[i]] <- NULL
        p.re.level.names[[i]] <- NULL
      }
    } else {
      stop(paste("error: det.formula for data source ", i, " is misspecified", sep = ''))
    }
  }
  x.p.names <- unlist(x.p.names)
  x.p.re.names <- unlist(x.p.re.names)

    # Get basic info from inputs ------------------------------------------
    # Total number of sites
    J.all <- nrow(X)
    if (length(X.p) != n.data | length(y) != n.data) {
      stop(paste("error: y and X.p must be lists of length ", n.data, ".", sep = ''))
    }
    # Number of occupancy parameters 
    p.occ <- ncol(X)
    # Number of occupancy random effect parameters
    p.occ.re <- ncol(X.re)
    # Number of detection random effect parameters for each data set
    p.det.re.by.data <- sapply(X.p.re, ncol)
    # Total number of detection random effects
    p.det.re <- sum(p.det.re.by.data)
    # Number of detection parameters for each data set
    p.det.long <- sapply(X.p, function(a) dim(a)[[2]])
    # Total number of detection parameters
    p.det <- sum(p.det.long)
    n.rep <- lapply(y, function(a1) apply(a1, 1, function(a2) sum(!is.na(a2))))
    # Number of latent occupancy random effect values
    n.occ.re <- length(unlist(apply(X.re, 2, unique)))
    n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
    # Number of levels for each detection random effect
    n.det.re.long <- unlist(sapply(X.p.re, function(a) apply(a, 2, function(b) length(unique(b)))))
    # Number of latent detection random effects for each data set
    n.det.re.by.data <- sapply(sapply(X.p.re, function(a) apply(a, 2, function(b) length(unique(b)))), sum)
    # Total number of detection random effect levels
    n.det.re <- sum(n.det.re.by.data)
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
    # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
    if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
      stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
    }
    if (!missing(k.fold)) {
      if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
        stop("error: k.fold must be a single integer value >= 2")  
      }
    }

    # Get indics to map z to y --------------------------------------------
    X.p.orig <- X.p
    y.big <- y
    names.long <- list()
    names.re.long <- list()
    # Remove missing observations when the covariate data are available but
    # there are missing detection-nondetection data
    for (i in 1:n.data) {
      if (nrow(X.p[[i]]) == length(y[[i]])) {
        X.p[[i]] <- X.p[[i]][!is.na(y[[i]]), , drop = FALSE]
      }
      if (nrow(X.p.re[[i]]) == length(y[[i]]) & p.det.re.by.data[i] > 0) {
        X.p.re[[i]] <- X.p.re[[i]][!is.na(y[[i]]), , drop = FALSE]
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
    # Get indices for WAIC calculation directly in C. 
    J.sum <- sum(J.long)
    waic.J.indx <- unlist(sites) - 1
    waic.n.obs.indx <- list()
    tmp.start <- 0
    for (i in 1:n.data) {
      tmp.vals <- rep(1:J.long[i], K.long.max[i])
      tmp.vals <- tmp.vals[!is.na(c(y[[i]]))]
      waic.n.obs.indx[[i]] <- tmp.vals + tmp.start
      tmp.start <- tmp.start + J.long[i]
    }
    waic.n.obs.indx <- unlist(waic.n.obs.indx) - 1
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

    # Get random effect matrices all set ----------------------------------
    # Make sure each level of each random effect has a different index value
    # for use when fitting the model. 
    # Occurrence REs ------------------
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
      }
    }
    # Detection REs -------------------
    # Need to give a different value for each level across different random
    # effects within a given data set and across a given data set.   
    # Total number of detection random effect observations. 
    n.obs.re <- sum(sapply(X.p.re, nrow))
    curr.max <- 0
    for (i in 1:n.data) {
      if (p.det.re.by.data[i] > 0) {
        for (j in 1:p.det.re.by.data[i]) {
          X.p.re[[i]][, j] <- X.p.re[[i]][, j] + curr.max
          curr.max <- max(X.p.re[[i]]) + 1
        }
      }
    }
    # Combine all detection REs into one group. 
    X.p.re.all <- matrix(NA, n.obs, max(p.det.re.by.data))
    indx <- 1
    for (i in 1:n.data) {
      if (p.det.re.by.data[i] > 0) {
        X.p.re.all[indx:(indx + nrow(X.p.re[[i]]) - 1), 1:p.det.re.by.data[i]] <- X.p.re[[i]]
      }
      indx <- indx + nrow(X.p[[i]])
    }
    # Number of random effects for each row of X.p.re.all
    alpha.n.re.indx <- apply(X.p.re.all, 1, function(a) sum(!is.na(a)))
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
      beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
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
      # Keep track of which detection random effect you're on. 
      alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
      # Index that indicates the column in X.p.re.all
      alpha.col.list <- list()
      indx <- 1
      for (i in 1:n.data) {
        if (p.det.re.by.data[i] > 0) { 
          for (j in 1:p.det.re.by.data[i]) {
            if (j > 1) {
              alpha.col.list[[i]] <- c(alpha.col.list[[i]], rep(j - 1, n.det.re.long[indx]))
	    } else {
              alpha.col.list[[i]] <- rep(j - 1, n.det.re.long[indx])
	    }
            indx <- indx + 1
	  }
	}
      }
      alpha.col.indx <- unlist(alpha.col.list)
      # Index that indicates the data source the random effect corresponds to. 
      alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
    } else {
      sigma.sq.p.inits <- 0
      alpha.star.indx <- 0
      alpha.star.inits <- 0
      alpha.col.indx <- 0
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
    consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, n.det.re, n.data, 
                fixed.sigma.sq, sigma.sq.ig)
    storage.mode(p.det.long) <- "integer"
    storage.mode(n.obs.long) <- "integer"
    storage.mode(J.long) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
    storage.mode(waic.J.indx) <- "integer"
    storage.mode(waic.n.obs.indx) <- "integer"
    # Calculate waic. Might make this an argument in the future. 
    waic.calc <- TRUE
    storage.mode(waic.calc) <- "integer" 
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
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    samples.info <- c(n.burn, n.thin, n.post.samples)
    storage.mode(samples.info) <- "integer"
    storage.mode(tuning.c) <- "double"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    # For detection random effects
    storage.mode(X.p.re.all) <- "integer"
    storage.mode(p.det.re.by.data) <- "integer"
    alpha.level.indx <- sort(unique(c(X.p.re.all)))
    storage.mode(alpha.level.indx) <- "integer"
    storage.mode(n.det.re.long) <- "integer"
    storage.mode(sigma.sq.p.inits) <- "double"
    storage.mode(sigma.sq.p.a) <- "double"
    storage.mode(sigma.sq.p.b) <- "double"
    storage.mode(alpha.star.inits) <- "double"
    storage.mode(alpha.star.indx) <- "integer"
    storage.mode(alpha.n.re.indx) <- "integer"
    storage.mode(alpha.col.indx) <- "integer"
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
          if (p.occ.re > 0) {
            sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
            beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
          }
          if (p.det.re > 0) {
            sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
            alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
          }
        }
        storage.mode(chain.info) <- "integer" 
        # Run the model in C
        out.tmp[[i]] <- .Call("spIntPGOcc", y, X, X.p.all, coords.D, X.re, X.p.re.all, 
                              consts, p.det.long, J.long, n.obs.long, n.occ.re.long, 
                              n.det.re.long, beta.inits, alpha.inits, 
                              sigma.sq.psi.inits, sigma.sq.p.inits, 
                              beta.star.inits, alpha.star.inits, z.inits, w.inits, 
                              phi.inits, sigma.sq.inits, nu.inits, 
                              z.long.indx.c, data.indx.c, alpha.indx.c, 
                              beta.star.indx, beta.level.indx, alpha.star.indx, 
                              alpha.level.indx, alpha.n.re.indx, alpha.col.indx, 
                              waic.J.indx, waic.n.obs.indx, waic.calc, mu.beta, mu.alpha, 
                              Sigma.beta, sigma.alpha, sigma.sq.psi.a, sigma.sq.psi.b, 
                              sigma.sq.p.a, sigma.sq.p.b, phi.a, phi.b, 
                              sigma.sq.a, sigma.sq.b, nu.a, nu.b, tuning.c, cov.model.indx, 
                              n.batch, batch.length, accept.rate,  
                              n.omp.threads, verbose, n.report, samples.info, chain.info)  
        chain.info[1] <- chain.info[1] + 1
      }
      # Calculate R-Hat ---------------
      out <- list()
      out$rhat <- list()
      if (n.chains > 1) {
        out$rhat$beta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
        out$rhat$alpha <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.samples)))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
          out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$theta.samples)))), 
          			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$sigma.sq.psi.samples)))), 
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        }
        if (p.det.re > 0) {
          out$rhat$sigma.sq.p <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$sigma.sq.p.samples)))), 
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
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
      # p.samples is returned as a list, where each element 
      out$p.samples <- do.call(rbind, lapply(out.tmp, function(a) a$p.samples))
      # corresponds to a different data set. 
      tmp <- list()
      indx <- 1
      for (q in 1:n.data) {
        tmp[[q]] <- array(NA, dim = c(J.long[q] * K.long.max[q], n.post.samples * n.chains))
        tmp[[q]][names.long[[q]], ] <- out$p.samples[indx:(indx + n.obs.long[q] - 1), ] 
        tmp[[q]] <- array(tmp[[q]], dim = c(J.long[q], K.long.max[q], n.post.samples * n.chains))
        tmp[[q]] <- aperm(tmp[[q]], c(3, 1, 2))
        indx <- indx + n.obs.long[q]
      }
      out$p.samples <- tmp
      # Likelihood samples for WAIC calculation. 
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
      if (p.occ.re > 0) {
        out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
      }
      if (p.det.re > 0) {
        out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
      }
      out$ESS$theta <- effectiveSize(out$theta.samples)
      out$X <- X
      out$X.p <- X.p.orig
      out$X.re <- X.re
      out$X.p.re <- X.p.re
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
      if (p.occ.re > 0) {
        out$psiRE <- TRUE
      } else {
        out$psiRE <- FALSE
      }
      out$pRELong <- ifelse(p.det.re.by.data > 0, TRUE, FALSE)
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
	coords.D.fit <- coords.D[curr.set.fit, curr.set.fit, drop = FALSE]
	coords.D.0 <- coords.D[curr.set.pred, curr.set.pred, drop = FALSE]
        X.p.fit <- X.p.all[y.indx, , drop = FALSE]
        X.p.0 <- X.p.all[!y.indx, , drop = FALSE]
        X.fit <- X[curr.set.fit, , drop = FALSE]
        X.0 <- X[curr.set.pred, , drop = FALSE]
        J.fit <- nrow(X.fit)
	data.indx.c.fit <- data.indx.c[y.indx]
	data.indx.0 <- data.indx.c[!y.indx] + 1
	data.indx.r.fit <- data.indx.c.fit + 1
	# Random Detection Effects
        X.p.re.all.fit <- X.p.re.all[y.indx, , drop = FALSE]
        X.p.re.all.0 <- X.p.re.all[!y.indx, , drop = FALSE]
        n.det.re.fit <- length(unique(c(X.p.re.all.fit)[!is.na(c(X.p.re.all.fit))]))
	# Number of random effect levels for each data set
	n.det.re.by.data.fit <- rep(NA, n.data)
	for (i in 1:n.data) {
          tmp <- c(X.p.re.all.fit[which(data.indx.r.fit == i, ), , drop = FALSE]) 
          n.det.re.by.data.fit[i] <- length(unique(tmp[!is.na(tmp)]))
	}
	# Number of unique levels in each fitted random effect. 
	n.det.re.long.fit <- rep(NA, p.det.re)
	curr.indx <- 1
        for (i in 1:n.data) {
          tmp <- X.p.re.all.fit[which(data.indx.r.fit == i), , drop = FALSE]
	  n.det.re.fit.curr <- apply(tmp, 2, function(a) length(unique(a[!is.na(a)])))
	  curr.res <- n.det.re.fit.curr[n.det.re.fit.curr > 0]
	  if (length(curr.res > 0)) {
	    n.det.re.long.fit[curr.indx:(curr.indx + length(curr.res) - 1)] <- curr.res 
	    curr.indx <- curr.indx + length(curr.res)
	  }
	}	
        if (p.det.re > 0) {	
          # Keep track of which detection random effect you're on. 
          alpha.star.indx.fit <- rep(0:(p.det.re - 1), n.det.re.long.fit)
          # Index that indicates the column in X.p.re.all
          alpha.col.fit.list <- list()
          indx <- 1
          for (i in 1:n.data) {
            if (p.det.re.by.data[i] > 0) { 
              for (j in 1:p.det.re.by.data[i]) {
                if (j > 1) {
                  alpha.col.fit.list[[i]] <- c(alpha.col.fit.list[[i]], rep(j - 1, n.det.re.long.fit[indx]))
                } else {
                  alpha.col.fit.list[[i]] <- rep(j - 1, n.det.re.long.fit[indx])
                }
                indx <- indx + 1
              }
            }
          }
          alpha.col.indx.fit <- unlist(alpha.col.fit.list)
          # Index that indicates the data source the random effect corresponds to. 
          alpha.star.inits.fit <- rnorm(n.det.re.fit, 0, sqrt(sigma.sq.p.inits[alpha.star.indx.fit + 1]))
          alpha.n.re.indx.fit <- apply(X.p.re.all.fit, 1, function(a) sum(!is.na(a)))
        } else {
          alpha.star.indx.fit <- alpha.star.indx
          alpha.level.indx.fit <- alpha.level.indx
          alpha.star.inits.fit <- alpha.star.inits
	  alpha.n.re.indx.fit <- alpha.n.re.indx
	  alpha.col.indx.fit <- alpha.col.indx
        }
	# Random occurrence effects
        X.re.fit <- X.re[curr.set.fit, , drop = FALSE]
        X.re.0 <- X.re[curr.set.pred, , drop = FALSE]
        n.occ.re.fit <- length(unique(c(X.re.fit)))
        n.occ.re.long.fit <- apply(X.re.fit, 2, function(a) length(unique(a)))
        if (p.occ.re > 0) {	
          beta.star.indx.fit <- rep(0:(p.occ.re - 1), n.occ.re.long.fit)
          beta.level.indx.fit <- sort(unique(c(X.re.fit)))
          beta.star.inits.fit <- rnorm(n.occ.re.fit, 0,
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
	# Site indices for fitted data
	sites.fit <- sapply(sites, function(a) a[a %in% curr.set.fit])
	vals <- sort(unique(unlist(sites.fit)))
	for (q in 1:n.data) {
          # This is needed to ensure that you don't pull the data from data source k.fold.data at sites where there 
	  # is another data source.
	  if (!is.null(k.fold.data)) {
            if (q == k.fold.data) {
              bad.indx <- sites.fit[[q]] %in% curr.set.pred
	      sites.fit[[q]] <- sites.fit[[q]][!bad.indx]
	    }
	  }
          for (j in 1:length(sites.fit[[q]])) {
            sites.fit[[q]][j] <- which(vals == sites.fit[[q]][[j]]) 
	  }
	}
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
	sites.0 <- sapply(sites, function(a) a[a %in% curr.set.pred])
	vals <- sort(unique(unlist(sites.0)))
	for (q in 1:n.data) {
          # Only get missing sites at the single hold out data set if that's the type you are using.
	  if (!is.null(k.fold.data)) {
            if (q != k.fold.data) {
              sites.0[[q]] <- NA
	    } else {
              for (j in 1:length(sites.0[[q]])) {
                sites.0[[q]][j] <- which(vals == sites.0[[q]][[j]]) 
	      }
	    }
	  } else {
            for (j in 1:length(sites.0[[q]])) {
              sites.0[[q]][j] <- which(vals == sites.0[[q]][[j]]) 
	    }
	  }
	}

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
	n.obs.long.fit <- as.vector(table(data.indx.c.fit))
	n.obs.long.0 <- n.obs.long - n.obs.long.fit
	J.long.fit <- as.vector(tapply(z.long.indx.fit, factor(data.indx.c.fit), 
			     FUN = function(a) length(unique(a))))


        storage.mode(y.fit) <- "double"
        storage.mode(z.inits.fit) <- "double"
        storage.mode(X.p.fit) <- "double"
        storage.mode(X.fit) <- "double"
        consts.fit <- c(J.fit, n.obs.fit, p.occ, p.occ.re, n.occ.re.fit, p.det, p.det.re, 
			n.det.re.fit, n.data, fixed.sigma.sq, sigma.sq.ig)
        storage.mode(consts.fit) <- 'integer'
        storage.mode(coords.D.fit) <- "double"
        storage.mode(n.obs.long.fit) <- "integer"
        storage.mode(J.long.fit) <- "integer"
        storage.mode(w.inits.fit) <- "double"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(data.indx.c.fit) <- "integer"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(n.report) <- "integer"
	chain.info[1] <- 1
	storage.mode(chain.info) <- 'integer'
	# Occurrence random effects
        storage.mode(X.re.fit) <- "integer"
        storage.mode(n.occ.re.long.fit) <- "integer"
        storage.mode(beta.star.inits.fit) <- "double"
        storage.mode(beta.star.indx.fit) <- "integer"
        storage.mode(beta.level.indx.fit) <- "integer"
        # For detection random effects
        storage.mode(X.p.re.all.fit) <- "integer"
        storage.mode(p.det.re.by.data) <- "integer"
        alpha.level.indx.fit <- sort(unique(c(X.p.re.all.fit)))
        storage.mode(alpha.level.indx.fit) <- "integer"
        storage.mode(n.det.re.long.fit) <- "integer"
        storage.mode(sigma.sq.p.inits) <- "double"
        storage.mode(sigma.sq.p.a) <- "double"
        storage.mode(sigma.sq.p.b) <- "double"
        storage.mode(alpha.star.inits.fit) <- "double"
        storage.mode(alpha.star.indx.fit) <- "integer"
        storage.mode(alpha.n.re.indx.fit) <- "integer"
        storage.mode(alpha.col.indx.fit) <- "integer"
	waic.calc.fit <- FALSE
	storage.mode(waic.calc.fit) <- "integer"

        out.fit <- .Call("spIntPGOcc", y.fit, X.fit, X.p.fit, coords.D.fit, X.re.fit, 
			 X.p.re.all.fit, consts.fit, p.det.long, 
		         J.long.fit, n.obs.long.fit, n.occ.re.long.fit, 
			 n.det.re.long.fit, beta.inits, alpha.inits, sigma.sq.psi.inits, 
			 sigma.sq.p.inits, beta.star.inits, alpha.star.inits, 
			 z.inits.fit, w.inits.fit, 
		         phi.inits, sigma.sq.inits, nu.inits, 
		         z.long.indx.fit, data.indx.c.fit, alpha.indx.c, 
			 beta.star.indx.fit, beta.level.indx.fit, alpha.star.indx.fit, 
			 alpha.level.indx.fit, alpha.n.re.indx.fit, alpha.col.indx.fit, 
			 waic.J.indx, waic.n.obs.indx, waic.calc.fit, mu.beta, mu.alpha, 
		         Sigma.beta, sigma.alpha, sigma.sq.psi.a, sigma.sq.psi.b, 
			 sigma.sq.p.a, sigma.sq.p.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		         nu.a, nu.b, tuning.c, cov.model.indx, 
		         n.batch, batch.length, accept.rate,  
		         n.omp.threads.fit, verbose.fit, n.report, 
			 samples.info, chain.info)

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
        if (p.occ.re > 0) {
          out.fit$psiRE <- TRUE
        } else {
          out.fit$psiRE <- FALSE	
        }
	class(out.fit) <- "spIntPGOcc"

	# Get RE levels correct for when they aren't supplied at values starting at 1.
	if (p.occ.re > 0) {
	  tmp <- unlist(re.level.names)
	  X.re.0 <- matrix(tmp[c(X.re.0 + 1)], nrow(X.re.0), ncol(X.re.0))
	  colnames(X.re.0) <- x.re.names
	}
        # Predict occurrence at new sites. 
	if (p.occ.re > 0) {X.0 <- cbind(X.0, X.re.0)}
        out.pred <- predict.spIntPGOcc(out.fit, X.0, coords.0, verbose = FALSE)
        # Detection 
        # Get full random effects if certain levels aren't in the fitted values
        if (p.det.re > 0) {
          if (n.det.re.fit != n.det.re) {
            tmp <- matrix(NA, n.det.re, n.post.samples)  
            tmp[alpha.level.indx.fit + 1, ] <- out.fit$alpha.star.samples
            out.fit$alpha.star.samples <- tmp
          }
          # Samples missing NA values
          tmp.indx <- which(apply(out.fit$alpha.star.samples, 1, function(a) sum(is.na(a))) == n.post.samples)
          for (l in tmp.indx) {
            out.fit$alpha.star.samples[l, ] <- rnorm(n.post.samples, 0, 
        					     sqrt(out.fit$sigma.sq.p.samples[alpha.star.indx[l] + 1, ]))
          }
        }

	p.0.samples <- matrix(NA, n.post.samples, nrow(X.p.0))
        like.samples <- rep(NA, nrow(X.p.0))
        for (j in 1:nrow(X.p.0)) {
          if (p.det.re > 0) {
            det.re.sum <- apply(out.fit$alpha.star.samples[which(alpha.level.indx %in% X.p.re.all.0[j, ]), , drop = FALSE], 2, sum)
            p.0.samples[, j] <- logit.inv(X.p.0[j, 1:sum(alpha.indx.r == data.indx.0[j])] %*% out.fit$alpha.samples[which(alpha.indx.r == data.indx.0[j]), ] + det.re.sum)
	  } else {
            p.0.samples[, j] <- logit.inv(X.p.0[j, 1:sum(alpha.indx.r == data.indx.0[j])] %*% out.fit$alpha.samples[which(alpha.indx.r == data.indx.0[j]), ])
	  }
          like.samples[j] <- mean(dbinom(y.0[j], 1, 
					 p.0.samples[, j] * out.pred$z.0.samples[, z.long.indx.0[j] + 1]))
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
    storage.mode(p.det.long) <- "integer"
    consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, n.det.re, n.data, 
                fixed.sigma.sq, sigma.sq.ig)
    storage.mode(n.obs.long) <- "integer"
    storage.mode(J.long) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
    storage.mode(waic.J.indx) <- "integer"
    storage.mode(waic.n.obs.indx) <- "integer"
    # Calculate waic. Might make this an argument in the future. 
    waic.calc <- TRUE
    storage.mode(waic.calc) <- "integer" 
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
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    samples.info <- c(n.burn, n.thin, n.post.samples)
    storage.mode(samples.info) <- "integer"
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
    # For detection random effects
    storage.mode(X.p.re.all) <- "integer"
    storage.mode(p.det.re.by.data) <- "integer"
    alpha.level.indx <- sort(unique(c(X.p.re.all)))
    storage.mode(alpha.level.indx) <- "integer"
    storage.mode(n.det.re.long) <- "integer"
    storage.mode(sigma.sq.p.inits) <- "double"
    storage.mode(sigma.sq.p.a) <- "double"
    storage.mode(sigma.sq.p.b) <- "double"
    storage.mode(alpha.star.inits) <- "double"
    storage.mode(alpha.star.indx) <- "integer"
    storage.mode(alpha.n.re.indx) <- "integer"
    storage.mode(alpha.col.indx) <- "integer"
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
          if (p.occ.re > 0) {
            sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
            beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
          }
          if (p.det.re > 0) {
            sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
            alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
          }
        }
        storage.mode(chain.info) <- "integer" 
        # Run the model in C
        out.tmp[[i]] <- .Call("spIntPGOccNNGP", y, X, X.p.all, coords, X.re, X.p.re.all, 
                              consts, p.det.long, J.long, n.obs.long, n.occ.re.long, n.det.re.long, 
                              n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
                              beta.inits, alpha.inits, sigma.sq.psi.inits, sigma.sq.p.inits, 
                              beta.star.inits, alpha.star.inits, z.inits, w.inits, 
                              phi.inits, sigma.sq.inits, nu.inits, 
                              z.long.indx.c, data.indx.c, alpha.indx.c, 
                              beta.star.indx, beta.level.indx, alpha.star.indx, 
                              alpha.level.indx, alpha.n.re.indx, alpha.col.indx, 
                              waic.J.indx, waic.n.obs.indx, waic.calc, mu.beta, mu.alpha, 
                              Sigma.beta, sigma.alpha, sigma.sq.psi.a, sigma.sq.psi.b, 
                              sigma.sq.p.a, sigma.sq.p.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
                              nu.a, nu.b, tuning.c, cov.model.indx,
                              n.batch, batch.length, accept.rate,  
                              n.omp.threads, verbose, n.report, samples.info, chain.info)
        chain.info[1] <- chain.info[1] + 1
      }
      # Calculate R-Hat ---------------
      out <- list()
      out$rhat <- list()
      if (n.chains > 1) {
        out$rhat$beta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
        out$rhat$alpha <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.samples)))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
        out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$theta.samples)))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$sigma.sq.psi.samples)))), 
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        }
        if (p.det.re > 0) {
          out$rhat$sigma.sq.p <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$sigma.sq.p.samples)))), 
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
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
      # p.samples is returned as a list, where each element 
      out$p.samples <- do.call(rbind, lapply(out.tmp, function(a) a$p.samples))
      # corresponds to a different data set. 
      tmp <- list()
      indx <- 1
      for (q in 1:n.data) {
        tmp[[q]] <- array(NA, dim = c(J.long[q] * K.long.max[q], n.post.samples * n.chains))
        tmp[[q]][names.long[[q]], ] <- out$p.samples[indx:(indx + n.obs.long[q] - 1), ] 
        tmp[[q]] <- array(tmp[[q]], dim = c(J.long[q], K.long.max[q], n.post.samples * n.chains))
        tmp[[q]] <- aperm(tmp[[q]], c(3, 1, 2))
        indx <- indx + n.obs.long[q]
      }
      out$p.samples <- tmp
      # Likelihood samples for WAIC calculation. 
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
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
      if (p.occ.re > 0) {
        out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
      }
      if (p.det.re > 0) {
        out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
      }
      out$ESS$theta <- effectiveSize(out$theta.samples)
      out$X.p <- X.p.orig
      out$X.re <- X.re
      out$X.p.re <- X.p.re
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
      if (p.occ.re > 0) {
        out$psiRE <- TRUE
      } else {
        out$psiRE <- FALSE
      }
      out$pRELong <- ifelse(p.det.re.by.data > 0, TRUE, FALSE)
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
	data.indx.c.fit <- data.indx.c[y.indx]
	data.indx.0 <- data.indx.c[!y.indx] + 1
	data.indx.r.fit <- data.indx.c.fit + 1
	# Random Detection Effects
        X.p.re.all.fit <- X.p.re.all[y.indx, , drop = FALSE]
        X.p.re.all.0 <- X.p.re.all[!y.indx, , drop = FALSE]
        n.det.re.fit <- length(unique(c(X.p.re.all.fit)[!is.na(c(X.p.re.all.fit))]))
	# Number of random effect levels for each data set
	n.det.re.by.data.fit <- rep(NA, n.data)
	for (i in 1:n.data) {
          tmp <- c(X.p.re.all.fit[which(data.indx.r.fit == i, ), , drop = FALSE]) 
          n.det.re.by.data.fit[i] <- length(unique(tmp[!is.na(tmp)]))
	}
	# Number of unique levels in each fitted random effect. 
	n.det.re.long.fit <- rep(NA, p.det.re)
	curr.indx <- 1
        for (i in 1:n.data) {
          tmp <- X.p.re.all.fit[which(data.indx.r.fit == i), , drop = FALSE]
	  n.det.re.fit.curr <- apply(tmp, 2, function(a) length(unique(a[!is.na(a)])))
	  curr.res <- n.det.re.fit.curr[n.det.re.fit.curr > 0]
	  if (length(curr.res > 0)) {
	    n.det.re.long.fit[curr.indx:(curr.indx + length(curr.res) - 1)] <- curr.res 
	    curr.indx <- curr.indx + length(curr.res)
	  }
	}	
        if (p.det.re > 0) {	
          # Keep track of which detection random effect you're on. 
          alpha.star.indx.fit <- rep(0:(p.det.re - 1), n.det.re.long.fit)
          # Index that indicates the column in X.p.re.all
          alpha.col.fit.list <- list()
          indx <- 1
          for (i in 1:n.data) {
            if (p.det.re.by.data[i] > 0) { 
              for (j in 1:p.det.re.by.data[i]) {
                if (j > 1) {
                  alpha.col.fit.list[[i]] <- c(alpha.col.fit.list[[i]], rep(j - 1, n.det.re.long.fit[indx]))
                } else {
                  alpha.col.fit.list[[i]] <- rep(j - 1, n.det.re.long.fit[indx])
                }
                indx <- indx + 1
              }
            }
          }
          alpha.col.indx.fit <- unlist(alpha.col.fit.list)
          # Index that indicates the data source the random effect corresponds to. 
          alpha.star.inits.fit <- rnorm(n.det.re.fit, 0, sqrt(sigma.sq.p.inits[alpha.star.indx.fit + 1]))
          alpha.n.re.indx.fit <- apply(X.p.re.all.fit, 1, function(a) sum(!is.na(a)))
        } else {
          alpha.star.indx.fit <- alpha.star.indx
          alpha.level.indx.fit <- alpha.level.indx
          alpha.star.inits.fit <- alpha.star.inits
	  alpha.n.re.indx.fit <- alpha.n.re.indx
	  alpha.col.indx.fit <- alpha.col.indx
        }
	# Random occurrence effects
        X.re.fit <- X.re[curr.set.fit, , drop = FALSE]
        X.re.0 <- X.re[curr.set.pred, , drop = FALSE]
        n.occ.re.fit <- length(unique(c(X.re.fit)))
        n.occ.re.long.fit <- apply(X.re.fit, 2, function(a) length(unique(a)))
        if (p.occ.re > 0) {	
          beta.star.indx.fit <- rep(0:(p.occ.re - 1), n.occ.re.long.fit)
          beta.level.indx.fit <- sort(unique(c(X.re.fit)))
          beta.star.inits.fit <- rnorm(n.occ.re.fit, 0,
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
	# Site indices for fitted data
	sites.fit <- sapply(sites, function(a) a[a %in% curr.set.fit])
	vals <- sort(unique(unlist(sites.fit)))
	for (q in 1:n.data) {
          # This is needed to ensure that you don't pull the data from data source k.fold.data at sites where there 
	  # is another data source.
	  if (!is.null(k.fold.data)) {
            if (q == k.fold.data) {
              bad.indx <- sites.fit[[q]] %in% curr.set.pred
	      sites.fit[[q]] <- sites.fit[[q]][!bad.indx]
	    }
	  }
          for (j in 1:length(sites.fit[[q]])) {
            sites.fit[[q]][j] <- which(vals == sites.fit[[q]][[j]]) 
	  }
	}
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
	sites.0 <- sapply(sites, function(a) a[a %in% curr.set.pred])
	vals <- sort(unique(unlist(sites.0)))
	for (q in 1:n.data) {
          # Only get missing sites at the single hold out data set if that's the type you are using.
	  if (!is.null(k.fold.data)) {
            if (q != k.fold.data) {
              sites.0[[q]] <- NA
	    } else {
              for (j in 1:length(sites.0[[q]])) {
                sites.0[[q]][j] <- which(vals == sites.0[[q]][[j]]) 
	      }
	    }
	  } else {
            for (j in 1:length(sites.0[[q]])) {
              sites.0[[q]][j] <- which(vals == sites.0[[q]][[j]]) 
	    }
	  }
	}

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
        consts.fit <- c(J.fit, n.obs.fit, p.occ, p.occ.re, n.occ.re.fit, p.det, p.det.re, 
			n.det.re.fit, n.data, fixed.sigma.sq, sigma.sq.ig)
        storage.mode(consts.fit) <- 'integer'
        storage.mode(coords.fit) <- "double"
        storage.mode(n.obs.long.fit) <- "integer"
        storage.mode(J.long.fit) <- "integer"
        storage.mode(w.inits.fit) <- "double"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(data.indx.c.fit) <- "integer"
        storage.mode(n.omp.threads.fit) <- "integer"
        storage.mode(verbose.fit) <- "integer"
        storage.mode(n.report) <- "integer"
	chain.info[1] <- 1
	storage.mode(chain.info) <- 'integer'
	# Occurrence random effects
        storage.mode(X.re.fit) <- "integer"
        storage.mode(n.occ.re.long.fit) <- "integer"
        storage.mode(beta.star.inits.fit) <- "double"
        storage.mode(beta.star.indx.fit) <- "integer"
        storage.mode(beta.level.indx.fit) <- "integer"
        # For detection random effects
        storage.mode(X.p.re.all.fit) <- "integer"
        storage.mode(p.det.re.by.data) <- "integer"
        alpha.level.indx.fit <- sort(unique(c(X.p.re.all.fit)))
        storage.mode(alpha.level.indx.fit) <- "integer"
        storage.mode(n.det.re.long.fit) <- "integer"
        storage.mode(sigma.sq.p.inits) <- "double"
        storage.mode(sigma.sq.p.a) <- "double"
        storage.mode(sigma.sq.p.b) <- "double"
        storage.mode(alpha.star.inits.fit) <- "double"
        storage.mode(alpha.star.indx.fit) <- "integer"
        storage.mode(alpha.n.re.indx.fit) <- "integer"
        storage.mode(alpha.col.indx.fit) <- "integer"
	waic.calc.fit <- FALSE
	storage.mode(waic.calc.fit) <- "integer"
        storage.mode(nn.indx.fit) <- "integer"
        storage.mode(nn.indx.lu.fit) <- "integer"
        storage.mode(u.indx.fit) <- "integer"
        storage.mode(u.indx.lu.fit) <- "integer"
        storage.mode(ui.indx.fit) <- "integer"
        storage.mode(n.neighbors) <- "integer"

        out.fit <- .Call("spIntPGOccNNGP", y.fit, X.fit, X.p.fit, coords.fit, 
			 X.re.fit, X.p.re.all.fit, 
			 consts.fit, p.det.long, J.long.fit, n.obs.long.fit, 
			 n.occ.re.long.fit, n.det.re.long.fit, 
          	         n.neighbors, nn.indx.fit, nn.indx.lu.fit, u.indx.fit, 
			 u.indx.lu.fit, ui.indx.fit, 
          	         beta.inits, alpha.inits, sigma.sq.psi.inits, sigma.sq.p.inits, 
			 beta.star.inits, alpha.star.inits, z.inits.fit, w.inits.fit, 
          	         phi.inits, sigma.sq.inits, nu.inits, 
          	         z.long.indx.fit, data.indx.c.fit, alpha.indx.c, 
			 beta.star.indx.fit, beta.level.indx.fit, alpha.star.indx.fit, 
			 alpha.level.indx.fit, alpha.n.re.indx.fit, alpha.col.indx.fit, 
			 waic.J.indx, waic.n.obs.indx, waic.calc.fit, mu.beta, mu.alpha, 
          	         Sigma.beta, sigma.alpha, sigma.sq.psi.a, sigma.sq.psi.b, 
			 sigma.sq.p.a, sigma.sq.p.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	         nu.a, nu.b, tuning.c, cov.model.indx,
          	         n.batch, batch.length, accept.rate,  
          	         n.omp.threads, verbose.fit, n.report, samples.info, chain.info)
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
        if (p.occ.re > 0) {
          out.fit$psiRE <- TRUE
        } else {
          out.fit$psiRE <- FALSE	
        }
	class(out.fit) <- "spIntPGOcc"

	# Get RE levels correct for when they aren't supplied at values starting at 1.
	if (p.occ.re > 0) {
	  tmp <- unlist(re.level.names)
	  X.re.0 <- matrix(tmp[c(X.re.0 + 1)], nrow(X.re.0), ncol(X.re.0))
	  colnames(X.re.0) <- x.re.names
	}
        # Predict occurrence at new sites. 
	if (p.occ.re > 0) {X.0 <- cbind(X.0, X.re.0)}
        out.pred <- predict.spIntPGOcc(out.fit, X.0, coords.0, verbose = FALSE)
        # Detection 
        # Get full random effects if certain levels aren't in the fitted values
        if (p.det.re > 0) {
          if (n.det.re.fit != n.det.re) {
            tmp <- matrix(NA, n.det.re, n.post.samples)  
            tmp[alpha.level.indx.fit + 1, ] <- out.fit$alpha.star.samples
            out.fit$alpha.star.samples <- tmp
          }
          # Samples missing NA values
          tmp.indx <- which(apply(out.fit$alpha.star.samples, 1, function(a) sum(is.na(a))) == n.post.samples)
          for (l in tmp.indx) {
            out.fit$alpha.star.samples[l, ] <- rnorm(n.post.samples, 0, 
        					     sqrt(out.fit$sigma.sq.p.samples[alpha.star.indx[l] + 1, ]))
          }
        }

	p.0.samples <- matrix(NA, n.post.samples, nrow(X.p.0))
        like.samples <- rep(NA, nrow(X.p.0))
        for (j in 1:nrow(X.p.0)) {
          if (p.det.re > 0) {
            det.re.sum <- apply(out.fit$alpha.star.samples[which(alpha.level.indx %in% X.p.re.all.0[j, ]), , drop = FALSE], 2, sum)
            p.0.samples[, j] <- logit.inv(X.p.0[j, 1:sum(alpha.indx.r == data.indx.0[j])] %*% out.fit$alpha.samples[which(alpha.indx.r == data.indx.0[j]), ] + det.re.sum)
	  } else {
            p.0.samples[, j] <- logit.inv(X.p.0[j, 1:sum(alpha.indx.r == data.indx.0[j])] %*% out.fit$alpha.samples[which(alpha.indx.r == data.indx.0[j]), ])
	  }
          like.samples[j] <- mean(dbinom(y.0[j], 1, 
					 p.0.samples[, j] * out.pred$z.0.samples[, z.long.indx.0[j] + 1]))
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
  out
}
