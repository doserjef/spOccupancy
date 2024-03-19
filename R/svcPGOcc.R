svcPGOcc <- function(occ.formula, det.formula, data, inits, priors, tuning, 
                     svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
		     n.neighbors = 15, search.type = 'cb', n.batch, 
		     batch.length, accept.rate = 0.43,
		     n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		     n.burn = round(.10 * n.batch * batch.length), 
		     n.thin = 1, n.chains = 1, k.fold, k.fold.threads = 1, 
		     k.fold.seed = 100, k.fold.only = FALSE, ...){

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
  # Check grid index
  if (!'grid.index' %in% names(data)) {
    if (nrow(data$coords) != nrow(data$y)) {
      stop("data$grid.index must be specified if nrow(data$coords) != nrow(data$y)")
    }
    grid.index <- 1:nrow(coords)
  } else {
    if (!is.atomic(data$grid.index) | !is.numeric(data$grid.index)) {
      stop("data$grid.index must be a numeric vector")
    }
    if (length(data$grid.index) < nrow(data$coords)) {
      stop("length(data$grid.index) must be greater than or equal to nrow(data$coords)")
    }
    grid.index <- data$grid.index
  }
  # Check if all spatial coordinates are unique. 
  unique.coords <- unique(data$coords)
  if (nrow(unique.coords) < nrow(data$coords)) {
    stop("coordinates provided in coords are not all unique. spOccupancy requires each site to have its own unique pair of spatial coordinates. This may be the result of an error in preparing the data list, or you will need to change what you consider a 'site' in order to meet this requirement. Alternatively, you can use data$grid.index to specify the SVC at a larger spatial level than the individual sites (e.g., a grid). See ?svcPGOcc for details.") 
  }
  if (!is.list(data$det.covs)) {
    stop("error: det.covs must be a list of matrices, data frames, and/or vectors")
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
    tmp <- lapply(ord, function (a) which(grid.index == a))
    tmp.2 <- sapply(tmp, length)
    grid.index.c <- unlist(lapply(1:length(tmp.2), function(a) rep(a, tmp.2[a]))) - 1
    grid.index.r <- grid.index.c + 1
    long.ord <- unlist(lapply(ord, function(a) which(grid.index == a)))
    # Reorder everything to align with NN ordering
    y <- y[long.ord, , drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Occupancy covariates
    data$occ.covs <- data$occ.covs[long.ord, , drop = FALSE]
    for (i in 1:length(data$det.covs)) {
      if (!is.null(dim(data$det.covs[[i]]))) {
        data$det.covs[[i]] <- data$det.covs[[i]][long.ord, , drop = FALSE]
      } else {
        data$det.covs[[i]] <- data$det.covs[[i]][long.ord]
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
    y.missing <- which(is.na(y))
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
  J <- nrow(y.big)
  # Number of coordinates in coordinate grid
  J.w <- nrow(coords)
  # Number of occupancy parameters 
  p.occ <- ncol(X)
  # Number of detection parameters
  p.det <- dim(X.p)[2]
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
  # Number of replicates at each site
  n.rep <- apply(y.big, 1, function(a) sum(!is.na(a)))
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- which(!is.na(y.big[j, ]))
  }
  # Max number of repeat visits
  K.max <- dim(y)[2] 
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
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
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
    phi.a <- rep(3 / max(coords.D), p.svc)
    phi.b <- rep(3 / sort(unique(c(coords.D)))[2], p.svc)
  }
  # sigma.sq --------------------------
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
    # Reorder the user supplied inits values
    if (NNGP) {
      z.inits <- z.inits[long.ord]
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
      	   p.svc, " x ", J.w, sep = ""))
    }
    if (nrow(w.inits) != p.svc | ncol(w.inits) != J.w) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   p.svc, " x ", J.w, sep = ""))
    }
    if (NNGP) {
      w.inits <- w.inits[, ord, drop = FALSE]
    }
  } else {
    w.inits <- matrix(0, p.svc, J.w)
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
    stop("error: svcPGOcc is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 
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
    storage.mode(J.w) <- "integer"

    if(verbose){
      cat("----------------------------------------\n");
      cat("Building the neighbors of neighbors list\n");
      cat("----------------------------------------\n");
    }
    
    indx <- mkUIndx(J.w, n.neighbors, nn.indx, nn.indx.lu, u.search.type)
    
    u.indx <- indx$u.indx
    u.indx.lu <- indx$u.indx.lu
    ui.indx <- indx$ui.indx
    u.indx.run.time <- indx$run.time

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(X.w) <- "double"
    consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, n.det.re, p.svc, J.w)
    storage.mode(consts) <- "integer"
    storage.mode(grid.index.c) <- "integer"
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
          if (p.det.re > 0) {
            sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
            alpha.star.inits <- rnorm(n.det.re, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
          }
          if (p.occ.re > 0) {
            sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
            beta.star.inits <- rnorm(n.occ.re, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
          }
        }
        storage.mode(chain.info) <- "integer"
        # Run the model in C    
        out.tmp[[i]] <- .Call("svcPGOccNNGP", y, X, X.w, X.p, coords, X.re, X.p.re, consts, 
        	                    K, n.occ.re.long, n.det.re.long, 
            	            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
                              beta.inits, alpha.inits, sigma.sq.psi.inits, sigma.sq.p.inits, 
        	                    beta.star.inits, alpha.star.inits, z.inits,
                              w.inits, phi.inits, sigma.sq.inits, nu.inits, z.long.indx, 
                              beta.star.indx, beta.level.indx, alpha.star.indx, 
          		    alpha.level.indx, mu.beta, mu.alpha, 
                              Sigma.beta, Sigma.alpha, phi.a, phi.b, 
                              sigma.sq.a, sigma.sq.b, nu.a, nu.b, 
          		    sigma.sq.psi.a, sigma.sq.psi.b, sigma.sq.p.a, sigma.sq.p.b, 
        	                    tuning.c, cov.model.indx,
                              n.batch, batch.length, 
                              accept.rate, n.omp.threads, verbose, n.report, 
                              samples.info, chain.info, fixed.params, sigma.sq.ig, grid.index.c)
        chain.info[1] <- chain.info[1] + 1
      }
      # Calculate R-Hat ---------------
      out$rhat <- list()
      if (n.chains > 1) {
        # as.vector removes the "Upper CI" when there is only 1 variable. 
        out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					         mcmc(t(a$beta.samples)))), 
        			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        out$rhat$alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$alpha.samples)))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
          out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					        mcmc(t(a$theta.samples)))), 
        			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
        if (p.det.re > 0) {
          out$rhat$sigma.sq.p <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$sigma.sq.p.samples)))), 
          			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        }
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
          					      mcmc(t(a$sigma.sq.psi.samples)))), 
          			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        }
      } else {
        out$rhat$beta <- rep(NA, p.occ)
        out$rhat$alpha <- rep(NA, p.det)
        out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3 * p.svc, 2 * p.svc))
        if (p.det.re > 0) {
          out$rhat$sigma.sq.p <- rep(NA, p.det.re)
        }
        if (p.occ.re > 0) {
          out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
        }
      }
      out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(do.call(rbind, 
        				lapply(out.tmp, function(a) t(a$alpha.samples))))
      colnames(out$alpha.samples) <- x.p.names
      out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
      if (cov.model != 'matern') {
        theta.names <- paste(rep(c('sigma.sq', 'phi'), each = p.svc), x.names[svc.cols], sep = '-')
      } else {
        theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = p.svc), x.names[svc.cols], sep = '-')
      } 
      colnames(out$theta.samples) <- theta.names
      # Get everything back in the original order
      out$coords <- coords[order(ord), ]
      out$z.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$z.samples))))
      out$z.samples <- mcmc(out$z.samples[, order(long.ord), drop = FALSE])
      out$X <- X[order(long.ord), , drop = FALSE]
      out$X.re <- X.re[order(long.ord), , drop = FALSE]
      out$X.w <- X.w[order(long.ord), , drop = FALSE]
      # Account for case when intercept only spatial model. 
      if (p.svc == 1) {
        tmp <- do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples)))
        tmp <- tmp[, order(ord), drop = FALSE]
        out$w.samples <- array(NA, dim = c(p.svc, J.w, n.post.samples * n.chains))
        out$w.samples[1, , ] <- t(tmp)
      } else {
        out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
          								dim = c(p.svc, J.w, n.post.samples))))
        out$w.samples <- out$w.samples[, order(ord), ]
      }
      out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
      out$psi.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$psi.samples))))
      out$psi.samples <- mcmc(out$psi.samples[, order(long.ord), drop = FALSE])
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
      out$like.samples <- mcmc(out$like.samples[, order(long.ord), drop = FALSE])
      # Get detection covariate stuff in right order. Method of doing this
      # depends on if there are observation level covariates or not. 
      if (!binom) {
        tmp <- matrix(NA, J * K.max, p.det)
        tmp[names.long, ] <- X.p
        tmp <- array(tmp, dim = c(J, K.max, p.det))
        tmp <- tmp[order(long.ord), , ]
        out$X.p <- matrix(tmp, J * K.max, p.det)
        out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
        colnames(out$X.p) <- x.p.names
        tmp <- matrix(NA, J * K.max, p.det.re)
        tmp[names.long, ] <- X.p.re
        tmp <- array(tmp, dim = c(J, K.max, p.det.re))
        tmp <- tmp[order(long.ord), , ]
        out$X.p.re <- matrix(tmp, J * K.max, p.det.re)
        out$X.p.re <- out$X.p.re[apply(out$X.p.re, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
        colnames(out$X.p.re) <- x.p.re.names
        tmp <- matrix(NA, J * K.max, n.det.re)
      } else {
        out$X.p <- X.p[order(long.ord), , drop = FALSE]
        out$X.p.re <- X.p.re[order(long.ord), , drop = FALSE]
      }
      out$y <- y.big[order(long.ord), , drop = FALSE]
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
      out$svc.cols <- svc.cols
      out$type <- "NNGP"
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
      sites.random <- sample(1:J.w)    
      sites.k.fold <- split(sites.random, sites.random %% k.fold)
      registerDoParallel(k.fold.threads)
      model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
        curr.set.small <- sort(sites.random[sites.k.fold[[i]]])
        curr.set <- which(grid.index.r %in% curr.set.small)
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
	X.w.fit <- X.w[-curr.set, , drop = FALSE]
	X.w.0 <- X.w[curr.set, , drop = FALSE]
	J.fit <- nrow(X.fit)
	J.0 <- nrow(X.0)
        coords.fit <- coords[-curr.set.small, , drop = FALSE]
        coords.0 <- coords[curr.set.small, , drop = FALSE]
	J.w.fit <- nrow(coords.fit)
        tmp.grid.fit <- grid.index.r[-curr.set]
        tmp.grid.0 <- grid.index.r[curr.set]
        # Reorder the grids to get their indices to work
        grid.fit <- vector('numeric', length = length(tmp.grid.fit))
        for (j in 1:nrow(coords.fit)) {
          grid.fit[which(tmp.grid.fit == (1:J.w)[-curr.set.small][j])] <- j
        }
        grid.0 <- vector('numeric', length = length(tmp.grid.0))
        for (j in 1:nrow(coords.0)) {
          grid.0[which(tmp.grid.0 == curr.set.small[j])] <- j
        }
	grid.fit.c <- grid.fit - 1
	K.fit <- K[-curr.set]
	K.0 <- K[curr.set]
	rep.indx.fit <- rep.indx[-curr.set]
        rep.indx.0 <- rep.indx[curr.set]
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
        if (!binom) {
          z.long.indx.fit <- rep(1:J.fit, dim(y.big.fit)[2])
          z.long.indx.fit <- z.long.indx.fit[!is.na(c(y.big.fit))]
          # Subtract 1 for indices in C
          z.long.indx.fit <- z.long.indx.fit - 1
          z.0.long.indx <- rep(1:J.0, dim(y.big.0)[2])
          z.0.long.indx <- z.0.long.indx[!is.na(c(y.big.0))]
	  # Don't subtract 1 for z.0.long.indx since its used in R only 
        } else {
          z.long.indx.fit <- 0:(J.fit - 1)
	  z.0.long.indx <- 1:J.0
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
        
        indx <- mkUIndx(J.w.fit, n.neighbors, nn.indx.fit, 
      		  nn.indx.lu.fit, u.search.type)
        
        u.indx.fit <- indx$u.indx
        u.indx.lu.fit <- indx$u.indx.lu
        ui.indx.fit <- indx$ui.indx

        storage.mode(y.fit) <- "double"
        storage.mode(z.inits.fit) <- "double"
        storage.mode(X.p.fit) <- "double"
        storage.mode(X.fit) <- "double"
	storage.mode(X.w.fit) <- "double"
        storage.mode(coords.fit) <- "double"
        storage.mode(K.fit) <- "double"
        consts.fit <- c(J.fit, n.obs.fit, p.occ, p.occ.re, n.occ.re.fit, 
                        p.det, p.det.re, n.det.re.fit, p.svc, J.w.fit)
        storage.mode(consts.fit) <- "integer"
        storage.mode(z.long.indx.fit) <- "integer"
        storage.mode(nn.indx.fit) <- "integer"
        storage.mode(nn.indx.lu.fit) <- "integer"
        storage.mode(u.indx.fit) <- "integer"
        storage.mode(u.indx.lu.fit) <- "integer"
        storage.mode(ui.indx.fit) <- "integer"
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
	storage.mode(grid.fit.c) <- 'integer'
        # Run the model in C
        out.fit <- .Call("svcPGOccNNGP", y.fit, X.fit, X.w.fit, X.p.fit, coords.fit, X.re.fit, X.p.re.fit,
			 consts.fit, K.fit, n.occ.re.long.fit, n.det.re.long.fit, 
        	         n.neighbors, nn.indx.fit, nn.indx.lu.fit, u.indx.fit, 
			 u.indx.lu.fit, ui.indx.fit, beta.inits, 
			 alpha.inits, sigma.sq.psi.inits, sigma.sq.p.inits, beta.star.inits.fit, 
			 alpha.star.inits.fit, z.inits.fit, w.inits, phi.inits, sigma.sq.inits, 
			 nu.inits, z.long.indx.fit, beta.star.indx.fit, 
			 beta.level.indx.fit, alpha.star.indx.fit, alpha.level.indx.fit, mu.beta, 
			 mu.alpha, Sigma.beta, Sigma.alpha, phi.a, phi.b, 
			 sigma.sq.a, sigma.sq.b, nu.a, nu.b, sigma.sq.psi.a, sigma.sq.psi.b, 
			 sigma.sq.p.a, sigma.sq.p.b, tuning.c, cov.model.indx, 
			 n.batch, batch.length, accept.rate, n.omp.threads.fit, verbose.fit, 
			 n.report, samples.info, chain.info, fixed.params, sigma.sq.ig, grid.fit.c)
        out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
        colnames(out.fit$beta.samples) <- x.names
        out.fit$alpha.samples <- mcmc(t(out.fit$alpha.samples))
        colnames(out.fit$alpha.samples) <- x.p.names
        out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
        # Account for case when intercept only spatial model. 
        if (p.svc == 1) {
          tmp <- t(out.fit$w.samples)
          out.fit$w.samples <- array(NA, dim = c(p.svc, J.w.fit, n.post.samples * n.chains))
          out.fit$w.samples[1, , ] <- tmp
        } else {
          out.fit$w.samples <- array(out.fit$w.samples, dim = c(p.svc, J.w.fit, n.post.samples))
        }
        out.fit$w.samples <- aperm(out.fit$w.samples, c(3, 1, 2))
        out.fit$X <- X.fit
        out.fit$y <- y.big.fit
        out.fit$X.p <- X.p.fit
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
        class(out.fit) <- "svcPGOcc"

	# Predict occurrence at new sites
	if (p.occ.re > 0) {
	  tmp <- unlist(re.level.names)
	  X.re.0 <- matrix(tmp[c(X.re.0 + 1)], nrow(X.re.0), ncol(X.re.0))
	  colnames(X.re.0) <- x.re.names
	}
        if (p.occ.re > 0) {
	  X.0 <- cbind(X.0, X.re.0)
	}
        out.pred <- predict.svcPGOcc(out.fit, X.0, coords.0, verbose = FALSE, grid.index.0 = grid.0)

	# Detection 
        # Generate detection values
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
            for (k in rep.indx.0[[j]]) {
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
    } # cross-validation
  } # NNGP or GP
  class(out) <- "svcPGOcc"
  out$run.time <- proc.time() - ptm
  return(out)
}

