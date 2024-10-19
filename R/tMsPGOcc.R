tMsPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
                     tuning, n.batch, batch.length, 
                     accept.rate = 0.43, n.omp.threads = 1, 
                     verbose = TRUE, ar1 = FALSE, n.report = 100, 
                     n.burn = round(.10 * n.batch * batch.length), 
                     n.thin = 1, n.chains = 1, ...) {

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
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  if (length(dim(data$y)) != 4) {
    stop("error: detection-nondetection data y must be a four-dimensional array with dimensions corresponding to species, sites, primary time periods, and replicates.")
  }
  y <- data$y
  sp.names <- attr(y, 'dimnames')[[1]]
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- list(int = array(1, dim = c(dim(y)[2], dim(y)[3])))
    } else {
      stop("error: occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!'det.covs' %in% names(data)) {
    if (det.formula == ~ 1) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
      }
      data$det.covs <- list(int = array(1, dim = dim(y)[-1]))
    } else {
      stop("error: det.covs must be specified in data for a detection model with covariates")
    }
  }

  # Reformat covariates ---------------------------------------------------
  # Get detection covariates in proper format
  # First subset detection covariates to only use those that are included in the analysis.
  data$det.covs <- data$det.covs[names(data$det.covs) %in% all.vars(det.formula)]
  # Null model support
  if (length(data$det.covs) == 0) {
    data$det.covs <- list(int = array(1, dim = dim(y)[-1]))
  }
  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  # Ordered by visit, year, then site. 
  data$det.covs <- data.frame(lapply(data$det.covs, function(a) unlist(c(a))))
  # Get detection covariates in site x year x replicate format
  if (nrow(data$det.covs) == dim(y)[2]) { # if only site-level covariates. 
    data$det.covs <- as.data.frame(mapply(rep, data$det.covs, dim(y)[3] * dim(y)[4]))
  } else if (nrow(data$det.covs) == dim(y)[2] * dim(y)[3]) { # if only site/year level covariates
    data$det.covs <- as.data.frame(mapply(rep, data$det.covs, dim(y)[4]))
  }
  y.big <- y
  # Get occurrence covariates in proper format
  # Subset covariates to only use those that are included in the analysis
  data$occ.covs <- data$occ.covs[names(data$occ.covs) %in% all.vars(occ.formula)]
  # Null model support
  if (length(data$occ.covs) == 0) {
    data$occ.covs <- list(int = matrix(1, nrow = dim(y)[2], ncol = dim(y)[3]))
  }
  # Ordered by year, then site within year. 
  data$occ.covs <- data.frame(lapply(data$occ.covs, function(a) unlist(c(a))))
  # Check if only site-level covariates are included
  if (nrow(data$occ.covs) == dim(y)[2]) {
    data$occ.covs <- as.data.frame(mapply(rep, data$occ.covs, dim(y)[3]))
  }

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
  for (i in 1:ncol(data$det.covs)) {
    if (sum(is.na(data$det.covs[, i])) > sum(is.na(y.big[1, , , ]))) {
      stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  if (det.formula != ~ 1) {
    # Misalignment between y and det.covs
    y.missing <- which(is.na(y[1, , , ]))
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
  # Number of sites
  J <- ncol(y)
  # Number of species 
  N <- dim(y)[1]
  # Total number of years
  n.years.max <- dim(y.big)[3]
  # Number of years for each site
  n.years <- rep(NA, J)
  # NOTE: assuming the same primary replicate history for each species. 
  for (j in 1:J) {
    n.years[j] <- sum(apply(y.big[1, j, , ], 1, function(a) sum(!is.na(a))) != 0)
  }
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
  # Number of repeat visits
  n.rep <- apply(y.big[1, , , , drop = FALSE], c(2, 3), function(a) sum(!is.na(a)))
  K.max <- dim(y.big)[4]
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

  # Get indices to map z to y -------------------------------------------
  z.site.indx <- rep(1:J, n.years.max) - 1
  z.year.indx <- rep(1:n.years.max, each = J) - 1
  z.dat.indx <- c(ifelse(K > 0, 1, 0))
  # This corresponds to the specific value in the n.years.max * J length 
  # vector of latent occurrence values, and corresponds with the z.site.indx
  # and z.year.indx
  z.long.indx <- rep(1:(J * n.years.max), K.max)
  # Order of c(y.big): visit, year within visit, site within year. 
  z.long.indx <- z.long.indx[!is.na(c(y.big[1, , , ]))]
  # Subtract 1 for indices in C
  z.long.indx <- z.long.indx - 1
  # Index that links observations to sites. 
  z.long.site.indx <- rep(rep(1:J, n.years.max), K.max)
  z.long.site.indx <- z.long.site.indx[!is.na(c(y.big[1, , , ]))]
  # Subtract 1 for indices in C
  z.long.site.indx <- z.long.site.indx - 1

  # y is order as follows: sorted by visit, year within visit, site within year, species within site
  # This will match up with z, which is sorted by year, then site within year, then species within site. 
  y <- c(y)
  # Assumes the missing data are constant across species, which seems likely, 
  # but may eventually need some updating. 
  # Removing missing observations when covariate data are available but 
  # there are missing detection-nondetection data. 
  names.long <- which(!is.na(c(y.big[1, , , ])))
  if (nrow(X.p) != length(names.long)) {
    X.p <- X.p[which(!is.na(c(y.big[1, , , ]))), , drop = FALSE]
  }
  if (p.det.re > 0) {
    if (nrow(X.p.re) != length(names.long)) {
      X.p.re <- X.p.re[which(!is.na(c(y.big[1, , , ]))), , drop = FALSE]
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

  # Independent beta parameters -----
  if ('independent.betas' %in% names(priors)) {
    if (priors$independent.betas == TRUE) {
      message("beta parameters will be estimated independently\n")
      ind.betas <- TRUE
    } else if (priors$independent.betas == FALSE) {
      ind.betas <- FALSE 
    }
  } else {
    ind.betas <- FALSE
  }
  # Independent alpha parameters -----
  if ('independent.alphas' %in% names(priors)) {
    if (priors$independent.alphas == TRUE) {
      message("alpha parameters will be estimated independently\n")
      ind.alphas <- TRUE
    } else if (priors$independent.alphas == FALSE) {
      ind.alphas <- FALSE 
    }
  } else {
    ind.alphas <- FALSE
  }

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
    if (verbose & !ind.betas) {
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
    if (verbose & !ind.alphas) {
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
    if (verbose & !ind.betas) {	    
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
    if (verbose & !ind.alphas) {	    
      message("No prior specified for tau.sq.alpha.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.alpha.a <- rep(0.1, p.det)
    tau.sq.alpha.b <- rep(0.1, p.det)
  }

  if (ar1) {
    # rho ---------------------------
    if ("rho.unif" %in% names(priors)) {
      if (!is.list(priors$rho.unif) | length(priors$rho.unif) != 2) {
        stop("error: rho.unif must be a list of length 2")
      }
      rho.a <- priors$rho.unif[[1]]
      rho.b <- priors$rho.unif[[2]]
      if (length(rho.a) != N & length(rho.a) != 1) {
        stop(paste("error: rho.unif[[1]] must be a vector of length ", 
        	   N, " or 1 with elements corresponding to rhos' lower bound for each species", sep = ""))
      }
      if (length(rho.b) != N & length(rho.b) != 1) {
        stop(paste("error: rho.unif[[2]] must be a vector of length ", 
        	   N, " or 1 with elements corresponding to rhos' upper bound for each species", sep = ""))
      }
      if (length(rho.a) != N) {
        rho.a <- rep(rho.a, N)
      }
      if (length(rho.b) != N) {
        rho.b <- rep(rho.b, N)
      }
    } else {
      if (verbose) {
      message("No prior specified for rho.unif.\nSetting uniform bounds to -1 and 1.\n")
      }
      rho.a <- rep(-1, N)
      rho.b <- rep(1, N)
    }
    # sigma.sq.t ---------------------- 
    if ("sigma.sq.t.ig" %in% names(priors)) { # inverse-gamma prior
      if (!is.list(priors$sigma.sq.t.ig) | length(priors$sigma.sq.t.ig) != 2) {
        stop("error: sigma.sq.t.ig must be a list of length 2")
      }
      sigma.sq.t.a <- priors$sigma.sq.t.ig[[1]]
      sigma.sq.t.b <- priors$sigma.sq.t.ig[[2]]
      if (length(sigma.sq.t.a) != N & length(sigma.sq.t.a) != 1) {
        stop(paste("error: sigma.sq.t.ig[[1]] must be a vector of length ", 
        	   N, " or 1 with elements corresponding to sigma.sq.ts' shape for each species", sep = ""))
      }
      if (length(sigma.sq.t.b) != N & length(sigma.sq.t.b) != 1) {
        stop(paste("error: sigma.sq.t.ig[[2]] must be a vector of length ", 
        	   N, " or 1 with elements corresponding to sigma.sq.ts' scale for each species", sep = ""))
      }
      if (length(sigma.sq.t.a) != N) {
        sigma.sq.t.a <- rep(sigma.sq.t.a, N)
      }
      if (length(sigma.sq.t.b) != N) {
        sigma.sq.t.b <- rep(sigma.sq.t.b, N)
      }
    } else {
      if (verbose) {
        message("No prior specified for sigma.sq.t.\nUsing an inverse-Gamma prior with the shape parameter to 2 and scale parameter to 0.5.\n")
      }
      sigma.sq.t.a <- rep(2, N)
      sigma.sq.t.b <- rep(0.5, N)
    }
  } else {
    rho.a <- rep(1, N)
    rho.b <- rep(1, N)
    sigma.sq.t.a <- rep(1, N)
    sigma.sq.t.b <- rep(1, N)
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
  # ORDER: an N x J x n.years array with values sorted by year, then site within year, 
  #        then species within site. 
  if ("z" %in% names(inits)) {
    z.inits <- inits$z
    if (length(dim(z.inits)) != 3) {
      stop(paste("error: initial values for z must be a 3-dimensional array with dimensions ", 
      	   N, " (species) x ", J, " (sites) x ", n.years.max, " (time periods)", sep = ""))
    }
    if (nrow(z.inits) != N | ncol(z.inits) != J | dim(z.inits)[3] != n.years.max) {
      stop(paste("error: initial values for z must be a 3-dimensional array with dimensions ", 
      	   N, " (species) x ", J, " (sites) x ", n.years.max, " (time periods)", sep = ""))
    }
    z.test <- apply(y.big, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    init.test <- sum(z.inits < z.test)
    if (init.test > 0) {
      stop("error: initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site.")
    }
  } else {
    z.inits <- apply(y.big, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    if (verbose) {
      message("z is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
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
  if (ar1) {
    # rho ------------------------
    if ("rho" %in% names(inits)) {
      rho.inits <- inits[["rho"]]
      if (length(rho.inits) != N & length(rho.inits) != 1) {
        stop(paste("error: initial values for rho must be of length ", N,  " or 1",
        	   sep = ""))
      }
      if (length(rho.inits) != N) {
        rho.inits <- rep(rho.inits, N)
      }
    } else {
      if (verbose) {
        message("rho is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      rho.inits <- runif(N, rho.a, rho.b)
    }
    # sigma.sq.t ----------------------
    if ("sigma.sq.t" %in% names(inits)) {
      sigma.sq.t.inits <- inits[["sigma.sq.t"]]
      if (length(sigma.sq.t.inits) != N & length(sigma.sq.t.inits) != 1) {
        stop(paste("error: initial values for sigma.sq.t must be of length ", N, 
      	     " or 1", sep = ""))
      }
      if (length(sigma.sq.t.inits) != N) {
        sigma.sq.t.inits <- rep(sigma.sq.t.inits, N)
      }
    } else {
      sigma.sq.t.inits <- runif(N, 0.5, 5)
      if (verbose) {
        message("sigma.sq.t is not specified in initial values.\nSetting initial values to random values between 0.5 and 5\n")
      }
    }
  } else {
    rho.inits <- rep(0, N)
    sigma.sq.t.inits <- rep(0, N)
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

  # Get tuning values ---------------------------------------------------
  rho.tuning <- rep(1, N)
  sigma.sq.t.tuning <- rep(1, N)
  if (ar1) {
    if (missing(tuning)) {
      rho.tuning <- rep(1, N)
      sigma.sq.t.tuning <- rep(1, N)
    } else {
      names(tuning) <- tolower(names(tuning))
      # rho ---------------------------
      if(!"rho" %in% names(tuning)) {
        stop("error: rho must be specified in tuning value list")
      }
      rho.tuning <- tuning$rho
      if (length(rho.tuning) == 1) {
        rho.tuning <- rep(tuning$rho, N)
      } else if (length(rho.tuning) != N) {
        stop(paste("error: rho tuning must be either a single value or a vector of length ",
        	   N, sep = ""))
      }
    }
  }
  tuning.c <- log(c(sigma.sq.t.tuning, rho.tuning))
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA
  curr.chain <- 1

  # Set storage for all variables ---------------------------------------
  storage.mode(y) <- "double"
  storage.mode(z.inits) <- "double"
  storage.mode(X.p) <- "double"
  consts <- c(N, J, n.obs, p.occ, p.occ.re, n.occ.re, 
      	p.det, p.det.re, n.det.re, n.years.max, ind.betas, ind.alphas)
  storage.mode(consts) <- "integer"
  storage.mode(beta.inits) <- "double"
  storage.mode(alpha.inits) <- "double"
  storage.mode(beta.comm.inits) <- "double"
  storage.mode(alpha.comm.inits) <- "double"
  storage.mode(tau.sq.beta.inits) <- "double"
  storage.mode(tau.sq.alpha.inits) <- "double"
  storage.mode(z.long.indx) <- "integer"
  storage.mode(z.year.indx) <- "integer"
  storage.mode(z.dat.indx) <- "integer"
  storage.mode(z.site.indx) <- "integer"
  storage.mode(mu.beta.comm) <- "double"
  storage.mode(Sigma.beta.comm) <- "double"
  storage.mode(mu.alpha.comm) <- "double"
  storage.mode(Sigma.alpha.comm) <- "double"
  storage.mode(tau.sq.beta.a) <- "double"
  storage.mode(tau.sq.beta.b) <- "double"
  storage.mode(tau.sq.alpha.a) <- "double"
  storage.mode(tau.sq.alpha.b) <- "double"
  storage.mode(tuning.c) <- "double"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(verbose) <- "integer"
  storage.mode(n.report) <- "integer"
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
  # AR1 parameters
  storage.mode(ar1) <- "integer"
  ar1.vals <- c(rho.a, rho.b, sigma.sq.t.a, sigma.sq.t.b, 
                rho.inits, sigma.sq.t.inits)
  storage.mode(ar1.vals) <- "double"

  # Fit the model -------------------------------------------------------
  out.tmp <- list()
  out <- list()
  for (i in 1:n.chains) {
    # Change initial values if i > 1
    if ((i > 1) & (!fix.inits)) {
      if (!ind.betas) {
        beta.comm.inits <- rnorm(p.occ, mu.beta.comm, sqrt(sigma.beta.comm))
        tau.sq.beta.inits <- runif(p.occ, 0.5, 10)
      }
      if (!ind.alphas) {
        alpha.comm.inits <- rnorm(p.det, mu.alpha.comm, sqrt(sigma.alpha.comm))
        tau.sq.alpha.inits <- runif(p.det, 0.5, 10)
      }
      beta.inits <- matrix(rnorm(N * p.occ, beta.comm.inits, 
                                 sqrt(tau.sq.beta.inits)), N, p.occ)
      beta.inits <- c(beta.inits)
      alpha.inits <- matrix(rnorm(N * p.det, alpha.comm.inits, 
                                  sqrt(tau.sq.alpha.inits)), N, p.det)
      alpha.inits <- c(alpha.inits)
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
      if (ar1) {
        rho.inits <- runif(N, rho.a, rho.b)
        sigma.sq.t.inits <- runif(N, 0.5, 5)
        ar1.vals <- c(rho.a, rho.b, sigma.sq.t.a, sigma.sq.t.b, 
                      rho.inits, sigma.sq.t.inits)
        storage.mode(ar1.vals) <- "double"
      }
    }

    storage.mode(chain.info) <- "integer"
    # Run the model in C
    out.tmp[[i]] <- .Call("tMsPGOcc", y, X, X.p, X.re, X.p.re, consts, 
                          n.occ.re.long, n.det.re.long, beta.inits, 
                          alpha.inits, z.inits, beta.comm.inits, 
                          alpha.comm.inits, tau.sq.beta.inits, 
                          tau.sq.alpha.inits, sigma.sq.psi.inits, sigma.sq.p.inits, 
                          beta.star.inits, alpha.star.inits, z.long.indx, z.year.indx, 
                          z.dat.indx, z.site.indx, beta.star.indx, 
                          beta.level.indx, alpha.star.indx, alpha.level.indx, mu.beta.comm, 
                          mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
                          tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
                          tau.sq.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
                          sigma.sq.p.a, sigma.sq.p.b, ar1, ar1.vals,
                          tuning.c, n.batch, batch.length, accept.rate, 
                          n.omp.threads, verbose, n.report, samples.info, chain.info)
    chain.info[1] <- chain.info[1] + 1
  }
  # Calculate R-Hat ---------------
  out$rhat <- list()
  if (n.chains > 1) {
    # as.vector removes the "Upper CI" when there is only 1 variable. 
    out$rhat$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
    					      mcmc(t(a$beta.comm.samples)))), 
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$alpha.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
    					      mcmc(t(a$alpha.comm.samples)))), 
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
    					      mcmc(t(a$tau.sq.beta.samples)))), 
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$tau.sq.alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
    					      mcmc(t(a$tau.sq.alpha.samples)))), 
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
    					         mcmc(t(a$beta.samples)))), 
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
    					      mcmc(t(a$alpha.samples)))), 
    			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    if (ar1) {
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$theta.samples)))), 
      			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
    }
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
    out$rhat$beta.comm <- rep(NA, p.occ)
    out$rhat$alpha.comm <- rep(NA, p.det)
    out$rhat$tau.sq.beta <- rep(NA, p.occ)
    out$rhat$tau.sq.alpha <- rep(NA, p.det)
    out$rhat$beta <- rep(NA, p.occ * N)
    out$rhat$alpha <- rep(NA, p.det * N)
    out$rhat$theta <- rep(NA, 2 * N)
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

  # Return things back in the original order. 
  out$z.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$z.samples, 
    								dim = c(N, J, n.years.max, n.post.samples))))
  out$z.samples <- aperm(out$z.samples, c(4, 1, 2, 3))

  out$psi.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$psi.samples, 
    								dim = c(N, J, n.years.max, n.post.samples))))
  out$psi.samples <- aperm(out$psi.samples, c(4, 1, 2, 3))
  out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples, 
    								dim = c(N, J, n.years.max, n.post.samples))))
  out$like.samples <- aperm(out$like.samples, c(4, 1, 2, 3))
  if (ar1) {
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    theta.names <- paste(rep(c('sigma.sq.t', 'rho'), each = N), sp.names, sep = '-')
    colnames(out$theta.samples) <- theta.names
    out$eta.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$eta.samples, 
      								dim = c(n.years.max, N, n.post.samples))))
    out$eta.samples <- aperm(out$eta.samples, c(3, 2, 1))
  }
  out$X.p <- X.p
  out$X.p.re <- X.p.re
  # Calculate effective sample sizes
  out$ESS <- list()
  out$ESS$beta.comm <- effectiveSize(out$beta.comm.samples)
  out$ESS$alpha.comm <- effectiveSize(out$alpha.comm.samples)
  out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
  out$ESS$tau.sq.alpha <- effectiveSize(out$tau.sq.alpha.samples)
  out$ESS$beta <- effectiveSize(out$beta.samples)
  out$ESS$alpha <- effectiveSize(out$alpha.samples)
  if (ar1) {
    out$ESS$theta <- effectiveSize(out$theta.samples)
  }
  if (p.det.re > 0) {
    out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
  }
  if (p.occ.re > 0) {
    out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
  }
  out$X <- array(X, dim = c(J, n.years.max, p.occ))
  dimnames(out$X)[[3]] <- x.names
  out$X.re <- array(X.re, dim = c(J, n.years.max, p.occ.re))
  dimnames(out$X.re)[[3]] <- x.re.names
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
  out$ar1 <- ar1
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
  
  
  class(out) <- "tMsPGOcc"

  out$run.time <- proc.time() - ptm
  return(out)
}
