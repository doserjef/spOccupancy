spMsPGOcc <- function(occ.formula, det.formula, data, starting, priors, 
		      tuning, cov.model = 'exponential', NNGP = TRUE, 
		      n.neighbors = 15, search.type = "cb", n.batch, 
		      batch.length, accept.rate = 0.43,
		      n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		      n.burn = round(.10 * n.batch * batch.length), 
		      n.thin = 1, k.fold, k.fold.threads = 1, 
		      k.fold.seed = 100, ...){

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
  if (length(dim(data$y)) != 3) {
    stop("error: detection-nondetection data y must be a three-dimensional array with dimensions corresponding to species, sites, and replicates.")
  }
  y <- data$y
  sp.names <- attr(y, 'dimnames')[[1]]
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("occupancy covariates (occ.covs) not specified in data. Assuming intercept only occupancy model.")
      }
      data$occ.covs <- matrix(1, dim(y)[2], 1)
    } else {
      stop("error: occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!'det.covs' %in% names(data)) {
    if (det.formula == ~ 1) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data. Assuming interept only detection model.")
      }
      data$det.covs <- list(int = matrix(1, dim(y)[2], dim(y)[3]))
    } else {
      stop("error: det.covs must be specified in data for a detection model with covariates")
    }
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  coords <- data$coords

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[, ord, , drop = FALSE]
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

  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  data$det.covs <- data.frame(lapply(data$det.covs, function(a) unlist(c(a))))
  # Replicate det.covs if only covariates are at the site level. 
  if (nrow(data$det.covs) == nrow(y)) {
    data$det.covs <- data.frame(sapply(data$det.covs, rep, times = dim(y)[2]))
  }
  data$occ.covs <- as.data.frame(data$occ.covs)

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
  # Number of detection parameters
  p.det <- ncol(X.p)
  # Number of detection random effect parameters
  p.det.re <- ncol(X.p.re)
  # Number of latent detection random effect values
  n.det.re <- length(unlist(apply(X.p.re, 2, unique)))
  n.det.re.long <- apply(X.p.re, 2, function(a) length(unique(a)))
  # Number of pseudoreplicates
  n.obs <- nrow(X.p)
  # Number of sites
  J <- nrow(X)
  # Number of repeat visits
  n.rep <- apply(y[1, , ], 1, function(a) sum(!is.na(a)))
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

  # Get indices to map z to y -------------------------------------------
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
  # Subtract 1 for indices in C
  z.long.indx <- z.long.indx - 1
  # y is stored in the following order: species, site, visit
  y.big <- y
  y <- c(y)
  # Assumes the missing data are constant across species, which seems likely, 
  # but may eventually need some updating. 
  # Removing missing observations when covariate data are available but 
  # there are missing detection-nondetection data. 
  names.long <- which(!is.na(c(y.big[1, , ])))
  if (nrow(X.p) == length(y) / N) {
    X.p <- X.p[!is.na(c(y.big[1, , ])), ]
    X.p.re <- X.p.re[!is.na(y), , drop = FALSE]
  }
  y <- y[!is.na(y)]

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

  # Starting values -----------------------------------------------------
  # Function to get starting values for species-level regression 
  # coefficients. 
  tmp.f <- function(a, x) {
    coefficients(glm((a) ~ x - 1, family = 'binomial'))
  }
  if (missing(starting)) {
    stop("error: starting value list for the parameters must be specified.")
  }
  names(starting) <- tolower(names(starting))
  # z -------------------------------
  if ("z" %in% names(starting)) {
    z.starting <- starting$z
    if (!is.matrix(z.starting)) {
      stop(paste("error: starting values for z must be a matrix with dimensions ", 
      	   N, " x ", J, sep = ""))
    }
    if (nrow(z.starting) != N | ncol(z.starting) != J) {
      stop(paste("error: starting values for z must be a matrix with dimensions ", 
      	   N, " x ", J, sep = ""))
    }
    # Reorder the user supplied starting values for NNGP models
    if (NNGP) {
      z.starting <- z.starting[, ord]
    }
  } else {
    # In correct order since you reordered y for NNGP. 
    z.starting <- apply(y.big, c(1, 2), max, na.rm = TRUE)
  }
  # beta ----------------------------
  if ("beta" %in% names(starting)) {
    beta.starting <- starting[["beta"]]
    if (!is.matrix(beta.starting)) {
      stop(paste("error: starting values for beta must be a matrix with dimensions ", 
      	   N, " x ", p.occ, sep = ""))
    }
    if (ncol(beta.starting) != p.occ | nrow(beta.starting) != N) {
      stop(paste("error: starting values for beta must be a matrix with dimensions ", 
      	   N, "x", p.occ, sep = ""))
    }
  } else {
    beta.starting <- t(apply(z.starting, 1, tmp.f, X))
  }
  # beta.comm -----------------------
  if ("beta.comm" %in% names(starting)) {
    beta.comm.starting <- starting[["beta.comm"]]
    if (length(beta.comm.starting) != p.occ) {
      stop(paste("error: starting values for beta.comm must be of length ", p.occ, 
      	   sep = ""))
    }
  } else {
    beta.comm.starting <- apply(beta.starting, 2, mean)
  }
  # alpha ----------------------------
  if ("alpha" %in% names(starting)) {
    alpha.starting <- starting[["alpha"]]
    if (ncol(alpha.starting) != p.det | nrow(alpha.starting) != N) {
      stop(paste("error: starting values for alpha must be a matrix with dimensions ", 
      	   N, "x", p.det, sep = ""))
    }
  } else {
    alpha.starting <- matrix(0, N, p.det)
  }
  # alpha.comm -----------------------
  if ("alpha.comm" %in% names(starting)) {
    alpha.comm.starting <- starting[["alpha.comm"]]
    if (length(alpha.comm.starting) != p.det) {
      stop(paste("error: starting values for alpha.comm must be of length ", p.det, 
      	   sep = ""))
    }
  } else {
    alpha.comm.starting <- apply(alpha.starting, 2, mean)
  }
  # tau.sq.beta ------------------------
  if ("tau.sq.beta" %in% names(starting)) {
    tau.sq.beta.starting <- starting[["tau.sq.beta"]]
    if (length(tau.sq.beta.starting) != p.occ) {
      stop(paste("error: starting values for tau.sq.beta must be of length ", p.occ, 
      	   sep = ""))
    }
  } else {
    tau.sq.beta.starting <- runif(p.occ, 0.1, 2)
  }
  # tau.sq.alpha -----------------------
  if ("tau.sq.alpha" %in% names(starting)) {
    tau.sq.alpha.starting <- starting[["tau.sq.alpha"]]
    if (length(tau.sq.alpha.starting) != p.det) {
      stop(paste("error: starting values for tau.sq.alpha must be of length ", p.det, 
      	   sep = ""))
    }
  } else {
    tau.sq.alpha.starting <- runif(p.det, 0.1, 2)
  }
  # phi -----------------------------
  if ("phi" %in% names(starting)) {
    phi.starting <- starting[["phi"]]
    if (length(phi.starting) != N) {
      stop(paste("error: starting values for phi must be of length ", N, 
      	   sep = ""))
    }
  } else {
    phi.starting <- rep(3/mean(range(coords)), N)
    if (verbose) {
      message("phi is not specified in starting values. Setting starting value to 3/mean(range(coords))\n")
    }
  }
  # sigma.sq ------------------------
  if ("sigma.sq" %in% names(starting)) {
    sigma.sq.starting <- starting[["sigma.sq"]]
    if (length(sigma.sq.starting) != N) {
      stop(paste("error: starting values for sigma.sq must be of length ", N, 
      	   sep = ""))
    }
  } else {
    sigma.sq.starting <- rep(2, N)
    if (verbose) {
      message("sigma.sq is not specified in starting values. Setting starting value to 2\n")
    }
  }
  # w -----------------------------00
  if ("w" %in% names(starting)) {
    w.starting <- starting[["w"]]
    if (!is.matrix(w.starting)) {
      stop(paste("error: starting values for w must be a matrix with dimensions ",
      	   N, " x ", J, sep = ""))
    }
    if (nrow(w.starting) != N | ncol(w.starting) != J) {
      stop(paste("error: starting values for w must be a matrix with dimensions ",
      	   N, " x ", J, sep = ""))
    }
  } else {
    w.starting <- matrix(0, N, J)
    if (verbose) {
      message("w is not specified in starting values. Setting starting value to 0\n")
    }
  }
  # nu ------------------------
  if ("nu" %in% names(starting)) {
    nu.starting <- starting[["nu"]]
    if (length(nu.starting) != N) {
      stop(paste("error: starting values for nu must be of length ", N, 
      	   sep = ""))
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in starting values. Setting starting value to 1\n")
      }
      nu.starting <- rep(1, N)
    } else {
      nu.starting <- rep(0, N)
    }
  }

  # sigma.sq.p ------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p" %in% names(starting)) {
      sigma.sq.p.starting <- starting[["sigma.sq.p"]]
      if (length(sigma.sq.p.starting) != p.det.re) {
        stop(paste("error: starting values for sigma.sq.p must be of length ", p.det.re, 
      	     sep = ""))
      }
    } else {
      sigma.sq.p.starting <- rep(1, p.det.re)
      if (verbose) {
        message("sigma.sq.p is not specified in starting values. Setting starting value to 1\n")
      }
    }
    alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
    alpha.star.starting <- rnorm(n.det.re, sqrt(sigma.sq.p.starting[alpha.star.indx + 1]))
    alpha.star.starting <- rep(alpha.star.starting, N)
  }

  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    stop("error: prior list for the parameters must be specified")
  }
  names(priors) <- tolower(names(priors))

  # beta.comm -----------------------
  if ("beta.comm.normal" %in% names(priors)) {
    mu.beta.comm <- priors$beta.comm.normal[[1]]
    sigma.beta.comm <- priors$beta.comm.normal[[2]]
    if (!is.list(priors$beta.comm.normal) | length(priors$beta.comm.normal) != 2) {
      stop("error: beta.comm.normal must be a list of length 2")
    }
    if (length(mu.beta.comm) != p.occ) {
      stop(paste("error: beta.comm.normal[[1]] must be a vector of length ", 
      	   p.occ, " with elements corresponding to beta.comms' mean", sep = ""))
    }
    if (length(sigma.beta.comm) != p.occ) {
      stop(paste("error: beta.comm.normal[[2]] must be a vector of length ", 
      	   p.occ, " with elements corresponding to beta.comms' variance", sep = ""))
    }
    Sigma.beta.comm <- sigma.beta.comm * diag(p.occ)
  } else {
    if (verbose) {
      message("No prior specified for beta.comm.normal. Setting prior mean to 0 and prior variance to 2.73\n")
    }
    mu.beta.comm <- rep(0, p.occ)
    Sigma.beta.comm <- diag(p.occ) * 2.73
  }

  # alpha.comm -----------------------
  if ("alpha.comm.normal" %in% names(priors)) {
    mu.alpha.comm <- priors$alpha.comm.normal[[1]]
    sigma.alpha.comm <- priors$alpha.comm.normal[[2]]
    if (!is.list(priors$alpha.comm.normal) | length(priors$alpha.comm.normal) != 2) {
      stop("error: alpha.comm.normal must be a list of length 2")
    }
    if (length(mu.alpha.comm) != p.det) {
      stop(paste("error: alpha.comm.normal[[1]] must be a vector of length ", 
      	   p.det, " with elements corresponding to alpha.comms' mean", sep = ""))
    }
    if (length(sigma.alpha.comm) != p.det) {
      stop(paste("error: alpha.comm.normal[[2]] must be a vector of length ", 
      	   p.det, " with elements corresponding to alphas.comms' variance", sep = ""))
    }
    Sigma.alpha.comm <- sigma.alpha.comm * diag(p.det)
  } else {
    if (verbose) {
      message("No prior specified for alpha.comm.normal. Setting prior mean to 0 and prior variance to 2.73\n")
    }
    mu.alpha.comm <- rep(0, p.det)
    Sigma.alpha.comm <- diag(p.det) * 2.73
  }

  # tau.sq.beta -----------------------
  if ("tau.sq.beta.ig" %in% names(priors)) {
    tau.sq.beta.a <- priors$tau.sq.beta.ig[[1]]
    tau.sq.beta.b <- priors$tau.sq.beta.ig[[2]]
    if (!is.list(priors$tau.sq.beta.ig) | length(priors$tau.sq.beta.ig) != 2) {
      stop("error: tau.sq.beta.ig must be a list of length 2")
    }
    if (length(tau.sq.beta.a) != p.occ) {
      stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ", 
      	   p.occ, " with elements corresponding to tau.sq.betas' shape", sep = ""))
    }
    if (length(tau.sq.beta.b) != p.occ) {
      stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ", 
      	   p.occ, " with elements corresponding to tau.sq.betas' scale", sep = ""))
    }
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.beta.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.beta.a <- rep(0.1, p.occ)
    tau.sq.beta.b <- rep(0.1, p.occ)
  }

  # tau.sq.alpha -----------------------
  if ("tau.sq.alpha.ig" %in% names(priors)) {
    tau.sq.alpha.a <- priors$tau.sq.alpha.ig[[1]]
    tau.sq.alpha.b <- priors$tau.sq.alpha.ig[[2]]
    if (!is.list(priors$tau.sq.alpha.ig) | length(priors$tau.sq.alpha.ig) != 2) {
      stop("error: tau.sq.alpha.ig must be a list of length 2")
    }
    if (length(tau.sq.alpha.a) != p.det) {
      stop(paste("error: tau.sq.alpha.ig[[1]] must be a vector of length ", 
      	   p.det, " with elements corresponding to tau.sq.alphas' shape", sep = ""))
    }
    if (length(tau.sq.alpha.b) != p.det) {
      stop(paste("error: tau.sq.alpha.ig[[2]] must be a vector of length ", 
      	   p.det, " with elements corresponding to tau.sq.alphas' scale", sep = ""))
    }
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.alpha.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.alpha.a <- rep(0.1, p.det)
    tau.sq.alpha.b <- rep(0.1, p.det)
  }

  # phi -----------------------------
  if (!"phi.unif" %in% names(priors)) {
    stop("error: phi.unif must be specified in priors value list")
  }
  phi.a <- priors$phi.unif[[1]]
  phi.b <- priors$phi.unif[[2]]
  if (!is.list(priors$phi.unif) | length(priors$phi.unif) != 2) {
    stop("error: phi.unif must be a list of length 2")
  }
  if (length(phi.a) != N) {
    stop(paste("error: phi.unif[[1]] must be a vector of length ", 
    	   N, " with elements corresponding to phis' lower bound for each species", sep = ""))
  }
  if (length(phi.b) != N) {
    stop(paste("error: phi.unif[[2]] must be a vector of length ", 
    	   N, " with elements corresponding to phis' upper bound for each species", sep = ""))
  }

  # sigma.sq -----------------------------
  if (!"sigma.sq.ig" %in% names(priors)) {
    stop("error: sigma.sq.ig must be specified in priors value list")
  }
  sigma.sq.a <- priors$sigma.sq.ig[[1]]
  sigma.sq.b <- priors$sigma.sq.ig[[2]]
  if (!is.list(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
    stop("error: sigma.sq.ig must be a list of length 2")
  }
  if (length(sigma.sq.a) != N) {
    stop(paste("error: sigma.sq.ig[[1]] must be a vector of length ", 
    	   N, " with elements corresponding to sigma.sqs' shape for each species", sep = ""))
  }
  if (length(sigma.sq.b) != N) {
    stop(paste("error: sigma.sq.ig[[2]] must be a vector of length ", 
    	   N, " with elements corresponding to sigma.sqs' scale for each species", sep = ""))
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
    if (length(nu.a) != N) {
      stop(paste("error: nu.unif[[1]] must be a vector of length ", 
      	   N, " with elements corresponding to nus' lower bound for each species", sep = ""))
    }
    if (length(nu.b) != N) {
      stop(paste("error: nu.unif[[2]] must be a vector of length ", 
      	   N, " with elements corresponding to nus' upper bound for each species", sep = ""))
    }
  } else {
    nu.a <- rep(0, N)
    nu.b <- rep(0, N)
  }

  # sigma.sq.p --------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p.ig" %in% names(priors)) {
      sigma.sq.p.a <- priors$sigma.sq.p.ig[[1]]
      sigma.sq.p.b <- priors$sigma.sq.p.ig[[2]]
      if (!is.list(priors$sigma.sq.p.ig) | length(priors$sigma.sq.p.ig) != 2) {
        stop("error: sigma.sq.p.ig must be a list of length 2")
      }
      if (length(sigma.sq.p.a) != p.det.re) {
        stop(paste("error: sigma.sq.p.ig[[1]] must be a vector of length ", 
        	   p.det.re, " with elements corresponding to sigma.sq.ps' shape", sep = ""))
      }
      if (length(sigma.sq.p.b) != p.det.re) {
        stop(paste("error: sigma.sq.p.ig[[2]] must be a vector of length ", 
        	   p.det.re, " with elements corresponding to sigma.sq.ps' scale", sep = ""))
      }
  }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.p.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.p.a <- rep(0.1, p.det.re)
      sigma.sq.p.b <- rep(0.1, p.det.re)
    }
  }

  # Covariance Model ----------------------------------------------------
  if (missing(cov.model)) {
    stop("error: cov.model must be specified")
  }
  # Order must match util.cpp spCor.
  cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
  if(! cov.model %in% cov.model.names){
    stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ", 
         paste(cov.model.names, collapse=", ", sep="") ,".")}
  # Obo for cov model lookup on c side
  cov.model.indx <- which(cov.model == cov.model.names) - 1

  # Get tuning values ---------------------------------------------------
  # Not accessed, but necessary to keep things in line. 
  sigma.sq.tuning <- rep(0, N)
  phi.tuning <- rep(0, N)
  nu.tuning <- rep(0, N)
  if (missing(tuning)) {
    phi.tuning <- rep(1, N)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, N)
    }
  } else {
    names(tuning) <- tolower(names(tuning))
    # phi ---------------------------
    if(!"phi" %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (length(phi.tuning) == 1) {
      phi.tuning <- rep(tuning$phi, N)
    } else if (length(phi.tuning) != N) {
      stop(paste("error: phi tuning must be either a single value or a vector of length ",
      	   N, sep = ""))
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) == 1) {
        nu.tuning <- rep(tuning$nu, N)
      } else if (length(nu.tuning) != N) {
        stop(paste("error: nu tuning must be either a single value or a vector of length ",
        	   N, sep = ""))
      }
    }
  }
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA

  if (!NNGP) {
    # Get distance matrix for full GP -------------------------------------
    coords.D <- iDist(coords)

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.starting) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(coords.D) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(N) <- "integer"
    storage.mode(beta.starting) <- "double"
    storage.mode(alpha.starting) <- "double"
    storage.mode(beta.comm.starting) <- "double"
    storage.mode(alpha.comm.starting) <- "double"
    storage.mode(tau.sq.beta.starting) <- "double"
    storage.mode(tau.sq.alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(w.starting) <- "double"
    storage.mode(nu.starting) <- "double"
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
    storage.mode(sigma.sq.a) <- "double"
    storage.mode(sigma.sq.b) <- "double"
    storage.mode(tuning.c) <- "double"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
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
      storage.mode(sigma.sq.p.starting) <- "double"
      storage.mode(sigma.sq.p.a) <- "double"
      storage.mode(sigma.sq.p.b) <- "double"
      storage.mode(alpha.star.starting) <- "double"
      storage.mode(alpha.star.indx) <- "integer"
      storage.mode(lambda.p) <- "double"

      out <- .Call("spMsPGOccRE", y, X, X.p, coords.D, X.p.re, lambda.p, 
		   p.occ, p.det, p.det.re, J, K, N, n.det.re, n.det.re.long,
          	   beta.starting, alpha.starting, z.starting,
          	   beta.comm.starting, 
          	   alpha.comm.starting, tau.sq.beta.starting, 
          	   tau.sq.alpha.starting, w.starting, phi.starting, 
          	   sigma.sq.starting, nu.starting, sigma.sq.p.starting, 
		   alpha.star.starting, z.long.indx, alpha.star.indx, mu.beta.comm, 
          	   mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	   tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	   tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	   nu.a, nu.b, sigma.sq.p.a, sigma.sq.p.b, tuning.c, cov.model.indx, 
          	   n.batch, batch.length, accept.rate, 
          	   n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	   n.post.samples)
      
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
      coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- coef.names.det
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
      } else {
        theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
      } 
      colnames(out$theta.samples) <- theta.names
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$w.samples <- array(out$w.samples, dim = c(N, J, n.post.samples))
      out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      tmp <- array(NA, dim = c(N, J * K.max, n.post.samples))
      tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.post.samples))
      out$y.rep.samples <- array(tmp, dim = c(N, J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
      out$X <- X
      out$X.p <- X.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$theta.names <- theta.names
      out$type <- "GP"
      out$coords <- coords
      out$cov.model.indx <- cov.model.indx
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
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
        # Currently implemented without parellization. 
	set.seed(k.fold.seed)
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
	  y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	  y.fit <- y.fit[!is.na(y.fit)]
	  y.0 <- c(y.big[, curr.set, , drop = FALSE])
	  y.0 <- y.0[!is.na(y.0)]
	  z.starting.fit <- z.starting[-curr.set]
	  w.starting.fit <- w.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  coords.D.fit <- coords.D[-curr.set, -curr.set, drop = FALSE]
	  coords.D.0 <- coords.D[curr.set, curr.set, drop = FALSE]
	  J.fit <- nrow(X.fit)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  lambda.p.fit <- lambda.p[y.indx, , drop = FALSE]
	  lambda.p.0 <- lambda.p[!y.indx, , drop = FALSE]
	  X.p.re.fit <- X.p.re[y.indx, , drop = FALSE]
	  X.p.re.0 <- X.p.re[!y.indx, , drop = FALSE]
	  # Gotta be a better way, but will do for now. 
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
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1

          storage.mode(y.fit) <- "double"
          storage.mode(z.starting.fit) <- "double"
          storage.mode(X.p.fit) <- "double"
          storage.mode(X.fit) <- "double"
          storage.mode(coords.D.fit) <- "double"
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(beta.comm.starting) <- "double"
          storage.mode(alpha.comm.starting) <- "double"
          storage.mode(tau.sq.beta.starting) <- "double"
          storage.mode(tau.sq.alpha.starting) <- "double"
          storage.mode(phi.starting) <- "double"
          storage.mode(sigma.sq.starting) <- "double"
          storage.mode(w.starting.fit) <- "double"
          storage.mode(nu.starting) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
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
          storage.mode(sigma.sq.a) <- "double"
          storage.mode(sigma.sq.b) <- "double"
          storage.mode(tuning.c) <- "double"
          storage.mode(n.batch) <- "integer"
          storage.mode(batch.length) <- "integer"
          storage.mode(accept.rate) <- "double"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
          storage.mode(cov.model.indx) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"
          storage.mode(p.det.re) <- "integer"
          storage.mode(X.p.re.fit) <- "integer"
          storage.mode(n.det.re) <- "integer"
          storage.mode(n.det.re.long) <- "integer"
          storage.mode(sigma.sq.p.starting) <- "double"
          storage.mode(sigma.sq.p.a) <- "double"
          storage.mode(sigma.sq.p.b) <- "double"
          storage.mode(alpha.star.starting) <- "double"
          storage.mode(alpha.star.indx) <- "integer"
          storage.mode(lambda.p.fit) <- "double"


          out.fit <- .Call("spMsPGOccRE", y.fit, X.fit, X.p.fit, coords.D.fit, X.p.re.fit, lambda.p.fit, 
		           p.occ, p.det, p.det.re, J.fit, K.fit, N, n.det.re, n.det.re.long,
          	           beta.starting, alpha.starting, z.starting.fit,
          	           beta.comm.starting, 
          	           alpha.comm.starting, tau.sq.beta.starting, 
          	           tau.sq.alpha.starting, w.starting.fit, phi.starting, 
          	           sigma.sq.starting, nu.starting, sigma.sq.p.starting, 
		           alpha.star.starting, z.long.indx.fit, alpha.star.indx, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	           nu.a, nu.b, sigma.sq.p.a, sigma.sq.p.b, tuning.c, cov.model.indx, 
          	           n.batch, batch.length, accept.rate, 
          	           n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
          	           n.post.samples)

          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
          } else {
            theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
          } 
          colnames(out.fit$theta.samples) <- theta.names
          out.fit$w.samples <- array(out.fit$w.samples, dim = c(N, J, n.post.samples))
          out.fit$w.samples <- aperm(out.fit$w.samples, c(3, 1, 2))
          out.fit$X <- X.fit
	  out.fit$y <- y.big
          out.fit$X.p <- X.p.fit
          out.fit$call <- cl
          out.fit$n.samples <- n.samples
          out.fit$type <- "GP"
          out.fit$coords <- coords.fit
          out.fit$cov.model.indx <- cov.model.indx
          out.fit$n.post <- n.post.samples
          out.fit$n.thin <- n.thin
          out.fit$n.burn <- n.burn
          out.fit$pRE <- TRUE
	  class(out.fit) <- "spMsPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.spMsPGOcc(out.fit, X.0, coords.0, verbose = FALSE)

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  like.samples <- matrix(NA, N, nrow(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
          sp.re.indx <- rep(1:N, each = nrow(out.fit$alpha.star.samples) / N)
	  for (q in 1:N) {
            p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] + 
					    lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
	    for (j in 1:nrow(X.p.0)) {
              like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
						p.0.samples[j, q, ] * 
					        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a)))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "spMsPGOcc"

    }

    if (p.det.re == 0) {

      # Run the model in C    
      out <- .Call("spMsPGOcc", y, X, X.p, coords.D, p.occ, p.det, J, K, N, 
          	 beta.starting, alpha.starting, z.starting,
          	 beta.comm.starting, 
          	 alpha.comm.starting, tau.sq.beta.starting, 
          	 tau.sq.alpha.starting, w.starting, phi.starting, 
          	 sigma.sq.starting, nu.starting, z.long.indx, mu.beta.comm, 
          	 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	 tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	 tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	 nu.a, nu.b, tuning.c, cov.model.indx, 
          	 n.batch, batch.length, accept.rate, 
          	 n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	 n.post.samples)

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
      coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- coef.names.det
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
      } else {
        theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
      } 
      colnames(out$theta.samples) <- theta.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$w.samples <- array(out$w.samples, dim = c(N, J, n.post.samples))
      out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      tmp <- array(NA, dim = c(N, J * K.max, n.post.samples))
      tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.post.samples))
      out$y.rep.samples <- array(tmp, dim = c(N, J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
      out$X <- X
      out$X.p <- X.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$theta.names <- theta.names
      out$type <- "GP"
      out$coords <- coords
      out$cov.model.indx <- cov.model.indx
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
        # Currently implemented without parellization. 
	set.seed(k.fold.seed)
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
	  y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	  y.fit <- y.fit[!is.na(y.fit)]
	  y.0 <- c(y.big[, curr.set, , drop = FALSE])
	  y.0 <- y.0[!is.na(y.0)]
	  z.starting.fit <- z.starting[-curr.set]
	  w.starting.fit <- w.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  coords.D.fit <- coords.D[-curr.set, -curr.set, drop = FALSE]
	  coords.D.0 <- coords.D[curr.set, curr.set, drop = FALSE]
	  J.fit <- nrow(X.fit)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  # Gotta be a better way, but will do for now. 
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
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1

          storage.mode(y.fit) <- "double"
          storage.mode(z.starting.fit) <- "double"
          storage.mode(X.p.fit) <- "double"
          storage.mode(X.fit) <- "double"
          storage.mode(coords.D.fit) <- "double"
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(beta.comm.starting) <- "double"
          storage.mode(alpha.comm.starting) <- "double"
          storage.mode(tau.sq.beta.starting) <- "double"
          storage.mode(tau.sq.alpha.starting) <- "double"
          storage.mode(phi.starting) <- "double"
          storage.mode(sigma.sq.starting) <- "double"
          storage.mode(w.starting.fit) <- "double"
          storage.mode(nu.starting) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
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
          storage.mode(sigma.sq.a) <- "double"
          storage.mode(sigma.sq.b) <- "double"
          storage.mode(tuning.c) <- "double"
          storage.mode(n.batch) <- "integer"
          storage.mode(batch.length) <- "integer"
          storage.mode(accept.rate) <- "double"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
          storage.mode(cov.model.indx) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"

          out.fit <- .Call("spMsPGOcc", y.fit, X.fit, X.p.fit, coords.D.fit, 
		           p.occ, p.det, J.fit, K.fit, N,
          	           beta.starting, alpha.starting, z.starting.fit,
          	           beta.comm.starting, 
          	           alpha.comm.starting, tau.sq.beta.starting, 
          	           tau.sq.alpha.starting, w.starting.fit, phi.starting, 
          	           sigma.sq.starting, nu.starting,
		           z.long.indx.fit, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	           nu.a, nu.b, tuning.c, cov.model.indx, 
          	           n.batch, batch.length, accept.rate, 
          	           n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
          	           n.post.samples)

          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
          } else {
            theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
          } 
          colnames(out.fit$theta.samples) <- theta.names
          out.fit$w.samples <- array(out.fit$w.samples, dim = c(N, J, n.post.samples))
          out.fit$w.samples <- aperm(out.fit$w.samples, c(3, 1, 2))
          out.fit$X <- X.fit
	  out.fit$y <- y.big
          out.fit$X.p <- X.p.fit
          out.fit$call <- cl
          out.fit$n.samples <- n.samples
          out.fit$type <- "GP"
          out.fit$coords <- coords.fit
          out.fit$cov.model.indx <- cov.model.indx
          out.fit$n.post <- n.post.samples
          out.fit$n.thin <- n.thin
          out.fit$n.burn <- n.burn
          out.fit$pRE <- FALSE
	  class(out.fit) <- "spMsPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.spMsPGOcc(out.fit, X.0, coords.0, verbose = FALSE)

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  like.samples <- matrix(NA, N, nrow(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
	  for (q in 1:N) {
            p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
	    for (j in 1:nrow(X.p.0)) {
              like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
						p.0.samples[j, q, ] * 
					        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a)))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "spMsPGOcc"

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
    
    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.starting) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(N) <- "integer"
    storage.mode(beta.starting) <- "double"
    storage.mode(alpha.starting) <- "double"
    storage.mode(beta.comm.starting) <- "double"
    storage.mode(alpha.comm.starting) <- "double"
    storage.mode(tau.sq.beta.starting) <- "double"
    storage.mode(tau.sq.alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    storage.mode(w.starting) <- "double"
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
    storage.mode(cov.model.indx) <- "integer"
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
      storage.mode(sigma.sq.p.starting) <- "double"
      storage.mode(sigma.sq.p.a) <- "double"
      storage.mode(sigma.sq.p.b) <- "double"
      storage.mode(alpha.star.starting) <- "double"
      storage.mode(alpha.star.indx) <- "integer"
      storage.mode(lambda.p) <- "double"

      # Run the model in C    
      out <- .Call("spMsPGOccNNGPRE", y, X, X.p, coords, X.p.re, lambda.p, 
		   p.occ, p.det, p.det.re, J, K, N, n.det.re, n.det.re.long, 
          	 n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
          	 beta.starting, alpha.starting, z.starting,
          	 beta.comm.starting, 
          	 alpha.comm.starting, tau.sq.beta.starting, 
          	 tau.sq.alpha.starting, w.starting, phi.starting, 
          	 sigma.sq.starting, nu.starting, sigma.sq.p.starting, 
		 alpha.star.starting, z.long.indx, alpha.star.indx, mu.beta.comm, 
          	 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	 tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	 tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	 nu.a, nu.b, sigma.sq.p.a, sigma.sq.p.b, tuning.c, cov.model.indx, n.batch, 
          	 batch.length, accept.rate, n.omp.threads, verbose, n.report, 
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
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
      } else {
        theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
      } 
      colnames(out$theta.samples) <- theta.names
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      # Return things back in the original order. 
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- out$z.samples[, order(ord), ]
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$w.samples <- array(out$w.samples, dim = c(N, J, n.post.samples))
      out$w.samples <- out$w.samples[, order(ord), ]
      out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- out$psi.samples[, order(ord), ]
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      tmp <- matrix(NA, J * K.max, p.det)
      tmp[names.long, ] <- X.p
      tmp <- array(tmp, dim = c(J, K.max, p.det))
      tmp <- tmp[order(ord), , ]
      out$X.p <- matrix(tmp, J * K.max, p.det)
      out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, ]
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$X <- X[order(ord), , drop = FALSE]
      out$y <- y.big[, order(ord), , drop = FALSE]
      out$call <- cl
      out$n.samples <- n.samples
      tmp <- array(NA, dim = c(N, J * K.max, n.post.samples))
      tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.post.samples))
      tmp <- array(tmp, dim = c(N, J, K.max, n.post.samples))
      out$y.rep.samples <- tmp[, order(ord), , ]
      out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$theta.names <- theta.names
      out$type <- "NNGP"
      out$coords <- coords[order(ord), ]
      out$cov.model.indx <- cov.model.indx
      out$n.neighbors <- n.neighbors
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
        # Currently implemented without parellization. 
	set.seed(k.fold.seed)
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
	  y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	  y.fit <- y.fit[!is.na(y.fit)]
	  y.0 <- c(y.big[, curr.set, , drop = FALSE])
	  y.0 <- y.0[!is.na(y.0)]
	  z.starting.fit <- z.starting[-curr.set]
	  w.starting.fit <- w.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  J.fit <- nrow(X.fit)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  lambda.p.fit <- lambda.p[y.indx, , drop = FALSE]
	  lambda.p.0 <- lambda.p[!y.indx, , drop = FALSE]
	  X.p.re.fit <- X.p.re[y.indx, , drop = FALSE]
	  X.p.re.0 <- X.p.re[!y.indx, , drop = FALSE]
	  # Gotta be a better way, but will do for now. 
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
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1
	  # Don't need to reorder things, since they are already sorted by 
	  # the horizontal location in the coordinates. 

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
          storage.mode(z.starting.fit) <- "double"
          storage.mode(X.p.fit) <- "double"
          storage.mode(X.fit) <- "double"
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(beta.comm.starting) <- "double"
          storage.mode(alpha.comm.starting) <- "double"
          storage.mode(tau.sq.beta.starting) <- "double"
          storage.mode(tau.sq.alpha.starting) <- "double"
          storage.mode(phi.starting) <- "double"
          storage.mode(sigma.sq.starting) <- "double"
          storage.mode(nu.starting) <- "double"
          storage.mode(w.starting.fit) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
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
          storage.mode(cov.model.indx) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"
          storage.mode(p.det.re) <- "integer"
          storage.mode(X.p.re.fit) <- "integer"
          storage.mode(n.det.re) <- "integer"
          storage.mode(n.det.re.long) <- "integer"
          storage.mode(sigma.sq.p.starting) <- "double"
          storage.mode(sigma.sq.p.a) <- "double"
          storage.mode(sigma.sq.p.b) <- "double"
          storage.mode(alpha.star.starting) <- "double"
          storage.mode(alpha.star.indx) <- "integer"
          storage.mode(lambda.p.fit) <- "double"

          out.fit <- .Call("spMsPGOccNNGPRE", y.fit, X.fit, X.p.fit, coords.fit, 
			   X.p.re.fit, lambda.p.fit, 
		           p.occ, p.det, p.det.re, J.fit, K.fit, N, n.det.re, n.det.re.long, 
          	           n.neighbors, nn.indx.fit, nn.indx.lu.fit, 
			   u.indx.fit, u.indx.lu.fit, ui.indx.fit,
          	           beta.starting, alpha.starting, z.starting.fit,
          	           beta.comm.starting, 
          	           alpha.comm.starting, tau.sq.beta.starting, 
          	           tau.sq.alpha.starting, w.starting.fit, phi.starting, 
          	           sigma.sq.starting, nu.starting, sigma.sq.p.starting, 
		           alpha.star.starting, z.long.indx.fit, alpha.star.indx, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	           nu.a, nu.b, sigma.sq.p.a, sigma.sq.p.b, tuning.c, cov.model.indx, n.batch, 
          	           batch.length, accept.rate, n.omp.threads.fit, verbose.fit, n.report, 
          	           n.burn, n.thin, n.post.samples)

          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
          } else {
            theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
          } 
          colnames(out.fit$theta.samples) <- theta.names
          out.fit$w.samples <- array(out.fit$w.samples, dim = c(N, J, n.post.samples))
          out.fit$w.samples <- aperm(out.fit$w.samples, c(3, 1, 2))
          out.fit$X <- X.fit
	  out.fit$y <- y.big
          out.fit$X.p <- X.p.fit
          out.fit$call <- cl
          out.fit$n.samples <- n.samples
          out.fit$type <- "NNGP"
	  out.fit$n.neighbors <- n.neighbors
          out.fit$coords <- coords.fit
          out.fit$cov.model.indx <- cov.model.indx
          out.fit$n.post <- n.post.samples
          out.fit$n.thin <- n.thin
          out.fit$n.burn <- n.burn
          out.fit$pRE <- TRUE
	  class(out.fit) <- "spMsPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.spMsPGOcc(out.fit, X.0, coords.0, verbose = FALSE)

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  like.samples <- matrix(NA, N, nrow(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
          sp.re.indx <- rep(1:N, each = nrow(out.fit$alpha.star.samples) / N)
	  for (q in 1:N) {
            p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] + 
					    lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
	    for (j in 1:nrow(X.p.0)) {
              like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
						p.0.samples[j, q, ] * 
					        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a)))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }
   
      class(out) <- "spMsPGOcc"
    }

    if (p.det.re == 0) {

      # Run the model in C    
      out <- .Call("spMsPGOccNNGP", y, X, X.p, coords, p.occ, p.det, J, K, N, 
          	 n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
          	 beta.starting, alpha.starting, z.starting,
          	 beta.comm.starting, 
          	 alpha.comm.starting, tau.sq.beta.starting, 
          	 tau.sq.alpha.starting, w.starting, phi.starting, 
          	 sigma.sq.starting, nu.starting, z.long.indx, mu.beta.comm, 
          	 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	 tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	 tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	 nu.a, nu.b, tuning.c, cov.model.indx, n.batch, 
          	 batch.length, accept.rate, n.omp.threads, verbose, n.report, 
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
      out$theta.samples <- mcmc(t(out$theta.samples))
      if (cov.model != 'matern') {
        theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
      } else {
        theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
      } 
      colnames(out$theta.samples) <- theta.names
      # Return things back in the original order. 
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- out$z.samples[, order(ord), ]
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$w.samples <- array(out$w.samples, dim = c(N, J, n.post.samples))
      out$w.samples <- out$w.samples[, order(ord), ]
      out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- out$psi.samples[, order(ord), ]
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      tmp <- matrix(NA, J * K.max, p.det)
      tmp[names.long, ] <- X.p
      tmp <- array(tmp, dim = c(J, K.max, p.det))
      tmp <- tmp[order(ord), , ]
      out$X.p <- matrix(tmp, J * K.max, p.det)
      out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, ]
      out$X <- X[order(ord), , drop = FALSE]
      out$y <- y.big[, order(ord), , drop = FALSE]
      out$call <- cl
      out$n.samples <- n.samples
      tmp <- array(NA, dim = c(N, J * K.max, n.post.samples))
      tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.post.samples))
      tmp <- array(tmp, dim = c(N, J, K.max, n.post.samples))
      out$y.rep.samples <- tmp[, order(ord), , ]
      out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$theta.names <- theta.names
      out$type <- "NNGP"
      out$coords <- coords[order(ord), ]
      out$cov.model.indx <- cov.model.indx
      out$n.neighbors <- n.neighbors
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
        # Currently implemented without parellization. 
	set.seed(k.fold.seed)
	# Number of sites in each hold out data set. 
	sites.random <- sample(1:J)    
        sites.k.fold <- split(sites.random, sites.random %% k.fold)
	registerDoParallel(k.fold.threads)
	model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
	  y.fit <- c(y.big[, -curr.set, , drop = FALSE])
	  y.fit <- y.fit[!is.na(y.fit)]
	  y.0 <- c(y.big[, curr.set, , drop = FALSE])
	  y.0 <- y.0[!is.na(y.0)]
	  z.starting.fit <- z.starting[-curr.set]
	  w.starting.fit <- w.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  coords.fit <- coords[-curr.set, , drop = FALSE]
	  coords.0 <- coords[curr.set, , drop = FALSE]
	  J.fit <- nrow(X.fit)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  # Gotta be a better way, but will do for now. 
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
	  verbose.fit <- FALSE
	  n.omp.threads.fit <- 1
	  # Don't need to reorder things, since they are already sorted by 
	  # the horizontal location in the coordinates. 

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
          storage.mode(z.starting.fit) <- "double"
          storage.mode(X.p.fit) <- "double"
          storage.mode(X.fit) <- "double"
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(N) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(beta.comm.starting) <- "double"
          storage.mode(alpha.comm.starting) <- "double"
          storage.mode(tau.sq.beta.starting) <- "double"
          storage.mode(tau.sq.alpha.starting) <- "double"
          storage.mode(phi.starting) <- "double"
          storage.mode(sigma.sq.starting) <- "double"
          storage.mode(nu.starting) <- "double"
          storage.mode(w.starting.fit) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
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
          storage.mode(cov.model.indx) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"

          out.fit <- .Call("spMsPGOccNNGP", y.fit, X.fit, X.p.fit, coords.fit, 
		           p.occ, p.det, J.fit, K.fit, N, 
          	           n.neighbors, nn.indx.fit, nn.indx.lu.fit, 
			   u.indx.fit, u.indx.lu.fit, ui.indx.fit,
          	           beta.starting, alpha.starting, z.starting.fit,
          	           beta.comm.starting, 
          	           alpha.comm.starting, tau.sq.beta.starting, 
          	           tau.sq.alpha.starting, w.starting.fit, phi.starting, 
          	           sigma.sq.starting, nu.starting, 
		           z.long.indx.fit, mu.beta.comm, 
          	           mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	           tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
          	           tau.sq.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
          	           nu.a, nu.b, tuning.c, cov.model.indx, n.batch, 
          	           batch.length, accept.rate, n.omp.threads.fit, verbose.fit, n.report, 
          	           n.burn, n.thin, n.post.samples)

          if (is.null(sp.names)) {
            sp.names <- paste('sp', 1:N, sep = '')
          }
          coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
          out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
          colnames(out.fit$beta.samples) <- coef.names
          out.fit$theta.samples <- mcmc(t(out.fit$theta.samples))
          if (cov.model != 'matern') {
            theta.names <- paste(rep(c('sigma.sq', 'phi'), each = N), sp.names, sep = '-')
          } else {
            theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = N), sp.names, sep = '-')
          } 
          colnames(out.fit$theta.samples) <- theta.names
          out.fit$w.samples <- array(out.fit$w.samples, dim = c(N, J, n.post.samples))
          out.fit$w.samples <- aperm(out.fit$w.samples, c(3, 1, 2))
          out.fit$X <- X.fit
	  out.fit$y <- y.big
          out.fit$X.p <- X.p.fit
          out.fit$call <- cl
          out.fit$n.samples <- n.samples
          out.fit$type <- "NNGP"
	  out.fit$n.neighbors <- n.neighbors
          out.fit$coords <- coords.fit
          out.fit$cov.model.indx <- cov.model.indx
          out.fit$n.post <- n.post.samples
          out.fit$n.thin <- n.thin
          out.fit$n.burn <- n.burn
          out.fit$pRE <- TRUE
	  class(out.fit) <- "spMsPGOcc"

	  # Predict occurrence at new sites. 
	  out.pred <- predict.spMsPGOcc(out.fit, X.0, coords.0, verbose = FALSE)

	  # Detection 
          sp.indx <- rep(1:N, ncol(X.p.0))
	  like.samples <- matrix(NA, N, nrow(X.p.0))
	  p.0.samples <- array(NA, dim = c(nrow(X.p.0), N, n.post.samples))
	  for (q in 1:N) {
            p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
	    for (j in 1:nrow(X.p.0)) {
              like.samples[q, j] <- mean(dbinom(y.0[N * (j - 1) + q], 1,  
						p.0.samples[j, q, ] * 
					        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
            }
          }
	  apply(like.samples, 1, function(a) sum(log(a)))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }
   
      class(out) <- "spMsPGOcc"
    }
  }

  out$run.time <- proc.time() - ptm
  return(out)
}
