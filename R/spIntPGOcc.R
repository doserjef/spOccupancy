spIntPGOcc <- function(occ.formula, det.formula, data, starting, n.batch, 
		       batch.length, accept.rate = 0.43, priors, 
		       cov.model = "exponential", tuning, 
		       n.omp.threads = 1, verbose = TRUE, NNGP = FALSE, 
		       n.neighbors = 15, search.type = "cb", 
		       n.report = 100, ...){

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
  if (!is.list(y)) {
    stop("error: y must be a list of detection-nondetection data sets")
  }
  y <- data$y
  n.data <- length(y)
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
        message("occupancy covariates (occ.covs) not specified in data. Assuming intercept only occupancy model.")
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
        message("detection covariates (det.covs) not specified in data. Assuming interept only detection model for each data source.")
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
  coords <- data$coords

  # Neighbors and Ordering ----------------------------------------------
  u.search.type <- 2 
  # Order by x column. Could potentailly allow this to be user defined.
  ord <- order(coords[,1])
  # Reorder everything to align with NN ordering. 
  coords <- coords[ord, ]
  X <- X[ord, , drop = FALSE]
  sites.orig <- sites
  # Don't need to actually reorder the data, can just reorder the site 
  # indices. 
  for (i in 1:n.data) {
    for (j in 1:length(sites[[i]])) {
    sites[[i]][j] <- which(ord == sites[[i]][j])
    }
  }

  # Make all covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  for (i in 1:n.data) {
    data$det.covs[[i]] <- data.frame(lapply(data$det.covs[[i]], function(a) unlist(c(a))))
    # Replicate det.covs if only covariates are at the site level. 
    if (nrow(data$det.covs[[i]]) == nrow(y[[i]])) {
      data$det.covs[[i]] <- data.frame(sapply(data$det.covs[[i]], rep, times = dim(y[[i]][2])))
    }
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
      if (nrow(X.p[[i]]) == length(y[[1]])) {
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

    # Starting values -----------------------------------------------------
    if (missing(starting)) {
      stop("error: starting value list for the parameters must be specified.")
    }
    names(starting) <- tolower(names(starting))
    # z -------------------------------
    if ("z" %in% names(starting)) {
      z.starting <- starting$z
      if (!is.vector(z.starting)) {
        stop(paste("error: starting values for z must be a vector of length ",
		   J, sep = ""))
      }
      if (length(z.starting) != J) {
        stop(paste("error: starting values for z must be a vector of length ",
		   J, sep = ""))
      }
    } else {
      z.starting <- tapply(y, z.long.indx.r, max, na.rm = TRUE)
    }
    # beta -----------------------
    if ("beta" %in% names(starting)) {
      beta.starting <- starting[["beta"]]
      if (length(beta.starting) != p.occ) {
        stop(paste("error: starting values for beta must be of length ", p.occ,
		   sep = ""))
      }
    } else {
      y.max <- tapply(y, z.long.indx.r, max, na.rm = TRUE)
      beta.starting <- coefficients(glm((y.max)~X-1, family="binomial"))
    }
    # alpha -----------------------
    if ("alpha" %in% names(starting)) {
      alpha.starting <- starting[["alpha"]]
      if (length(alpha.starting) != n.data | !is.list(alpha.starting)) {
        stop(paste("error: starting values for alpha must be a list of length ", n.data,
		   sep = ""))
      }
      for (q in 1:n.data) {
        if (length(alpha.starting[[q]]) != p.det.long[q] | !is.vector(alpha.starting[[q]])) {
          stop(paste("error: starting values for alpha[[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
        }
      }
      alpha.starting <- unlist(alpha.starting)
    } else {
      if (verbose) {
        message("alpha is not specified in starting values. Setting starting value to 0\n")
      }
      alpha.starting <- rep(0, p.det)
    }

    alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
    alpha.indx.c <- alpha.indx.r - 1

    # phi -----------------------------
    if ("phi" %in% names(starting)) {
      phi.starting <- starting[["phi"]]
      if (length(phi.starting) != 1) {
        stop("error: starting values for phi must be of length 1")
      }
    } else {
      phi.starting <- 3 / mean(range(coords))
      if (verbose) {
        message("phi is not specified in starting values. Setting starting value to 3/mean(range(coords))\n")
      }
    }
    # sigma.sq ------------------------
    if ("sigma.sq" %in% names(starting)) {
      sigma.sq.starting <- starting[["sigma.sq"]]
      if (length(sigma.sq.starting) != 1) {
        stop("error: starting values for sigma.sq must be of length 1")
      }
    } else {
      sigma.sq.starting <- 2
      if (verbose) {
        message("sigma.sq is not specified in starting values. Setting starting value to 2\n")
      }
    }
    # w -----------------------------00
    if ("w" %in% names(starting)) {
      w.starting <- starting[["w"]]
      if (!is.vector(w.starting)) {
        stop(paste("error: starting values for w must be a vector of length ",
		   J, sep = ""))
      }
      if (length(w.starting) != J) {
        stop(paste("error: starting values for w must be a vector of length ",
		   J, sep = ""))
      }
    } else {
      w.starting <- rep(0, J)
      if (verbose) {
        message("w is not specified in starting values. Setting starting value to 0\n")
      }
    }
    # nu ------------------------
    if ("nu" %in% names(starting)) {
      nu.starting <- starting[["nu"]]
      if (length(nu.starting) != 1) {
        stop("error: starting values for nu must be of length 1")
      }
    } else {
      if (cov.model == 'matern') {
        if (verbose) {
          message("nu is not specified in starting values. Setting starting value to 1\n")
        }
        nu.starting <- 1
      } else {
        nu.starting <- 0
      }
    }

    # Priors --------------------------------------------------------------
    if (missing(priors)) {
      stop("error: prior list for the parameters must be specified")
    }
    names(priors) <- tolower(names(priors))
    # beta -----------------------
    if ("beta.normal" %in% names(priors)) {
      mu.beta <- priors$beta.normal[[1]]
      sigma.beta <- priors$beta.normal[[2]]
      if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
        stop("error: beta.normal must be a list of length 2")
      }
      if (length(mu.beta) != p.occ) {
        stop(paste("error: beta.normal[[1]] must be a vector of length ", 
		   p.occ, " with elements corresponding to betas' mean", sep = ""))
      }
      if (length(sigma.beta) != p.occ) {
        stop(paste("error: beta.normal[[2]] must be a vector of length ", 
		   p.occ, " with elements corresponding to betas' variance", sep = ""))
      }
      Sigma.beta <- sigma.beta * diag(p.occ)
    } else {
      if (verbose) {
        message("No prior specified for beta.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      }
      mu.beta <- rep(0, p.occ)
      Sigma.beta <- diag(p.occ) * 2.73
    }

    # alpha -----------------------
    if ("alpha.normal" %in% names(priors)) {
      mu.alpha <- priors$alpha.normal[[1]]
      sigma.alpha <- priors$alpha.normal[[2]]
      if (!is.list(priors$alpha.normal) | length(priors$alpha.normal) != 2) {
        stop("error: alpha.normal must be a list of length 2")
      }
      if (length(mu.alpha) != n.data | !is.list(mu.alpha)) {
        stop(paste("error: alpha.normal[[1]] must be a list of length ", 
		   n.data, " with elements corresponding to alphas' mean for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(mu.alpha[[q]]) != p.det.long[q] | !is.vector(mu.alpha[[q]])) {
          stop(paste("error: prior means for alpha.normal[['mean']][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
        }
      }
      mu.alpha <- unlist(mu.alpha)
      if (length(sigma.alpha) != n.data | !is.list(sigma.alpha)) {
        stop(paste("error: alpha.normal[[2]] must be a list of length ", 
		   n.data, " with elements corresponding to alphas' variance for each data set", sep = ""))
      }
      for (q in 1:n.data) {
        if (length(sigma.alpha[[q]]) != p.det.long[q] | !is.vector(sigma.alpha[[q]])) {
          stop(paste("error: prior variances for alpha.normal[['var']][[", q, "]] must be a vector of length ", 
		     p.det.long[q], sep = ""))
        }
      }
      sigma.alpha <- unlist(sigma.alpha)
    } else {
      if (verbose) {
        message("No prior specified for alpha.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      }
      mu.alpha <- rep(0, p.det)
      sigma.alpha <- rep(2.73, p.det) 
    }

    # phi -----------------------------
    if (!"phi.unif" %in% names(priors)) {
      stop("error: phi.unif must be specified in priors value list")
    }
    if (!is.vector(priors$phi.unif) | !is.atomic(priors$phi.unif) | length(priors$phi.unif) != 2) {
      stop("error: phi.unif must be a vector of length 2 with elements corresponding to phi's lower and upper bounds")
    }
    phi.a <- priors$phi.unif[1]
    phi.b <- priors$phi.unif[2]
    # sigma.sq -----------------------------
    if (!"sigma.sq.ig" %in% names(priors)) {
      stop("error: sigma.sq.ig must be specified in priors value list")
    }
    if (!is.vector(priors$sigma.sq.ig) | !is.atomic(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
      stop("error: sigma.sq.ig must be a vector of length 2 with elements corresponding to sigma.sq's shape and scale parameters")
    }
    sigma.sq.a <- priors$sigma.sq.ig[1]
    sigma.sq.b <- priors$sigma.sq.ig[2]
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
    sigma.sq.tuning <- 0
    phi.tuning <- 0
    nu.tuning <- 0
    if (missing(tuning)) {
      phi.tuning <- 1
      if (cov.model == 'matern') {
        nu.tuning <- 1
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
    }
    tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))

  if (!NNGP) { 

    # Get distance matrix for full GP -------------------------------------
    coords.D <- iDist(coords)

    # Specify storage modes -----------------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.starting) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
    storage.mode(coords.D) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.det.long) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(n.obs) <- "integer"
    storage.mode(n.obs.long) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(J.long) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(n.data) <- "integer"
    storage.mode(beta.starting) <- "double"
    storage.mode(alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    storage.mode(w.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
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
    storage.mode(tuning.c) <- "double"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"

    ptm <- proc.time()

    # Run the model in C    
    out <- .Call("spIntPGOcc", y, X, X.p.all, coords.D, p.occ, p.det, p.det.long, 
		 J, J.long, K, n.obs, n.obs.long, n.data, 
		 beta.starting, alpha.starting, z.starting, w.starting, 
		 phi.starting, sigma.sq.starting, nu.starting, 
		 z.long.indx.c, data.indx.c, alpha.indx.c, mu.beta, mu.alpha, 
		 Sigma.beta, sigma.alpha, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		 nu.a, nu.b, tuning.c, cov.model.indx, 
		 n.batch, batch.length, accept.rate,  
		 n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- x.names
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    colnames(out$alpha.samples) <- x.p.names
    out$theta.samples <- mcmc(t(out$theta.samples))
    if (cov.model != 'matern') {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi')
    } else {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
    }
    if (cov.model != 'matern') {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi')
    } else {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
    }
    out$z.samples <- mcmc(t(out$z.samples))
    out$w.samples <- mcmc(t(out$w.samples))
    out$psi.samples <- mcmc(t(out$psi.samples))
    # y.rep.samples is returned as a list, where each element 
    # corresponds to a different data set. 
    tmp <- list()
    indx <- 1
    n.samples <- batch.length * n.batch
    for (q in 1:n.data) {
      tmp[[q]] <- array(NA, dim = c(J.long[q] * K.long.max[q], n.samples))
      tmp[[q]][names.long[[q]], ] <- out$y.rep.samples[indx:(indx + n.obs.long[q] - 1), ] 
      tmp[[q]] <- array(tmp[[q]], dim = c(J.long[q], K.long.max[q], n.samples))
      tmp[[q]] <- aperm(tmp[[q]], c(3, 1, 2))
      indx <- indx + n.obs.long[q]
    }
    out$y.rep.samples <- tmp
    out$X <- X
    out$X.p <- X.p.orig
    out$y <- y.big
    out$call <- cl
    out$n.samples <- n.samples
    out$sites <- sites
    out$cov.model.indx <- cov.model.indx
    out$type <- "GP"
    out$coords <- coords

    class(out) <- "spIntPGOcc"
    
    out
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
    storage.mode(z.starting) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.det.long) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(n.obs) <- "integer"
    storage.mode(n.obs.long) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(J.long) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(n.data) <- "integer"
    storage.mode(beta.starting) <- "double"
    storage.mode(alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    storage.mode(w.starting) <- "double"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
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


    ptm <- proc.time()

    # Run the model in C    
    out <- .Call("spIntPGOccNNGP", y, X, X.p.all, coords, p.occ, p.det, p.det.long, 
		 J, J.long, K, n.obs, n.obs.long, n.data, 
		 n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
		 beta.starting, alpha.starting, z.starting, w.starting, 
		 phi.starting, sigma.sq.starting, nu.starting, 
		 z.long.indx.c, data.indx.c, alpha.indx.c, mu.beta, mu.alpha, 
		 Sigma.beta, sigma.alpha, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		 nu.a, nu.b, tuning.c, cov.model.indx,
		 n.batch, batch.length, accept.rate,  
		 n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    # Get everything back in the original order
    out$coords <- coords[order(ord), ]
    out$z.samples <- mcmc(t(out$z.samples[order(ord), , drop = FALSE]))
    out$X <- X[order(ord), , drop = FALSE]
    out$w.samples <- mcmc(t(out$w.samples[order(ord), , drop = FALSE]))
    out$psi.samples <- mcmc(t(out$psi.samples[order(ord), , drop = FALSE]))

    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- x.names
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    colnames(out$alpha.samples) <- x.p.names
    out$theta.samples <- mcmc(t(out$theta.samples))
    if (cov.model != 'matern') {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi')
    } else {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
    }
    # Get y in a useful format. 
    tmp <- list()
    indx <- 1
    n.samples <- batch.length * n.batch
    for (q in 1:n.data) {
      tmp[[q]] <- array(NA, dim = c(J.long[q] * K.long.max[q], n.samples))
      tmp[[q]][names.long[[q]], ] <- out$y.rep.samples[indx:(indx + n.obs.long[q] - 1), ] 
      tmp[[q]] <- array(tmp[[q]], dim = c(J.long[q], K.long.max[q], n.samples))
      tmp[[q]] <- aperm(tmp[[q]], c(3, 1, 2))
      indx <- indx + n.obs.long[q]
    }
    out$y.rep.samples <- tmp
    out$X.p <- X.p.orig
    out$y <- y.big
    out$call <- cl
    out$n.samples <- batch.length * n.batch
    out$sites <- sites.orig
    out$n.neighbors <- n.neighbors
    out$cov.model.indx <- cov.model.indx
    out$type <- "NNGP"

    class(out) <- "spIntPGOcc"
    
    out

  }
}
