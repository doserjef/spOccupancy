spMsPGOcc <- function(y, X, X.p, starting, coords, n.rep, n.batch, 
		      batch.length, accept.rate = 0.43, priors, 
		      cov.model = 'exponential', tuning,
		      n.omp.threads = 1, verbose = TRUE, NNGP = FALSE, 
		      n.neighbors = 15, search.type = "cb",
		      n.report = 100, ...){

  if (!NNGP) {
    # Make it look nice
    cat("----------------------------------------\n");
    cat("\tPreparing the data\n");
    cat("----------------------------------------\n");
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
    if (missing(y)) {
      stop("error: data (y) must be specified") 
    }
    if (missing(X)) {
      message("X is not specified. Assuming intercept only occupancy model.\n")
      X <- matrix(1, J, 1)
    }
    if (missing(X.p)) { 
      message("X.p is not specified. Assuming intercept only detection model.\n")
      X.p <- array(1, dim = c(J, dim(y)[3], 1))
    }
    if (missing(coords)) {
      stop("error: must specify coords for spatial model. For non spatial model, use msPGOcc")
    }

    # Get basic info from inputs ------------------------------------------
    # Number of species 
    N <- dim(y)[1]
    # Number of sites
    J <- dim(y)[2]
    # Number of occupancy parameters 
    p.occ <- ncol(X)
    # Number of detection parameters
    p.det <- dim(X.p)[3]
    # Get names for occupancy covariates
    if (length(colnames(X)) == 0) {
      x.names <- paste('x', 1:p.occ, sep = '')
    } else {
      x.names <- colnames(X)
    }
    # Get names for detection parameters and reorganize
    if (length(dimnames(X.p)[3]) == 0) {
      x.p.names <- paste('x.p', 1:p.det, sep = '')
    } else {
      x.p.names <- dimnames(X.p)[3]
    }
    X.p.orig <- X.p
    X.p <- matrix(X.p, nrow = dim(X.p)[1] * dim(X.p)[2], 
		  ncol = dim(X.p)[3])
    rownames(X.p) <- 1:nrow(X.p)
    # Remove missing values in design matrix
    # Get missing values from the first species, which should always be there
    tmp <- c(apply(y[1, , ], c(1, 2), function (a) sum(is.na(a))))
    X.p <- X.p[tmp == 0, , drop = FALSE]
    # Number of pseudoreplicates
    n.obs <- nrow(X.p)
    # Need these rownames later on for reordering
    names.long <- as.numeric(rownames(X.p))
    if (missing(n.rep)) {
      stop("error: number of replicates (n.rep) must be specified")
    }
    if (length(n.rep) == 1) {
      n.rep <- rep(n.rep, J) 
    } else if (length(n.rep) != J) {
      stop("error: n.rep must be of length 1 or J")
    }
    # Max number of repeat visits
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
    y <- y[!is.na(y)]

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
    } else {
      # In correct order since you reordered y. 
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
    # tau.beta ------------------------
    if ("tau.beta" %in% names(starting)) {
      tau.beta.starting <- starting[["tau.beta"]]
      if (length(tau.beta.starting) != p.occ) {
        stop(paste("error: starting values for tau.beta must be of length ", p.occ, 
		   sep = ""))
      }
    } else {
      tau.beta.starting <- runif(p.occ, 0.1, 2)
    }
    # tau.alpha -----------------------
    if ("tau.alpha" %in% names(starting)) {
      tau.alpha.starting <- starting[["tau.alpha"]]
      if (length(tau.alpha.starting) != p.det) {
        stop(paste("error: starting values for tau.alpha must be of length ", p.det, 
		   sep = ""))
      }
    } else {
      tau.alpha.starting <- runif(p.det, 0.1, 2)
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
      message("phi is not specified in starting values. Setting starting value to 3/mean(range(coords))\n")
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
      message("sigma.sq is not specified in starting values. Setting starting value to 2\n")
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
      message("w is not specified in starting values. Setting starting value to 0\n")
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
        message("nu is not specified in starting values. Setting starting value to 1\n")
        nu.starting <- rep(1, N)
      } else {
        nu.starting <- rep(0, N)
      }
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
      message("No prior specified for beta.comm.normal. Setting prior mean to 0 and prior variance to 2.73\n")
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
      message("No prior specified for alpha.comm.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      mu.alpha.comm <- rep(0, p.det)
      Sigma.alpha.comm <- diag(p.det) * 2.73
    }


    # tau.beta -----------------------
    if ("tau.beta.ig" %in% names(priors)) {
      tau.beta.a <- priors$tau.beta.ig[[1]]
      tau.beta.b <- priors$tau.beta.ig[[2]]
      if (!is.list(priors$tau.beta.ig) | length(priors$tau.beta.ig) != 2) {
        stop("error: tau.beta.ig must be a list of length 2")
      }
      if (length(tau.beta.a) != p.occ) {
        stop(paste("error: tau.beta.ig[[1]] must be a vector of length ", 
		   p.occ, " with elements corresponding to tau.betas' shape", sep = ""))
      }
      if (length(tau.beta.b) != p.occ) {
        stop(paste("error: tau.beta.ig[[2]] must be a vector of length ", 
		   p.occ, " with elements corresponding to tau.betas' scale", sep = ""))
      }
    } else {
      message("No prior specified for tau.beta.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
      tau.beta.a <- rep(0.1, p.occ)
      tau.beta.b <- rep(0.1, p.occ)
    }

    # tau.alpha -----------------------
    if ("tau.alpha.ig" %in% names(priors)) {
      tau.alpha.a <- priors$tau.alpha.ig[[1]]
      tau.alpha.b <- priors$tau.alpha.ig[[2]]
      if (!is.list(priors$tau.alpha.ig) | length(priors$tau.alpha.ig) != 2) {
        stop("error: tau.alpha.ig must be a list of length 2")
      }
      if (length(tau.alpha.a) != p.det) {
        stop(paste("error: tau.alpha.ig[[1]] must be a vector of length ", 
		   p.det, " with elements corresponding to tau.alphas' shape", sep = ""))
      }
      if (length(tau.alpha.b) != p.det) {
        stop(paste("error: tau.alpha.ig[[2]] must be a vector of length ", 
		   p.det, " with elements corresponding to tau.alphas' scale", sep = ""))
      }
    } else {
      message("No prior specified for tau.alpha.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
      tau.alpha.a <- rep(0.1, p.det)
      tau.alpha.b <- rep(0.1, p.det)
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
      phi.tuning <- rep(exp(1), N)
      if (cov.model == 'matern') {
        nu.tuning <- rep(exp(1), N)
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
    tuning.c <- c(sigma.sq.tuning, phi.tuning, nu.tuning)

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
    storage.mode(tau.beta.starting) <- "double"
    storage.mode(tau.alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(w.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(mu.alpha.comm) <- "double"
    storage.mode(Sigma.alpha.comm) <- "double"
    storage.mode(tau.beta.a) <- "double"
    storage.mode(tau.beta.b) <- "double"
    storage.mode(tau.alpha.a) <- "double"
    storage.mode(tau.alpha.b) <- "double"
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

    ####################################################
    ##Other stuff
    ####################################################

    ptm <- proc.time()

    # Run the model in C    
    out <- .Call("spMsPGOcc", y, X, X.p, coords.D, p.occ, p.det, J, K, N, 
		 beta.starting, alpha.starting, z.starting,
		 beta.comm.starting, 
		 alpha.comm.starting, tau.beta.starting, 
		 tau.alpha.starting, w.starting, phi.starting, 
		 sigma.sq.starting, nu.starting, z.long.indx, mu.beta.comm, 
		 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
		 tau.beta.a, tau.beta.b, tau.alpha.a, 
		 tau.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		 nu.a, nu.b, tuning.c, cov.model.indx, 
		 n.batch, batch.length, accept.rate, 
		 n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    out$beta.comm.samples <- mcmc(t(out$beta.comm.samples))
    colnames(out$beta.comm.samples) <- x.names
    out$alpha.comm.samples <- mcmc(t(out$alpha.comm.samples))
    colnames(out$alpha.comm.samples) <- x.p.names
    out$tau.beta.samples <- mcmc(t(out$tau.beta.samples))
    colnames(out$tau.beta.samples) <- x.names
    out$tau.alpha.samples <- mcmc(t(out$tau.alpha.samples))
    colnames(out$tau.alpha.samples) <- x.p.names
    names.1 <- rep(paste("sp", 1:N, sep = ''), p.occ)
    names.2 <- rep(paste('x', 1:p.occ, sep = ''), each = N)
    coef.names <- paste(names.1, names.2, sep = '-')
    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- coef.names
    names.1.det <- rep(paste("sp", 1:N, sep = ''), p.det)
    names.2.det <- rep(paste('x', 1:p.det, sep = ''), each = N)
    coef.names.det <- paste(names.1.det, names.2.det, sep = '-')
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    colnames(out$alpha.samples) <- coef.names.det
    out$theta.samples <- mcmc(t(out$theta.samples))
    out$z.samples <- array(out$z.samples, dim = c(N, J, n.samples))
    out$w.samples <- array(out$w.samples, dim = c(N, J, n.samples))
    out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.samples))
    tmp <- array(NA, dim = c(N, J * K.max, n.samples))
    tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.samples))
    out$y.rep.samples <- array(tmp, dim = c(N, J, K.max, n.samples))
    out$X.occ <- X
    out$X.p <- X.p
    out$y <- y.big
    out$call <- cl
    out$n.samples <- n.samples

    # class(out) <- "PGOcc"
    
    out
  } else {

    # Make it look nice
    cat("----------------------------------------\n");
    cat("\tPreparing the data\n");
    cat("----------------------------------------\n");
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
    if (missing(y)) {
      stop("error: data (y) must be specified") 
    }
    if (missing(X)) {
      message("X is not specified. Assuming intercept only occupancy model.\n")
      X <- matrix(1, J, 1)
    }
    if (missing(X.p)) { 
      message("X.p is not specified. Assuming intercept only detection model.\n")
      X.p <- array(1, dim = c(J, dim(y)[3], 1))
    }
    if (missing(coords)) {
      stop("error: must specify coords for spatial model. For non spatial model, use msPGOcc")
    }

    # Neighbors and Ordering ----------------------------------------------
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y.big <- y
    y <- y[, ord, , drop = FALSE]
    X <- X[ord, , drop = FALSE]
    X.p <- X.p[ord, , , drop = FALSE]

    # Get basic info from inputs ------------------------------------------
    # Number of species 
    N <- dim(y)[1]
    # Number of sites
    J <- dim(y)[2]
    # Number of occupancy parameters 
    p.occ <- ncol(X)
    # Number of detection parameters
    p.det <- dim(X.p)[3]
    # Get names for occupancy covariates
    if (length(colnames(X)) == 0) {
      x.names <- paste('x', 1:p.occ, sep = '')
    } else {
      x.names <- colnames(X)
    }
    # Get names for detection parameters and reorganize
    if (length(dimnames(X.p)[3]) == 0) {
      x.p.names <- paste('x.p', 1:p.det, sep = '')
    } else {
      x.p.names <- dimnames(X.p)[3]
    }
    X.p.orig <- X.p
    X.p <- matrix(X.p, nrow = dim(X.p)[1] * dim(X.p)[2], 
		  ncol = dim(X.p)[3])
    rownames(X.p) <- 1:nrow(X.p)
    # Remove missing values in design matrix
    # Get missing values from the first species, which should always be there
    tmp <- c(apply(y[1, , ], c(1, 2), function (a) sum(is.na(a))))
    X.p <- X.p[tmp == 0, , drop = FALSE]
    # Number of pseudoreplicates
    n.obs <- nrow(X.p)
    # Need these rownames later on for reordering
    names.long <- as.numeric(rownames(X.p))
    if (missing(n.rep)) {
      stop("error: number of replicates (n.rep) must be specified")
    }
    if (length(n.rep) == 1) {
      n.rep <- rep(n.rep, J) 
    } else if (length(n.rep) != J) {
      stop("error: n.rep must be of length 1 or J")
    }
    # Reorder n.rep after the fact. 
    n.rep <- n.rep[ord]
    # Max number of repeat visits
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
    y <- c(y)
    y <- y[!is.na(y)]

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
      # Reorder the user supply starting values
      z.starting <- z.starting[, ord]
    } else {
      # In correct order since you reordered y. 
      z.starting <- apply(y.big[, ord, , drop = FALSE], c(1, 2), max, na.rm = TRUE)
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
    # tau.beta ------------------------
    if ("tau.beta" %in% names(starting)) {
      tau.beta.starting <- starting[["tau.beta"]]
      if (length(tau.beta.starting) != p.occ) {
        stop(paste("error: starting values for tau.beta must be of length ", p.occ, 
		   sep = ""))
      }
    } else {
      tau.beta.starting <- runif(p.occ, 0.1, 2)
    }
    # tau.alpha -----------------------
    if ("tau.alpha" %in% names(starting)) {
      tau.alpha.starting <- starting[["tau.alpha"]]
      if (length(tau.alpha.starting) != p.det) {
        stop(paste("error: starting values for tau.alpha must be of length ", p.det, 
		   sep = ""))
      }
    } else {
      tau.alpha.starting <- runif(p.det, 0.1, 2)
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
      message("phi is not specified in starting values. Setting starting value to 3/mean(range(coords))\n")
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
      message("sigma.sq is not specified in starting values. Setting starting value to 2\n")
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
      message("w is not specified in starting values. Setting starting value to 0\n")
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
        message("nu is not specified in starting values. Setting starting value to 1\n")
        nu.starting <- rep(1, N)
      } else {
        nu.starting <- rep(0, N)
      }
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
      message("No prior specified for beta.comm.normal. Setting prior mean to 0 and prior variance to 2.73\n")
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
      message("No prior specified for alpha.comm.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      mu.alpha.comm <- rep(0, p.det)
      Sigma.alpha.comm <- diag(p.det) * 2.73
    }


    # tau.beta -----------------------
    if ("tau.beta.ig" %in% names(priors)) {
      tau.beta.a <- priors$tau.beta.ig[[1]]
      tau.beta.b <- priors$tau.beta.ig[[2]]
      if (!is.list(priors$tau.beta.ig) | length(priors$tau.beta.ig) != 2) {
        stop("error: tau.beta.ig must be a list of length 2")
      }
      if (length(tau.beta.a) != p.occ) {
        stop(paste("error: tau.beta.ig[[1]] must be a vector of length ", 
		   p.occ, " with elements corresponding to tau.betas' shape", sep = ""))
      }
      if (length(tau.beta.b) != p.occ) {
        stop(paste("error: tau.beta.ig[[2]] must be a vector of length ", 
		   p.occ, " with elements corresponding to tau.betas' scale", sep = ""))
      }
    } else {
      message("No prior specified for tau.beta.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
      tau.beta.a <- rep(0.1, p.occ)
      tau.beta.b <- rep(0.1, p.occ)
    }

    # tau.alpha -----------------------
    if ("tau.alpha.ig" %in% names(priors)) {
      tau.alpha.a <- priors$tau.alpha.ig[[1]]
      tau.alpha.b <- priors$tau.alpha.ig[[2]]
      if (!is.list(priors$tau.alpha.ig) | length(priors$tau.alpha.ig) != 2) {
        stop("error: tau.alpha.ig must be a list of length 2")
      }
      if (length(tau.alpha.a) != p.det) {
        stop(paste("error: tau.alpha.ig[[1]] must be a vector of length ", 
		   p.det, " with elements corresponding to tau.alphas' shape", sep = ""))
      }
      if (length(tau.alpha.b) != p.det) {
        stop(paste("error: tau.alpha.ig[[2]] must be a vector of length ", 
		   p.det, " with elements corresponding to tau.alphas' scale", sep = ""))
      }
    } else {
      message("No prior specified for tau.alpha.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
      tau.alpha.a <- rep(0.1, p.det)
      tau.alpha.b <- rep(0.1, p.det)
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
    storage.mode(cov.model.indx) <- "integer"

    # Get tuning values ---------------------------------------------------
    # Not accessed, but necessary to keep things in line. 
    sigma.sq.tuning <- rep(0, N)
    phi.tuning <- rep(0, N)
    nu.tuning <- rep(0, N)
    if (missing(tuning)) {
      phi.tuning <- rep(exp(1), N)
      if (cov.model == 'matern') {
        nu.tuning <- rep(exp(1), N)
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
    tuning.c <- c(sigma.sq.tuning, phi.tuning, nu.tuning)

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
    storage.mode(tau.beta.starting) <- "double"
    storage.mode(tau.alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    storage.mode(w.starting) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(mu.alpha.comm) <- "double"
    storage.mode(Sigma.alpha.comm) <- "double"
    storage.mode(tau.beta.a) <- "double"
    storage.mode(tau.beta.b) <- "double"
    storage.mode(tau.alpha.a) <- "double"
    storage.mode(tau.alpha.b) <- "double"
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

    ptm <- proc.time()

    # Run the model in C    
    out <- .Call("spMsPGOccNNGP", y, X, X.p, coords, p.occ, p.det, J, K, N, 
		 n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
		 beta.starting, alpha.starting, z.starting,
		 beta.comm.starting, 
		 alpha.comm.starting, tau.beta.starting, 
		 tau.alpha.starting, w.starting, phi.starting, 
		 sigma.sq.starting, nu.starting, z.long.indx, mu.beta.comm, 
		 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
		 tau.beta.a, tau.beta.b, tau.alpha.a, 
		 tau.alpha.b, phi.a, phi.b, sigma.sq.a, sigma.sq.b, 
		 nu.a, nu.b, tuning.c, cov.model.indx, n.batch, 
		 batch.length, accept.rate, n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    out$beta.comm.samples <- mcmc(t(out$beta.comm.samples))
    colnames(out$beta.comm.samples) <- x.names
    out$alpha.comm.samples <- mcmc(t(out$alpha.comm.samples))
    colnames(out$alpha.comm.samples) <- x.p.names
    out$tau.beta.samples <- mcmc(t(out$tau.beta.samples))
    colnames(out$tau.beta.samples) <- x.names
    out$tau.alpha.samples <- mcmc(t(out$tau.alpha.samples))
    colnames(out$tau.alpha.samples) <- x.p.names
    names.1 <- rep(paste("sp", 1:N, sep = ''), p.occ)
    names.2 <- rep(paste('x', 1:p.occ, sep = ''), each = N)
    coef.names <- paste(names.1, names.2, sep = '-')
    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- coef.names
    names.1.det <- rep(paste("sp", 1:N, sep = ''), p.det)
    names.2.det <- rep(paste('x', 1:p.det, sep = ''), each = N)
    coef.names.det <- paste(names.1.det, names.2.det, sep = '-')
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    colnames(out$alpha.samples) <- coef.names.det
    out$theta.samples <- mcmc(t(out$theta.samples))
    # Return things back in the original order. 
    out$z.samples <- array(out$z.samples, dim = c(N, J, n.samples))
    out$z.samples <- out$z.samples[, order(ord), ]
    out$w.samples <- array(out$w.samples, dim = c(N, J, n.samples))
    out$w.samples <- out$w.samples[, order(ord), ]
    out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.samples))
    out$psi.samples <- out$psi.samples[, order(ord), ]
    out$X.occ <- X[order(ord), , drop = FALSE]
    out$X.p <- X.p.orig[order(ord), , , drop = FALSE]
    out$y <- y.big
    out$call <- cl
    out$n.samples <- n.samples
    tmp <- array(NA, dim = c(N, J * K.max, n.samples))
    tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.samples))
    tmp <- array(tmp, dim = c(N, J, K.max, n.samples))
    out$y.rep.samples <- tmp[, order(ord), , ]
   
    # class(out) <- "PGOcc"
    
    out
  }
}
