spPGOcc <- function(y, X, X.p, starting, coords, n.rep, n.batch, 
		    batch.length, accept.rate = 0.43, priors, 
		    cov.model = "exponential", tuning, 
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
    # Number of sites
    J <- nrow(y)
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
    tmp <- c(apply(y, c(1, 2), function (a) sum(is.na(a))))
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
    z.long.indx <- z.long.indx[!is.na(c(y))]
    # Subtract 1 for indices in C
    z.long.indx <- z.long.indx - 1
    # y is stored in the following order: species, site, visit
    y.big <- y
    y <- c(y)
    y <- y[!is.na(y)]

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
      # In correct order since you reordered y.
      z.starting <- apply(y.big, 1, max, na.rm = TRUE)
    }
    # beta -----------------------
    if ("beta" %in% names(starting)) {
      beta.starting <- starting[["beta"]]
      if (length(beta.starting) != p.occ) {
        stop(paste("error: starting values for beta must be of length ", p.occ,
		   sep = ""))
      }
    } else {
      beta.starting <- coefficients(glm((z.starting)~X-1, family="binomial"))
    }
    # alpha -----------------------
    if ("alpha" %in% names(starting)) {
      alpha.starting <- starting[["alpha"]]
      if (length(alpha.starting) != p.det) {
        stop(paste("error: starting values for alpha must be of length ", p.det,
		   sep = ""))
      }
    } else {
      alpha.starting <- coefficients(glm((y)~X.p-1, family = 'binomial'))
    }
    # phi -----------------------------
    if ("phi" %in% names(starting)) {
      phi.starting <- starting[["phi"]]
      if (length(phi.starting) != 1) {
        stop("error: starting values for phi must be of length 1")
      }
    } else {
      phi.starting <- 3 / mean(range(coords))
      message("phi is not specified in starting values. Setting starting value to 3/mean(range(coords))\n")
    }
    # sigma.sq ------------------------
    if ("sigma.sq" %in% names(starting)) {
      sigma.sq.starting <- starting[["sigma.sq"]]
      if (length(sigma.sq.starting) != 1) {
        stop("error: starting values for sigma.sq must be of length 1")
      }
    } else {
      sigma.sq.starting <- 2
      message("sigma.sq is not specified in starting values. Setting starting value to 2\n")
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
      message("w is not specified in starting values. Setting starting value to 0\n")
    }
    # nu ------------------------
    if ("nu" %in% names(starting)) {
      nu.starting <- starting[["nu"]]
      if (length(nu.starting) != 1) {
        stop("error: starting values for nu must be of length 1")
      }
    } else {
      if (cov.model == 'matern') {
        message("nu is not specified in starting values. Setting starting value to 1\n")
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
      message("No prior specified for beta.normal. Setting prior mean to 0 and prior variance to 2.73\n")
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
      if (length(mu.alpha) != p.det) {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ", 
		   p.det, " with elements corresponding to alphas' mean", sep = ""))
      }
      if (length(sigma.alpha) != p.det) {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ", 
		   p.det, " with elements corresponding to alphas.comms' variance", sep = ""))
      }
      Sigma.alpha <- sigma.alpha * diag(p.det)
    } else {
      message("No prior specified for alpha.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      mu.alpha <- rep(0, p.det)
      Sigma.alpha <- diag(p.det) * 2.73
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
    storage.mode(cov.model.indx) <- "integer"

    # Get tuning values ---------------------------------------------------
    # Not accessed, but necessary to keep things in line. 
    sigma.sq.tuning <- 0
    phi.tuning <- 0
    nu.tuning <- 0
    if (missing(tuning)) {
      phi.tuning <- exp(1)
      if (cov.model == 'matern') {
        nu.tuning <- exp(1)
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
    tuning.c <- c(sigma.sq.tuning, phi.tuning, nu.tuning)


    # Get distance matrix for full GP -------------------------------------
    coords.D <- iDist(coords)

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.starting) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(coords.D) <- "double"
    storage.mode(beta.starting) <- "double"
    storage.mode(alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    storage.mode(w.starting) <- "double"
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
    storage.mode(tuning.c) <- "double"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"

    ptm <- proc.time()

    # Run the model in C    
    out <- .Call("spPGOcc", y, X, X.p, coords.D, p.occ, p.det, J, K, 
   	         beta.starting, alpha.starting, z.starting,
   	         w.starting, phi.starting, sigma.sq.starting, nu.starting, z.long.indx, 
   	         mu.beta, mu.alpha, 
   	         Sigma.beta, Sigma.alpha, phi.a, phi.b, 
   	         sigma.sq.a, sigma.sq.b, nu.a, nu.b, tuning.c, cov.model.indx,
   	         n.batch, batch.length, 
   	         accept.rate, n.omp.threads, verbose, n.report)

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
    out$z.samples <- mcmc(t(out$z.samples))
    out$psi.samples <- mcmc(t(out$psi.samples))
    tmp <- array(NA, dim = c(J * K.max, n.samples))
    tmp[names.long, ] <- out$y.rep.samples
    out$y.rep.samples <- array(tmp, dim = c(J, K.max, n.samples))
    out$w.samples <- mcmc(t(out$w.samples))
    out$X.occ <- X
    out$X.p <- X.p.orig
    out$y <- y.big
    out$call <- cl
    out$n.samples <- batch.length * n.batch

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
    y <- y[ord, , drop = FALSE]
    X <- X[ord, , drop = FALSE]
    X.p <- X.p[ord, , , drop = FALSE]

    # Get basic info from inputs ------------------------------------------
    # Number of sites
    J <- nrow(y)
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
    tmp <- c(apply(y, c(1, 2), function (a) sum(is.na(a))))
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
    z.long.indx <- z.long.indx[!is.na(c(y))]
    # Subtract 1 for indices in C
    z.long.indx <- z.long.indx - 1
    # y is stored in the following order: species, site, visit
    y <- c(y)
    y <- y[!is.na(y)]

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
      # Reorder the user supplied starting values
      z.starting <- z.starting[ord]
    } else {
      # In correct order since you reordered y.
      z.starting <- apply(y.big, 1, max, na.rm = TRUE)
    }
    # beta -----------------------
    if ("beta" %in% names(starting)) {
      beta.starting <- starting[["beta"]]
      if (length(beta.starting) != p.occ) {
        stop(paste("error: starting values for beta must be of length ", p.occ,
		   sep = ""))
      }
    } else {
      beta.starting <- coefficients(glm((z.starting)~X-1, family="binomial"))
    }
    # alpha -----------------------
    if ("alpha" %in% names(starting)) {
      alpha.starting <- starting[["alpha"]]
      if (length(alpha.starting) != p.det) {
        stop(paste("error: starting values for alpha must be of length ", p.det,
		   sep = ""))
      }
    } else {
      alpha.starting <- coefficients(glm((y)~X.p-1, family = 'binomial'))
    }
    # phi -----------------------------
    if ("phi" %in% names(starting)) {
      phi.starting <- starting[["phi"]]
      if (length(phi.starting) != 1) {
        stop("error: starting values for phi must be of length 1")
      }
    } else {
      phi.starting <- 3 / mean(range(coords))
      message("phi is not specified in starting values. Setting starting value to 3/mean(range(coords))\n")
    }
    # sigma.sq ------------------------
    if ("sigma.sq" %in% names(starting)) {
      sigma.sq.starting <- starting[["sigma.sq"]]
      if (length(sigma.sq.starting) != 1) {
        stop("error: starting values for sigma.sq must be of length 1")
      }
    } else {
      sigma.sq.starting <- 2
      message("sigma.sq is not specified in starting values. Setting starting value to 2\n")
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
      message("w is not specified in starting values. Setting starting value to 0\n")
    }
    # nu ------------------------
    if ("nu" %in% names(starting)) {
      nu.starting <- starting[["nu"]]
      if (length(nu.starting) != 1) {
        stop("error: starting values for nu must be of length 1")
      }
    } else {
      if (cov.model == 'matern') {
        message("nu is not specified in starting values. Setting starting value to 1\n")
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
      message("No prior specified for beta.normal. Setting prior mean to 0 and prior variance to 2.73\n")
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
      if (length(mu.alpha) != p.det) {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ", 
		   p.det, " with elements corresponding to alphas' mean", sep = ""))
      }
      if (length(sigma.alpha) != p.det) {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ", 
		   p.det, " with elements corresponding to alphas.comms' variance", sep = ""))
      }
      Sigma.alpha <- sigma.alpha * diag(p.det)
    } else {
      message("No prior specified for alpha.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      mu.alpha <- rep(0, p.det)
      Sigma.alpha <- diag(p.det) * 2.73
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
    if (!"nu.unif" %in% names(priors)) {
      stop("error: nu.unif must be specified in priors value list")
    }
    if (!is.vector(priors$nu.unif) | !is.atomic(priors$nu.unif) | length(priors$nu.unif) != 2) {
      stop("error: nu.unif must be a vector of length 2 with elements corresponding to nu's lower and upper bounds")
    }
    nu.a <- priors$nu.unif[1]
    nu.b <- priors$nu.unif[2]

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
    sigma.sq.tuning <- 0
    phi.tuning <- 0
    nu.tuning <- 0
    if (missing(tuning)) {
      phi.tuning <- exp(1)
      if (cov.model == 'matern') {
        nu.tuning <- exp(1)
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
    storage.mode(z.starting) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p.det) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(K) <- "integer"
    storage.mode(coords) <- "double"
    storage.mode(beta.starting) <- "double"
    storage.mode(alpha.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(sigma.sq.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    storage.mode(w.starting) <- "double"
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

    # Run the model in C --------------------------------------------------
    ptm <- proc.time()

    out <- .Call("spPGOccNNGP", y, X, X.p, coords, p.occ, p.det, J, K, 
		 n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
    	         beta.starting, alpha.starting, z.starting,
    	         phi.starting, sigma.sq.starting, nu.starting, 
		 z.long.indx, mu.beta, mu.alpha, 
    	         Sigma.beta, Sigma.alpha, phi.a, phi.b, 
    	         sigma.sq.a, sigma.sq.b, nu.a, nu.b, tuning.c, 
		 cov.model.indx, n.batch, batch.length, 
    	         accept.rate, n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    # Get everything back in the original order
    out$coords <- coords[order(ord), ]
    out$z.samples <- mcmc(t(out$z.samples[order(ord), , drop = FALSE]))
    out$X.occ <- X[order(ord), , drop = FALSE]
    out$w.samples <- mcmc(t(out$w.samples[order(ord), , drop = FALSE]))
    out$psi.samples <- mcmc(t(out$psi.samples[order(ord), , drop = FALSE]))
    out$X.p <- X.p.orig[order(ord), , , drop = FALSE]
    out$y <- y.big
    tmp <- array(NA, dim = c(J * K.max, n.samples))
    tmp[names.long, ] <- out$y.rep.samples
    tmp <- array(tmp, dim = c(J, K.max, n.samples))
    out$y.rep.samples <- tmp[order(ord), , ]

    # Return all the other good stuff. 
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
    out$call <- cl
    out$n.samples <- batch.length * n.batch

    # class(out) <- "PGOcc"
    
    out
  }
}
