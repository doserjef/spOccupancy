msPGOcc <- function(occ.formula, det.formula, data, starting, n.samples, 
		    priors, n.omp.threads = 1, verbose = TRUE, n.report = 100, ...){
 
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
        message("occupancy covariates (occ.covs) not specified in data. Assuming intercept only occupancy model.")
        data$occ.covs <- matrix(1, dim(y)[2], 1)
      } else {
        stop("error: occ.covs must be specified in data for an occupancy model with covariates")
      }
    }
    if (!'det.covs' %in% names(data)) {
      if (det.formula == ~ 1) {
        message("detection covariates (det.covs) not specified in data. Assuming interept only detection model.")
        data$det.covs <- list(int = matrix(1, dim(y)[2], dim(y)[3]))
      } else {
        stop("error: det.covs must be specified in data for a detection model with covariates")
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
    # Number of pseudoreplicates
    n.obs <- nrow(X.p)
    # Number of sites
    J <- nrow(X)
    # Number of repeat visits
    n.rep <- apply(y[1, , ], 1, function(a) sum(!is.na(a)))
    K.max <- max(n.rep)
    # Because I like K better than n.rep
    K <- n.rep
    if (missing(n.samples)) {
      stop("error: must specify number of MCMC samples")
    }
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
    names.long <- which(!is.na(c(y.big[1, , ])))
    if (nrow(X.p) == length(y) / N) {
      X.p <- X.p[!is.na(c(y.big[1, , ])), ]
    }
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

    # Separate out priors -------------------------------------------------
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
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(mu.alpha.comm) <- "double"
    storage.mode(Sigma.alpha.comm) <- "double"
    storage.mode(tau.beta.a) <- "double"
    storage.mode(tau.beta.b) <- "double"
    storage.mode(tau.alpha.a) <- "double"
    storage.mode(tau.alpha.b) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"

    ptm <- proc.time()

    # Run the model in C    
    out <- .Call("msPGOcc", y, X, X.p, p.occ, p.det, J, K, N, 
		 beta.starting, alpha.starting, z.starting,
		 beta.comm.starting, 
		 alpha.comm.starting, tau.beta.starting, 
		 tau.alpha.starting, z.long.indx, mu.beta.comm, 
		 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
		 tau.beta.a, tau.beta.b, tau.alpha.a, 
		 tau.alpha.b, n.samples, n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    out$beta.comm.samples <- mcmc(t(out$beta.comm.samples))
    colnames(out$beta.comm.samples) <- x.names
    out$alpha.comm.samples <- mcmc(t(out$alpha.comm.samples))
    colnames(out$alpha.comm.samples) <- x.p.names
    out$tau.beta.samples <- mcmc(t(out$tau.beta.samples))
    colnames(out$tau.beta.samples) <- x.names
    out$tau.alpha.samples <- mcmc(t(out$tau.alpha.samples))
    colnames(out$tau.alpha.samples) <- x.p.names
    if (is.null(sp.names)) {
      sp.names <- paste('sp', 1:N, sep = '')
    }
    coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- coef.names
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    coef.names.det <- paste(rep(x.p.names, each = N), sp.names, sep = '-')
    colnames(out$alpha.samples) <- coef.names.det
    out$z.samples <- array(out$z.samples, dim = c(N, J, n.samples))
    out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
    out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.samples))
    out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
    tmp <- array(NA, dim = c(N, J * K.max, n.samples))
    tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.samples))
    out$y.rep.samples <- array(tmp, dim = c(N, J, K.max, n.samples))
    out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
    out$X <- X
    out$X.p <- X.p
    out$y <- y.big
    out$call <- cl
    out$n.samples <- n.samples
    out$x.names <- x.names
    out$sp.names <- sp.names
    out$x.p.names <- x.p.names

    class(out) <- "msPGOcc"
    
    out
}
