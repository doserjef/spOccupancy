msPGOcc <- function(occ.formula, det.formula, data, starting, n.samples, 
		    priors, n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		    n.burn = round(.10 * n.samples), n.thin = 1, ...){
 
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
      X.re <- as.matrix(tmp[[4]])
      x.re.names <- colnames(X.re)
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
      X.p.re <- X.p.re[!is.na(y), , drop = FALSE]
    }
    y <- y[!is.na(y)]

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
    lambda.psi <- matrix(0, J, n.occ.re)
    lambda.p <- matrix(0, n.obs, n.det.re)
    if (p.occ.re > 0) {
      for (i in 1:n.occ.re) {
        lambda.psi[which(X.re == (i - 1), arr.ind = TRUE)[, 1], i] <- 1
      }
    }
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

    # sigma.sq.psi -------------------
    if (p.occ.re > 0) {
      if ("sigma.sq.psi" %in% names(starting)) {
        sigma.sq.psi.starting <- starting[["sigma.sq.psi"]]
        if (length(sigma.sq.psi.starting) != p.occ.re) {
          stop(paste("error: starting values for sigma.sq.psi must be of length ", p.occ.re, 
		     sep = ""))
        }
      } else {
        sigma.sq.psi.starting <- rep(1, p.occ.re)
        if (verbose) {
          message("sigma.sq.psi is not specified in starting values. Setting starting value to 1\n")
        }
      }
      beta.star.indx <- rep(0:(p.occ.re - 1), n.occ.re.long)
      beta.star.starting <- rnorm(n.occ.re, sqrt(sigma.sq.psi.starting[beta.star.indx + 1]))
      # Starting values for all species 
      beta.star.starting <- rep(beta.star.starting, N)
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
      if (verbose) {	    
        message("No prior specified for tau.beta.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
      }
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
      if (verbose) {
        message("No prior specified for tau.alpha.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
      }
      tau.alpha.a <- rep(0.1, p.det)
      tau.alpha.b <- rep(0.1, p.det)
    }

    # sigma.sq.psi --------------------
    if (p.occ.re > 0) {
      if ("sigma.sq.psi.ig" %in% names(priors)) {
        sigma.sq.psi.a <- priors$sigma.sq.psi.ig[[1]]
        sigma.sq.psi.b <- priors$sigma.sq.psi.ig[[2]]
        if (!is.list(priors$sigma.sq.psi.ig) | length(priors$sigma.sq.psi.ig) != 2) {
          stop("error: sigma.sq.psi.ig must be a list of length 2")
        }
        if (length(sigma.sq.psi.a) != p.occ.re) {
          stop(paste("error: sigma.sq.psi.ig[[1]] must be a vector of length ", 
          	   p.occ.re, " with elements corresponding to sigma.sq.psis' shape", sep = ""))
        }
        if (length(sigma.sq.psi.b) != p.occ.re) {
          stop(paste("error: sigma.sq.psi.ig[[2]] must be a vector of length ", 
          	   p.occ.re, " with elements corresponding to sigma.sq.psis' scale", sep = ""))
        }
    }   else {
        if (verbose) {	    
          message("No prior specified for sigma.sq.psi.ig. Setting prior shape to 0.1 and prior scale to 0.1\n")
        }
        sigma.sq.psi.a <- rep(0.1, p.occ.re)
        sigma.sq.psi.b <- rep(0.1, p.occ.re)
      }
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
    storage.mode(n.burn) <- "integer"
    storage.mode(n.thin) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"

    if (p.occ.re > 0 & p.det.re == 0) {
      storage.mode(p.occ.re) <- "integer"
      storage.mode(X.re) <- "integer"
      storage.mode(n.occ.re) <- "integer"
      storage.mode(n.occ.re.long) <- "integer"
      storage.mode(sigma.sq.psi.starting) <- "double"
      storage.mode(sigma.sq.psi.a) <- "double"
      storage.mode(sigma.sq.psi.b) <- "double"
      storage.mode(beta.star.starting) <- "double"
      storage.mode(beta.star.indx) <- "integer"
      storage.mode(lambda.psi) <- "double"

      ptm <- proc.time()

      out <- .Call("msPGOccREOcc", y, X, X.p, X.re, 
		   lambda.psi, p.occ, p.det, p.occ.re, 
		   J, K, N, n.occ.re, n.occ.re.long, 
          	   beta.starting, alpha.starting, z.starting,
          	   beta.comm.starting, 
          	   alpha.comm.starting, tau.beta.starting, 
          	   tau.alpha.starting, sigma.sq.psi.starting, 
		   beta.star.starting, 
		   z.long.indx, beta.star.indx, mu.beta.comm, 
          	   mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	   tau.beta.a, tau.beta.b, tau.alpha.a, 
          	   tau.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
		   n.samples, n.omp.threads, 
		   verbose, n.report, n.burn, n.thin, n.post.samples)

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
      out$sigma.sq.psi.samples <- mcmc(t(out$sigma.sq.psi.samples))
      colnames(out$sigma.sq.psi.samples) <- x.re.names
      out$beta.star.samples <- mcmc(t(out$beta.star.samples))
      tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
      beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
      beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      tmp <- array(NA, dim = c(N, J * K.max, n.post.samples))
      tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.post.samples))
      out$y.rep.samples <- array(tmp, dim = c(N, J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
      out$X <- X
      out$X.p <- X.p
      out$X.re <- X.re
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE
      out$psiRE <- TRUE

      class(out) <- "msPGOcc"

    }

    if (p.occ.re == 0 & p.det.re > 0) {
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

      ptm <- proc.time()

      out <- .Call("msPGOccREDet", y, X, X.p, X.p.re, 
		   lambda.p, p.occ, p.det, p.det.re, 
		   J, K, N, n.det.re, n.det.re.long,
          	   beta.starting, alpha.starting, z.starting,
          	   beta.comm.starting, 
          	   alpha.comm.starting, tau.beta.starting, 
          	   tau.alpha.starting,  
		   sigma.sq.p.starting, alpha.star.starting, 
		   z.long.indx,
		   alpha.star.indx, mu.beta.comm, 
          	   mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	   tau.beta.a, tau.beta.b, tau.alpha.a, 
          	   tau.alpha.b, sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads, 
		   verbose, n.report, n.burn, n.thin, n.post.samples)

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
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      tmp <- array(NA, dim = c(N, J * K.max, n.post.samples))
      tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.post.samples))
      out$y.rep.samples <- array(tmp, dim = c(N, J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
      out$X <- X
      out$X.p <- X.p
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE
      out$psiRE <- FALSE

      class(out) <- "msPGOcc"

    }

    if (p.occ.re > 0 & p.det.re > 0) {

      # Testing
      #alpha.star.starting <- c(dat$alpha.star)
      #sigma.sq.p.starting <- p.RE$sigma.sq.p 
      #alpha.starting <- alpha.true

      storage.mode(p.occ.re) <- "integer"
      storage.mode(p.det.re) <- "integer"
      storage.mode(X.re) <- "integer"
      storage.mode(X.p.re) <- "integer"
      storage.mode(n.occ.re) <- "integer"
      storage.mode(n.det.re) <- "integer"
      storage.mode(n.occ.re.long) <- "integer"
      storage.mode(n.det.re.long) <- "integer"
      storage.mode(sigma.sq.psi.starting) <- "double"
      storage.mode(sigma.sq.p.starting) <- "double"
      storage.mode(sigma.sq.psi.a) <- "double"
      storage.mode(sigma.sq.psi.b) <- "double"
      storage.mode(sigma.sq.p.a) <- "double"
      storage.mode(sigma.sq.p.b) <- "double"
      storage.mode(beta.star.starting) <- "double"
      storage.mode(beta.star.indx) <- "integer"
      storage.mode(alpha.star.starting) <- "double"
      storage.mode(alpha.star.indx) <- "integer"
      storage.mode(lambda.psi) <- "double"
      storage.mode(lambda.p) <- "double"

      ptm <- proc.time()

      out <- .Call("msPGOccREBoth", y, X, X.p, X.re, X.p.re, 
		   lambda.psi, lambda.p, p.occ, p.det, p.occ.re, p.det.re, 
		   J, K, N, n.occ.re, n.det.re, n.occ.re.long, n.det.re.long,
          	   beta.starting, alpha.starting, z.starting,
          	   beta.comm.starting, 
          	   alpha.comm.starting, tau.beta.starting, 
          	   tau.alpha.starting, sigma.sq.psi.starting, 
		   sigma.sq.p.starting, beta.star.starting, 
		   alpha.star.starting, z.long.indx, beta.star.indx, 
		   alpha.star.indx, mu.beta.comm, 
          	   mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	   tau.beta.a, tau.beta.b, tau.alpha.a, 
          	   tau.alpha.b, sigma.sq.psi.a, sigma.sq.psi.b, 
		   sigma.sq.p.a, sigma.sq.p.b, n.samples, n.omp.threads, 
		   verbose, n.report, n.burn, n.thin, n.post.samples)

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
      out$sigma.sq.psi.samples <- mcmc(t(out$sigma.sq.psi.samples))
      colnames(out$sigma.sq.psi.samples) <- x.re.names
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$beta.star.samples <- mcmc(t(out$beta.star.samples))
      tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
      beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
      beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.occ.re), sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
      out$psi.samples <- array(out$psi.samples, dim = c(N, J, n.post.samples))
      out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
      tmp <- array(NA, dim = c(N, J * K.max, n.post.samples))
      tmp[, names.long, ] <- array(out$y.rep.samples, dim = c(N, n.obs, n.post.samples))
      out$y.rep.samples <- array(tmp, dim = c(N, J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(4, 1, 2, 3))
      out$X <- X
      out$X.p <- X.p
      out$X.re <- X.re
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$y <- y.big
      out$call <- cl
      out$n.samples <- n.samples
      out$x.names <- x.names
      out$sp.names <- sp.names
      out$x.p.names <- x.p.names
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE
      out$psiRE <- TRUE

      class(out) <- "msPGOcc"

    }

    if (p.occ.re == 0 & p.det.re == 0) {
      ptm <- proc.time()

      # Run the model in C    
      out <- .Call("msPGOcc", y, X, X.p, p.occ, p.det, J, K, N, 
          	 beta.starting, alpha.starting, z.starting,
          	 beta.comm.starting, 
          	 alpha.comm.starting, tau.beta.starting, 
          	 tau.alpha.starting, z.long.indx, mu.beta.comm, 
          	 mu.alpha.comm, Sigma.beta.comm, Sigma.alpha.comm, 
          	 tau.beta.a, tau.beta.b, tau.alpha.a, 
          	 tau.alpha.b, n.samples, n.omp.threads, verbose, n.report, 
          	 n.burn, n.thin, n.post.samples)

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
      out$z.samples <- array(out$z.samples, dim = c(N, J, n.post.samples))
      out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
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
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE
      out$psiRE <- FALSE

      class(out) <- "msPGOcc"

    }

    
    out
}
