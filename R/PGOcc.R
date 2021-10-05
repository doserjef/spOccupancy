PGOcc <- function(occ.formula, det.formula, data, starting, priors, 
		  n.samples, n.omp.threads = 1, verbose = TRUE,
		  n.report = 100, n.burn = round(.10 * n.samples), n.thin = 1, 
		  k.fold, k.fold.threads = 1, k.fold.seed = 100, ...){

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
    y <- as.matrix(data$y)
    if (!'occ.covs' %in% names(data)) {
      if (occ.formula == ~ 1) {
        if (verbose) {
          message("occupancy covariates (occ.covs) not specified in data. Assuming intercept only occupancy model.")
        }
        data$occ.covs <- matrix(1, dim(y)[1], 1)
      } else {
        stop("error: occ.covs must be specified in data for an occupancy model with covariates")
      }
    }
    if (!'det.covs' %in% names(data)) {
      if (det.formula == ~ 1) {
        if (verbose) {
          message("detection covariates (det.covs) not specified in data. Assuming interept only detection model.")
	}
        data$det.covs <- list(int = matrix(1, dim(y)[1], dim(y)[2]))
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

    # Get basic info from inputs ------------------------------------------
    # Number of sites
    J <- nrow(y)
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
    # Number of replicates at each site
    n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    # Max number of repeat visits
    K.max <- max(n.rep)
    # Because I like K better than n.rep
    K <- n.rep

    # Get indices to map z to y -------------------------------------------
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y))]
    # Subtract 1 for indices in C
    z.long.indx <- z.long.indx - 1
    # y is stored in the following order: species, site, visit
    y.big <- y
    y <- c(y)
    names.long <- which(!is.na(y))
    # Remove missing observations when the covariate data are available but
    # there are missing detection-nondetection data. 
    if (nrow(X.p) == length(y)) {
      X.p <- X.p[!is.na(y), , drop = FALSE]  
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
      if (verbose) {
        message("No prior specified for alpha.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      }
      mu.alpha <- rep(0, p.det)
      Sigma.alpha <- diag(p.det) * 2.73
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
    storage.mode(beta.starting) <- "double"
    storage.mode(alpha.starting) <- "double"
    storage.mode(z.long.indx) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(Sigma.alpha) <- "double"
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
    # Set model.deviance to NA for returning when no cross-validation
    model.deviance <- NA

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
      

      out <- .Call("PGOccREOcc", y, X, X.p, X.re,  
		   lambda.psi, p.occ, p.det, 
		   p.occ.re, J, K, n.occ.re,  
		   n.occ.re.long, beta.starting, 
		   alpha.starting, sigma.sq.psi.starting, beta.star.starting, 
		   z.starting, z.long.indx, beta.star.indx, 
		   mu.beta, mu.alpha, Sigma.beta, Sigma.alpha, 
		   sigma.sq.psi.a, sigma.sq.psi.b, n.samples,
          	   n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	   n.post.samples)


      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$z.samples <- mcmc(t(out$z.samples))
      out$psi.samples <- mcmc(t(out$psi.samples))
      tmp <- array(NA, dim = c(J * K.max, n.post.samples))
      tmp[names.long, ] <- out$y.rep.samples
      out$y.rep.samples <- array(tmp, dim = c(J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(3, 1, 2))
      out$sigma.sq.psi.samples <- mcmc(t(out$sigma.sq.psi.samples))
      colnames(out$sigma.sq.psi.samples) <- x.re.names
      out$beta.star.samples <- mcmc(t(out$beta.star.samples))
      tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
      beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$X <- X
      out$X.p <- X.p
      out$X.re <- X.re
      out$y <- y.big
      out$n.samples <- n.samples
      out$call <- cl
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE
      out$psiRE <- TRUE

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
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  z.starting.fit <- z.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  J.fit <- nrow(X.fit)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  lambda.psi.fit <- lambda.psi[-curr.set, , drop = FALSE]
	  lambda.psi.0 <- lambda.psi[curr.set, , drop = FALSE]
	  X.re.fit <- X.re[-curr.set, , drop = FALSE]
	  X.re.0 <- X.re[curr.set, , drop = FALSE]
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
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
          storage.mode(n.samples) <- "integer"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"
          storage.mode(p.occ.re) <- "integer"
          storage.mode(X.re.fit) <- "integer"
          storage.mode(n.occ.re) <- "integer"
          storage.mode(n.occ.re.long) <- "integer"
          storage.mode(sigma.sq.psi.starting) <- "double"
          storage.mode(sigma.sq.psi.a) <- "double"
          storage.mode(sigma.sq.psi.b) <- "double"
          storage.mode(beta.star.starting) <- "double"
          storage.mode(beta.star.indx) <- "integer"
          storage.mode(lambda.psi.fit) <- "double"

          # Run the model in C
          out.fit <- .Call("PGOccREOcc", y.fit, X.fit, X.p.fit, X.re.fit,  
		           lambda.psi.fit, p.occ, p.det, 
		           p.occ.re, J.fit, K.fit, n.occ.re,  
		           n.occ.re.long, beta.starting, 
		           alpha.starting, sigma.sq.psi.starting, beta.star.starting, 
		           z.starting.fit, z.long.indx.fit, beta.star.indx, 
		           mu.beta, mu.alpha, Sigma.beta, Sigma.alpha, 
		           sigma.sq.psi.a, sigma.sq.psi.b, n.samples,
          	           n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
          	           n.post.samples)

	  # Predict occurrence at new sites. 
          # Now can get the output
          psi.0.samples <- mcmc(logit.inv(t(X.0 %*% out.fit$beta.samples + 
					  lambda.psi.0 %*% out.fit$beta.star.samples)))
	  z.0.samples <- matrix(NA, n.post.samples, nrow(X.0))
	  for (j in 1:nrow(X.0)) {
            z.0.samples[, j] <- rbinom(n.post.samples, 1, psi.0.samples[, j])
          }
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples))
	  like.samples <- rep(NA, nrow(X.p.0))
	  for (j in 1:nrow(X.p.0)) {
            like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * z.0.samples[, z.0.long.indx[j]]))
          }
	  sum(log(like.samples))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }
      class(out) <- "PGOcc"
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
      
      out <- .Call("PGOccREDet", y, X, X.p, X.p.re, lambda.p, p.occ, p.det, 
		   p.det.re, J, K, n.det.re, n.det.re.long, beta.starting, 
		   alpha.starting, sigma.sq.p.starting, 
		   alpha.star.starting, z.starting, z.long.indx, alpha.star.indx, 
		   mu.beta, mu.alpha, Sigma.beta, Sigma.alpha, 
		   sigma.sq.p.a, sigma.sq.p.b, n.samples,
          	   n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	   n.post.samples)

      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$z.samples <- mcmc(t(out$z.samples))
      out$psi.samples <- mcmc(t(out$psi.samples))
      tmp <- array(NA, dim = c(J * K.max, n.post.samples))
      tmp[names.long, ] <- out$y.rep.samples
      out$y.rep.samples <- array(tmp, dim = c(J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(3, 1, 2))
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$X <- X
      out$X.p <- X.p
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$y <- y.big
      out$n.samples <- n.samples
      out$call <- cl
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE
      out$psiRE <- FALSE

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
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  z.starting.fit <- z.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
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
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
          storage.mode(n.samples) <- "integer"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
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

          # Run the model in C
          out.fit <- .Call("PGOccREDet", y.fit, X.fit, X.p.fit, X.p.re.fit,  
		           lambda.p.fit, p.occ, p.det, 
		           p.det.re, J.fit, K.fit, n.det.re,  
		           n.det.re.long, beta.starting, 
		           alpha.starting, sigma.sq.p.starting, alpha.star.starting, 
		           z.starting.fit, z.long.indx.fit, alpha.star.indx, 
		           mu.beta, mu.alpha, Sigma.beta, Sigma.alpha, 
		           sigma.sq.p.a, sigma.sq.p.b, n.samples,
          	           n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
          	           n.post.samples)

	  # Predict occurrence at new sites. 
          # Now can get the output
          psi.0.samples <- mcmc(logit.inv(t(X.0 %*% out.fit$beta.samples)))
	  z.0.samples <- matrix(NA, n.post.samples, nrow(X.0))
	  for (j in 1:nrow(X.0)) {
            z.0.samples[, j] <- rbinom(n.post.samples, 1, psi.0.samples[, j])
          }
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples + 
				     lambda.p.0 %*% out.fit$alpha.star.samples))
	  like.samples <- rep(NA, nrow(X.p.0))
	  for (j in 1:nrow(X.p.0)) {
            like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * z.0.samples[, z.0.long.indx[j]]))
          }
	  sum(log(like.samples))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "PGOcc"

    }

    if (p.occ.re > 0 & p.det.re > 0) {

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
      

      out <- .Call("PGOccREBoth", y, X, X.p, X.re, X.p.re, 
		   lambda.psi, lambda.p, p.occ, p.det, 
		   p.occ.re, p.det.re, J, K, n.occ.re, n.det.re, 
		   n.occ.re.long, n.det.re.long, beta.starting, 
		   alpha.starting, sigma.sq.psi.starting, 
		   sigma.sq.p.starting, beta.star.starting, 
		   alpha.star.starting, z.starting, z.long.indx, 
		   beta.star.indx, alpha.star.indx, 
		   mu.beta, mu.alpha, Sigma.beta, Sigma.alpha, 
		   sigma.sq.psi.a, sigma.sq.psi.b, 
		   sigma.sq.p.a, sigma.sq.p.b, n.samples,
          	   n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	   n.post.samples)

      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$z.samples <- mcmc(t(out$z.samples))
      out$psi.samples <- mcmc(t(out$psi.samples))
      tmp <- array(NA, dim = c(J * K.max, n.post.samples))
      tmp[names.long, ] <- out$y.rep.samples
      out$y.rep.samples <- array(tmp, dim = c(J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(3, 1, 2))
      out$sigma.sq.psi.samples <- mcmc(t(out$sigma.sq.psi.samples))
      colnames(out$sigma.sq.psi.samples) <- x.re.names
      out$sigma.sq.p.samples <- mcmc(t(out$sigma.sq.p.samples))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$beta.star.samples <- mcmc(t(out$beta.star.samples))
      tmp.names <- unlist(sapply(n.occ.re.long, function(a) 1:a))
      beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$alpha.star.samples <- mcmc(t(out$alpha.star.samples))
      tmp.names <- unlist(sapply(n.det.re.long, function(a) 1:a))
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$X <- X
      out$X.p <- X.p
      out$X.re <- X.re
      out$X.p.re <- X.p.re
      out$lambda.p <- lambda.p
      out$y <- y.big
      out$n.samples <- n.samples
      out$call <- cl
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- TRUE
      out$psiRE <- TRUE


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
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  # y.indx is for the actual model-fitting component. 
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  z.starting.fit <- z.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
	  J.fit <- nrow(X.fit)
	  K.fit <- K[-curr.set]
	  K.0 <- K[curr.set]
	  lambda.psi.fit <- lambda.psi[-curr.set, , drop = FALSE]
	  lambda.psi.0 <- lambda.psi[curr.set, , drop = FALSE]
	  X.re.fit <- X.re[-curr.set, , drop = FALSE]
	  X.re.0 <- X.re[curr.set, , drop = FALSE]
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
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
          storage.mode(n.samples) <- "integer"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
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
          storage.mode(sigma.sq.psi.starting) <- "double"
          storage.mode(sigma.sq.psi.a) <- "double"
          storage.mode(sigma.sq.psi.b) <- "double"
          storage.mode(beta.star.starting) <- "double"
          storage.mode(beta.star.indx) <- "integer"
          storage.mode(lambda.psi.fit) <- "double"

          # Run the model in C
          out.fit <- .Call("PGOccREBoth", y.fit, X.fit, X.p.fit, X.re.fit, X.p.re.fit, 
		           lambda.psi.fit, lambda.p.fit, p.occ, p.det, 
		           p.occ.re, p.det.re, J.fit, K.fit, n.occ.re, n.det.re, 
		           n.occ.re.long, n.det.re.long, beta.starting, 
		           alpha.starting, sigma.sq.psi.starting, 
		           sigma.sq.p.starting, beta.star.starting, 
		           alpha.star.starting, z.starting.fit, z.long.indx.fit, 
		           beta.star.indx, alpha.star.indx, 
		           mu.beta, mu.alpha, Sigma.beta, Sigma.alpha, 
		           sigma.sq.psi.a, sigma.sq.psi.b, 
		           sigma.sq.p.a, sigma.sq.p.b, n.samples,
          	           n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
          	           n.post.samples)

	  # Predict occurrence at new sites. 
          # Now can get the output
          psi.0.samples <- mcmc(logit.inv(t(X.0 %*% out.fit$beta.samples + 
					  lambda.psi.0 %*% out.fit$beta.star.samples)))
	  z.0.samples <- matrix(NA, n.post.samples, nrow(X.0))
	  for (j in 1:nrow(X.0)) {
            z.0.samples[, j] <- rbinom(n.post.samples, 1, psi.0.samples[, j])
          }
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples + 
				     lambda.p.0 %*% out.fit$alpha.star.samples))
	  like.samples <- rep(NA, nrow(X.p.0))
	  for (j in 1:nrow(X.p.0)) {
            like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * z.0.samples[, z.0.long.indx[j]]))
          }
	  sum(log(like.samples))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "PGOcc"
      
    }

    if (p.occ.re == 0 & p.det.re == 0) {

      # Run the model in C
      out <- .Call("PGOcc", y, X, X.p, p.occ, p.det, J, K,
          	 beta.starting, alpha.starting, z.starting,
          	 z.long.indx, mu.beta, mu.alpha,
          	 Sigma.beta, Sigma.alpha, n.samples,
          	 n.omp.threads, verbose, n.report, n.burn, n.thin, 
          	 n.post.samples)

      out$beta.samples <- mcmc(t(out$beta.samples))
      colnames(out$beta.samples) <- x.names
      out$alpha.samples <- mcmc(t(out$alpha.samples))
      colnames(out$alpha.samples) <- x.p.names
      out$z.samples <- mcmc(t(out$z.samples))
      out$psi.samples <- mcmc(t(out$psi.samples))
      tmp <- array(NA, dim = c(J * K.max, n.post.samples))
      tmp[names.long, ] <- out$y.rep.samples
      out$y.rep.samples <- array(tmp, dim = c(J, K.max, n.post.samples))
      out$y.rep.samples <- aperm(out$y.rep.samples, c(3, 1, 2))
      out$X <- X
      out$X.p <- X.p
      out$y <- y.big
      out$n.samples <- n.samples
      out$call <- cl
      out$n.post <- n.post.samples
      out$n.thin <- n.thin
      out$n.burn <- n.burn
      out$pRE <- FALSE
      out$psiRE <- FALSE

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
	model.deviance <- foreach (i = 1:k.fold, .combine = sum) %dopar% {
          curr.set <- sort(sites.k.fold[[i]])
	  y.indx <- !((z.long.indx + 1) %in% curr.set)
          y.fit <- y[y.indx]
	  y.0 <- y[!y.indx]
	  z.starting.fit <- z.starting[-curr.set]
	  X.p.fit <- X.p[y.indx, , drop = FALSE]
	  X.p.0 <- X.p[!y.indx, , drop = FALSE]
	  X.fit <- X[-curr.set, , drop = FALSE]
	  X.0 <- X[curr.set, , drop = FALSE]
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
          storage.mode(p.det) <- "integer"
          storage.mode(p.occ) <- "integer"
          storage.mode(J.fit) <- "integer"
          storage.mode(K.fit) <- "integer"
          storage.mode(beta.starting) <- "double"
          storage.mode(alpha.starting) <- "double"
          storage.mode(z.long.indx.fit) <- "integer"
          storage.mode(mu.beta) <- "double"
          storage.mode(Sigma.beta) <- "double"
          storage.mode(mu.alpha) <- "double"
          storage.mode(Sigma.alpha) <- "double"
          storage.mode(n.samples) <- "integer"
          storage.mode(n.omp.threads.fit) <- "integer"
          storage.mode(verbose.fit) <- "integer"
          storage.mode(n.report) <- "integer"
          storage.mode(n.burn) <- "integer"
          storage.mode(n.thin) <- "integer"

          # Run the model in C
          out.fit <- .Call("PGOcc", y.fit, X.fit, X.p.fit, p.occ, p.det, J.fit, K.fit,
              	           beta.starting, alpha.starting, z.starting.fit,
              	           z.long.indx.fit, mu.beta, mu.alpha,
              	           Sigma.beta, Sigma.alpha, n.samples,
              	           n.omp.threads.fit, verbose.fit, n.report, n.burn, n.thin, 
              	           n.post.samples)

	  # Predict occurrence at new sites. 
	  psi.0.samples <- logit.inv(t(X.0 %*% out.fit$beta.samples))
	  # Check this if you encounter a bug. 
	  z.0.samples <- matrix(NA, n.post.samples, nrow(X.0))
	  for (j in 1:nrow(X.0)) {
            z.0.samples[, j] <- rbinom(n.post.samples, 1, psi.0.samples[, j])
          }
	  # Detection 
	  p.0.samples <- logit.inv(t(X.p.0 %*% out.fit$alpha.samples))
	  like.samples <- rep(NA, nrow(X.p.0))
	  for (j in 1:nrow(X.p.0)) {
            like.samples[j] <- mean(dbinom(y.0[j], 1, p.0.samples[, j] * z.0.samples[, z.0.long.indx[j]]))
          }
	  sum(log(like.samples))
        }
	model.deviance <- -2 * model.deviance
	# Return objects from cross-validation
	out$k.fold.deviance <- model.deviance
	stopImplicitCluster()
      }

      class(out) <- "PGOcc"
    }
    out$run.time <- proc.time() - ptm
    out
}
