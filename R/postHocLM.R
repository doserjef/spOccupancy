postHocLM <- function(formula, data, inits, priors, verbose = FALSE, 
		      n.report = 100, n.samples, n.chains = 1, ...){

    ptm <- proc.time()

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
    if (missing(formula)) {
      stop("error: formula must be specified")
    }
    if (!'y' %in% names(data)) {
      stop("error: data y must be specified in data")
    }
    if (!'covs' %in% names(data)) {
      stop("error: covs must be specified in data")
    }
    if (!is.matrix(data$covs) & !is.data.frame(data$covs)) {
      stop("error: covs must be a matrix or data frame")
    }
    data$covs <- as.data.frame(data$covs)

    # Check whether random effects are sent in as numeric, and
    # return error if they are. 
    if (!is.null(findbars(formula))) {
      re.names <- sapply(findbars(formula), all.vars)
      for (i in 1:length(re.names)) {
        if (is(data$covs[, re.names[i]], 'factor')) {
          stop(paste("error: random effect variable ", re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
        } 
        if (is(data$covs[, re.names[i]], 'character')) {
          stop(paste("error: random effect variable ", re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
        }
      }
    }

    # Checking missing values ---------------------------------------------
    # y -------------------------------
    if (sum(is.na(data$y)) > 0) {
      stop("error: missing values in data$y. Missing response values are not allowed in postHocLM")
    }
    # covs ------------------------
    if (sum(is.na(data$covs)) > 0) {
      stop("error: missing covariate values in data$covs. Missing covariate values are not allowed in postHocLM.")
    }

    # Formula -------------------------------------------------------------
    if (is(formula, 'formula')) {
      tmp <- parseFormula(formula, data$covs)
      X <- as.matrix(tmp[[1]])
      X.re <- as.matrix(tmp[[4]])
      x.re.names <- colnames(X.re)
      x.names <- tmp[[2]]
    } else {
      stop("error: formula is misspecified")
    }
    # Get RE level names
    re.level.names <- lapply(data$covs[, x.re.names, drop = FALSE],
			     function (a) sort(unique(a)))

    # Get basic info from inputs ------------------------------------------
    # Number of observations
    N <- ncol(data$y)
    # Number of samples first stage model is fit.
    n.samples.y <- nrow(data$y)
    if (!missing(n.samples)) {
      if (n.samples %% n.samples.y != 0) {
        stop("if specified, n.samples must be divisible by the number of samples the first stage model was fit.")
      } else {
        n.times <- n.samples / n.samples.y
	y <- matrix(NA, n.samples, N)
	for (i in 1:n.times) {
          row.indx <- ((i - 1) * n.samples.y + 1):(i * n.samples.y)
          y[row.indx, ] <- data$y
	}
      }
    } else {
      n.samples <- n.samples.y
      y <- data$y
    }
    # Number of occupancy parameters
    p <- ncol(X)
    # Number of occupancy random effect parameters
    p.re <- ncol(X.re)
    # Number of latent occupancy random effect values
    n.re <- length(unlist(apply(X.re, 2, unique)))
    n.re.long <- apply(X.re, 2, function(a) length(unique(a)))
    # y ORDER: iteration, then observation. 
    y <- c(t(y))
    n.obs <- n.samples * N

    # Get random effect matrices all set ----------------------------------
    if (p.re > 1) {
      for (j in 2:p.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
      }
    }
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
      if (length(mu.beta) != p & length(mu.beta) != 1) {
        if (p == 1) {
          stop(paste("error: beta.normal[[1]] must be a vector of length ",
	  	     p, " with elements corresponding to betas' mean", sep = ""))
	} else {
          stop(paste("error: beta.normal[[1]] must be a vector of length ",
	  	     p, " or 1 with elements corresponding to betas' mean", sep = ""))
        }
      }
      if (length(sigma.beta) != p & length(sigma.beta) != 1) {
        if (p == 1) {
          stop(paste("error: beta.normal[[2]] must be a vector of length ",
		   p, " with elements corresponding to betas' variance", sep = ""))
        } else {
          stop(paste("error: beta.normal[[2]] must be a vector of length ",
		   p, " or 1 with elements corresponding to betas' variance", sep = ""))
        }
      }
      if (length(sigma.beta) != p) {
        sigma.beta <- rep(sigma.beta, p)
      }
      if (length(mu.beta) != p) {
        mu.beta <- rep(mu.beta, p)
      }
      Sigma.beta <- sigma.beta * diag(p)
    } else {
      if (verbose) {
        message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 100\n")
      }
      mu.beta <- rep(0, p)
      sigma.beta <- rep(100, p)
      Sigma.beta <- diag(p) * 100
    }
    # tau.sq --------------------------
    if ("tau.sq.ig" %in% names(priors)) {
      if (!is.vector(priors$tau.sq.ig) | !is.atomic(priors$tau.sq.ig) | length(priors$tau.sq.ig) != 2) {
        stop("error: tau.sq.ig must be a vector of length 2 with elements corresponding to tau.sq's shape and scale parameters")
      }
      tau.sq.a <- priors$tau.sq.ig[1]
      tau.sq.b <- priors$tau.sq.ig[2]
    } else {
      if (verbose) {
        message("No prior specified for tau.sq.\nSetting the inverse-Gamma shape and scale parameter to 0.001.\n")
      }
      tau.sq.a <- 0.001
      tau.sq.b <- 0.001
    }
    # sigma.sq ------------------------
    if (p.re > 0) {
      if ("sigma.sq.ig" %in% names(priors)) {
        if (!is.list(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
          stop("error: sigma.sq.ig must be a list of length 2")
        }
        sigma.sq.a <- priors$sigma.sq.ig[[1]]
        sigma.sq.b <- priors$sigma.sq.ig[[2]]
        if (length(sigma.sq.a) != p.re & length(sigma.sq.a) != 1) {
          if (p.re == 1) {
          stop(paste("error: sigma.sq.ig[[1]] must be a vector of length ", 
          	   p.re, " with elements corresponding to sigma.sqs' shape", sep = ""))
	  } else {
          stop(paste("error: sigma.sq.ig[[1]] must be a vector of length ", 
          	   p.re, " or 1 with elements corresponding to sigma.sqs' shape", sep = ""))
          }
        }
        if (length(sigma.sq.b) != p.re & length(sigma.sq.b) != 1) {
          if (p.re == 1) {
            stop(paste("error: sigma.sq.ig[[2]] must be a vector of length ", 
          	   p.re, " with elements corresponding to sigma.sqs' scale", sep = ""))
	  } else {
            stop(paste("error: sigma.sq.ig[[2]] must be a vector of length ", 
          	   p.re, " or 1 with elements corresponding to sigma.sqs' scale", sep = ""))
          }
        }
	if (length(sigma.sq.a) != p.re) {
          sigma.sq.a <- rep(sigma.sq.a, p.re)
        }
	if (length(sigma.sq.b) != p.re) {
          sigma.sq.b <- rep(sigma.sq.b, p.re)
        }
    }   else {
        if (verbose) {	    
          message("No prior specified for sigma.sq.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
        }
        sigma.sq.a <- rep(0.1, p.re)
        sigma.sq.b <- rep(0.1, p.re)
      }
    } else {
      sigma.sq.a <- 0
      sigma.sq.b <- 0
    }

    # Starting values -----------------------------------------------------
    if (missing(inits)) {
      inits <- list()
    }
    names(inits) <- tolower(names(inits))
    # beta -----------------------
    if ("beta" %in% names(inits)) {
      beta.inits <- inits[["beta"]]
      if (length(beta.inits) != p & length(beta.inits) != 1) {
        if (p == 1) {
          stop(paste("error: initial values for beta must be of length ", p,
		     sep = ""))

        } else {
          stop(paste("error: initial values for beta must be of length ", p, " or 1",
	  	     sep = ""))
        }
      }
      if (length(beta.inits) != p) {
        beta.inits <- rep(beta.inits, p)
      }
    } else {
      beta.inits <- rnorm(p)
      if (verbose) {
        message('beta is not specified in initial values.\nSetting initial values to random standard normal values.\n')
      }
    }
    # tau.sq ------------------------
    if ("tau.sq" %in% names(inits)) {
      tau.sq.inits <- inits[["tau.sq"]]
      if (length(tau.sq.inits) != 1) {
        stop("error: initial values for tau.sq must be of length 1")
      }
    } else {
      tau.sq.inits <- runif(1, 0.5, 10)
      if (verbose) {
        message("tau.sq is not specified in initial values.\nSetting initial value to random value between 0.5 and 10.\n")
      }
    }
    # sigma.sq -------------------
    if (p.re > 0) {
      if ("sigma.sq" %in% names(inits)) {
        sigma.sq.inits <- inits[["sigma.sq"]]
        if (length(sigma.sq.inits) != p.re & length(sigma.sq.inits) != 1) {
          if (p.re == 1) {
            stop(paste("error: initial values for sigma.sq must be of length ", p.re, 
		       sep = ""))
	  } else {
            stop(paste("error: initial values for sigma.sq must be of length ", p.re, 
		       " or 1", sep = ""))
          }
        }
	if (length(sigma.sq.inits) != p.re) {
          sigma.sq.inits <- rep(sigma.sq.inits, p.re)
        }
      } else {
        sigma.sq.inits <- runif(p.re, 0.5, 10)
        if (verbose) {
          message("sigma.sq is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
        }
      }
      beta.star.indx <- rep(0:(p.re - 1), n.re.long)
      beta.star.inits <- rnorm(n.re, sqrt(sigma.sq.inits[beta.star.indx + 1]))
    } else {
      sigma.sq.inits <- 0
      beta.star.indx <- 0
      beta.star.inits <- 0
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

    curr.chain <- 1
    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    consts <- c(N, p, p.re, n.re)
    storage.mode(consts) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(tau.sq.inits) <- "double"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(tau.sq.a) <- "double"
    storage.mode(tau.sq.b) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    storage.mode(n.samples) <- "integer"
    # For occurrence random effects
    storage.mode(X.re) <- "integer"
    beta.level.indx <- sort(unique(c(X.re)))
    storage.mode(beta.level.indx) <- "integer"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(sigma.sq.a) <- "double"
    storage.mode(sigma.sq.b) <- "double"
    storage.mode(n.re.long) <- "integer"
    storage.mode(beta.star.inits) <- "double"
    storage.mode(beta.star.indx) <- "integer"

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        beta.inits <- rnorm(p)
        tau.sq.inits <- runif(1, 0.5, 10)
	if (p.re > 0) {
          sigma.sq.inits <- runif(p.re, 0.5, 10)
          beta.star.inits <- rnorm(n.re, sqrt(sigma.sq.inits[beta.star.indx + 1]))
	}
      }
      storage.mode(chain.info) <- "integer"
      # Run the model in C
      out.tmp[[i]] <- .Call("postHocLM", y, X, X.re, consts, n.re.long, 
			    beta.inits, tau.sq.inits, sigma.sq.inits, beta.star.inits, 
			    beta.star.indx, beta.level.indx, mu.beta, Sigma.beta, 
			    tau.sq.a, tau.sq.b, sigma.sq.a, sigma.sq.b, 
			    n.samples, verbose, n.report, chain.info)
      chain.info[1] <- chain.info[1] + 1
    } # i   
    # Calculate R-Hat ---------------
    out <- list()
    out$rhat <- list()
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable. 
      out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$beta.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      out$rhat$tau.sq <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$tau.sq.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      if (p.re > 0) {
      out$rhat$sigma.sq <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$sigma.sq.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      }
    } else {
      out$rhat$beta <- rep(NA, p)
      out$rhat$tau.sq.beta <- NA
      if (p.re > 0) {
        out$rhat$sigma.sq <- rep(NA, p.re)
      }
    }
    # Put everything into MCMC objects
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- x.names
    out$tau.sq.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$tau.sq.samples))))
    out$y.hat.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.hat.samples))))
    if (p.re > 0) {
      out$sigma.sq.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.samples))))
      colnames(out$sigma.sq.samples) <- x.re.names
      out$beta.star.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
      tmp.names <- unlist(re.level.names)
      beta.star.names <- paste(rep(x.re.names, n.re.long), tmp.names, sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$re.level.names <- re.level.names
    }
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$tau.sq <- effectiveSize(out$tau.sq.samples)
    if (p.re > 0) {
      out$ESS$sigma.sq <- effectiveSize(out$sigma.sq.samples)
    }
    out$X <- X
    out$X.re <- X.re
    out$y <- y
    out$n.samples <- n.samples
    out$call <- cl
    out$n.chains <- n.chains
    if (p.re > 0) {
      out$RE <- TRUE
    } else {
      out$RE <- FALSE
    }
    # Calculate Bayesian R2 -----------
    # Note that this is the conditional Bayes R2, not marginal.
    var.y.hat <- apply(out$y.hat.samples, 1, var) 
    if (p.re > 0) {
      var.resid <- out$tau.sq.samples
    } else {
      var.resid <- out$tau.sq.samples
    }
    out$bayes.R2 <- mcmc(var.y.hat / (var.y.hat + var.resid))
    
    class(out) <- "postHocLM"
    out$run.time <- proc.time() - ptm
    out
  }
