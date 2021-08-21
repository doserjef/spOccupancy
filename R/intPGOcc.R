intPGOcc <- function(y, X, X.p, sites, n.data, starting, n.rep, n.samples, 
		     priors, n.omp.threads = 1, verbose = TRUE, 
		     n.report = 1000, ...){
 
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
    if (missing(sites)) {
      stop("error: site ids must be specified for each data source.\n")
    }
    # Number of sites with at least one data source. 
    J <- length(unique(unlist(sites)))
    # Number of sites for each data set
    J.long <- sapply(y, function(a) dim(a)[[1]])
    if (missing(X)) {
      message("X is not specified. Assuming intercept only occupancy model.\n")
      X <- matrix(1, J, 1)
    }
    # Total number of sites
    J.all <- nrow(X)
    if (missing(n.data)) {
      message(paste("n.data is not specified. Assuming a total of ", length(y), " datasets.\n", 
	      sep = ''))
      n.data <- length(y)
    }
    if (missing(X.p)) { 
      message("X.p is not specified. Assuming intercept only detection model for all data sources.\n")
      X.p <- list()
      for (q in 1:n.data) {
        X.p[[q]] <- array(1, dim = c(J.long[[q]], dim(y[[q]])[2], 1))
      }
    }
    if (length(X.p) != n.data | length(y) != n.data) {
      stop(paste("error: y and X.p must be lists of length ", n.data, ".", sep = ''))
    }

    # Extract data from inputs --------------------------------------------
    # Number of occupancy parameters 
    p.occ <- ncol(X)
    # Number of detection parameters for each data set
    p.det.long <- sapply(X.p, function(a) dim(a)[[3]])
    # Total number of detection parameters
    p.det <- sum(p.det.long)
    if (missing(n.rep)) {
      stop("error: number of replicates (n.rep) must be specified")
    }
    if (length(n.rep) != n.data | !is.list(n.rep)) {
      stop(paste("error: n.rep must be a list of length ", n.data, ".", sep = ""))
    }
    for (q in 1:n.data) {
      if (length(n.rep[[q]]) == 1) {
        n.rep[[q]] <- rep(n.rep[[q]], J.long[[q]])
      } else if (length(n.rep[[q]]) != J.long[[q]]) {
        stop(paste("error: n.rep[[", q, "]] must be of length ", J.long[[q]], 
      	     " or 1", sep = ''))
      }
    } # q
    # Number of repeat visits for each data set
    K.long.max <- sapply(n.rep, max)
    # Number of repeat visits for each data set site. 
    K <- unlist(n.rep)

    # Get covariate names for later ---------------------------------------
    if (length(colnames(X)) == 0) {
      x.names <- paste('x', 1:p.occ, sep = '')
    } else {
      x.names <- colnames(X)
    }
    x.p.names <- rep(NA, p.det)
    curr.indx <- 1
    for (i in 1:n.data) {
      if (length(colnames(X.p[[i]])) == 0) {
        x.p.names[curr.indx:(curr.indx + p.det.long[i] - 1)] <- 
	  paste('x.p', 1:p.det.long[i], '-D', i, sep = '')
      } else {
        x.p.names[curr.indx:(curr.indx + p.det.long[i] - 1)] <- 
	  paste(colnames(X.p), '-D', i, sep = '')
      }
      curr.indx <- curr.indx + p.det.long[i]
    } # i

    # Get detection covariates into n.obs x p.det matrices ----------------
    X.p.orig <- X.p
    y.big <- y
    X.p <- lapply(X.p, function(a) matrix(a, dim(a)[1] * dim(a)[2], dim(a)[3]))
    for (i in 1:n.data) {
      rownames(X.p[[i]]) <- 1:nrow(X.p[[i]])
      tmp <- c(apply(y[[i]], c(1, 2), function(a) sum(is.na(a))))
      X.p[[i]] <- X.p[[i]][tmp == 0, , drop = FALSE]
    }
    names.long <- lapply(X.p, function(a) as.numeric(rownames(a)))

    # Get data all organized for C ----------------------------------------
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
      message("alpha is not specified in starting values. Setting starting value to 0\n")
      alpha.starting <- rep(0, p.det)
    }

    alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
    alpha.indx.c <- alpha.indx.r - 1

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
      message("No prior specified for alpha.normal. Setting prior mean to 0 and prior variance to 2.73\n")
      mu.alpha <- rep(0, p.det)
      sigma.alpha <- rep(2.73, p.det) 
    }

    # Specify storage modes -----------------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.starting) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
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
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(data.indx.c) <- "integer"
    storage.mode(alpha.indx.c) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(sigma.alpha) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"

    ptm <- proc.time()

    # Run the model in C    
    out <- .Call("intPGOcc", y, X, X.p.all, p.occ, p.det, p.det.long, 
		 J, J.long, K, n.obs, n.obs.long, n.data, 
		 beta.starting, alpha.starting, z.starting,
		 z.long.indx.c, data.indx.c, alpha.indx.c, mu.beta, mu.alpha, 
		 Sigma.beta, sigma.alpha, n.samples, 
		 n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- x.names
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    colnames(out$alpha.samples) <- x.p.names
    out$z.samples <- mcmc(t(out$z.samples))
    out$psi.samples <- mcmc(t(out$psi.samples))
    # y.rep.samples is returned as a list, where each element 
    # corresponds to a different data set. 
    tmp <- list()
    indx <- 1
    for (q in 1:n.data) {
      tmp[[q]] <- array(NA, dim = c(J.long[q] * K.long.max[q], n.samples))
      tmp[[q]][names.long[[q]], ] <- out$y.rep.samples[indx:(indx + n.obs.long[q] - 1), ] 
      tmp[[q]] <- array(tmp[[q]], dim = c(J.long[q], K.long.max[q], n.samples))
      indx <- indx + n.obs.long[q]
    }
    out$y.rep.samples <- tmp
    out$X.occ <- X
    out$X.p <- X.p.orig
    out$y <- y.big
    out$call <- cl

    # class(out) <- "PGOcc"
    
    out
}
