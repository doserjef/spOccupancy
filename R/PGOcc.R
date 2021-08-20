#' Function for Fitting Single Species Occupancy Models Using Polya-Gamma Latent Variables
#'
#' @param y Data
#' @param X Occupancy covariates
#' @param X.p Detection covariates
#' @param starting Starting values
#' @param n.rep Temporal replicates
#' @param n.samples MCMC samples
#' @param priors Priors
#' @param n.omp.threads Number of threads
#' @param verbose Be wordy
#' @param n.report Value to report
#' @param ... currently no additional arguments
#'
#' @return stuff
#' @export
#'
#' @examples
PGOcc <- function(y, X, X.p, starting, n.rep, n.samples,
		  priors, n.omp.threads = 1, verbose = TRUE,
		  n.report = 100, ...){

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
    # Number of sites
    J <- nrow(y)
    if (missing(X)) {
      message("X is not specified. Assuming intercept only occupancy model.\n")
      X <- matrix(1, J, 1)
    }
    if (missing(X.p)) {
      message("X.p is not specified. Assuming intercept only detection model.\n")
      X.p <- array(1, dim = c(J, dim(y)[2], 1))
    }

    # Get basic info from inputs ------------------------------------------
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

    ptm <- proc.time()

    # Run the model in C
    out <- .Call("PGOcc", y, X, X.p, p.occ, p.det, J, K,
		 beta.starting, alpha.starting, z.starting,
		 z.long.indx, mu.beta, mu.alpha,
		 Sigma.beta, Sigma.alpha, n.samples,
		 n.omp.threads, verbose, n.report)

    out$run.time <- proc.time() - ptm

    out$beta.samples <- mcmc(t(out$beta.samples))
    colnames(out$beta.samples) <- x.names
    out$alpha.samples <- mcmc(t(out$alpha.samples))
    colnames(out$alpha.samples) <- x.p.names
    out$z.samples <- mcmc(t(out$z.samples))
    out$psi.samples <- mcmc(t(out$psi.samples))
    tmp <- array(NA, dim = c(J * K.max, n.samples))
    tmp[names.long, ] <- out$y.rep.samples
    out$y.rep.samples <- array(tmp, dim = c(J, K.max, n.samples))
    out$X.occ <- X
    out$X.p <- X.p.orig
    out$y <- y.big
    out$call <- cl

    # class(out) <- "PGOcc"

    out
}
