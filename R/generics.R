# PGOcc -------------------------------------------------------------------
predict.PGOcc <- function(object, X.0, sub.sample, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (class(object) != "PGOcc") {
    stop("error: requires an output object of class PGOcc\n")
  }

  # Get samples for composition sampling ----------------------------------
  n.samples <- object$n.samples
  if (missing(sub.sample)) {
    message("sub.sample is not specified. Using all posterior samples for prediction.")
    s.indx <- 1:n.samples
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
    if (!is.numeric(start) || start >= n.samples){ 
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){ 
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){ 
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  p.occ <- ncol(object$X)
  if (ncol(X.0) != p.occ) {
    stop(paste("error: X.0 must have ", p.occ, " columns\n", sep = ''))
  }

  # Composition sampling --------------------------------------------------
  beta.samples <- as.matrix(object$beta.samples[s.indx, , drop = FALSE])
  out <- list()
  out$psi.0.samples <- mcmc(logit.inv(t(as.matrix(X.0) %*% t(beta.samples))))
  out$z.0.samples <- mcmc(matrix(rbinom(length(out$psi.0.samples), 1, c(out$psi.0.samples)), 
		      nrow = n.samples, ncol = nrow(X.0)))
  out$sub.sample <- sub.sample
  out$s.indx <- s.indx
  out$call <- cl

  class(out) <- "predict.PGOcc"
  out
}

print.PGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

fitted.PGOcc <- function(object, ...) {
  return(object$y.rep.samples)
}

summary.PGOcc <- function(object, sub.sample, 
			  quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
			  digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.samples <- nrow(object$beta.samples)

  if (missing(sub.sample)) {
    s.indx <- 1:n.samples
    start <- 1
    end <- n.samples
    thin <- 1
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
    if (!is.numeric(start) || start >= n.samples){ 
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){ 
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){ 
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }

  cat("Chain sub.sample:\n")
  cat(paste("start = ",start,"\n", sep=""))
  cat(paste("end = ",end,"\n", sep=""))
  cat(paste("thin = ",thin,"\n", sep=""))
  cat(paste("sample size = ",length(s.indx),"\n\n", sep=""))
  
  # Occupancy
  cat("Occupancy: \n")
  print(noquote(round(t(apply(object$beta.samples[s.indx,, drop = FALSE], 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))
  cat("\n")
  # Detection
  cat("Detection: \n")
  print(noquote(round(t(apply(object$alpha.samples[s.indx,, drop = FALSE], 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))
}

summary.ppcOcc <- function(object, level, 
			   digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:", deparse(object$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")

  cat("Chain sub.sample:\n")
  cat(paste("start = ",object$start,"\n", sep=""))
  cat(paste("end = ",object$end,"\n", sep=""))
  cat(paste("thin = ",object$thin,"\n", sep=""))
  cat(paste("sample size = ",object$sample.size,"\n\n", sep=""))

  if (object$class %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc')) {
    cat("Bayesian p-value: ", mean(object$fit.y.rep > object$fit.y), "\n")
    cat("Fit statistic: ", object$fit.stat, "\n")
  } else {

    if (missing(level)) {
      stop("error: must specify level of parameters to display. Valid values are 'species', 'community', or 'both'")
    }
    if (tolower(level) == 'community') {
      cat("----------------------------------------\n");
      cat("\tCommunity Level\n");
      cat("----------------------------------------\n");
      cat("Bayesian p-value: ", mean(object$fit.y.rep > object$fit.y), "\n")
      cat("Fit statistic: ", object$fit.stat, "\n")
    }

    if (tolower(level) == 'species') {
      cat("----------------------------------------\n");
      cat("\tSpecies Level\n");
      cat("----------------------------------------\n");
      N <- ncol(object$fit.y)
      for (i in 1:N) {
        cat(paste(object$sp.names[i], " Bayesian p-value: ", 
		  mean(object$fit.y.rep[, i] > object$fit.y[, i]), "\n", sep = '')) 
      }
      cat("Fit statistic: ", object$fit.stat, "\n")
    }

    if (tolower(level) == 'both') {
      cat("----------------------------------------\n");
      cat("\tCommunity Level\n");
      cat("----------------------------------------\n");
      cat("Bayesian p-value: ", mean(object$fit.y.rep > object$fit.y), "\n")
      cat("\n")
      cat("----------------------------------------\n");
      cat("\tSpecies Level\n");
      cat("----------------------------------------\n");
      N <- ncol(object$fit.y)
      for (i in 1:N) {
        cat(paste(object$sp.names[i], " Bayesian p-value: ", 
		  mean(object$fit.y.rep[, i] > object$fit.y[, i]), "\n", sep = '')) 
      }
      cat("Fit statistic: ", object$fit.stat, "\n")
    }
  }

}

# spPGOcc -----------------------------------------------------------------

predict.spPGOcc <- function(object, X.0, coords.0, sub.sample, n.omp.threads = 1, 
			    verbose = TRUE, n.report = 100, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (class(object) != "spPGOcc") {
    stop("error: requires an output object of class spPGOcc\n")
  }

  # Get samples for predictive sampler ------------------------------------
  n.samples <- object$n.samples
  if (missing(sub.sample)) {
    message("sub.sample is not specified. Using all posterior samples for prediction.")
    s.indx <- 1:n.samples
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
    if (!is.numeric(start) || start >= n.samples){ 
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){ 
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){ 
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }


  # Data prep -------------------------------------------------------------
  X <- object$X
  coords <- object$coords 
  J <- nrow(X)
  p.occ <- ncol(X)
  theta.samples <- object$theta.samples
  beta.samples <- object$beta.samples
  w.samples <- object$w.samples
  n.neighbors <- object$n.neighbors
  cov.model.indx <- object$cov.model.indx
  type <- object$type

  # Sub-sample previous 
  theta.samples <- t(theta.samples[s.indx, , drop = FALSE])
  beta.samples <- t(beta.samples[s.indx, , drop = FALSE])
  w.samples <- t(w.samples[s.indx, , drop = FALSE])

  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))){
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  if (ncol(X.0) != ncol(X)){
    stop(paste("error: X.0 must have ", p.occ," columns\n"))
  }
  X.0 <- as.matrix(X.0)
  
  if (missing(coords.0)) {
    stop("error: coords.0 must be specified\n")
  }
  if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
    stop("error: coords.0 must be a data.frame or matrix\n")
  }
  if (!ncol(coords.0) == 2){
    stop("error: coords.0 must have two columns\n")
  }
  coords.0 <- as.matrix(coords.0)
  
  q <- nrow(X.0)

  if (type == 'GP') {
  
    obs.pred.D <- iDist(coords, coords.0)
    obs.D <- iDist(coords)
    
    storage.mode(obs.pred.D) <- "double"
    storage.mode(obs.D) <- "double"
    storage.mode(J) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(X.0) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(beta.samples) <- "double"
    storage.mode(theta.samples) <- "double"
    storage.mode(w.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    
    ptm <- proc.time()

    out <- .Call("spPGOccPredict", J, p.occ, X.0, q, obs.D, 
		 obs.pred.D, beta.samples, theta.samples, 
		 w.samples, n.samples, cov.model.indx, verbose, 
		 n.omp.threads, n.report)

    out$z.0.samples <- mcmc(t(out$z.0.samples))
    out$w.0.samples <- mcmc(t(out$w.0.samples))  
    out$psi.0.samples <- mcmc(t(out$psi.0.samples))
    out$run.time <- proc.time() - ptm
    out$call <- cl
    out$object.class <- class(object)

    class(out) <- "predict.spPGOcc"

    out

  } else { 
    # Get nearest neighbors 
    # nn2 is a function from RANN. 
    nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1 

    storage.mode(coords) <- "double"
    storage.mode(J) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(X.0) <- "double"
    storage.mode(coords.0) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(beta.samples) <- "double"
    storage.mode(theta.samples) <- "double"
    storage.mode(w.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(nn.indx.0) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    ptm <- proc.time()

    out <- .Call("spPGOccNNGPPredict", coords, J, p.occ, n.neighbors, 
                 X.0, coords.0, q, nn.indx.0, beta.samples, 
                 theta.samples, w.samples, n.samples, 
                 cov.model.indx, n.omp.threads, verbose, n.report)

    out$z.0.samples <- mcmc(t(out$z.0.samples))
    out$w.0.samples <- mcmc(t(out$w.0.samples))  
    out$psi.0.samples <- mcmc(t(out$psi.0.samples))
    out$run.time <- proc.time() - ptm
    out$call <- cl
    out$object.class <- class(object)

    class(out) <- "predict.spPGOcc"

    out

  }

}

print.spPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.spPGOcc <- function(object, sub.sample, 
			    quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
			    digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.samples <- nrow(object$beta.samples)

  if (missing(sub.sample)) {
    s.indx <- 1:n.samples
    start <- 1
    end <- n.samples
    thin <- 1
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
    if (!is.numeric(start) || start >= n.samples){ 
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){ 
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){ 
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }

  cat("Chain sub.sample:\n")
  cat(paste("start = ",start,"\n", sep=""))
  cat(paste("end = ",end,"\n", sep=""))
  cat(paste("thin = ",thin,"\n", sep=""))
  cat(paste("sample size = ",length(s.indx),"\n\n", sep=""))
  
  # Occupancy
  cat("Occupancy: \n")
  print(noquote(round(t(apply(object$beta.samples[s.indx,, drop = FALSE], 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))
  cat("\n")
  # Detection
  cat("Detection: \n")
  print(noquote(round(t(apply(object$alpha.samples[s.indx,, drop = FALSE], 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))

  cat("\n")
  # Covariance
  cat("Covariance: \n")
  print(noquote(round(t(apply(object$theta.samples[s.indx, , drop = FALSE], 2, 
			      function(x) quantile(x, prob = quantiles))), digits)))
}


fitted.spPGOcc <- function(object, ...) {
  return(object$y.rep.samples)
}

# msPGOcc -----------------------------------------------------------------

predict.msPGOcc <- function(object, X.0, sub.sample, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (class(object) != "msPGOcc") {
    stop("error: requires an output object of class msPGOcc\n")
  }

  # Get samples for composition sampling ----------------------------------
  n.samples <- object$n.samples
  if (missing(sub.sample)) {
    message("sub.sample is not specified. Using all posterior samples for prediction.")
    s.indx <- 1:n.samples
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
    if (!is.numeric(start) || start >= n.samples){ 
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){ 
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){ 
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  p.occ <- ncol(object$X)
  if (ncol(X.0) != p.occ) {
    stop(paste("error: X.0 must have ", p.occ, " columns\n", sep = ''))
  }

  # Composition sampling --------------------------------------------------
  N <- dim(object$y)[1]
  sp.indx <- rep(1:N, p.occ)
  beta.samples <- as.matrix(object$beta.samples[s.indx, , drop = FALSE])
  out <- list()
  out$psi.0.samples <- array(NA, dim = c(n.samples, N, nrow(X.0)))
  out$z.0.samples <- array(NA, dim = c(n.samples, N, nrow(X.0)))
  for (i in 1:N) {
    out$psi.0.samples[, i, ] <- logit.inv(t(as.matrix(X.0) %*% t(beta.samples[, sp.indx == i])))
    out$z.0.samples[, i, ] <- matrix(rbinom(n.samples * nrow(X.0), 1, 
					    c(out$psi.0.samples[, i, ])), 
				     nrow = n.samples, ncol = nrow(X.0))
  }
  out$sub.sample <- sub.sample
  out$s.indx <- s.indx
  out$call <- cl

  class(out) <- "predict.msPGOcc"
  out
}

print.msPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.msPGOcc <- function(object, sub.sample,
			    quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
			    level,
			    digits = max(3L, getOption("digits") - 3L), ...) {

  if (missing(level)) {
    stop("error: must specify level of parameters to display. Valid values are 'species', 'community', or 'both'")
  }

  print(object)

  n.samples <- nrow(object$beta.comm.samples)

  if (missing(sub.sample)) {
    s.indx <- 1:n.samples
    start <- 1
    end <- n.samples
    thin <- 1
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
    if (!is.numeric(start) || start >= n.samples){
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }

  cat("Chain sub.sample:\n")
  cat(paste("start = ",start,"\n", sep=""))
  cat(paste("end = ",end,"\n", sep=""))
  cat(paste("thin = ",thin,"\n", sep=""))
  cat(paste("sample size = ",length(s.indx),"\n\n", sep=""))

  if (tolower(level) == 'community') {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occupancy
    cat("Occupancy Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccupancy Variances: \n")
    print(noquote(round(t(apply(object$tau.beta.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.alpha.samples[s.indx,, drop = FALSE], 2,
			      function(x) quantile(x, prob=quantiles))), digits)))

  }

  if (tolower(level) == 'species') {
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occupancy: \n")
    print(noquote(round(t(apply(object$beta.samples[s.indx,, drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples[s.indx,, drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))

  }

  if (tolower(level) == 'both') {
    
    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occupancy
    cat("Occupancy Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccupancy Variances: \n")
    print(noquote(round(t(apply(object$tau.beta.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.alpha.samples[s.indx,, drop = FALSE], 2,
			      function(x) quantile(x, prob=quantiles))), digits)))

    cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occupancy: \n")
    print(noquote(round(t(apply(object$beta.samples[s.indx,, drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples[s.indx,, drop = FALSE], 2,
  			        function(x) quantile(x, prob=quantiles))), digits)))
  }

}

fitted.msPGOcc <- function(object, ...) {
  return(object$y.rep.samples)
}

# spMsPGOcc ---------------------------------------------------------------
summary.spMsPGOcc <- function(object, sub.sample,
			      quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
			      level,
			      digits = max(3L, getOption("digits") - 3L), ...) {

  if (missing(level)) {
    stop("error: must specify level of parameters to display. Valid values are 'species', 'community', or 'both'")
  }

  print(object)

  n.samples <- nrow(object$beta.comm.samples)

  if (missing(sub.sample)) {
    s.indx <- 1:n.samples
    start <- 1
    end <- n.samples
    thin <- 1
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
    if (!is.numeric(start) || start >= n.samples){
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }

  cat("Chain sub.sample:\n")
  cat(paste("start = ",start,"\n", sep=""))
  cat(paste("end = ",end,"\n", sep=""))
  cat(paste("thin = ",thin,"\n", sep=""))
  cat(paste("sample size = ",length(s.indx),"\n\n", sep=""))

  if (tolower(level) == 'community') {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occupancy
    cat("Occupancy Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccupancy Variances: \n")
    print(noquote(round(t(apply(object$tau.beta.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.alpha.samples[s.indx,, drop = FALSE], 2,
			      function(x) quantile(x, prob=quantiles))), digits)))

  }

  if (tolower(level) == 'species') {
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occupancy: \n")
    print(noquote(round(t(apply(object$beta.samples[s.indx,, drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples[s.indx,, drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))

    cat("\n")
    # Covariance
    cat("Covariance: \n")
    print(noquote(round(t(apply(object$theta.samples[s.indx, , drop = FALSE], 2, 
  			        function(x) quantile(x, prob = quantiles))), digits)))

  }

  if (tolower(level) == 'both') {
    
    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occupancy
    cat("Occupancy Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccupancy Variances: \n")
    print(noquote(round(t(apply(object$tau.beta.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples[s.indx,, drop = FALSE], 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.alpha.samples[s.indx,, drop = FALSE], 2,
			      function(x) quantile(x, prob=quantiles))), digits)))

    cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occupancy: \n")
    print(noquote(round(t(apply(object$beta.samples[s.indx,, drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples[s.indx,, drop = FALSE], 2,
  			        function(x) quantile(x, prob=quantiles))), digits)))

    cat("\n")
    # Covariance
    cat("Covariance: \n")
    print(noquote(round(t(apply(object$theta.samples[s.indx, , drop = FALSE], 2, 
  			        function(x) quantile(x, prob = quantiles))), digits)))
  }

}

fitted.spMsPGOcc <- function(object, ...) {
  return(object$y.rep.samples)
}

print.spMsPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}


predict.spMsPGOcc <- function(object, X.0, coords.0, sub.sample, n.omp.threads = 1, 
			      verbose = TRUE, n.report = 100, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (class(object) != "spMsPGOcc") {
    stop("error: requires an output object of class spMsPGOcc\n")
  }

  # Get samples for predictive sampler ------------------------------------
  n.samples <- object$n.samples
  if (missing(sub.sample)) {
    message("sub.sample is not specified. Using all posterior samples for prediction.")
    s.indx <- 1:n.samples
  } else {
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)   
    if (!is.numeric(start) || start >= n.samples){ 
      stop("invalid start")
    }
    if (!is.numeric(end) || end > n.samples){ 
      stop("invalid end")
    }
    if (!is.numeric(thin) || thin >= n.samples){ 
      stop("invalid thin")
    }
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    n.samples <- length(s.indx)
  }


  # Data prep -------------------------------------------------------------
  X <- object$X
  y <- object$y
  coords <- object$coords 
  J <- nrow(X)
  N <- dim(y)[1]
  p.occ <- ncol(X)
  theta.samples <- object$theta.samples
  beta.samples <- object$beta.samples
  w.samples <- object$w.samples
  n.neighbors <- object$n.neighbors
  cov.model.indx <- object$cov.model.indx
  type <- object$type

  # Sub-sample previous 
  theta.samples <- t(theta.samples[s.indx, , drop = FALSE])
  beta.samples <- t(beta.samples[s.indx, , drop = FALSE])
  w.samples <- aperm(w.samples[s.indx, , , drop = FALSE], c(2, 3, 1))

  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))){
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  if (ncol(X.0) != ncol(X)){
    stop(paste("error: X.0 must have ", p.occ," columns\n"))
  }
  X.0 <- as.matrix(X.0)
  
  if (missing(coords.0)) {
    stop("error: coords.0 must be specified\n")
  }
  if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
    stop("error: coords.0 must be a data.frame or matrix\n")
  }
  if (!ncol(coords.0) == 2){
    stop("error: coords.0 must have two columns\n")
  }
  coords.0 <- as.matrix(coords.0)
  
  q <- nrow(X.0)

  if (type == 'GP') {
  
    obs.pred.D <- iDist(coords, coords.0)
    obs.D <- iDist(coords)
    
    storage.mode(obs.pred.D) <- "double"
    storage.mode(obs.D) <- "double"
    storage.mode(J) <- "integer"
    storage.mode(N) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(X.0) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(beta.samples) <- "double"
    storage.mode(theta.samples) <- "double"
    storage.mode(w.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    
    ptm <- proc.time()

    out <- .Call("spMsPGOccPredict", J, N, p.occ, X.0, q, obs.D, 
		 obs.pred.D, beta.samples, theta.samples, 
		 w.samples, n.samples, cov.model.indx, verbose, 
		 n.omp.threads, n.report)

    out$z.0.samples <- array(out$z.0.samples, dim = c(N, q, n.samples))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$w.0.samples <- array(out$w.0.samples, dim = c(N, q, n.samples))
    out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, q, n.samples))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
    out$run.time <- proc.time() - ptm
    out$call <- cl
    out$object.class <- class(object)

    class(out) <- "predict.spMsPGOcc"

    out

  } else { 
    # Get nearest neighbors 
    # nn2 is a function from RANN. 
    nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1 

    storage.mode(coords) <- "double"
    storage.mode(N) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(X.0) <- "double"
    storage.mode(coords.0) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(beta.samples) <- "double"
    storage.mode(theta.samples) <- "double"
    storage.mode(w.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(nn.indx.0) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    ptm <- proc.time()

    out <- .Call("spMsPGOccNNGPPredict", coords, J, N, p.occ, n.neighbors, 
                 X.0, coords.0, q, nn.indx.0, beta.samples, 
                 theta.samples, w.samples, n.samples, 
                 cov.model.indx, n.omp.threads, verbose, n.report)

    out$z.0.samples <- array(out$z.0.samples, dim = c(N, q, n.samples))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$w.0.samples <- array(out$w.0.samples, dim = c(N, q, n.samples))
    out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, q, n.samples))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
    out$run.time <- proc.time() - ptm
    out$call <- cl
    out$object.class <- class(object)

    class(out) <- "predict.spPGOcc"

    out

  }

}
