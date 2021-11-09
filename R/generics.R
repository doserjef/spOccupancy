# PGOcc -------------------------------------------------------------------
predict.PGOcc <- function(object, X.0, ...) {
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

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  p.occ <- ncol(object$X)
  p.design <- p.occ
  if (object$psiRE) {
    p.design <- p.occ + ncol(object$sigma.sq.psi.samples)
  }
  if (ncol(X.0) != p.design) {
    stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
  }

  # Composition sampling --------------------------------------------------
  beta.samples <- as.matrix(object$beta.samples)
  n.post <- object$n.post
  out <- list()
  if (object$psiRE) {
    beta.star.samples <- as.matrix(object$beta.star.samples)
    # Get columns in design matrix with random effects. 
    x.re.names <- colnames(object$X.re)
    indx <- which(colnames(X.0) %in% x.re.names)
    X.re <- as.matrix(X.0[, indx, drop = FALSE])
    X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
    n.occ.re <- ncol(object$beta.star.samples)
    # Form design matrix for random effects
    lambda.psi <- matrix(0, nrow(X.0), n.occ.re)
    for (i in 1:n.occ.re) {
      lambda.psi[which(X.re == (i - 1), arr.ind = TRUE)[, 1], i] <- 1
    }
    # Now can get the output
    out$psi.0.samples <- mcmc(logit.inv(t(X.fix %*% t(beta.samples) + 
					  lambda.psi %*% t(beta.star.samples))))
  } else {
    out$psi.0.samples <- mcmc(logit.inv(t(X.0 %*% t(beta.samples))))
  }
  out$z.0.samples <- mcmc(matrix(rbinom(length(out$psi.0.samples), 1, c(out$psi.0.samples)), 
  		                 nrow = n.post, ncol = nrow(X.0)))
  out$call <- cl

  class(out) <- "predict.PGOcc"
  out
}

print.PGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

fitted.PGOcc <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (class(object) != "PGOcc") {
    stop("error: object must be one of class PGOcc\n")
  }
  n.post <- object$n.post
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  J <- nrow(y)
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y))]
  if (nrow(X.p) == nrow(y)) {
    X.p <- do.call(rbind, replicate(ncol(y), X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y)), , drop = FALSE]
    if (object$pRE) {
      lambda.p <- do.call(rbind, replicate(ncol(y), object$lambda.p, simplify = FALSE))
      lambda.p <- lambda.p[!is.na(c(y)), , drop = FALSE]
    }
  } else {
    if (object$pRE) {
      lambda.p <- object$lambda.p
    }
  }
  y <- c(y)
  y <- y[!is.na(y)]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  if (object$pRE) {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) +
      			      lambda.p %*% t(object$alpha.star.samples)))
  } else {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
  }
  det.prob.samples <- det.prob.samples * z.samples[, z.long.indx]
  y.rep.samples <- t(apply(det.prob.samples, 2, function(a) rbinom(n.post, 1, a)))
  tmp <- array(NA, dim = c(J * K.max, n.post))
  names.long <- which(!is.na(c(object$y)))
  tmp[names.long, ] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(J, K.max, n.post))
  y.rep.samples <- aperm(y.rep.samples, c(3, 1, 2))
  return(y.rep.samples)
}

summary.PGOcc <- function(object,
			  quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
			  digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin

  cat("Chain Information:\n")
  cat(paste("Total samples: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thin: ",n.thin,"\n", sep=""))
  cat(paste("Total Posterior Samples: ",n.post,"\n\n", sep=""))
  
  # Occurrence
  cat("Occurrence: \n")
  print(noquote(round(t(apply(object$beta.samples, 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))
  if (object$psiRE) {
    cat("\n")
    cat("Occurrence Random Effect Variances: \n")
    print(noquote(round(t(apply(object$sigma.sq.psi.samples, 2, 
			        function(x) quantile(x, prob=quantiles))), digits)))
  }
  cat("\n")
  # Detection
  cat("Detection: \n")
  print(noquote(round(t(apply(object$alpha.samples, 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))
  if (object$pRE) {
    cat("\n")
    cat("Detection Random Effect Variances: \n")
    print(noquote(round(t(apply(object$sigma.sq.p.samples, 2, 
			        function(x) quantile(x, prob=quantiles))), digits)))
  }
}

# ppcOcc ------------------------------------------------------------------ 
summary.ppcOcc <- function(object, level, 
			   digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:", deparse(object$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin

  cat("Chain Information:\n")
  cat(paste("Total samples: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thin: ",n.thin,"\n", sep=""))
  cat(paste("Total Posterior Samples: ",n.post,"\n\n", sep=""))

  if (object$class %in% c('PGOcc', 'spPGOcc')) {
    cat("Bayesian p-value: ", mean(object$fit.y.rep > object$fit.y), "\n")
    cat("Fit statistic: ", object$fit.stat, "\n")
  }

  if (object$class %in% c('msPGOcc', 'spMsPGOcc')) {

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

  if (object$class %in% c('intPGOcc', 'spIntPGOcc')) {
    n.data <- length(object$fit.y.rep)
    for (q in 1:n.data) {
      cat("Data Source", q, "\n\n")	    
      cat("Bayesian p-value:", mean(object$fit.y.rep[[q]] > object$fit.y[[q]]), "\n")
      cat("Fit statistic:", object$fit.stat, "\n\n")
    }
  }

}

# spPGOcc -----------------------------------------------------------------

predict.spPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1, 
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

  # Data prep -------------------------------------------------------------
  n.samples <- object$n.post
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
  theta.samples <- t(theta.samples)
  beta.samples <- t(beta.samples)
  w.samples <- t(w.samples)

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
		 w.samples, n.samples, cov.model.indx, n.omp.threads, 
		 verbose, n.report)

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

summary.spPGOcc <- function(object,
			    quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
			    digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin

  cat("Chain Information:\n")
  cat(paste("Total samples: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thin: ",n.thin,"\n", sep=""))
  cat(paste("Total Posterior Samples: ",n.post,"\n\n", sep=""))
  
  # Occurrence
  cat("Occurrence: \n")
  print(noquote(round(t(apply(object$beta.samples, 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))
  cat("\n")
  # Detection
  cat("Detection: \n")
  print(noquote(round(t(apply(object$alpha.samples, 2, 
			      function(x) quantile(x, prob=quantiles))), digits)))

  if (object$pRE) {
    cat("\n")
    cat("Detection Random Effect Variances: \n")
    print(noquote(round(t(apply(object$sigma.sq.p.samples, 2, 
			        function(x) quantile(x, prob=quantiles))), digits)))
  }

  cat("\n")
  # Covariance
  cat("Covariance: \n")
  print(noquote(round(t(apply(object$theta.samples, 2, 
			      function(x) quantile(x, prob = quantiles))), digits)))
}


fitted.spPGOcc <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (class(object) != "spPGOcc") {
    stop("error: object must be one of class spPGOcc\n")
  }
  n.post <- object$n.post
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  J <- nrow(y)
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y))]
  if (nrow(X.p) == nrow(y)) {
    X.p <- do.call(rbind, replicate(ncol(y), X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y)), , drop = FALSE]
    if (object$pRE) {
      lambda.p <- do.call(rbind, replicate(ncol(y), object$lambda.p, simplify = FALSE))
      lambda.p <- lambda.p[!is.na(c(y)), , drop = FALSE]
    }
  } else {
    if (object$pRE) {
      lambda.p <- object$lambda.p
    }
  }
  y <- c(y)
  y <- y[!is.na(y)]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  if (object$pRE) {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) +
      			      lambda.p %*% t(object$alpha.star.samples)))
  } else {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
  }
  det.prob.samples <- det.prob.samples * z.samples[, z.long.indx]
  y.rep.samples <- t(apply(det.prob.samples, 2, function(a) rbinom(n.post, 1, a)))
  tmp <- array(NA, dim = c(J * K.max, n.post))
  names.long <- which(!is.na(c(object$y)))
  tmp[names.long, ] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(J, K.max, n.post))
  y.rep.samples <- aperm(y.rep.samples, c(3, 1, 2))
  return(y.rep.samples)
}

# msPGOcc -----------------------------------------------------------------

predict.msPGOcc <- function(object, X.0, ...) {
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

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  p.occ <- ncol(object$X)
  p.design <- p.occ
  if (object$psiRE) {
    p.design <- p.occ + ncol(object$sigma.sq.psi.samples)
  }
  if (ncol(X.0) != p.design) {
    stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
  }
  # Composition sampling --------------------------------------------------
  N <- dim(object$y)[1]
  sp.indx <- rep(1:N, p.occ)
  n.post <- object$n.post
  beta.samples <- as.matrix(object$beta.samples)
  out <- list()
  out$psi.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
  out$z.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
  if (object$psiRE) {
    beta.star.samples <- as.matrix(object$beta.star.samples)
    # Get columns in design matrix with random effects. 
    x.re.names <- colnames(object$X.re)
    indx <- which(colnames(X.0) %in% x.re.names)
    X.re <- as.matrix(X.0[, indx, drop = FALSE])
    X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
    n.occ.re <- ncol(object$beta.star.samples) / N
    # Form design matrix for random effects
    lambda.psi <- matrix(0, nrow(X.0), n.occ.re)
    for (i in 1:n.occ.re) {
      lambda.psi[which(X.re == (i - 1), arr.ind = TRUE)[, 1], i] <- 1
    }
    sp.re.indx <- rep(1:N, each = n.occ.re)
    for (i in 1:N) {
      out$psi.0.samples[, i, ] <- logit.inv(t(X.fix %*% t(beta.samples[, sp.indx == i]) + 
          				  lambda.psi %*% t(beta.star.samples[, sp.re.indx == i])))
      out$z.0.samples[, i, ] <- matrix(rbinom(n.post * nrow(X.0), 1, 
          				    c(out$psi.0.samples[, i, ])), 
          			     nrow = n.post, ncol = nrow(X.0))
    }
  } else {
    for (i in 1:N) {
      out$psi.0.samples[, i, ] <- logit.inv(t(as.matrix(X.0) %*% t(beta.samples[, sp.indx == i])))
      out$z.0.samples[, i, ] <- matrix(rbinom(n.post * nrow(X.0), 1, 
  					    c(out$psi.0.samples[, i, ])), 
  				     nrow = n.post, ncol = nrow(X.0))
    }
  }
  out$call <- cl

  class(out) <- "predict.msPGOcc"
  out
}

print.msPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.msPGOcc <- function(object,
			    level,
			    quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {

  if (missing(level)) {
    stop("error: must specify level of parameters to display. Valid values are 'species', 'community', or 'both'")
  }

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin

  cat("Chain Information:\n")
  cat(paste("Total samples: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thin: ",n.thin,"\n", sep=""))
  cat(paste("Total Posterior Samples: ",n.post,"\n\n", sep=""))

  if (tolower(level) == 'community') {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccurrence Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.beta.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    if (object$psiRE) {
      cat("\n")
      cat("Occurrence Random Effect Variances: \n")
      print(noquote(round(t(apply(object$sigma.sq.psi.samples, 2, 
          		        function(x) quantile(x, prob=quantiles))), digits)))
    }
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.alpha.samples, 2,
			      function(x) quantile(x, prob=quantiles))), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances: \n")
      print(noquote(round(t(apply(object$sigma.sq.p.samples, 2, 
          		        function(x) quantile(x, prob=quantiles))), digits)))
    }

  }

  if (tolower(level) == 'species') {
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence: \n")
    print(noquote(round(t(apply(object$beta.samples, 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples, 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))

  }

  if (tolower(level) == 'both') {
    
    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccurrence Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.beta.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    if (object$psiRE) {
      cat("\n")
      cat("Occurrence Random Effect Variances: \n")
      print(noquote(round(t(apply(object$sigma.sq.psi.samples, 2, 
          		        function(x) quantile(x, prob=quantiles))), digits)))
    }
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.alpha.samples, 2,
			      function(x) quantile(x, prob=quantiles))), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances: \n")
      print(noquote(round(t(apply(object$sigma.sq.p.samples, 2, 
          		        function(x) quantile(x, prob=quantiles))), digits)))
    }

    cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence: \n")
    print(noquote(round(t(apply(object$beta.samples, 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples, 2,
  			        function(x) quantile(x, prob=quantiles))), digits)))
  }

}

fitted.msPGOcc <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (class(object) != "msPGOcc") {
    stop("error: object must be of class msPGOcc\n")
  }
  n.post <- object$n.post
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  J <- dim(y)[2]
  N <- dim(y)[1]
  if (nrow(X.p) == dim(y)[2]) {
    X.p <- do.call(rbind, replicate(dim(y)[3], X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y[1, , ])), , drop = FALSE]
    if (object$pRE) {
      lambda.p <- do.call(rbind, replicate(dim(y)[3], object$lambda.p, simplify = FALSE))
      lambda.p <- lambda.p[!is.na(c(y[1, , ])), , drop = FALSE]
    }
  } else {
    if (object$pRE) {
      lambda.p <- object$lambda.p
    }
  }
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  n.obs <- nrow(X.p)
  det.prob.samples <- array(NA, dim = c(n.obs, N, n.post))
  sp.indx <- rep(1:N, ncol(X.p))
  y <- matrix(y, N, J * K.max)
  y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
  if (object$pRE) {
    sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
    for (i in 1:N) {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + 
      				   lambda.p %*% t(object$alpha.star.samples[, sp.re.indx == i]))
    }
  } else {
    for (i in 1:N) {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }
  }
  # Need to be careful here that all arrays line up. 
  det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
  det.prob.samples <- det.prob.samples * z.samples[, , z.long.indx]
  y.rep.samples <- array(NA, dim = dim(det.prob.samples))
  for (i in 1:N) {
    y.rep.samples[, i, ] <- apply(det.prob.samples[, i, ], 2, function(a) rbinom(n.post, 1, a))
  }
  tmp <- array(NA, dim = c(n.post, N, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(n.post, N, J, K.max))
  return(y.rep.samples)
}

# spMsPGOcc ---------------------------------------------------------------
summary.spMsPGOcc <- function(object, 
			      level,
			      quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {

  if (missing(level)) {
    stop("error: must specify level of parameters to display. Valid values are 'species', 'community', or 'both'")
  }

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin

  cat("Chain Information:\n")
  cat(paste("Total samples: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thin: ",n.thin,"\n", sep=""))
  cat(paste("Total Posterior Samples: ",n.post,"\n\n", sep=""))

  if (tolower(level) == 'community') {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccurrence Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.beta.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.alpha.samples, 2,
			      function(x) quantile(x, prob=quantiles))), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances: \n")
      print(noquote(round(t(apply(object$sigma.sq.p.samples, 2, 
          		        function(x) quantile(x, prob=quantiles))), digits)))
    }

  }

  if (tolower(level) == 'species') {
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence: \n")
    print(noquote(round(t(apply(object$beta.samples, 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples, 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))

    cat("\n")
    # Covariance
    cat("Covariance: \n")
    print(noquote(round(t(apply(object$theta.samples, 2, 
  			        function(x) quantile(x, prob = quantiles))), digits)))

  }

  if (tolower(level) == 'both') {
    
    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means: \n")
    print(noquote(round(t(apply(object$beta.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nOccurrence Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.beta.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection Means: \n")
    print(noquote(round(t(apply(object$alpha.comm.samples, 2,
          		      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\nDetection Variances: \n")
    print(noquote(round(t(apply(object$tau.sq.alpha.samples, 2,
			      function(x) quantile(x, prob=quantiles))), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances: \n")
      print(noquote(round(t(apply(object$sigma.sq.p.samples, 2, 
          		        function(x) quantile(x, prob=quantiles))), digits)))
    }
    cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence: \n")
    print(noquote(round(t(apply(object$beta.samples, 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    cat("\n")
    # Detection
    cat("Detection: \n")
    print(noquote(round(t(apply(object$alpha.samples, 2,
  			        function(x) quantile(x, prob=quantiles))), digits)))

    cat("\n")
    # Covariance
    cat("Covariance: \n")
    print(noquote(round(t(apply(object$theta.samples, 2, 
  			        function(x) quantile(x, prob = quantiles))), digits)))
  }

}

fitted.spMsPGOcc <- function(object, ...) {
  # NOTE: this is identical to fitted.msPGOcc
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (class(object) != "spMsPGOcc") {
    stop("error: object must be of class spMsPGOcc\n")
  }
  n.post <- object$n.post
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  J <- dim(y)[2]
  N <- dim(y)[1]
  if (nrow(X.p) == dim(y)[2]) {
    X.p <- do.call(rbind, replicate(dim(y)[3], X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y[1, , ])), , drop = FALSE]
    if (object$pRE) {
      lambda.p <- do.call(rbind, replicate(dim(y)[3], object$lambda.p, simplify = FALSE))
      lambda.p <- lambda.p[!is.na(c(y[1, , ])), , drop = FALSE]
    }
  } else {
    if (object$pRE) {
      lambda.p <- object$lambda.p
    }
  }
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  n.obs <- nrow(X.p)
  det.prob.samples <- array(NA, dim = c(n.obs, N, n.post))
  sp.indx <- rep(1:N, ncol(X.p))
  y <- matrix(y, N, J * K.max)
  y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
  if (object$pRE) {
    sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
    for (i in 1:N) {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + 
      				   lambda.p %*% t(object$alpha.star.samples[, sp.re.indx == i]))
    }
  } else {
    for (i in 1:N) {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }
  }
  # Need to be careful here that all arrays line up. 
  det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
  det.prob.samples <- det.prob.samples * z.samples[, , z.long.indx]
  y.rep.samples <- array(NA, dim = dim(det.prob.samples))
  for (i in 1:N) {
    y.rep.samples[, i, ] <- apply(det.prob.samples[, i, ], 2, function(a) rbinom(n.post, 1, a))
  }
  tmp <- array(NA, dim = c(n.post, N, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(n.post, N, J, K.max))
  return(y.rep.samples)
}

print.spMsPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}


predict.spMsPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1, 
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

  # Data prep -------------------------------------------------------------
  n.samples <- object$n.post
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
  theta.samples <- t(theta.samples)
  beta.samples <- t(beta.samples)
  w.samples <- aperm(w.samples, c(2, 3, 1))

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
		 w.samples, n.samples, cov.model.indx, n.omp.threads, 
		 verbose, n.report)

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

# intPGOcc ----------------------------------------------------------------
predict.intPGOcc <- function(object, X.0, ...) {
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
  if (class(object) != "intPGOcc") {
    stop("error: requires an output object of class intPGOcc\n")
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
  n.post <- object$n.post
  beta.samples <- as.matrix(object$beta.samples)
  out <- list()
  out$psi.0.samples <- mcmc(logit.inv(t(as.matrix(X.0) %*% t(beta.samples))))
  out$z.0.samples <- mcmc(matrix(rbinom(length(out$psi.0.samples), 1, c(out$psi.0.samples)),
		      nrow = n.post, ncol = nrow(X.0)))
  out$call <- cl

  class(out) <- "predict.intPGOcc"
  out
}

print.intPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.intPGOcc <- function(object, ...) {
  return(object$y.rep.samples)
}

summary.intPGOcc <- function(object,
			     quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
			     digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin

  cat("Chain Information:\n")
  cat(paste("Total samples: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thin: ",n.thin,"\n", sep=""))
  cat(paste("Total Posterior Samples: ",n.post,"\n\n", sep=""))

  n.data <- length(object$y)
  p.det.long <- sapply(object$X.p, function(a) dim(a)[[2]])

  # Occurrence
  cat("Occurrence: \n")
  print(noquote(round(t(apply(object$beta.samples, 2,
			      function(x) quantile(x, prob=quantiles))), digits)))
  cat("\n")
  # Detection
  indx <- 1
  for (i in 1:n.data) {
    cat(paste("Data source ", i, " Detection: \n", sep = ""))
    print(noquote(round(t(apply(object$alpha.samples[,indx:(indx+p.det.long[i] - 1), drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    indx <- indx + p.det.long[i]
    cat("\n")
  }
}

# spIntPGOcc --------------------------------------------------------------
print.spIntPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.spIntPGOcc <- function(object, ...) {
  return(object$y.rep.samples)
}

summary.spIntPGOcc <- function(object,
			       quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
			       digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin

  cat("Chain Information:\n")
  cat(paste("Total samples: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thin: ",n.thin,"\n", sep=""))
  cat(paste("Total Posterior Samples: ",n.post,"\n\n", sep=""))

  n.data <- length(object$y)
  p.det.long <- sapply(object$X.p, function(a) dim(a)[[2]])

  # Occurrence
  cat("Occurrence: \n")
  print(noquote(round(t(apply(object$beta.samples, 2,
			      function(x) quantile(x, prob=quantiles))), digits)))
  cat("\n")
  # Detection
  indx <- 1
  for (i in 1:n.data) {
    cat(paste("Data source ", i, " Detection: \n", sep = ""))
    print(noquote(round(t(apply(object$alpha.samples[,indx:(indx+p.det.long[i] - 1), drop = FALSE], 2,
  			      function(x) quantile(x, prob=quantiles))), digits)))
    indx <- indx + p.det.long[i]
    cat("\n")
  }
  # Covariance
  cat("Covariance: \n")
  print(noquote(round(t(apply(object$theta.samples, 2, 
			      function(x) quantile(x, prob = quantiles))), digits)))

}

predict.spIntPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1,
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
  if (class(object) != "spIntPGOcc") {
    stop("error: requires an output object of class spIntPGOcc\n")
  }

  # Data prep -------------------------------------------------------------
  n.samples <- object$n.post
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
  theta.samples <- t(theta.samples)
  beta.samples <- t(beta.samples)
  w.samples <- t(w.samples)

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

