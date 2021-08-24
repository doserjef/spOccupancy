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
  p.occ <- ncol(object$X.occ)
  if (ncol(X.0) != p.occ) {
    stop(paste("error: X.0 must have ", p.occ, " columns\n", sep = ''))
  }

  # Composition sampling --------------------------------------------------
  beta.samples <- as.matrix(object$beta.samples[s.indx, , drop = FALSE])
  out <- list()
  out$psi.hat <- logit.inv(t(as.matrix(X.0) %*% t(beta.samples)))
  out$z.hat <- matrix(rbinom(length(out$psi.hat), 1, c(out$psi.hat)), 
		      nrow = n.samples, ncol = ncol(X.0))
  out$sub.sample <- sub.sample
  out$s.indx <- s.indx

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
  print(noquote(apply(t(apply(object$beta.samples[s.indx,], 2, 
			      function(x) quantile(x, prob=quantiles))), 
		      2, function(x) formatC(x, format = "f", digits = digits))))
  cat("\n")
  # Detection
  cat("Detection: \n")
  print(noquote(apply(t(apply(object$alpha.samples[s.indx,], 2, 
			      function(x) quantile(x, prob=quantiles))), 
		      2, function(x) formatC(x, format = "f", digits = digits))))
  

}
