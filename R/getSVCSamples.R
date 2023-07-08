getSVCSamples <- function(object, pred.object, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (!(class(object) %in% c('svcPGOcc', 'svcPGBinom', 
			     'svcTPGOcc', 'svcTPGBinom', 
			     'svcMsPGOcc', 'svcTMsPGOcc'))) {
    stop("error: object must be of class svcPGOcc, svcPGBinom, svcTPGOcc, svcTPGBinom, svcMsPGOcc, svcTMsPGOcc\n")
  }

  n.post <- object$n.post * object$n.chains
  if (!missing(pred.object)) {
    J <- dim(pred.object$w.0.samples)[3]
  } else {
    J <- dim(object$w.samples)[3]
  }
  if (length(dim(object$X)) == 3) {
    svc.names <- dimnames(object$X)[[3]][object$svc.cols]
  } else {
    svc.names <- colnames(object$X)[object$svc.cols]
  }
  svc.cols <- object$svc.cols
  p.svc <- length(svc.cols)
  # Single-species models -------------------------------------------------
  if (class(object) %in% c('svcPGOcc', 'svcPGBinom', 
			   'svcTPGOcc', 'svcTPGBinom')) {
    if (!missing(pred.object)) {
      svc.samples <- lapply(svc.cols, function(a) mcmc(object$beta.samples[, a] + pred.object$w.0.samples[, which(svc.cols == a), ]))
    } else {
      svc.samples <- lapply(svc.cols, function(a) mcmc(object$beta.samples[, a] + object$w.samples[, which(svc.cols == a), ]))
    }
  }
  # Multi-species models --------------------------------------------------
  if (class(object) %in% c('svcMsPGOcc', 'svcTMsPGOcc')) {
    N <- nrow(object$y)
    if (!missing(pred.object)) {
      J <- dim(pred.object$z.0.samples)[3]
    } else {
      J <- ncol(object$y)
    }
    q <- object$q
    svc.samples <- list()
    for (i in 1:p.svc) {
      svc.samples[[i]] <- array(NA, dim = c(N, J, n.post)) 
    }
    lambda.samples <- array(object$lambda.samples, dim = c(n.post, N, q, p.svc))
    beta.samples <- array(object$beta.samples, dim = c(n.post, N, ncol(object$X)))
    for (i in 1:n.post) {
        for (j in 1:p.svc) {
          tmp <- matrix(lambda.samples[i, , , j], N, q)
          if (!missing(pred.object)) {
            tmp.2 <- matrix(pred.object$w.0.samples[i, , , j], q, J)
	  } else {
            tmp.2 <- matrix(object$w.samples[i, , , j], q, J)
	  }
          svc.samples[[j]][, , i] <- tmp %*% tmp.2 + beta.samples[i, , svc.cols[j]]
        }
    }
    svc.samples <- lapply(svc.samples, aperm, c(3, 1, 2))
  }
  names(svc.samples) <- svc.names
  return(svc.samples)
}
