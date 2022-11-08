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
			     'svcTPGOcc', 'svcTPGBinom'))) {
    stop("error: object must be of class svcPGOcc, svcPGBinom, svcTPGOcc, svcTPGBinom\n")
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
  if (!missing(pred.object)) {
    svc.samples <- lapply(svc.cols, function(a) mcmc(object$beta.samples[, a] + pred.object$w.0.samples[, a, ]))
  } else {
    svc.samples <- lapply(svc.cols, function(a) mcmc(object$beta.samples[, a] + object$w.samples[, a, ]))
  }
  names(svc.samples) <- svc.names
  return(svc.samples)
}
