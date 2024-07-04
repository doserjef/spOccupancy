waicOcc <- function(object, by.sp = FALSE, ...) {

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
  if (!(class(object) %in% c('PGOcc', 'spPGOcc', 'msPGOcc', 
                             'spMsPGOcc', 'intPGOcc', 'spIntPGOcc', 
                             'lfMsPGOcc', 'sfMsPGOcc', 'lfJSDM', 'sfJSDM', 
                             'tPGOcc', 'stPGOcc', 'svcPGBinom', 'svcPGOcc', 
                             'svcTPGBinom', 'svcTPGOcc', 'tMsPGOcc', 'intMsPGOcc', 
                             'svcMsPGOcc', 'stMsPGOcc', 'svcTMsPGOcc', 
                             'tIntPGOcc', 'stIntPGOcc', 'svcTIntPGOcc'))) {
    stop("error: object must be one of the following classes: PGOcc, spPGOcc, msPGOcc, spMsPGOcc, intPGOcc, spIntPGOcc, lfMsPGOcc, sfMsPGOcc, lfJSDM, sfJSDM, svcPGOcc, tPGOcc, stPGOcc, svcPGBinom, svcPGOcc, svcTPGBinom, svcTPGOcc, tMsPGOcc, intMsPGOcc, svcMsPGOcc, stMsPGOcc, svcTMsPGOcc, tIntPGOcc, stIntPGOcc, svcTIntPGOcc\n")
  }

  n.post <- object$n.post * object$n.chains
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  if (class(object) %in% c('PGOcc', 'spPGOcc', 'svcPGBinom', 'svcPGOcc')) {
    elpd <- sum(apply(object$like.samples, 2, function(a) log(mean(a))), na.rm = TRUE)
    pD <- sum(apply(object$like.samples, 2, function(a) var(log(a))), na.rm = TRUE)
    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }
  if (class(object) %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 
			   'sfMsPGOcc', 'lfJSDM', 'sfJSDM', 
			   'tPGOcc', 'stPGOcc', 'svcTPGBinom', 'svcTPGOcc', 'intMsPGOcc', 
			   'svcMsPGOcc')) {
    if (!by.sp) {
      elpd <- sum(apply(object$like.samples, c(2, 3), function(a) log(mean(a))), na.rm = TRUE)
      pD <- sum(apply(object$like.samples, c(2, 3), function(a) var(log(a))), na.rm = TRUE)
      out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
    } else {
      elpd <- apply(apply(object$like.samples, c(2, 3), function(a) log(mean(a))), 
                    1, sum, na.rm = TRUE) 
      pD <- apply(apply(object$like.samples, c(2, 3), function(a) var(log(a))), 
                  1, sum, na.rm = TRUE)
      out <- data.frame(elpd = elpd, 
			pD = pD, 
			WAIC = -2 * (elpd - pD))
    }
  }

  if (class(object) %in% c('tMsPGOcc', 'svcTMsPGOcc', 'stMsPGOcc')) {
    if (!by.sp) {
      elpd <- sum(apply(object$like.samples, c(2, 3, 4), function(a) log(mean(a))), na.rm = TRUE)
      pD <- sum(apply(object$like.samples, c(2, 3, 4), function(a) var(log(a))), na.rm = TRUE)
      out <- c(elpd, pD, -2 * (elpd - pD))
      names(out) <- c('elpd', 'pD', 'WAIC')
    } else {
      elpd <- apply(apply(object$like.samples, c(2, 3, 4), function(a) log(mean(a))), 
                    1, sum, na.rm = TRUE) 
      pD <- apply(apply(object$like.samples, c(2, 3, 4), function(a) var(log(a))), 
                  1, sum, na.rm = TRUE)
      out <- data.frame(elpd = elpd, 
			pD = pD, 
			WAIC = -2 * (elpd - pD))
    }
  }

  if (class(object) %in% c('intPGOcc', 'spIntPGOcc')) {
    n.data <- length(object$sites)
    J.long <- sapply(object$sites, length)
    indx <- 1
    elpd <- rep(0, n.data)
    pD <- rep(0, n.data)
    for (i in 1:n.data) {
      elpd[i] <- sum(apply(object$like.samples[, indx:(indx + J.long[i] - 1), drop = FALSE],
    		       2, function(a) log(mean(a))), na.rm = TRUE)
      pD[i] <- sum(apply(object$like.samples[, indx:(indx + J.long[i] - 1), drop = FALSE], 
    		     2, function(a) var(log(a))), na.rm = TRUE)
      indx <- indx + J.long[i]
    }
    out <- data.frame(elpd = elpd, 
		      pD = pD, 
		      WAIC = -2 * (elpd - pD))
  }
  if (class(object) %in% c('tIntPGOcc', 'stIntPGOcc', 'svcTIntPGOcc')) {
    n.data <- length(object$sites)
    n.waic.long <- rep(0, n.data)
    for (q in 1:n.data) {
      n.waic.long[q] <- length(object$sites[[q]]) * length(object$seasons[[q]])
    }
    indx <- 1
    elpd <- rep(0, n.data)
    pD <- rep(0, n.data)
    for (i in 1:n.data) {
      elpd[i] <- sum(apply(object$like.samples[, indx:(indx + n.waic.long[i] - 1), drop = FALSE],
                           2, function(a) log(mean(a))), na.rm = TRUE)
      pD[i] <- sum(apply(object$like.samples[, indx:(indx + n.waic.long[i] - 1), drop = FALSE], 
                         2, function(a) var(log(a))), na.rm = TRUE)
      indx <- indx + n.waic.long[i]
    }
    out <- data.frame(elpd = elpd, 
		      pD = pD, 
		      WAIC = -2 * (elpd - pD))
  }

  return(out)

}
