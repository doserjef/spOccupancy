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
			     'svcMsPGOcc', 'stMsPGOcc', 'svcTMsPGOcc'))) {
    stop("error: object must be one of the following classes: PGOcc, spPGOcc, msPGOcc, spMsPGOcc, intPGOcc, spIntPGOcc, lfMsPGOcc, sfMsPGOcc, lfJSDM, sfJSDM, svcPGOcc, tPGOcc, stPGOcc, svcPGBinom, svcPGOcc, svcTPGBinom, svcTPGOcc, tMsPGOcc, intMsPGOcc, svcMsPGOcc, stMsPGOcc, svcTMsPGOcc\n")
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
    X.p <- object$X.p
    y <- object$y
    y.sum <- lapply(y, function(q) apply(q, 1, sum, na.rm = TRUE))
    y.ind <- lapply(y.sum, function(q) as.numeric(q == 0))
    n.data <- length(y)
    sites <- object$sites
    p.det.long <- sapply(X.p, function(a) dim(a)[2])
    n.rep <- sapply(y, function(a1) apply(a1, 1, function(a2) sum(!is.na(a2))))
    J.long <- sapply(y, nrow)
    alpha.indx.r <- unlist(lapply(1:n.data, function(a) rep(a, p.det.long[a])))

    elpd <- rep(0, n.data)
    pD <- rep(0, n.data)

    for (q in 1:n.data) {

      n.rep <- apply(y[[q]], 1, function(a) sum(!is.na(a)))
      J <- nrow(y[[q]])
      z.long.indx <- rep(sites[[q]], dim(y[[q]])[2])
      z.long.indx <- z.long.indx[!is.na(c(y[[q]]))]
      y[[q]] <- c(y[[q]])
      y[[q]] <- y[[q]][!is.na(y[[q]])]
      psi.samples <- object$psi.samples
      alpha.samples <- object$alpha.samples[, alpha.indx.r == q, drop = FALSE]
      det.prob.samples <- t(logit.inv(X.p[[q]] %*% t(alpha.samples)))
      n.obs <- length(y[[q]])
      curr.sites <- object$sites[[q]]

      # Convert det.prob.samples to likelihood form. 
      det.prob.samples[, y[[q]] == 0] <- 1 - det.prob.samples[, y[[q]] == 0]


      for (j in 1:J) {
        long.indx <- which(z.long.indx == curr.sites[j])
        L <- apply(det.prob.samples[, long.indx, drop = FALSE], 1, prod)  * 
                psi.samples[, curr.sites[j]] + y.ind[[q]][j] * (1 - psi.samples[, curr.sites[j]])
        elpd[q] <- elpd[q] + log(mean(L))
        if (is.na(elpd[q])) print(j)
          pD[q] <- pD[q] + var(log(L))
      }
    }

    out <- data.frame(elpd = elpd, 
		      pD = pD, 
		      WAIC = -2 * (elpd - pD))
  }

  return(out)

}
