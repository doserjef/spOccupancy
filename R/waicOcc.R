waicOcc <- function(object, ...) {

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
  if (!class(object) %in% c('PGOcc', 'spPGOcc', 'msPGOcc', 
			    'spMsPGOcc', 'intPGOcc', 'spIntPGOcc')) {
    stop("error: object must be one of the following classes: PGOcc, spPGOcc, msPGOcc, spMsPGOcc, intPGOcc, spIntPGOcc\n")
  }

  n.post <- object$n.post * object$n.chains
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  if (class(object) %in% c('PGOcc', 'spPGOcc')) {
    X.p <- object$X.p
    y <- object$y
    # Necessary for when X.p is only at site level
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
    n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    y.sum <- apply(y, 1, sum, na.rm = TRUE)
    y.ind <- as.numeric(y.sum == 0)
    K.max <- max(n.rep)
    J <- nrow(y)
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y))]
    y <- c(y)
    y <- y[!is.na(y)]
    psi.samples <- object$psi.samples
    alpha.samples <- object$alpha.samples
    if (object$pRE) {
      det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) + 
				      lambda.p %*% t(object$alpha.star.samples)))
    } else {
      det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
    }
    # Convert det.prob.samples to likelihood form. 
    det.prob.samples[, y == 0] <- 1 - det.prob.samples[, y == 0]
    n.obs <- length(y)

    elpd <- 0
    pD <- 0

    for (i in 1:J) {
      long.indx <- which(z.long.indx == i)
      L <- apply(det.prob.samples[, long.indx, drop = FALSE], 1, prod)  * 
	      psi.samples[, i] + y.ind[i] * (1 - psi.samples[, i])
      elpd <- elpd + log(mean(L))
      pD <- pD + var(log(L))
    }

    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }
  if (class(object) %in% c('msPGOcc', 'spMsPGOcc')) {
    X.p <- object$X.p
    y <- object$y
    y.ind <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) == 0))
    n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
    # Necessary for when X.p is only at site level
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
    K.max <- max(n.rep)
    J <- dim(y)[2]
    N <- dim(y)[1]
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
    y <- matrix(y, N, J * K.max)
    y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
    psi.samples <- object$psi.samples
    alpha.samples <- object$alpha.samples
    # Get detection probability
    n.obs <- nrow(X.p)
    det.prob.samples <- array(NA, dim = c(n.obs, N, n.post))
    sp.indx <- rep(1:N, ncol(X.p))
    if (object$pRE) {
      sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
      for (i in 1:N) {
        det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + 
  					   lambda.p %*% t(object$alpha.star.samples[, sp.re.indx == i]))
        # Convert for likelihood
        det.prob.samples[y[i, ] == 0, i, ] <- 1 - det.prob.samples[y[i, ] == 0, i, ]
      }
    } else {
      for (i in 1:N) {
        det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
        # Convert for likelihood
        det.prob.samples[y[i, ] == 0, i, ] <- 1 - det.prob.samples[y[i, ] == 0, i, ]
      }
    }
    det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
    n.obs <- length(y)

    elpd <- 0
    pD <- 0

    for (j in 1:J) {
      long.indx <- which(z.long.indx == j)
      for (i in 1:N) {
        L <- apply(det.prob.samples[, i, long.indx, drop = FALSE], 1, prod)  * 
	      psi.samples[, i, j] + y.ind[i, j] * (1 - psi.samples[, i, j])
        elpd <- elpd + log(mean(L))
        pD <- pD + var(log(L))
      }
    }

    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
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
      K.max <- max(n.rep)
      J <- nrow(y[[q]])
      z.long.indx <- rep(sites[[q]], K.max)
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
        long.indx <- which(z.long.indx == j)
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
