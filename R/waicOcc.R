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

  n.post <- object$n.post
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  if (class(object) %in% c('PGOcc', 'spPGOcc')) {
    X.p <- object$X.p
    y <- object$y
    n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    K.max <- max(n.rep)
    J <- nrow(y)
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y))]
    y <- c(y)
    y <- y[!is.na(y)]
    psi.samples <- object$psi.samples
    alpha.samples <- object$alpha.samples
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
    p.psi.samples <- det.prob.samples * psi.samples[, z.long.indx]
    n.obs <- length(y)

    elpd <- 0
    pD <- 0

    for (i in 1:n.obs) {
      prob <- p.psi.samples[, i] 
      L <- dbinom(y[i], 1, prob, log = FALSE)
      elpd <- elpd + log(mean(L))
      pD <- pD + var(log(L))
    }

    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }
  if (class(object) %in% c('msPGOcc', 'spMsPGOcc')) {
    X.p <- object$X.p
    y <- object$y
    n.rep <- apply(y[1, , ], 1, function(a) sum(!is.na(a)))
    K.max <- max(n.rep)
    J <- dim(y)[2]
    N <- dim(y)[1]
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y))]
    z.long.indx <- rep(1:J, K.max)
    z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
    y <- matrix(y, N, J * K.max)
    y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
    psi.samples <- object$psi.samples
    alpha.samples <- object$alpha.samples
    # Get detection probability
    n.obs <- nrow(X.p)
    det.prob.samples <- array(NA, dim = c(n.post, N, n.obs))
    sp.indx <- rep(1:N, ncol(X.p))
    for (i in 1:N) {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }

    p.psi.samples <- det.prob.samples * psi.samples[, , z.long.indx]

    elpd <- 0
    pD <- 0

    for (i in 1:n.obs) {
      prob <- p.psi.samples[, , i] 
      L <- c(dbinom(y[, i], 1, prob, log = FALSE))
      elpd <- elpd + log(mean(L))
      pD <- pD + var(log(L))
    }

    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }

  if (class(object) %in% c('intPGOcc', 'spIntPGOcc')) {
    X.p <- object$X.p
    y <- object$y
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
      p.psi.samples <- det.prob.samples * psi.samples[, z.long.indx]
      n.obs <- length(y[[q]])

      for (i in 1:n.obs) {
        prob <- p.psi.samples[, i] 
        L <- dbinom(y[[q]][i], 1, prob, log = FALSE)
        elpd[q] <- elpd[q] + log(mean(L))
        pD[q] <- pD[q] + var(log(L))
      }
    }

    out <- data.frame(elpd = elpd, 
		      pD = pD, 
		      WAIC = -2 * (elpd - pD))
  }

  return(out)

}
