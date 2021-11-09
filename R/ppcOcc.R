ppcOcc <- function(object, fit.stat, group, ...) {

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
  # Fit statistic ---------------------
  if (missing(fit.stat)) {
    stop("error: fit.stat must be specified")
  }
  if (!tolower(fit.stat) %in% c('chi-square', 'freeman-tukey')) {
    stop("error: fit.stat must be either 'chi-square' or 'freeman-tukey'")
  }
  fit.stat <- tolower(fit.stat)
  # Group -----------------------------
  if (missing(group)) {
    stop("error: group must be specified")
  }
  if (!(group %in% c(1, 2)) & class(object) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc')) {
    stop("error: group must be 1 (row) or 2 (columns) for objects of class PGOcc, spPGOcc, intPGOcc, spIntPGOcc")
  }
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  out <- list()
  # For single species models
  if (class(object) %in% c('PGOcc', 'spPGOcc')) {
    y <- object$y
    X.p <- object$X.p
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
    p.det <- dim(X.p)[2]
    n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    J <- nrow(y)
    if (class(object) == 'PGOcc') {
      y.rep.samples <- fitted.PGOcc(object)
    } else {
      y.rep.samples <- fitted.spPGOcc(object)
    }
    z.samples <- object$z.samples
    alpha.samples <- object$alpha.samples
    # Get detection probability
    if (object$pRE) {
      det.prob <- logit.inv(X.p %*% t(alpha.samples) + 
		            lambda.p %*% t(object$alpha.star.samples))
    } else {
      det.prob <- logit.inv(X.p %*% t(alpha.samples))
    }
    det.prob <- array(det.prob, dim(y.rep.samples))
    n.samples <- object$n.post
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    e <- 0.0001
    # Do the stuff 
    if (group == 1) {
      y.grouped <- apply(y, 1, sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2), sum, na.rm = TRUE)
      fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
      fit.big.y <- matrix(NA, length(y.grouped), n.samples)
      if (fit.stat == 'chi-square') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 2, sum, na.rm = TRUE)
          fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (y.rep.grouped[i,] - E.grouped)^2 / (E.grouped + e)
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      } else if (fit.stat == 'freeman-tukey') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 2, sum, na.rm = TRUE)
          fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[i,]) - sqrt(E.grouped))^2 
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      }
    } else if (group == 2) {
      y.grouped <- apply(y, 2, sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 3), sum, na.rm = TRUE)
      fit.big.y <- matrix(NA, length(y.grouped), n.samples)
      fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
      if (fit.stat == 'chi-square') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 3, sum, na.rm = TRUE)
          fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (y.rep.grouped[i,] - E.grouped)^2 / (E.grouped + e)
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      } else if (fit.stat == 'freeman-tukey') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 3, sum, na.rm = TRUE)
          fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[i,]) - sqrt(E.grouped))^2 
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      }
    }
    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    out$fit.y.group.quants <- apply(fit.big.y, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    out$fit.y.rep.group.quants <- apply(fit.big.y.rep, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
  } 
  # Multispecies models
  if (class(object) %in% c('msPGOcc', 'spMsPGOcc')) {
    y <- object$y
    X.p <- object$X.p
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
    p.det <- dim(X.p)[2]
    n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
    J <- dim(y)[2]
    N <- dim(y)[1]
    if (class(object) == 'msPGOcc') {
      y.rep.samples <- fitted.msPGOcc(object)
    } else {
      y.rep.samples <- fitted.spMsPGOcc(object)
    }
    z.samples <- object$z.samples
    alpha.samples <- object$alpha.samples
    # Get detection probability
    n.samples <- object$n.post
    det.prob <- array(NA, dim = c(nrow(X.p), N, n.samples))
    sp.indx <- rep(1:N, ncol(X.p))
    if (object$pRE) {
      sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
      for (i in 1:N) {
        det.prob[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + 
  					   lambda.p %*% t(object$alpha.star.samples[, sp.re.indx == i]))
      }
    } else {
      for (i in 1:N) {
        det.prob[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
      }
    }
    det.prob <- aperm(det.prob, c(3, 2, 1))
    det.prob.full <- array(NA, dim(y.rep.samples))
    # Null model support
    if (sum(dim(y.rep.samples)[1:3] == dim(det.prob)) == 3) {
      det.prob.full[, , , ] <- det.prob
    } else {
      det.prob.full[!is.na(y.rep.samples)] <- det.prob
    }
    fit.y <- matrix(NA, n.samples, N)
    fit.y.rep <- matrix(NA, n.samples, N)
    e <- 0.0001
    # Do the stuff 
    if (group == 1) {
      y.grouped <- apply(y, c(1, 2), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2, 3), sum, na.rm = TRUE)
      fit.big.y.rep <- array(NA, dim = c(n.samples, N, J))
      fit.big.y <- array(NA, dim = c(n.samples, N, J))
      for (i in 1:N) {
        print(paste("Currently on species ", i, " out of ", N, sep = ''))
        if (fit.stat == 'chi-square') {
            for (j in 1:n.samples) {
              E.grouped <- apply(det.prob.full[j, i, , , drop = FALSE] * z.samples[j, i, ], 3, sum, na.rm = TRUE)
              fit.big.y[j, i, ] <- (y.grouped[i, ] - E.grouped)^2 / (E.grouped + e)
              fit.y[j, i] <- sum(fit.big.y[j, i, ])
              fit.big.y.rep[j, i, ] <- (y.rep.grouped[j, i, ] - E.grouped)^2 / (E.grouped + e)
              fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
            }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob.full[j, i, , , drop = FALSE] * z.samples[j, i, ], 3, sum, na.rm = TRUE)
            fit.big.y[j, i, ] <- (sqrt(y.grouped[i, ]) - sqrt(E.grouped))^2 
            fit.y[j, i] <- sum(fit.big.y[j, i, ])
            fit.big.y.rep[j, i, ] <- (sqrt(y.rep.grouped[j, i, ]) - sqrt(E.grouped))^2 
            fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
          }
        }
      }
    } else if (group == 2) {
      y.grouped <- apply(y, c(1, 3), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2, 4), sum, na.rm = TRUE)
      fit.big.y <- array(NA, dim = c(n.samples, N, max(n.rep)))
      fit.big.y.rep <- array(NA, dim = c(n.samples, N, max(n.rep)))
      for (i in 1:N) {
        print(paste("Currently on species ", i, " out of ", N, sep = ''))
        if (fit.stat == 'chi-square') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob.full[j, i, , , drop = FALSE] * z.samples[j, i, ], 4, sum, na.rm = TRUE)
            fit.big.y[j, i, ] <- (y.grouped[i, ] - E.grouped)^2 / (E.grouped + e)
            fit.y[j, i] <- sum(fit.big.y[j, i, ])
            fit.big.y.rep[j, i, ] <- (y.rep.grouped[j, i, ] - E.grouped)^2 / (E.grouped + e)
            fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
          }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob.full[j, i, , , drop = FALSE] * z.samples[j, i, ], 4, sum, na.rm = TRUE)
            fit.big.y[j, i, ] <- (sqrt(y.grouped[i, ]) - sqrt(E.grouped))^2 
            fit.y[j, i] <- sum(fit.big.y[j, i, ])
            fit.big.y.rep[j, i, ] <- (sqrt(y.rep.grouped[j, i, ]) - sqrt(E.grouped))^2 
            fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
          }
        }
      }
    }
    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    out$fit.y.group.quants <- apply(fit.big.y, c(2, 3), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    out$fit.y.rep.group.quants <- apply(fit.big.y.rep, c(2, 3), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$sp.names <- object$sp.names
  }
  # For integrated models
  if (class(object) %in% c('intPGOcc', 'spIntPGOcc')) {
    y <- object$y
    n.data <- length(y)
    sites <- object$sites
    X.p <- object$X.p
    p.det.long <- sapply(X.p, function(a) dim(a)[2])
    n.rep <- sapply(y, function(a1) apply(a1, 1, function(a2) sum(!is.na(a2))))
    J.long <- sapply(y, nrow)
    fit.y.list <- list()
    fit.y.rep.list <- list()
    fit.y.group.quants.list <- list()
    fit.y.rep.group.quants.list <- list()

    for (q in 1:n.data) {
      y.rep.samples <- object$y.rep.samples[[q]]
      z.samples <- object$z.samples[, sites[[q]], drop = FALSE]
      alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
      alpha.samples <- object$alpha.samples[, alpha.indx.r == q, drop = FALSE]
      # Get detection probability
      det.prob <- logit.inv(X.p[[q]] %*% t(alpha.samples))
      det.prob <- array(det.prob, dim(y.rep.samples))
      n.samples <- object$n.post
      fit.y <- rep(NA, n.samples)
      fit.y.rep <- rep(NA, n.samples)
      e <- 0.0001
      # Do the stuff 
      if (group == 1) {
        y.grouped <- apply(y[[q]], 1, sum, na.rm = TRUE)
        y.rep.grouped <- apply(y.rep.samples, c(1, 2), sum, na.rm = TRUE)
        fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
        fit.big.y <- matrix(NA, length(y.grouped), n.samples)
        if (fit.stat == 'chi-square') {
          for (i in 1:n.samples) {
            E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 2, sum, na.rm = TRUE)
            fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
            fit.y[i] <- sum(fit.big.y[, i])
            fit.big.y.rep[, i] <- (y.rep.grouped[i,] - E.grouped)^2 / (E.grouped + e)
            fit.y.rep[i] <- sum(fit.big.y.rep[, i])
          }
        } else if (fit.stat == 'freeman-tukey') {
          for (i in 1:n.samples) {
            E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 2, sum, na.rm = TRUE)
            fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
            fit.y[i] <- sum(fit.big.y[, i])
            fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[i,]) - sqrt(E.grouped))^2 
            fit.y.rep[i] <- sum(fit.big.y.rep[, i])
          }
        }
      } else if (group == 2) {
        y.grouped <- apply(y[[q]], 2, sum, na.rm = TRUE)
        y.rep.grouped <- apply(y.rep.samples, c(1, 3), sum, na.rm = TRUE)
        fit.big.y <- matrix(NA, length(y.grouped), n.samples)
        fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
        if (fit.stat == 'chi-square') {
          for (i in 1:n.samples) {
            E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 3, sum, na.rm = TRUE)
            fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
            fit.y[i] <- sum(fit.big.y[, i])
            fit.big.y.rep[, i] <- (y.rep.grouped[i,] - E.grouped)^2 / (E.grouped + e)
            fit.y.rep[i] <- sum(fit.big.y.rep[, i])
          }
        } else if (fit.stat == 'freeman-tukey') {
          for (i in 1:n.samples) {
            E.grouped <- apply(det.prob[i, , , drop = FALSE] * z.samples[i, ], 3, sum, na.rm = TRUE)
            fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
            fit.y[i] <- sum(fit.big.y[, i])
            fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[i,]) - sqrt(E.grouped))^2 
            fit.y.rep[i] <- sum(fit.big.y.rep[, i])
          }
        }
      }
      fit.y.list[[q]] <- fit.y
      fit.y.rep.list[[q]] <- fit.y.rep
      fit.y.group.quants.list[[q]] <- apply(fit.big.y, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
      fit.y.rep.group.quants.list[[q]] <- apply(fit.big.y.rep, 1, 
					      quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    }
    out$fit.y <- fit.y.list
    out$fit.y.rep <- fit.y.rep.list
    out$fit.y.group.quants <- fit.y.group.quants.list
    out$fit.y.rep.group.quants <- fit.y.rep.group.quants.list
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
  }

  class(out) <- 'ppcOcc'

  return(out)

}
