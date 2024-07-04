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
  if (class(object) %in% c('lfJSDM', 'sfJSDM', 'svcPGBinom', 'svcTPGBinom')) {
    stop("error: ppcOcc is not implemented for lfJSDM, sfJSDM, svcPGBinom, svcTPGBinom")
  }
  if (!(class(object) %in% c('PGOcc', 'spPGOcc', 'msPGOcc', 
                             'spMsPGOcc', 'intPGOcc', 'spIntPGOcc', 
                             'lfMsPGOcc', 'sfMsPGOcc', 'tPGOcc', 'stPGOcc', 
                             'svcPGOcc', 'svcTPGOcc', 'tMsPGOcc', 'svcMsPGOcc', 
                             'svcTMsPGOcc', 'stMsPGOcc', 'tIntPGOcc', 
                             'stIntPGOcc', 'svcTIntPGOcc'))) {
    stop("error: object must be one of the following classes: PGOcc, spPGOcc, msPGOcc, spMsPGOcc, intPGOcc, spIntPGOcc, lfMsPGOcc, sfMsPGOcc, tPGOcc, stPGOcc, svcPGOcc, svcTPGOcc, tMsPGOcc, svcMsPGOcc, svcTMsPGOcc, stMsPGOcc, tIntPGOcc, stIntPGOcc, svcTIntPGOcc\n")
  }
  # Fit statistic ---------------------
  if (missing(fit.stat)) {
    stop("error: fit.stat must be specified")
  }
  if (!tolower(fit.stat) %in% c('chi-squared', 'freeman-tukey', 'chi-square')) {
    stop("error: fit.stat must be either 'chi-squared' or 'freeman-tukey'")
  }
  fit.stat <- tolower(fit.stat)
  # Group -----------------------------
  if (missing(group)) {
    stop("error: group must be specified")
  }
  if (!(group %in% c(1, 2))) {
    stop("error: group must be 1 (sites) or 2 (replicates)")
  }
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  out <- list()
  # Single-species models -------------------------------------------------
  #if (is(object, c('PGOcc', 'spPGOcc'))) {
  if (class(object) %in% c('PGOcc', 'spPGOcc', 'svcPGOcc')) {
    y <- object$y
    J <- nrow(y)
    if (is(object, 'PGOcc')) {
      fitted.out <- fitted.PGOcc(object)
    } else {
      fitted.out <- fitted.spPGOcc(object)
    }
    z.samples <- object$z.samples
    y.rep.samples <- fitted.out$y.rep.samples
    det.prob <- fitted.out$p.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    e <- 0.0001
    # Do the stuff 
    if (group == 1) {
      y.grouped <- apply(y, 1, sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2), sum, na.rm = TRUE)
      fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
      fit.big.y <- matrix(NA, length(y.grouped), n.samples)
      if (fit.stat %in% c('chi-squared', 'chi-square')) {
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
      if (fit.stat %in% c('chi-squared', 'chi-square')) {
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
    out$n.chains <- object$n.chains
  } 
  # Multispecies models ---------------------------------------------------
  # if (is(object, c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc'))) {
  if (class(object) %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc', 'svcMsPGOcc')) {
    y <- object$y
    J <- dim(y)[2]
    N <- dim(y)[1]
    # Fitted function is the same for all multispecies object types. 
    fitted.out <- fitted.msPGOcc(object)
    y.rep.samples <- fitted.out$y.rep.samples
    det.prob <- fitted.out$p.samples
    z.samples <- object$z.samples
    n.samples <- object$n.post * object$n.chains
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
        message(noquote(paste("Currently on species ", i, " out of ", N, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
            for (j in 1:n.samples) {
              E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * z.samples[j, i, ], 3, sum, na.rm = TRUE)
              fit.big.y[j, i, ] <- (y.grouped[i, ] - E.grouped)^2 / (E.grouped + e)
              fit.y[j, i] <- sum(fit.big.y[j, i, ])
              fit.big.y.rep[j, i, ] <- (y.rep.grouped[j, i, ] - E.grouped)^2 / (E.grouped + e)
              fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
            }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * z.samples[j, i, ], 3, sum, na.rm = TRUE)
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
      fit.big.y <- array(NA, dim = c(n.samples, N, dim(y)[3]))
      fit.big.y.rep <- array(NA, dim = c(n.samples, N, dim(y)[3]))
      for (i in 1:N) {
        message(noquote(paste("Currently on species ", i, " out of ", N, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * z.samples[j, i, ], 4, sum, na.rm = TRUE)
            fit.big.y[j, i, ] <- (y.grouped[i, ] - E.grouped)^2 / (E.grouped + e)
            fit.y[j, i] <- sum(fit.big.y[j, i, ])
            fit.big.y.rep[j, i, ] <- (y.rep.grouped[j, i, ] - E.grouped)^2 / (E.grouped + e)
            fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
          }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * z.samples[j, i, ], 4, sum, na.rm = TRUE)
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
    out$n.chains <- object$n.chains
    out$sp.names <- object$sp.names
  }
  # Integrated models -----------------------------------------------------
  if (class(object) %in% c('intPGOcc', 'spIntPGOcc')) {
    y <- object$y
    n.data <- length(y)
    sites <- object$sites
    X.p <- object$X.p
    p.det.long <- sapply(X.p, function(a) dim(a)[2])
    J.long <- sapply(y, nrow)
    fitted.out <- fitted.intPGOcc(object)
    y.rep.all <- fitted.out$y.rep.samples
    det.prob.all <- fitted.out$p.samples
    fit.y.list <- list()
    fit.y.rep.list <- list()
    fit.y.group.quants.list <- list()
    fit.y.rep.group.quants.list <- list()

    for (q in 1:n.data) {
      y.rep.samples <- y.rep.all[[q]] 
      z.samples <- object$z.samples[, sites[[q]], drop = FALSE]
      alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
      alpha.samples <- object$alpha.samples[, alpha.indx.r == q, drop = FALSE]
      # Get detection probability
      det.prob <- det.prob.all[[q]]
      n.samples <- object$n.post * object$n.chains
      fit.y <- rep(NA, n.samples)
      fit.y.rep <- rep(NA, n.samples)
      e <- 0.0001
      # Do the stuff 
      if (group == 1) {
        y.grouped <- apply(y[[q]], 1, sum, na.rm = TRUE)
        y.rep.grouped <- apply(y.rep.samples, c(1, 2), sum, na.rm = TRUE)
        fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
        fit.big.y <- matrix(NA, length(y.grouped), n.samples)
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
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
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
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
    out$n.chains <- object$n.chains
  }
  # Multi-season models ---------------------------------------------------
  # if (is(object, c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc'))) {
  if (class(object) %in% c('tPGOcc', 'stPGOcc', 'svcTPGOcc')) {
    y <- object$y
    J <- dim(y)[1]
    n.years.max <- dim(y)[2]
    fitted.out <- fitted.stPGOcc(object)
    y.rep.samples <- fitted.out$y.rep.samples
    det.prob <- fitted.out$p.samples
    z.samples <- object$z.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- matrix(NA, n.samples, n.years.max)
    fit.y.rep <- matrix(NA, n.samples, n.years.max)
    e <- 0.0001
    # Do the stuff 
    if (group == 1) { # Group by site
      y.grouped <- apply(y, c(1, 2), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2, 3), sum, na.rm = TRUE)
      fit.big.y.rep <- array(NA, dim = c(n.samples, J, n.years.max))
      fit.big.y <- array(NA, dim = c(n.samples, J, n.years.max))
      for (t in 1:n.years.max) {
        message(noquote(paste("Currently on time period ", t, " out of ", n.years.max, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
            for (j in 1:n.samples) {
              E.grouped <- apply(det.prob[j, , t, ] * z.samples[j, , t], 1, sum, na.rm = TRUE)
              fit.big.y[j, , t] <- (y.grouped[, t] - E.grouped)^2 / (E.grouped + e)
              fit.y[j, t] <- sum(fit.big.y[j, , t])
              fit.big.y.rep[j, , t] <- (y.rep.grouped[j, , t] - E.grouped)^2 / (E.grouped + e)
              fit.y.rep[j, t] <- sum(fit.big.y.rep[j, , t])
            }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, , t, ] * z.samples[j, , t], 1, sum, na.rm = TRUE)
            fit.big.y[j, , t] <- (sqrt(y.grouped[, t]) - sqrt(E.grouped))^2 
            fit.y[j, t] <- sum(fit.big.y[j, , t])
            fit.big.y.rep[j, , t] <- (sqrt(y.rep.grouped[j, , t]) - sqrt(E.grouped))^2 
            fit.y.rep[j, t] <- sum(fit.big.y.rep[j, , t])
          }
        }
      }
    } else if (group == 2) { # Group by visit
      y.grouped <- apply(y, c(2, 3), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 3, 4), sum, na.rm = TRUE)
      fit.big.y <- array(NA, dim = c(n.samples, n.years.max, dim(y)[3]))
      fit.big.y.rep <- array(NA, dim = c(n.samples, n.years.max, dim(y)[3]))
      for (t in 1:n.years.max) {
        message(noquote(paste("Currently on time period ", t, " out of ", n.years.max, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, , t, ] * z.samples[j, , t], 2, sum, na.rm = TRUE)
            fit.big.y[j, t, ] <- (y.grouped[t, ] - E.grouped)^2 / (E.grouped + e)
            fit.y[j, t] <- sum(fit.big.y[j, t, ])
            fit.big.y.rep[j, t, ] <- (y.rep.grouped[j, t, ] - E.grouped)^2 / (E.grouped + e)
            fit.y.rep[j, t] <- sum(fit.big.y.rep[j, t, ])
          }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, , t, ] * z.samples[j, , t], 2, sum, na.rm = TRUE)
            fit.big.y[j, t, ] <- (sqrt(y.grouped[t, ]) - sqrt(E.grouped))^2 
            fit.y[j, t] <- sum(fit.big.y[j, t, ])
            fit.big.y.rep[j, t, ] <- (sqrt(y.rep.grouped[j, t, ]) - sqrt(E.grouped))^2 
            fit.y.rep[j, t] <- sum(fit.big.y.rep[j, t, ])
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
    out$n.chains <- object$n.chains
  }
  # Multi-species, multi-season models ------------------------------------
  if (class(object) %in% c('tMsPGOcc', 'svcTMsPGOcc', 'stMsPGOcc')) {
    y <- object$y
    J <- dim(y)[2]
    N <- dim(y)[1]
    n.years.max <- dim(y)[3]
    # Fitted function is the same for all multi-species, multi-season object types. 
    fitted.out <- fitted.tMsPGOcc(object)
    y.rep.samples <- fitted.out$y.rep.samples
    det.prob <- fitted.out$p.samples
    z.samples <- object$z.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- array(NA, dim = c(n.samples, N, n.years.max))
    fit.y.rep <- array(NA, dim = c(n.samples, N, n.years.max))
    e <- 0.0001
    # Do the stuff 
    if (group == 1) {
      y.grouped <- apply(y, c(1, 2, 3), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2, 3, 4), sum, na.rm = TRUE)
      fit.big.y.rep <- array(NA, dim = c(n.samples, N, J, n.years.max))
      fit.big.y <- array(NA, dim = c(n.samples, N, J, n.years.max))
      for (t in 1:n.years.max) {
        for (i in 1:N) {
          message(noquote(paste("Currently on species ", i, " out of ", N, 
        			" and time period ", t, " out of ", n.years.max, sep = '')))
          if (fit.stat %in% c('chi-squared', 'chi-square')) {
              for (j in 1:n.samples) {
                E.grouped <- apply(det.prob[j, i, , t, , drop = FALSE] * z.samples[j, i, , t], 3, sum, na.rm = TRUE)
                fit.big.y[j, i, , t] <- (y.grouped[i, , t] - E.grouped)^2 / (E.grouped + e)
                fit.y[j, i, t] <- sum(fit.big.y[j, i, , t])
                fit.big.y.rep[j, i, , t] <- (y.rep.grouped[j, i, , t] - E.grouped)^2 / (E.grouped + e)
                fit.y.rep[j, i, t] <- sum(fit.big.y.rep[j, i, , t])
              }
          } else if (fit.stat == 'freeman-tukey') {
            for (j in 1:n.samples) {
              E.grouped <- apply(det.prob[j, i, , t, , drop = FALSE] * z.samples[j, i, , t], 3, sum, na.rm = TRUE)
              fit.big.y[j, i, , t] <- (sqrt(y.grouped[i, , t]) - sqrt(E.grouped))^2 
              fit.y[j, i, t] <- sum(fit.big.y[j, i, , t])
              fit.big.y.rep[j, i, , t] <- (sqrt(y.rep.grouped[j, i, , t]) - sqrt(E.grouped))^2 
              fit.y.rep[j, i, t] <- sum(fit.big.y.rep[j, i, , t])
            }
          }
        }
      }
    } else if (group == 2) {
      y.grouped <- apply(y, c(1, 3, 4), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2, 4, 5), sum, na.rm = TRUE)
      fit.big.y <- array(NA, dim = c(n.samples, N, n.years.max, dim(y)[4]))
      fit.big.y.rep <- array(NA, dim = c(n.samples, N, n.years.max, dim(y)[4]))
      for (t in 1:n.years.max) {
        for (i in 1:N) {
          message(noquote(paste("Currently on species ", i, " out of ", N, 
        			" and time period ", t, " out of ", n.years.max, sep = '')))
          if (fit.stat %in% c('chi-squared', 'chi-square')) {
            for (j in 1:n.samples) {
              E.grouped <- apply(det.prob[j, i, , t, , drop = FALSE] * z.samples[j, i, , t], 5, sum, na.rm = TRUE)
              fit.big.y[j, i, t, ] <- (y.grouped[i, t, ] - E.grouped)^2 / (E.grouped + e)
              fit.y[j, i, t] <- sum(fit.big.y[j, i, t, ])
              fit.big.y.rep[j, i, t, ] <- (y.rep.grouped[j, i, t, ] - E.grouped)^2 / (E.grouped + e)
              fit.y.rep[j, i, t] <- sum(fit.big.y.rep[j, i, t, ])
            }
          } else if (fit.stat == 'freeman-tukey') {
            for (j in 1:n.samples) {
              E.grouped <- apply(det.prob[j, i, , t, , drop = FALSE] * z.samples[j, i, , t], 5, sum, na.rm = TRUE)
              fit.big.y[j, i, t, ] <- (sqrt(y.grouped[i, t, ]) - sqrt(E.grouped))^2 
              fit.y[j, i, t] <- sum(fit.big.y[j, i, t, ])
              fit.big.y.rep[j, i, t, ] <- (sqrt(y.rep.grouped[j, i, t, ]) - sqrt(E.grouped))^2 
              fit.y.rep[j, i, t] <- sum(fit.big.y.rep[j, i, t, ])
            }
          }
        }
      }
    }
    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    out$fit.y.group.quants <- apply(fit.big.y, c(2, 3, 4), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    out$fit.y.rep.group.quants <- apply(fit.big.y.rep, c(2, 3, 4), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$n.chains <- object$n.chains
    out$sp.names <- object$sp.names
  }
  # Integrated multi-season models ----------------------------------------
  if (class(object) %in% c('tIntPGOcc', 'stIntPGOcc', 'svcTIntPGOcc')) {
    y <- object$y
    n.data <- length(y)
    sites <- object$sites
    seasons <- object$seasons
    X.p <- object$X.p
    p.det.long <- sapply(X.p, function(a) dim(a)[2])
    J.long <- sapply(y, nrow)
    n.years.long <- sapply(y, ncol)
    fitted.out <- fitted.tIntPGOcc(object)
    y.rep.all <- fitted.out$y.rep.samples
    det.prob.all <- fitted.out$p.samples
    fit.y.list <- list()
    fit.y.rep.list <- list()
    fit.y.group.quants.list <- list()
    fit.y.rep.group.quants.list <- list()

    for (q in 1:n.data) {
      y.rep.samples <- y.rep.all[[q]] 
      z.samples <- object$z.samples[, sites[[q]], seasons[[q]], drop = FALSE]
      alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
      alpha.samples <- object$alpha.samples[, alpha.indx.r == q, drop = FALSE]
      # Get detection probability
      det.prob <- det.prob.all[[q]]
      n.samples <- object$n.post * object$n.chains
      fit.y <- matrix(NA, n.samples, n.years.long[q])
      fit.y.rep <- matrix(NA, n.samples, n.years.long[q])
      e <- 0.0001
      # Do the stuff 
      if (group == 1) {
        y.grouped <- apply(y[[q]], c(1, 2), sum, na.rm = TRUE)
        y.rep.grouped <- apply(y.rep.samples, c(1, 2, 3), sum, na.rm = TRUE)
        fit.big.y.rep <- array(NA, dim = c(n.samples, J.long[q], n.years.long[q]))
        fit.big.y <- array(NA, dim = c(n.samples, J.long[q], n.years.long[q]))
        for (t in 1:n.years.long[q]) {
          message(noquote(paste("Currently on time period ", t, " out of ", n.years.long[q], 
                                ' for data set ', q, sep = '')))
          if (fit.stat %in% c('chi-squared', 'chi-square')) {
            for (j in 1:n.samples) {
              p.curr <- matrix(det.prob[j, , t, ], nrow = J.long[q])
              z.curr <- z.samples[j, , t]
              E.grouped <- apply(p.curr * z.curr, 1, sum, na.rm = TRUE)
              fit.big.y[j, , t] <- (y.grouped[, t] - E.grouped)^2 / (E.grouped + e)
              fit.y[j, t] <- sum(fit.big.y[j, , t])
              fit.big.y.rep[j, , t] <- (y.rep.grouped[j, , t] - E.grouped)^2 / (E.grouped + e)
              fit.y.rep[j, t] <- sum(fit.big.y.rep[j, , t])
            }
          } else if (fit.stat == 'freeman-tukey') {
            for (j in 1:n.samples) {
              p.curr <- matrix(det.prob[j, , t, ], nrow = J.long[q])
              z.curr <- z.samples[j, , t]
              E.grouped <- apply(p.curr * z.curr, 1, sum, na.rm = TRUE)
              fit.big.y[j, , t] <- (sqrt(y.grouped[, t]) - sqrt(E.grouped))^2 
              fit.y[j, t] <- sum(fit.big.y[j, , t])
              fit.big.y.rep[j, , t] <- (sqrt(y.rep.grouped[j, , t]) - sqrt(E.grouped))^2 
              fit.y.rep[j, t] <- sum(fit.big.y.rep[j, , t])
            }
          }
        }
      } else if (group == 2) { # Group by visit
        y.grouped <- apply(y[[q]], c(2, 3), sum, na.rm = TRUE)
        y.rep.grouped <- apply(y.rep.samples, c(1, 3, 4), sum, na.rm = TRUE)
        fit.big.y <- array(NA, dim = c(n.samples, n.years.long[q], dim(y[[q]])[3]))
        fit.big.y.rep <- array(NA, dim = c(n.samples, n.years.long[q], dim(y[[q]])[3]))
        for (t in 1:n.years.long[q]) {
          message(noquote(paste("Currently on time period ", t, " out of ", n.years.long[q], 
                                ' for data set ', q, sep = '')))
          if (fit.stat %in% c('chi-squared', 'chi-square')) {
            for (j in 1:n.samples) {
              p.curr <- matrix(det.prob[j, , t, ], nrow = J.long[q])
              z.curr <- z.samples[j, , t]
              E.grouped <- apply(p.curr * z.curr, 2, sum, na.rm = TRUE)
              fit.big.y[j, t, ] <- (y.grouped[t, ] - E.grouped)^2 / (E.grouped + e)
              fit.y[j, t] <- sum(fit.big.y[j, t, ])
              fit.big.y.rep[j, t, ] <- (y.rep.grouped[j, t, ] - E.grouped)^2 / (E.grouped + e)
              fit.y.rep[j, t] <- sum(fit.big.y.rep[j, t, ])
            }
          } else if (fit.stat == 'freeman-tukey') {
            for (j in 1:n.samples) {
              p.curr <- matrix(det.prob[j, , t, ], nrow = J.long[q])
              z.curr <- z.samples[j, , t]
              E.grouped <- apply(p.curr * z.curr, 2, sum, na.rm = TRUE)
              fit.big.y[j, t, ] <- (sqrt(y.grouped[t, ]) - sqrt(E.grouped))^2 
              fit.y[j, t] <- sum(fit.big.y[j, t, ])
              fit.big.y.rep[j, t, ] <- (sqrt(y.rep.grouped[j, t, ]) - sqrt(E.grouped))^2 
              fit.y.rep[j, t] <- sum(fit.big.y.rep[j, t, ])
            }
          }
        }
      }
      fit.y.list[[q]] <- fit.y
      fit.y.rep.list[[q]] <- fit.y.rep
      fit.y.group.quants.list[[q]] <- apply(fit.big.y, c(2, 3), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
      fit.y.rep.group.quants.list[[q]] <- apply(fit.big.y.rep, c(2, 3), 
                                                quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    }
    out$fit.y <- fit.y.list
    out$fit.y.rep <- fit.y.rep.list
    out$fit.y.group.quants <- fit.y.group.quants.list
    out$fit.y.rep.group.quants <- fit.y.rep.group.quants.list
    # For summaries
    out$group <- group
    out$seasons <- object$seasons
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$n.chains <- object$n.chains
  }

  class(out) <- 'ppcOcc'

  return(out)

}
