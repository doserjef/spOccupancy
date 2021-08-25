ppcOcc <- function(object, fit.stat, sub.sample, group, ...) {

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
  # Sub samples -----------------------
  n.samples <- object$n.samples
  if (missing(sub.sample)) {
    message("sub.sample is not specified. Using all posterior samples for posterior predictive check.")
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
  # Group -----------------------------
  if (missing(group)) {
    stop("error: group must be specified")
  }
  if (!(group %in% c(1, 2)) & class(object) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc')) {
    stop("error: group must be 1 (row) or 2 (columns) for objects of class PGOcc, spPGOcc, intPGOcc, spIntPGOcc")
  }

  out <- list()
  # For single species models
  if (class(object) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc')) {
    # Functions -------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
    y <- object$y
    X.p <- object$X.p
    p.det <- dim(X.p)[3]
    n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    J <- nrow(y)
    y.rep.samples <- object$y.rep.samples[, , s.indx, drop = FALSE]
    z.samples <- object$z.samples[s.indx, , drop = FALSE]
    X.p.long <- matrix(X.p, J * max(n.rep), p.det)
    alpha.samples <- object$alpha.samples[s.indx, , drop = FALSE]
    # Get detection probability
    det.prob <- logit.inv(X.p.long %*% t(alpha.samples))
    det.prob <- array(det.prob, dim(y.rep.samples))
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    e <- 0.0001
    # Do the stuff 
    if (group == 1) {
      y.grouped <- apply(y, 1, sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 3), sum, na.rm = TRUE)
      fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
      fit.big.y <- matrix(NA, length(y.grouped), n.samples)
      if (fit.stat == 'chi-square') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[, , i] * z.samples[i, ], 1, sum, na.rm = TRUE)
          fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (y.rep.grouped[, i] - E.grouped)^2 / (E.grouped + e)
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      } else if (fit.stat == 'freeman-tukey') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[, , i] * z.samples[i, ], 1, sum, na.rm = TRUE)
          fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[, i]) - sqrt(E.grouped))^2 
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      }
    } else if (group == 2) {
      y.grouped <- apply(y, 2, sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(2, 3), sum, na.rm = TRUE)
      fit.big.y <- matrix(NA, length(y.grouped), n.samples)
      fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
      if (fit.stat == 'chi-square') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[, , i] * z.samples[i, ], 2, sum, na.rm = TRUE)
          fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (y.rep.grouped[, i] - E.grouped)^2 / (E.grouped + e)
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      } else if (fit.stat == 'freeman-tukey') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[, , i] * z.samples[i, ], 2, sum, na.rm = TRUE)
          fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[, i]) - sqrt(E.grouped))^2 
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
    out$start <- start
    out$end <- end
    out$thin <- thin
    out$sample.size <- length(s.indx)
  }

  class(out) <- 'ppcOcc'

  return(out)

}
