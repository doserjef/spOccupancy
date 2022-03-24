# PGOcc -------------------------------------------------------------------
predict.PGOcc <- function(object, X.0, ignore.RE = FALSE, 
			  type = 'occupancy', ...) {
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
  if (!is(object, "PGOcc")) {
  # if (class(object) != "PGOcc") {
    stop("error: requires an output object of class PGOcc\n")
  }
  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  # Occurrence predictions ------------------------------------------------
  if (tolower(type) == 'occupancy') {  
    p.occ <- ncol(object$X)
    p.design <- p.occ
    if (object$psiRE & !ignore.RE) {
      p.design <- p.occ + ncol(object$sigma.sq.psi.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }

    # Composition sampling --------------------------------------------------
    beta.samples <- as.matrix(object$beta.samples)
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$psiRE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }

    if (object$psiRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0) %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.occ.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.occ.re)
      for (i in 1:p.occ.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.occ.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            beta.star.sites.0.samples[, j] <- 
              beta.star.samples[, X.re.ind[j, t]] + 
              beta.star.sites.0.samples[, j]
          } else {
            beta.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.psi.samples[, t])) + 
              beta.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.occ.re <- 0
    }
    J.str <- nrow(X.0)
    out$psi.0.samples <- mcmc(logit.inv(t(X.fix %*% t(beta.samples) + 
          				t(beta.star.sites.0.samples))))

    out$z.0.samples <- mcmc(matrix(rbinom(length(out$psi.0.samples), 1, c(out$psi.0.samples)), 
    		                 nrow = n.post, ncol = nrow(X.0)))

  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE & !ignore.RE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }

    # Composition sampling --------------------------------------------------
    alpha.samples <- as.matrix(object$alpha.samples)
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }

    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      indx <- which(colnames(X.0) %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(p.re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      alpha.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            alpha.star.sites.0.samples[, j] <- 
              alpha.star.samples[, X.re.ind[j, t]] + 
              alpha.star.sites.0.samples[, j]
          } else {
            alpha.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
              alpha.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    out$p.0.samples <- mcmc(logit.inv(t(X.fix %*% t(alpha.samples) + 
          				t(alpha.star.sites.0.samples))))
  }
  out$call <- cl

  class(out) <- "predict.PGOcc"
  out
}

print.PGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

fitted.PGOcc <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (!is(object, "PGOcc")) {
  # if (class(object) != "PGOcc") {
    stop("error: object must be one of class PGOcc\n")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  J <- nrow(y)
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y))]
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
  y <- c(y)
  y <- y[!is.na(y)]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  if (object$pRE) {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) +
      			      lambda.p %*% t(object$alpha.star.samples)))
  } else {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
  }
  y.rep.samples <- t(apply(det.prob.samples * z.samples[, z.long.indx], 
			   2, function(a) rbinom(n.post, 1, a)))
  tmp <- array(NA, dim = c(J * K.max, n.post))
  names.long <- which(!is.na(c(object$y)))
  tmp[names.long, ] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(J, K.max, n.post))
  y.rep.samples <- aperm(y.rep.samples, c(3, 1, 2))
  tmp <- array(NA, dim = c(J * K.max, n.post))
  tmp[names.long, ] <- det.prob.samples
  det.prob.samples <- array(tmp, dim = c(J, K.max, n.post))
  det.prob.samples <- aperm(det.prob.samples, c(3, 1, 2))
  out <- list()
  out$y.rep.samples <- y.rep.samples
  out$p.samples <- det.prob.samples
  return(out)
}

summary.PGOcc <- function(object,
			  quantiles = c(0.025, 0.5, 0.975), 
			  digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))
  
  # Occurrence
  cat("Occurrence (logit scale): \n")
  tmp.1 <- t(apply(object$beta.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$psiRE) {
    cat("\n")
    cat("Occurrence Random Effect Variances (logit scale): \n")
    tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")
  # Detection
  cat("Detection (logit scale): \n")
  tmp.1 <- t(apply(object$alpha.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$alpha.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$pRE) {
    cat("\n")
    cat("Detection Random Effect Variances (logit scale): \n")
    tmp.1 <- t(apply(object$sigma.sq.p.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.p.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

# ppcOcc ------------------------------------------------------------------ 
summary.ppcOcc <- function(object, level = 'both', 
			   digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:", deparse(object$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n\n", sep=""))

  if (object$class %in% c('PGOcc', 'spPGOcc')) {
    cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
    cat("Fit statistic: ", object$fit.stat, "\n")
  }

  if (object$class %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc')) {

    if (tolower(level) == 'community') {
      cat("----------------------------------------\n");
      cat("\tCommunity Level\n");
      cat("----------------------------------------\n");
      cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
      cat("Fit statistic: ", object$fit.stat, "\n")
    }

    if (tolower(level) == 'species') {
      cat("----------------------------------------\n");
      cat("\tSpecies Level\n");
      cat("----------------------------------------\n");
      N <- ncol(object$fit.y)
      for (i in 1:N) {
        cat(paste(object$sp.names[i], " Bayesian p-value: ", 
		  round(mean(object$fit.y.rep[, i] > object$fit.y[, i]), digits), "\n", sep = '')) 
      }
      cat("Fit statistic: ", object$fit.stat, "\n")
    }

    if (tolower(level) == 'both') {
      cat("----------------------------------------\n");
      cat("\tCommunity Level\n");
      cat("----------------------------------------\n");
      cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
      cat("\n")
      cat("----------------------------------------\n");
      cat("\tSpecies Level\n");
      cat("----------------------------------------\n");
      N <- ncol(object$fit.y)
      for (i in 1:N) {
        cat(paste(object$sp.names[i], " Bayesian p-value: ", 
		  round(mean(object$fit.y.rep[, i] > object$fit.y[, i]), digits), "\n", sep = '')) 
      }
      cat("Fit statistic: ", object$fit.stat, "\n")
    }
  }

  if (object$class %in% c('intPGOcc', 'spIntPGOcc')) {
    n.data <- length(object$fit.y.rep)
    for (q in 1:n.data) {
      cat("Data Source", q, "\n\n")	    
      cat("Bayesian p-value:", round(mean(object$fit.y.rep[[q]] > object$fit.y[[q]])), "\n")
      cat("Fit statistic:", object$fit.stat, "\n\n")
    }
  }

}

# spPGOcc -----------------------------------------------------------------

predict.spPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1, 
			    verbose = TRUE, n.report = 100, 
			    ignore.RE = FALSE, type = 'occupancy', ...) {

  ptm <- proc.time()
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
  #if (!is(object, c('spPGOcc', 'spIntPGOcc'))) {
  if (!(class(object) %in% c('spPGOcc', 'spIntPGOcc'))) {
    stop("error: requires an output object of class spPGOcc or spIntPGOcc\n")
  }

  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }
  
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))){
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  X.0 <- as.matrix(X.0)

  # Occurrence predictions ------------------------------------------------
  if (tolower(type) == 'occupancy') {
    if (missing(coords.0)) {
      stop("error: coords.0 must be specified\n")
    }
    if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
      stop("error: coords.0 must be a data.frame or matrix\n")
    }
    if (!ncol(coords.0) == 2){
      stop("error: coords.0 must have two columns\n")
    }
    coords.0 <- as.matrix(coords.0)
    n.post <- object$n.post * object$n.chains
    X <- object$X
    coords <- object$coords 
    J <- nrow(X)
    p.occ <- ncol(X)
    theta.samples <- object$theta.samples
    beta.samples <- object$beta.samples
    w.samples <- object$w.samples
    n.neighbors <- object$n.neighbors
    cov.model.indx <- object$cov.model.indx
    sp.type <- object$type
    if (object$psiRE & !ignore.RE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }
    if (ncol(X.0) != p.occ + p.occ.re){
      stop(paste("error: X.0 must have ", p.occ + p.occ.re," columns\n", sep = ''))
    }
    # Eliminate prediction sites that have already sampled been for now
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))
    coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
    X.0.new <- X.0[coords.0.indx, , drop = FALSE]

    if (length(coords.indx) == nrow(X.0)) {
      stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
    }

    if (object$psiRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0.new) %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
      X.re <- as.matrix(X.0.new[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0.new[, -indx, drop = FALSE])
      n.occ.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.occ.re)
      for (i in 1:p.occ.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.occ.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            beta.star.sites.0.samples[, j] <- 
              beta.star.samples[, X.re.ind[j, t]] + 
              beta.star.sites.0.samples[, j]
          } else {
            beta.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.psi.samples[, t])) + 
              beta.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0.new
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0.new))
      p.occ.re <- 0
    }

    # Sub-sample previous 
    theta.samples <- t(theta.samples)
    beta.samples <- t(beta.samples)
    w.samples <- t(w.samples)
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

    q <- nrow(X.0.new)

    if (sp.type == 'GP') {
    
      obs.pred.D <- iDist(coords, coords.0.new)
      obs.D <- iDist(coords)
      
      storage.mode(obs.pred.D) <- "double"
      storage.mode(obs.D) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(X.0.new) <- "double"
      storage.mode(q) <- "integer"
      storage.mode(beta.samples) <- "double"
      storage.mode(theta.samples) <- "double"
      storage.mode(w.samples) <- "double"
      storage.mode(beta.star.sites.0.samples) <- "double"
      storage.mode(n.post) <- "integer"
      storage.mode(cov.model.indx) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      storage.mode(verbose) <- "integer"
      storage.mode(n.report) <- "integer"
      storage.mode(n.omp.threads) <- "integer"

      out <- .Call("spPGOccPredict", J, p.occ, X.0.new, q, obs.D, 
          	 obs.pred.D, beta.samples, theta.samples, 
          	 w.samples, beta.star.sites.0.samples, 
          	 n.post, cov.model.indx, n.omp.threads, 
          	 verbose, n.report)
    } else { 
      # Get nearest neighbors 
      # nn2 is a function from RANN. 
      nn.indx.0 <- nn2(coords, coords.0.new, k=n.neighbors)$nn.idx-1 

      storage.mode(coords) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.0.new) <- "double"
      storage.mode(coords.0.new) <- "double"
      storage.mode(q) <- "integer"
      storage.mode(beta.samples) <- "double"
      storage.mode(theta.samples) <- "double"
      storage.mode(w.samples) <- "double"
      storage.mode(beta.star.sites.0.samples) <- "double"
      storage.mode(n.post) <- "integer"
      storage.mode(cov.model.indx) <- "integer"
      storage.mode(nn.indx.0) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      storage.mode(verbose) <- "integer"
      storage.mode(n.report) <- "integer"
      
      ptm <- proc.time()

      out <- .Call("spPGOccNNGPPredict", coords, J, p.occ, n.neighbors, 
                   X.0.new, coords.0.new, q, nn.indx.0, beta.samples, 
                   theta.samples, w.samples, beta.star.sites.0.samples, n.post, 
                   cov.model.indx, n.omp.threads, verbose, n.report)
    }

    if (nrow(X.0) == q) {
      out$z.0.samples <- mcmc(t(out$z.0.samples))
      out$psi.0.samples <- mcmc(t(out$psi.0.samples))
      out$w.0.samples <- mcmc(t(out$w.0.samples))
    } else {
      tmp <- matrix(NA, n.post, nrow(X.0))
      tmp[, coords.0.indx] <- t(out$z.0.samples)
      tmp[, coords.place.indx] <- object$z.samples[, coords.indx]
      out$z.0.samples <- mcmc(tmp)
      tmp <- matrix(NA, n.post, nrow(X.0))
      tmp[, coords.0.indx] <- t(out$psi.0.samples)
      tmp[, coords.place.indx] <- object$psi.samples[, coords.indx]
      out$psi.0.samples <- mcmc(tmp)
      tmp <- matrix(NA, n.post, nrow(X.0))
      tmp[, coords.0.indx] <- t(out$w.0.samples)
      tmp[, coords.place.indx] <- object$w.samples[, coords.indx]
      out$w.0.samples <- mcmc(tmp)
    }
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE & !ignore.RE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }

    # Composition sampling --------------------------------------------------
    alpha.samples <- as.matrix(object$alpha.samples)
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }

    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      indx <- which(colnames(X.0) %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(p.re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      alpha.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            alpha.star.sites.0.samples[, j] <- 
              alpha.star.samples[, X.re.ind[j, t]] + 
              alpha.star.sites.0.samples[, j]
          } else {
            alpha.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
              alpha.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    out$p.0.samples <- mcmc(logit.inv(t(X.fix %*% t(alpha.samples) + 
          				t(alpha.star.sites.0.samples))))
  }
  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)
  class(out) <- "predict.spPGOcc"
  out
}

print.spPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.spPGOcc <- function(object,
			    quantiles = c(0.025, 0.5, 0.975), 
			    digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))
  
  # Occurrence ------------------------
  cat("Occurrence (logit scale): \n")
  tmp.1 <- t(apply(object$beta.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$psiRE) {
    cat("\n")
    cat("Occurrence Random Effect Variances (logit scale): \n")
    tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")

  # Detection -------------------------
  cat("Detection (logit scale): \n")
  tmp.1 <- t(apply(object$alpha.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$alpha.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$pRE) {
    cat("\n")
    cat("Detection Random Effect Variances (logit scale): \n")
    tmp.1 <- t(apply(object$sigma.sq.p.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.p.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")
  # Covariance ------------------------
  cat("Spatial Covariance: \n")
  tmp.1 <- t(apply(object$theta.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$theta.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
}


fitted.spPGOcc <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (!is(object, 'spPGOcc')) {
  # if (class(object) != 'spPGOcc') {
    stop("error: object must be one of class spPGOcc\n")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  J <- nrow(y)
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y))]
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
  y <- c(y)
  y <- y[!is.na(y)]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  if (object$pRE) {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) +
      			      lambda.p %*% t(object$alpha.star.samples)))
  } else {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
  }
  y.rep.samples <- t(apply(det.prob.samples * z.samples[, z.long.indx], 
			   2, function(a) rbinom(n.post, 1, a)))
  tmp <- array(NA, dim = c(J * K.max, n.post))
  names.long <- which(!is.na(c(object$y)))
  tmp[names.long, ] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(J, K.max, n.post))
  y.rep.samples <- aperm(y.rep.samples, c(3, 1, 2))
  tmp <- array(NA, dim = c(J * K.max, n.post))
  tmp[names.long, ] <- det.prob.samples
  det.prob.samples <- array(tmp, dim = c(J, K.max, n.post))
  det.prob.samples <- aperm(det.prob.samples, c(3, 1, 2))
  out <- list()
  out$y.rep.samples <- y.rep.samples
  out$p.samples <- det.prob.samples
  return(out)
}

# msPGOcc -----------------------------------------------------------------

predict.msPGOcc <- function(object, X.0, ignore.RE = FALSE, 
			    type = 'occupancy', ...) {
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
  if (!is(object, 'msPGOcc')) {
  # if (class(object) != 'msPGOcc') {
    stop("error: requires an output object of class msPGOcc\n")
  }

  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }

  # Occurrence predictions ------------------------------------------------
  if (tolower(type) == 'occupancy') {
    p.occ <- ncol(object$X)
    p.design <- p.occ
    if (object$psiRE) {
      p.design <- p.occ + ncol(object$sigma.sq.psi.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    sp.indx <- rep(1:N, p.occ)
    n.post <- object$n.post * object$n.chains
    beta.samples <- as.matrix(object$beta.samples)
    out <- list()
    out$psi.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
    out$z.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
    if (object$psiRE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }
    if (object$psiRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0) %in% x.re.names)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.occ.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.occ.re)
      for (i in 1:p.occ.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      for (i in 1:N) {
        for (t in 1:p.occ.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                beta.star.samples[, (i - 1) * n.occ.re + X.re.ind[j, t]] + 
                beta.star.sites.0.samples[, (j - 1) * N + i]
            } else {
              beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.psi.samples[, t])) + 
                beta.star.sites.0.samples[, (j - 1) * N + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0))
      p.occ.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:N) {
      for (j in 1:J.str) {
        out$psi.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				     t(beta.samples[, sp.indx == i]) + 
                                               beta.star.sites.0.samples[, (j - 1) * N + i])
        out$z.0.samples[, i, j] <- rbinom(n.post, 1, out$psi.0.samples[, i, j])
      } # j
    } # i
  } # occurrence predictions
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    sp.indx <- rep(1:N, p.det)
    n.post <- object$n.post * object$n.chains
    alpha.samples <- as.matrix(object$alpha.samples)
    out <- list()
    out$p.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }
    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      indx <- which(colnames(X.0) %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(p.re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      for (i in 1:N) {
        for (t in 1:p.det.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                alpha.star.samples[, (i - 1) * n.det.re + X.re.ind[j, t]] + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            } else {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:N) {
      for (j in 1:J.str) {
        out$p.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				     t(alpha.samples[, sp.indx == i]) + 
                                               alpha.star.sites.0.samples[, (j - 1) * N + i])
      } # j
    } # i

  }
  out$call <- cl

  class(out) <- "predict.msPGOcc"
  out
}

print.msPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.msPGOcc <- function(object,
			    level = 'both',
			    quantiles = c(0.025, 0.5, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means (logit scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nOccurrence Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$psiRE) {
      cat("\n")
      cat("Occurrence Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
    cat("\n")
    # Detection
    cat("Detection Means (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\nDetection Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.alpha, round(object$ESS$tau.sq.alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.p.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.p.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence (logit scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\n")
    # Detection
    cat("Detection (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  }
}

fitted.msPGOcc <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  # if (!is(object, c("msPGOcc", "spMsPGOcc", "lfMsPGOcc", "sfMsPGOcc"))) {
  if (!(class(object) %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc'))) {
    stop("error: object must be of class msPGOcc, spMsPGOcc, lfMsPGOcc, or sfMsPGOcc\n")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  J <- dim(y)[2]
  N <- dim(y)[1]
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
  z.long.indx <- rep(1:J, K.max)
  z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  n.obs <- nrow(X.p)
  det.prob.samples <- array(NA, dim = c(n.obs, N, n.post))
  sp.indx <- rep(1:N, ncol(X.p))
  y <- matrix(y, N, J * K.max)
  y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
  if (object$pRE) {
    sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
    for (i in 1:N) {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + 
      				   lambda.p %*% t(object$alpha.star.samples[, sp.re.indx == i]))
    }
  } else {
    for (i in 1:N) {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }
  }
  out <- list()
  tmp <- array(NA, dim = c(n.post, N, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- det.prob.samples
  p.samples <- array(tmp, dim = c(n.post, N, J, K.max))
  out$p.samples <- p.samples
  # Need to be careful here that all arrays line up. 
  det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
  det.prob.samples <- det.prob.samples * z.samples[, , z.long.indx]
  y.rep.samples <- array(NA, dim = dim(det.prob.samples))
  for (i in 1:N) {
    y.rep.samples[, i, ] <- apply(det.prob.samples[, i, ], 2, function(a) rbinom(n.post, 1, a))
  }
  tmp <- array(NA, dim = c(n.post, N, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(n.post, N, J, K.max))
  out$y.rep.samples <- y.rep.samples
  return(out)
}

# spMsPGOcc ---------------------------------------------------------------
summary.spMsPGOcc <- function(object, 
			      level = 'both',
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means (logit scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nOccurrence Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$psiRE) {
      cat("\n")
      cat("Occurrence Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
    cat("\n")

    # Detection
    cat("Detection Means (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\nDetection Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.alpha, round(object$ESS$tau.sq.alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.p.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.p.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence (logit scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\n")
    # Detection
    cat("Detection (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\n")
    # Covariance
    cat("Spatial Covariance: \n")
    tmp.1 <- t(apply(object$theta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$theta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

fitted.spMsPGOcc <- function(object, ...) {
  fitted.msPGOcc(object)
}

print.spMsPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}


predict.spMsPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1, 
			      verbose = TRUE, n.report = 100, 
			      ignore.RE = FALSE, type = 'occupancy', ...) {

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
  if (!is(object, 'spMsPGOcc')) {
  # if (object != 'spMsPGOcc') {
    stop("error: requires an output object of class spMsPGOcc\n")
  }

  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }

  ptm <- proc.time()

  # Occurrence predictions ------------------------------------------------
  if (tolower(type == 'occupancy')) {
    n.post <- object$n.post * object$n.chains
    X <- object$X
    y <- object$y
    coords <- object$coords 
    J <- nrow(X)
    N <- dim(y)[1]
    p.occ <- ncol(X)
    theta.samples <- object$theta.samples
    beta.samples <- object$beta.samples
    w.samples <- object$w.samples
    n.neighbors <- object$n.neighbors
    cov.model.indx <- object$cov.model.indx
    sp.type <- object$type
    if (object$psiRE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }

    if (ncol(X.0) != p.occ + p.occ.re){
      stop(paste("error: X.0 must have ", p.occ + p.occ.re," columns\n", sep = ''))
    }
    X.0 <- as.matrix(X.0)

    if (missing(coords.0)) {
      stop("error: coords.0 must be specified\n")
    }
    if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
      stop("error: coords.0 must be a data.frame or matrix\n")
    }
    if (!ncol(coords.0) == 2){
      stop("error: coords.0 must have two columns\n")
    }
    coords.0 <- as.matrix(coords.0)

    # Eliminate prediction sites that have already been sampled for now
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))
    coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
    X.0.new <- X.0[coords.0.indx, , drop = FALSE]

    if (length(coords.indx) == nrow(X.0)) {
      stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
    }

    if (object$psiRE & !ignore.RE) {
      # Random effects of fitted values
      beta.star.samples <- object$beta.star.samples
      # Level names of the fitted random effects
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0.new) %in% x.re.names)
      # Random effect columns in predicted values
      X.re <- as.matrix(X.0.new[, indx, drop = FALSE])
      # Fixed effects in predicted values
      X.fix <- as.matrix(X.0.new[, -indx, drop = FALSE])
      # Number of random effect levels in the fitted values. 
      n.occ.re <- length(unlist(re.level.names))
      # Matrix that indicates which random effect level in 
      # beta.star.samples corresponds to the current 
      # random effect in the predicted values. 
      X.re.ind <- matrix(NA, nrow(X.re), p.occ.re)
      if (!ignore.RE) {
        for (i in 1:p.occ.re) {
          for (j in 1:nrow(X.re)) {
            # Which level in the fitted data equals the current level, and 
            # where is it located in beta.star.samples 	
            tmp <- which(re.level.names[[i]] == X.re[j, i])
            if (length(tmp) > 0) {
              if (i > 1) {
                X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
              } else {
                X.re.ind[j, i] <- tmp 
              }
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site. 
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      for (i in 1:N) {
        for (t in 1:p.occ.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                beta.star.samples[, (i - 1) * n.occ.re + X.re.ind[j, t]] + 
                beta.star.sites.0.samples[, (j - 1) * N + i]
            } else {
              beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.psi.samples[, t])) + 
                beta.star.sites.0.samples[, (j - 1) * N + i]
            }
          } # j
        } # t
      } # i
    } else {
      X.fix <- X.0.new
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0.new))
      p.occ.re <- 0
    }

    # Sub-sample previous
    theta.samples <- t(theta.samples)
    beta.samples <- t(beta.samples)
    w.samples <- aperm(w.samples, c(2, 3, 1))
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)
    
    q <- nrow(X.0.new)

    if (sp.type == 'GP') {
    
      obs.pred.D <- iDist(coords, coords.0.new)
      obs.D <- iDist(coords)
      
      storage.mode(obs.pred.D) <- "double"
      storage.mode(obs.D) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(N) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(X.0.new) <- "double"
      storage.mode(q) <- "integer"
      storage.mode(beta.samples) <- "double"
      storage.mode(theta.samples) <- "double"
      storage.mode(w.samples) <- "double"
      storage.mode(beta.star.sites.0.samples) <- "double"
      storage.mode(n.post) <- "integer"
      storage.mode(cov.model.indx) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      storage.mode(verbose) <- "integer"
      storage.mode(n.report) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      

      out <- .Call("spMsPGOccPredict", J, N, p.occ, X.0.new, q, obs.D, 
          	 obs.pred.D, beta.samples, theta.samples, 
          	 w.samples, beta.star.sites.0.samples, 
          	 n.post, cov.model.indx, n.omp.threads, 
          	 verbose, n.report)

    } else { 
      # Get nearest neighbors 
      # nn2 is a function from RANN. 
      nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1 

      storage.mode(coords) <- "double"
      storage.mode(N) <- "integer"
      storage.mode(J) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.fix) <- "double"
      storage.mode(coords.0.new) <- "double"
      storage.mode(q) <- "integer"
      storage.mode(beta.samples) <- "double"
      storage.mode(theta.samples) <- "double"
      storage.mode(w.samples) <- "double"
      storage.mode(beta.star.sites.0.samples) <- "double"
      storage.mode(n.post) <- "integer"
      storage.mode(cov.model.indx) <- "integer"
      storage.mode(nn.indx.0) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      storage.mode(verbose) <- "integer"
      storage.mode(n.report) <- "integer"
      
      ptm <- proc.time()

      out <- .Call("spMsPGOccNNGPPredict", coords, J, N, p.occ, n.neighbors, 
                   X.0.new, coords.0.new, q, nn.indx.0, beta.samples, 
                   theta.samples, w.samples, beta.star.sites.0.samples, n.post, 
                   cov.model.indx, n.omp.threads, verbose, n.report)

    }
    out$z.0.samples <- array(out$z.0.samples, dim = c(N, q, n.post))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$w.0.samples <- array(out$w.0.samples, dim = c(N, q, n.post))
    out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, q, n.post))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))

    # If some of the sites are sampled
    if (nrow(X.0) != q) {
      tmp <- array(NA, dim = c(n.post, N, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$z.0.samples
      tmp[, , coords.place.indx] <- object$z.samples[, , coords.indx]
      out$z.0.samples <- tmp
      tmp <- array(NA, dim = c(n.post, N, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$psi.0.samples
      tmp[, , coords.place.indx] <- object$psi.samples[, , coords.indx]
      out$psi.0.samples <- tmp
      tmp <- array(NA, dim = c(n.post, N, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$w.0.samples
      tmp[, , coords.place.indx] <- object$w.samples[, , coords.indx]
      out$w.0.samples <- tmp
    }
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    sp.indx <- rep(1:N, p.det)
    n.post <- object$n.post * object$n.chains
    alpha.samples <- as.matrix(object$alpha.samples)
    out <- list()
    out$p.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }
    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      indx <- which(colnames(X.0) %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(p.re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      for (i in 1:N) {
        for (t in 1:p.det.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                alpha.star.samples[, (i - 1) * n.det.re + X.re.ind[j, t]] + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            } else {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:N) {
      for (j in 1:J.str) {
        out$p.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				     t(alpha.samples[, sp.indx == i]) + 
                                               alpha.star.sites.0.samples[, (j - 1) * N + i])
      } # j
    } # i

  }

  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)

  class(out) <- "predict.spMsPGOcc"

  out

}

# intPGOcc ----------------------------------------------------------------
predict.intPGOcc <- function(object, X.0, ...) {
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
  if (!is(object, 'intPGOcc')) {
  # if (object != 'intPGOcc') {
    stop("error: requires an output object of class intPGOcc\n")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  p.occ <- ncol(object$X)
  if (ncol(X.0) != p.occ) {
    stop(paste("error: X.0 must have ", p.occ, " columns\n", sep = ''))
  }

  # Composition sampling --------------------------------------------------
  n.post <- object$n.post * object$n.chains
  beta.samples <- as.matrix(object$beta.samples)
  out <- list()
  out$psi.0.samples <- mcmc(logit.inv(t(as.matrix(X.0) %*% t(beta.samples))))
  out$z.0.samples <- mcmc(matrix(rbinom(length(out$psi.0.samples), 1, c(out$psi.0.samples)),
		      nrow = n.post, ncol = nrow(X.0)))
  out$call <- cl

  class(out) <- "predict.intPGOcc"
  out
}

print.intPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.intPGOcc <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  #if (!is(object, c('intPGOcc', 'spIntPGOcc'))) {
  if (!(class(object) %in% c('intPGOcc', 'spIntPGOcc'))) {
    stop("error: object must be one of class intPGOcc or spIntPGOcc\n")
  }

  y <- object$y
  n.data <- length(y)
  sites <- object$sites
  X.p <- object$X.p
  p.det.long <- sapply(X.p, function(a) dim(a)[2])
  n.rep <- sapply(y, function(a1) apply(a1, 1, function(a2) sum(!is.na(a2))))
  J.long <- sapply(y, nrow)
  det.prob <- list()

  for (q in 1:n.data) {
    y.rep.samples <- object$y.rep.samples[[q]]
    z.samples <- object$z.samples[, sites[[q]], drop = FALSE]
    alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
    alpha.samples <- object$alpha.samples[, alpha.indx.r == q, drop = FALSE]
    # Get detection probability
    det.prob[[q]] <- logit.inv(X.p[[q]] %*% t(alpha.samples))
    det.prob[[q]] <- array(det.prob[[q]], dim(y.rep.samples))
  }
  out <- list()
  out$y.rep.samples <- object$y.rep.samples
  out$p.samples <- det.prob
  return(out)
}

summary.intPGOcc <- function(object,
			     quantiles = c(0.025, 0.5, 0.975),
			     digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  n.data <- length(object$y)
  p.det.long <- sapply(object$X.p, function(a) dim(a)[[2]])

  # Occurrence ------------------------
  cat("Occurrence (logit scale): \n")
  tmp.1 <- t(apply(object$beta.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  cat("\n")
  # Detection -------------------------


  indx <- 1
  for (i in 1:n.data) {
    cat(paste("Data source ", i, " Detection (logit scale): \n", sep = ""))
    tmp.1 <- t(apply(object$alpha.samples[,indx:(indx+p.det.long[i] - 1), drop = FALSE], 2, 
		     function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples[,indx:(indx+p.det.long[i] - 1), drop = FALSE], 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha[indx:(indx+p.det.long[i] - 1)], 
		      round(object$ESS$alpha[indx:(indx+p.det.long[i] - 1)], 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    indx <- indx + p.det.long[i]
    cat("\n")
  }
}

# spIntPGOcc --------------------------------------------------------------
print.spIntPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.spIntPGOcc <- function(object, ...) {
  out <- fitted.intPGOcc(object, ...)
}

summary.spIntPGOcc <- function(object,
			       quantiles = c(0.025, 0.5, 0.975),
			       digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  n.data <- length(object$y)
  p.det.long <- sapply(object$X.p, function(a) dim(a)[[2]])

  # Occurrence ------------------------
  cat("Occurrence (logit scale): \n")
  tmp.1 <- t(apply(object$beta.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  cat("\n")
  # Detection -------------------------
  indx <- 1
  for (i in 1:n.data) {
    cat(paste("Data source ", i, " Detection (logit scale): \n", sep = ""))
    tmp.1 <- t(apply(object$alpha.samples[,indx:(indx+p.det.long[i] - 1), drop = FALSE], 2, 
		     function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples[,indx:(indx+p.det.long[i] - 1), drop = FALSE], 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha[indx:(indx+p.det.long[i] - 1)], 
		      round(object$ESS$alpha[indx:(indx+p.det.long[i] - 1)], 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    indx <- indx + p.det.long[i]
    cat("\n")
  }
  # Covariance ------------------------
  cat("Spatial Covariance: \n")
  tmp.1 <- t(apply(object$theta.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$theta.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
}

predict.spIntPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1,
			       verbose = TRUE, n.report = 100, ...) {
  out <- predict.spPGOcc(object, X.0, coords.0, n.omp.threads, 
			 verbose, n.report)
  class(out) <- "predict.spIntPGOcc"
  out
}

# lfMsPGOcc ---------------------------------------------------------------
print.lfMsPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.lfMsPGOcc <- function(object,
			      level = 'both',
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means (logit scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nOccurrence Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$psiRE) {
      cat("\n")
      cat("Occurrence Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
    cat("\n")

    # Detection
    cat("Detection Means (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\nDetection Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.alpha, round(object$ESS$tau.sq.alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.p.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.p.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence (logit scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\n")
    # Detection
    cat("Detection (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  }
}

fitted.lfMsPGOcc <- function(object, ...) {
  fitted.msPGOcc(object)
}

predict.lfMsPGOcc <- function(object, X.0, coords.0, ignore.RE = FALSE, 
			      type = 'occupancy', ...) {
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
  # if (!is(object, c('lfMsPGOcc', 'lfJSDM'))) {
  if (!(class(object) %in% c('lfMsPGOcc', 'lfJSDM'))) {
    stop("error: requires an output object of class lfMsPGOcc or lfJSDM\n")
  }

  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  
  # Occurrence predictions ------------------------------------------------
  if (tolower(type == 'occupancy')) {
    p.occ <- ncol(object$X)
    p.design <- p.occ
    if (object$psiRE) {
      p.design <- p.occ + ncol(object$sigma.sq.psi.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    J.0 <- nrow(X.0)
    q <- object$q
    coords <- object$coords
    sp.indx <- rep(1:N, p.occ)
    n.post <- object$n.post * object$n.chains
    beta.samples <- as.matrix(object$beta.samples)
    lambda.samples <- array(object$lambda.samples, dim = c(n.post, N, q))
    w.samples <- object$w.samples
    if (object$psiRE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }

    # Eliminate prediction sites that have already been sampled for now
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), 
          	      do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))
    coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
    X.0.new <- X.0[coords.0.indx, , drop = FALSE]

    if (length(coords.indx) == nrow(X.0)) {
      stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
    }

    if (object$psiRE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0.new) %in% x.re.names)
      X.re <- as.matrix(X.0.new[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0.new[, -indx, drop = FALSE])
      n.occ.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.occ.re)
      # TODO: need to double check this works with multiple random effects. 
      for (i in 1:p.occ.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      if (!ignore.RE) {
        for (i in 1:N) {
          for (t in 1:p.occ.re) {
            for (j in 1:nrow(X.re)) {
              if (!is.na(X.re.ind[j, t])) {
                beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                  beta.star.samples[, (i - 1) * n.occ.re + X.re.ind[j, t]] + 
                  beta.star.sites.0.samples[, (j - 1) * N + i]
              } else {
                beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                  rnorm(n.post, 0, sqrt(object$sigma.sq.psi.samples[, t])) + 
                  beta.star.sites.0.samples[, (j - 1) * N + i]
              }
            } # j
          } # t
        } # i 
      }
    } else {
      X.fix <- X.0.new
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0.new))
      p.occ.re <- 0
    }
    J.str <- nrow(X.0.new)
    # Create new random normal latent factors at unobserved sites. 
    w.0.samples <- array(rnorm(n.post * q * J.str), dim = c(n.post, q, J.str))
    w.star.0.samples <- array(NA, dim = c(n.post, N, J.str))

    for (i in 1:n.post) {
      w.star.0.samples[i, , ] <- matrix(lambda.samples[i, , ], N, q) %*%
                               matrix(w.0.samples[i, , ], q, J.str)
    }
    out <- list()
    out$psi.0.samples <- array(NA, dim = c(n.post, N, nrow(X.fix)))
    out$z.0.samples <- array(NA, dim = c(n.post, N, nrow(X.fix)))
    # Make predictions
    for (i in 1:N) {
      for (j in 1:J.str) {
      out$psi.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				   t(beta.samples[, sp.indx == i]) + 
          				   w.star.0.samples[, i, j] + 
                                             beta.star.sites.0.samples[, (j - 1) * N + i])
      out$z.0.samples[, i, j] <- rbinom(n.post, 1, out$psi.0.samples[, i, j])
    				     
      } # j
    } # i

    # If some of the sites are sampled
    if (nrow(X.0) != J.str) {
      tmp <- array(NA, dim = c(n.post, N, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$z.0.samples
      tmp[, , coords.place.indx] <- object$z.samples[, , coords.indx]
      out$z.0.samples <- tmp
      tmp <- array(NA, dim = c(n.post, N, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$psi.0.samples
      tmp[, , coords.place.indx] <- object$psi.samples[, , coords.indx]
      out$psi.0.samples <- tmp
    }
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    sp.indx <- rep(1:N, p.det)
    n.post <- object$n.post * object$n.chains
    alpha.samples <- as.matrix(object$alpha.samples)
    out <- list()
    out$p.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }
    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      indx <- which(colnames(X.0) %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(p.re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      for (i in 1:N) {
        for (t in 1:p.det.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                alpha.star.samples[, (i - 1) * n.det.re + X.re.ind[j, t]] + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            } else {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:N) {
      for (j in 1:J.str) {
        out$p.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				     t(alpha.samples[, sp.indx == i]) + 
                                               alpha.star.sites.0.samples[, (j - 1) * N + i])
      } # j
    } # i
  }

  out$call <- cl
  class(out) <- "predict.lfMsPGOcc"
  out
}

# sfMsPGOcc ---------------------------------------------------------------
print.sfMsPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.sfMsPGOcc <- function(object,
			      level = 'both',
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Occurrence Means (logit scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nOccurrence Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$psiRE) {
      cat("\n")
      cat("Occurrence Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
    cat("\n")
    # Detection
    cat("Detection Means (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\nDetection Variances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.alpha, round(object$ESS$tau.sq.alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$pRE) {
      cat("\n")
      cat("Detection Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.p.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.p.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Occurrence (logit scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\n")
    # Detection
    cat("Detection (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  }
    # Covariance
    cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpatial Covariance\n");
    cat("----------------------------------------\n");
    tmp.1 <- t(apply(object$theta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$theta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
}


fitted.sfMsPGOcc <- function(object, ...) {
  fitted.msPGOcc(object)
}

predict.sfMsPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1,
			      verbose = TRUE, n.report = 100, 
			      ignore.RE = FALSE, type = 'occupancy', ...) {

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
  # if (!is(object, c('sfMsPGOcc', 'sfJSDM'))) {
  if (!(class(object) %in% c('sfMsPGOcc', 'sfJSDM'))) {
    stop("error: requires an output object of class sfMsPGOcc or sfJSDM\n")
  }
  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }

  ptm <- proc.time()

  # Occurrence predictions ------------------------------------------------
  if (tolower(type == 'occupancy')) {
    n.post <- object$n.post * object$n.chains
    X <- object$X
    y <- object$y
    coords <- object$coords
    J <- nrow(X)
    N <- dim(y)[1]
    q <- object$q
    p.occ <- ncol(X)
    theta.samples <- object$theta.samples
    beta.samples <- object$beta.samples
    lambda.samples <- object$lambda.samples
    w.samples <- object$w.samples
    n.neighbors <- object$n.neighbors
    cov.model.indx <- object$cov.model.indx
    sp.type <- object$type
    if (object$psiRE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }

    if (ncol(X.0) != p.occ + p.occ.re){
      stop(paste("error: X.0 must have ", p.occ + p.occ.re," columns\n"))
    }
    X.0 <- as.matrix(X.0)

    if (missing(coords.0)) {
      stop("error: coords.0 must be specified\n")
    }
    if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
      stop("error: coords.0 must be a data.frame or matrix\n")
    }
    if (!ncol(coords.0) == 2){
      stop("error: coords.0 must have two columns\n")
    }
    coords.0 <- as.matrix(coords.0)

    # Eliminate prediction sites that have already been sampled for now
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))
    coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
    X.0.new <- X.0[coords.0.indx, , drop = FALSE]
    
    if (length(coords.indx) == nrow(X.0)) {
      stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
    }

    if (object$psiRE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0.new) %in% x.re.names)
      X.re <- as.matrix(X.0.new[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0.new[, -indx, drop = FALSE])
      n.occ.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.occ.re)
      # TODO: need to double check this works with multiple random effects. 
      if (!ignore.RE) {
        for (i in 1:p.occ.re) {
          for (j in 1:nrow(X.re)) {
            tmp <- which(re.level.names[[i]] == X.re[j, i])
            if (length(tmp) > 0) {
              if (i > 1) {
                X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
              } else {
                X.re.ind[j, i] <- tmp 
              }
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site. 
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      if (!ignore.RE) {
        for (i in 1:N) {
          for (t in 1:p.occ.re) {
            for (j in 1:nrow(X.re)) {
              if (!is.na(X.re.ind[j, t])) {
                beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                  beta.star.samples[, (i - 1) * n.occ.re + X.re.ind[j, t]] + 
                  beta.star.sites.0.samples[, (j - 1) * N + i]
              } else {
                beta.star.sites.0.samples[, (j - 1) * N + i] <- 
                  rnorm(n.post, 0, sqrt(object$sigma.sq.psi.samples[, t])) + 
                  beta.star.sites.0.samples[, (j - 1) * N + i]
              }
            } # j
          } # t
        } # i
      } 
    } else {
      X.fix <- X.0.new
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0.new))
      p.occ.re <- 0
    }

    # Sub-sample previous
    theta.samples <- t(theta.samples)
    lambda.samples <- t(lambda.samples)
    beta.samples <- t(beta.samples)
    w.samples <- aperm(w.samples, c(2, 3, 1))
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

    J.str <- nrow(X.0.new)

    if (sp.type == 'GP') {
      # Not currently implemented or accessed. 
    } else {
      # Get nearest neighbors
      # nn2 is a function from RANN.
      nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

      storage.mode(coords) <- "double"
      storage.mode(N) <- "integer"
      storage.mode(J) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.fix) <- "double"
      storage.mode(coords.0.new) <- "double"
      storage.mode(J.str) <- "integer"
      storage.mode(q) <- "integer"
      storage.mode(beta.samples) <- "double"
      storage.mode(theta.samples) <- "double"
      storage.mode(lambda.samples) <- "double"
      storage.mode(beta.star.sites.0.samples) <- "double"
      storage.mode(w.samples) <- "double"
      storage.mode(n.post) <- "integer"
      storage.mode(cov.model.indx) <- "integer"
      storage.mode(nn.indx.0) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      storage.mode(verbose) <- "integer"
      storage.mode(n.report) <- "integer"

      out <- .Call("sfMsPGOccNNGPPredict", coords, J, N, q, p.occ, n.neighbors,
                   X.0.new, coords.0.new, J.str, nn.indx.0, beta.samples,
                   theta.samples, lambda.samples, w.samples, 
          	 beta.star.sites.0.samples, n.post,
                   cov.model.indx, n.omp.threads, verbose, n.report)

    }
    out$z.0.samples <- array(out$z.0.samples, dim = c(N, J.str, n.post))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$w.0.samples <- array(out$w.0.samples, dim = c(q, J.str, n.post))
    out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, J.str, n.post))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))

    # If some of the sites are sampled
    if (nrow(X.0) != J.str) {
      tmp <- array(NA, dim = c(n.post, N, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$z.0.samples
      tmp[, , coords.place.indx] <- object$z.samples[, , coords.indx]
      out$z.0.samples <- tmp
      tmp <- array(NA, dim = c(n.post, N, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$psi.0.samples
      tmp[, , coords.place.indx] <- object$psi.samples[, , coords.indx]
      out$psi.0.samples <- tmp
      tmp <- array(NA, dim = c(n.post, q, nrow(X.0)))
      tmp[, , coords.0.indx] <- out$w.0.samples
      tmp[, , coords.place.indx] <- object$w.samples[, , coords.indx]
      out$w.0.samples <- tmp
    }
  } # occurrence predictions
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    sp.indx <- rep(1:N, p.det)
    n.post <- object$n.post * object$n.chains
    alpha.samples <- as.matrix(object$alpha.samples)
    out <- list()
    out$p.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0)))
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }
    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      indx <- which(colnames(X.0) %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(p.re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.re))
      for (i in 1:N) {
        for (t in 1:p.det.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                alpha.star.samples[, (i - 1) * n.det.re + X.re.ind[j, t]] + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            } else {
              alpha.star.sites.0.samples[, (j - 1) * N + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
                alpha.star.sites.0.samples[, (j - 1) * N + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:N) {
      for (j in 1:J.str) {
        out$p.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				     t(alpha.samples[, sp.indx == i]) + 
                                               alpha.star.sites.0.samples[, (j - 1) * N + i])
      } # j
    } # i
  }

  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)

  class(out) <- "predict.sfMsPGOcc"

  out

}

# lfJSDM ------------------------------------------------------------------
predict.lfJSDM <- function(object, X.0, coords.0, ignore.RE = FALSE, ...) {
  out <- predict.lfMsPGOcc(object, X.0, coords.0, ignore.RE = ignore.RE, ...)
  class(out) <- "predict.lfJSDM"
  out
}

print.lfJSDM <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

fitted.lfJSDM <- function(object, ...) {
  out <- list()
  out$z.samples <- object$z.samples
  out$psi.samples <- object$psi.samples
  return(out)
}

summary.lfJSDM <- function(object,
			   level = 'both',
			   quantiles = c(0.025, 0.5, 0.975),
			   digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Means (logit scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nVariances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$psiRE) {
      cat("\n")
      cat("Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Estimates (logit scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

# sfJSDM ------------------------------------------------------------------
predict.sfJSDM <- function(object, X.0, coords.0, n.omp.threads = 1, 
			   verbose = TRUE, n.report = 100, 
			   ignore.RE = FALSE, ...) {

  out <- predict.sfMsPGOcc(object, X.0, coords.0, n.omp.threads = n.omp.threads, 
			   verbose = verbose, n.report = n.report, 
			   ignore.RE = ignore.RE, ...)
  class(out) <- "predict.sfJSDM"
  out
}

fitted.sfJSDM <- function(object, ...) {
  return(fitted.lfJSDM(object, ...))
}


print.sfJSDM <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

summary.sfJSDM <- function(object,
			      level = 'both',
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Occurrence
    cat("Means (logit scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nVariances (logit scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$psiRE) {
      cat("\n")
      cat("Random Effect Variances (logit scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.psi.samples, 2, 
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.psi.samples, 2, 
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.psi, round(object$ESS$sigma.sq.psi, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Estimates (logit scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2, 
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2, 
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  # Covariance
  cat("\n")
  cat("----------------------------------------\n");
  cat("\tSpatial Covariance\n");
  cat("----------------------------------------\n");
  tmp.1 <- t(apply(object$theta.samples, 2, 
        	   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$theta.samples, 2, 
        	 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
}
