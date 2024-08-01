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
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  K.max <- dim(y)[2]
  J <- nrow(y)
  z.long.indx <- rep(1:J, dim(y)[2])
  z.long.indx <- z.long.indx[!is.na(c(y))]
  if (object$pRE) {
    if (nrow(object$X.p.re) == nrow(y)) {
      X.p.re <- do.call(rbind, replicate(ncol(y), object$X.p.re, simplify = FALSE))
      X.p.re <- X.p.re[!is.na(c(y)), , drop = FALSE]
      # Add 1 to get it to R indexing. 
      X.p.re <- X.p.re + 1
    } else {
      X.p.re <- object$X.p.re + 1
    }
  }

  if (nrow(X.p) == nrow(y)) {
    X.p <- do.call(rbind, replicate(ncol(y), X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y)), , drop = FALSE]
  } 
  y <- c(y)
  y <- y[!is.na(y)]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  tmp.samples <- matrix(0, n.post, length(y))
  if (object$pRE) {
    for (i in 1:ncol(X.p.re)) {
      tmp.samples <- tmp.samples + object$alpha.star.samples[, X.p.re[, i]]
    }
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) + t(tmp.samples)))
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
  tmp[names.long, ] <- t(det.prob.samples)
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

  # if (class(object) == 'tPGOcc') {
  if (is(object, 'tPGOcc')) {
    if (object$ar1) {
      cat("\n")
      cat("Occurrence AR(1) Temporal Covariance: \n")
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
}

residuals.PGOcc <- function(object, n.post.samples = 100, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  
  if (!(class(object) %in% c('PGOcc', 'spPGOcc', 'svcPGOcc'))) {
    stop('object must be of class PGOcc, spPGOcc, or svcPGOcc')
  }
  n.post <- object$n.post
  # Generate sub-sample of MCMC samples if relevant
  indx <- sample(1:n.post, n.post.samples, replace = FALSE) 
  # Occupancy residuals
  occ.resids <- object$z.samples[indx, , drop = FALSE] - 
                object$psi.samples[indx, , drop = FALSE]
  # Detection probability samples
  p.samples <- fitted.PGOcc(object)$p.samples[indx, , , drop = FALSE]
  # Form detection residuals
  det.resids <- array(NA, dim = c(n.post.samples, dim(object$y)[1], dim(object$y)[2]))
  for (l in 1:n.post.samples) {
    z.indx <- which(object$z.samples[indx[l], ] == 1) 
    # Don't need drop = FALSE b/c things will simplify if only 1 rep.
    det.resids[l, z.indx, ] <- object$y[z.indx, ] - p.samples[l, z.indx, ]
  }
  out <- list(occ.resids = occ.resids, 
              det.resids = det.resids) 
  out
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

  if (object$class %in% c('PGOcc', 'spPGOcc', 'svcPGOcc')) {
    cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
    cat("Fit statistic: ", object$fit.stat, "\n")
  }

  if (object$class %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc', 
			  'svcMsPGOcc')) {

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
      cat("Bayesian p-value:", round(mean(object$fit.y.rep[[q]] > object$fit.y[[q]]), digits), "\n")
      cat("Fit statistic:", object$fit.stat, "\n\n")
    }
  }

  if (object$class %in% c('tPGOcc', 'stPGOcc', 'svcTPGOcc')) {

    cat("----------------------------------------\n");
    cat("\tAll time periods combined\n");
    cat("----------------------------------------\n");
    cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
    cat("\n")
    cat("----------------------------------------\n");
    cat("\tIndividual time periods\n");
    cat("----------------------------------------\n");
    n.years.max <- ncol(object$fit.y)
    for (i in 1:n.years.max) {
      cat(paste("Time Period ", i,  " Bayesian p-value: ",
      	  round(mean(object$fit.y.rep[, i] > object$fit.y[, i]), digits), "\n", sep = ''))
    }
    cat("Fit statistic: ", object$fit.stat, "\n")
  }
  
  if (object$class %in% c('tIntPGOcc', 'stIntPGOcc', 'svcTIntPGOcc')) {
    n.data <- length(object$fit.y.rep)
    seasons <- object$seasons
    for (q in 1:n.data) {
      cat("\n----------------------------------------\n");
      cat("\tData source", q, "\n");
      cat("----------------------------------------\n");
      cat("----------------------------------------\n");
      cat("\tAll time periods combined\n");
      cat("----------------------------------------\n");
      cat("Bayesian p-value: ", round(mean(object$fit.y.rep[[q]] > object$fit.y[[q]]), digits), "\n")
      cat("\n")
      cat("----------------------------------------\n");
      cat("\tIndividual time periods\n");
      cat("----------------------------------------\n");
      n.years.max <- ncol(object$fit.y[[q]])
      for (i in 1:n.years.max) {
        cat(paste("Time Period ", seasons[[q]][i],  " Bayesian p-value: ",
        	  round(mean(object$fit.y.rep[[q]][, i] > object$fit.y[[q]][, i]), digits), "\n", sep = ''))
      }
    }
  }
  
  if (object$class %in% c('tMsPGOcc', 'svcTMsPGOcc', 'stMsPGOcc')) {
    n.years.max <- dim(object$fit.y)[3]

    if (tolower(level) == 'community' | tolower(level) == 'both') {
      cat("----------------------------------------\n");
      cat("\tCommunity Level\n");
      cat("----------------------------------------\n");
      for (i in 1:n.years.max) {
        cat(paste('Time Period ', i, " Bayesian p-value: ", round(mean(object$fit.y.rep[, , i] > object$fit.y[, , i]), digits), "\n", sep = ''))
      }
      if (tolower(level) == 'community') {
        cat("Fit statistic: ", object$fit.stat, "\n")
      }
    }

    if (tolower(level) == 'species' | tolower(level) == 'both') {
      cat("----------------------------------------\n");
      cat("\tSpecies Level\n");
      cat("----------------------------------------\n");
      N <- ncol(object$fit.y)
      for (i in 1:N) {
        cat("----------------------------------------\n");
        cat(paste('\tSpecies ', object$sp.names[i], '\n', sep = ''))
        cat("----------------------------------------\n");
        for (t in 1:n.years.max) {
          cat(paste('Time Period ', t, " Bayesian p-value: ", 
		    round(mean(object$fit.y.rep[, i, t] > object$fit.y[, i, t]), digits), "\n", sep = '')) 
	}
	cat(paste('Overall Bayesian p-value: ', 
		  round(mean(object$fit.y.rep[, i, ] > object$fit.y[, i, ]), digits), '\n',
		  sep = ''))
      }
      cat("Fit statistic: ", object$fit.stat, "\n")
    }
  }

}

# spPGOcc -----------------------------------------------------------------

predict.spPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1, 
			    verbose = TRUE, n.report = 100, 
			    ignore.RE = FALSE, type = 'occupancy',
                            grid.index.0, ...) {

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
  if (missing(grid.index.0)) {
    grid.index.0 <- 1:nrow(X.0)
  }
  grid.index.0.c <- grid.index.0 - 1

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
    J.w.0 <- nrow(coords.0)
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
    # Determine coordinates that have already been sampled
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))

    X.0.new <- X.0
    coords.0.new <- coords.0

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
      
      # Indicates whether a site has been sampled. 1 = sampled
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, J.w.0)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1
      
      storage.mode(obs.pred.D) <- "double"
      storage.mode(obs.D) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(X.fix) <- "double"
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
      storage.mode(sites.link) <- 'integer'
      storage.mode(sites.0.sampled) <- 'integer'
      
      out <- .Call("spPGOccPredict", J, p.occ, X.fix, q, obs.D, 
          	 obs.pred.D, beta.samples, theta.samples, 
          	 w.samples, beta.star.sites.0.samples, 
          	 n.post, cov.model.indx, n.omp.threads, 
          	 verbose, n.report, sites.link, sites.0.sampled)

      out$z.0.samples <- mcmc(t(out$z.0.samples))
      out$psi.0.samples <- mcmc(t(out$psi.0.samples))
      out$w.0.samples <- mcmc(t(out$w.0.samples))
    } else { 
      # Get nearest neighbors 
      # nn2 is a function from RANN. 
      nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1 

      # Indicates whether a site has been sampled. 1 = sampled
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, J.w.0)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1
      # Number of spatial REs for model fit
      J.w <- nrow(coords)

      storage.mode(coords) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(J.w) <- 'integer'
      storage.mode(J.w.0) <- 'integer'
      storage.mode(p.occ) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.fix) <- "double"
      storage.mode(coords.0) <- "double"
      storage.mode(grid.index.0.c) <- 'integer'
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
      storage.mode(sites.link) <- 'integer'
      storage.mode(sites.0.sampled) <- 'integer'
      
      ptm <- proc.time()

      out <- .Call("spPGOccNNGPPredict", coords, J, p.occ, n.neighbors, 
                   X.fix, coords.0, q, nn.indx.0, beta.samples, 
                   theta.samples, w.samples, beta.star.sites.0.samples, n.post, 
                   cov.model.indx, n.omp.threads, verbose, n.report, J.w.0, J.w, grid.index.0.c, 
                   sites.link, sites.0.sampled)
    
      out$z.0.samples <- mcmc(t(out$z.0.samples))
      out$psi.0.samples <- mcmc(t(out$psi.0.samples))
      out$w.0.samples <- mcmc(t(out$w.0.samples))
    }

  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.PGOcc(object, X.0, ignore.RE, type = 'detection')
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

  if (!(class(object) %in% c('svcPGBinom', 'svcTPGBinom'))) {
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
  }
  # Covariance ------------------------
  # if (class(object) == 'stPGOcc') {
  if (class(object) %in% c('stPGOcc', 'svcTPGOcc', 'svcTPGBinom')) {
    if (object$ar1) {
      cat("Spatio-temporal Covariance: \n")
    } else {
      cat("Spatial Covariance: \n")
    }
  } else {
    cat("Spatial Covariance: \n")
  }
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
  fitted.PGOcc(object)
}

residuals.spPGOcc <- function(object, n.post.samples = 100, ...) {
  residuals.PGOcc(object, n.post.samples)  
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
    if (is(object, 'intMsPGOcc')) {
      N <- object$N
    } else {
      N <- dim(object$y)[1]
    } 
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
    if (object$pRE & !ignore.RE) {
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

  if (is(object, 'tMsPGOcc')) {
    if (object$ar1) {
      cat("\n")
      cat("Occurrence AR(1) Temporal Covariance: \n")
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
  if (!(class(object) %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc', 
			     'svcMsPGOcc'))) {
    stop("error: object must be of class msPGOcc, spMsPGOcc, lfMsPGOcc, sfMsPGOcc, or svcMsPGOcc\n")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- dim(y)[3]
  J <- dim(y)[2]
  N <- dim(y)[1]
  if (object$pRE) {
    if (nrow(object$X.p.re) == dim(y)[2]) {
      X.p.re <- do.call(rbind, replicate(dim(y)[3], object$X.p.re, simplify = FALSE))
      X.p.re <- X.p.re[!is.na(c(y[1, , ])), , drop = FALSE]
      # Add 1 to get it to R indexing. 
      X.p.re <- X.p.re + 1
    } else {
      # Add 1 to get it to R indexing.
      X.p.re <- object$X.p.re + 1
    }
  }
  if (nrow(X.p) == dim(y)[2]) {
    X.p <- do.call(rbind, replicate(dim(y)[3], X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y[1, , ])), , drop = FALSE]
  }
  z.long.indx <- rep(1:J, dim(y)[3])
  z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  n.obs <- nrow(X.p)
  det.prob.samples <- array(NA, dim = c(n.obs, N, n.post))
  sp.indx <- rep(1:N, ncol(X.p))
  y <- matrix(y, N, J * K.max)
  y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
  for (i in 1:N) {
    if (object$pRE) {
      sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
      tmp.samples <- matrix(0, n.post, n.obs)
      tmp.alpha.star <- object$alpha.star.samples[, sp.re.indx == i]
      for (j in 1:ncol(X.p.re)) {
        tmp.samples <- tmp.samples + tmp.alpha.star[, X.p.re[, j]]
      }
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + t(tmp.samples))
    } else {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }
  }

  out <- list()
  # Get detection probability
  # Need to be careful here that all arrays line up. 
  det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
  tmp <- array(NA, dim = c(n.post, N, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- det.prob.samples
  p.samples <- array(tmp, dim = c(n.post, N, J, K.max))
  out$p.samples <- p.samples
  # Get fitted values
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
    if (object$psiRE & !ignore.RE) {
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

    # Determine coordinates that have already been sampled
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))
    
    X.0.new <- X.0
    coords.0.new <- coords.0
    
    if (object$psiRE & !ignore.RE) {
      # Random effects of fitted values
      beta.star.samples <- object$beta.star.samples
      # Level names of the fitted random effects
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0.new) %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
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
      
      # Indicates whether a site has been sampled. 1 = sampled
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, q)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1
      
      storage.mode(obs.pred.D) <- "double"
      storage.mode(obs.D) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(N) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(X.fix) <- "double"
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
      storage.mode(sites.link) <- 'integer'
      storage.mode(sites.0.sampled) <- 'integer'

      out <- .Call("spMsPGOccPredict", J, N, p.occ, X.fix, q, obs.D, 
          	 obs.pred.D, beta.samples, theta.samples, 
          	 w.samples, beta.star.sites.0.samples, 
          	 n.post, cov.model.indx, n.omp.threads, 
          	 verbose, n.report, sites.link, sites.0.sampled)
      
      out$z.0.samples <- array(out$z.0.samples, dim = c(N, q, n.post))
      out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
      out$w.0.samples <- array(out$w.0.samples, dim = c(N, q, n.post))
      out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
      out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, q, n.post))
      out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
      
    } else { 
      # Get nearest neighbors 
      # nn2 is a function from RANN. 
      nn.indx.0 <- nn2(coords, coords.0.new, k=n.neighbors)$nn.idx-1 
      
      # Indicates whether a site has been sampled. 1 = sampled
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, q)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1
      
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
      storage.mode(sites.link) <- 'integer'
      storage.mode(sites.0.sampled) <- 'integer'
      
      ptm <- proc.time()

      out <- .Call("spMsPGOccNNGPPredict", coords, J, N, p.occ, n.neighbors, 
                   X.fix, coords.0.new, q, nn.indx.0, beta.samples, 
                   theta.samples, w.samples, beta.star.sites.0.samples, n.post, 
                   cov.model.indx, n.omp.threads, verbose, n.report, sites.link,
                   sites.0.sampled)

    }
    out$z.0.samples <- array(out$z.0.samples, dim = c(N, q, n.post))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$w.0.samples <- array(out$w.0.samples, dim = c(N, q, n.post))
    out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, q, n.post))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.msPGOcc(object, X.0, ignore.RE, type)
  }

  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)

  class(out) <- "predict.spMsPGOcc"

  out

}

# intPGOcc ----------------------------------------------------------------
predict.intPGOcc <- function(object, X.0, ignore.RE = FALSE, 
			     type = 'occupancy', ...) {
  # Occupancy predictions -------------------------------------------------
  if (tolower(type == 'occupancy')) {	
    out <- predict.PGOcc(object, X.0, ignore.RE, type)
  }

  # Detection predictions -------------------------------------------------
  if (tolower(type == 'detection')) {
    stop("detection prediction is not currently implemented.")
    out <- list()
  }

  class(out) <- "predict.intPGOcc"
  out
}

print.intPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

# TODO: should check this a bit more.
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

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (!(class(object) %in% c('intPGOcc', 'spIntPGOcc'))) {
    stop("error: object must be one of class intPGOcc or spIntPGOcc\n")
  }

  y <- object$y
  n.data <- length(y)
  sites <- object$sites
  X.p <- object$X.p
  y.rep.samples <- list()
  for (i in 1:n.data) {
    y.rep.samples[[i]] <- array(NA, dim = dim(object$p.samples[[i]]))
  }
  n.post <- object$n.post * object$n.chains
  z.long.indx.r <- list()
  for (i in 1:n.data) {
    z.long.indx.r[[i]] <- rep(sites[[i]], dim(y[[i]])[2])
    z.long.indx.r[[i]] <- z.long.indx.r[[i]][!is.na(c(y[[i]]))]
  }

  for (q in 1:n.data) {
    for (j in 1:ncol(y.rep.samples[[q]])) {
      for (k in 1:dim(y.rep.samples[[q]])[3]) {
        if (sum(!is.na(object$p.samples[[q]][, j, k])) != 0) {
          y.rep.samples[[q]][, j, k] <- rbinom(n.post, 1, 
	  				     object$z.samples[, sites[[q]][j]] * 
                                               object$p.samples[[q]][, j, k])
	}
      }
    }
  }

  out <- list()
  out$y.rep.samples <- y.rep.samples
  out$p.samples <- object$p.samples
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
  p.det.re.long <- sapply(object$X.p.re, function(a) dim(a)[[2]])

  # Occurrence ------------------------
  cat("----------------------------------------\n")
  cat("Occurrence\n")
  cat("----------------------------------------\n")
  cat("Fixed Effects (logit scale):\n")
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
    cat("Random Effect Variances (logit scale):\n")
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
  indx <- 1
  indx.re <- 1
  for (i in 1:n.data) {
    cat("----------------------------------------\n")
    cat(paste("Data source ", i, " Detection\n", sep = ""))
    cat("----------------------------------------\n")
    cat("Fixed Effects (logit scale):\n")
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
    if (object$pRELong[i]) {
      tmp.samples <- object$sigma.sq.p.samples[, indx.re:(indx.re+p.det.re.long[i] - 1), 
					       drop = FALSE]
      cat("Random Effect Variances (logit scale):\n")
      tmp.1 <- t(apply(tmp.samples, 2, function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(tmp.samples, 2, function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p[indx.re:(indx.re+p.det.re.long[i] - 1)], 
			round(object$ESS$sigma.sq.p[indx.re:(indx.re+p.det.re.long[i] - 1)],
			      0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
      cat("\n")
      indx.re <- indx.re + p.det.re.long[i]
    }
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
  p.det.re.long <- sapply(object$X.p.re, function(a) dim(a)[[2]])

  # Occurrence ------------------------
  cat("----------------------------------------\n")
  cat("Occurrence\n")
  cat("----------------------------------------\n")
  cat("Fixed Effects (logit scale):\n")
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
    cat("Random Effect Variances (logit scale):\n")
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
  indx <- 1
  indx.re <- 1
  for (i in 1:n.data) {
    cat("----------------------------------------\n")
    cat(paste("Data source ", i, " Detection\n", sep = ""))
    cat("----------------------------------------\n")
    cat("Fixed Effects (logit scale):\n")
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
    if (object$pRELong[i]) {
      tmp.samples <- object$sigma.sq.p.samples[, indx.re:(indx.re+p.det.re.long[i] - 1), 
					       drop = FALSE]
      cat("Random Effect Variances (logit scale):\n")
      tmp.1 <- t(apply(tmp.samples, 2, function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(tmp.samples, 2, function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p[indx.re:(indx.re+p.det.re.long[i] - 1)], 
			round(object$ESS$sigma.sq.p[indx.re:(indx.re+p.det.re.long[i] - 1)],
			      0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
      cat("\n")
      indx.re <- indx.re + p.det.re.long[i]
    }
  }
  # Covariance ------------------------
  cat("----------------------------------------\n")
  cat("Spatial Covariance\n")
  cat("----------------------------------------\n")
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
			       verbose = TRUE, n.report = 100, 
			       ignore.RE = FALSE, type = 'occupancy', ...) {
  # Occupancy predictions -------------------------------------------------
  if (tolower(type == 'occupancy')) {	
    out <- predict.spPGOcc(object, X.0, coords.0, n.omp.threads, 
			   verbose, n.report, ignore.RE, type)
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type == 'detection')) {
    stop("detection prediction is not currently implemented.")
    out <- list()
  }

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
			      type = 'occupancy', include.w = TRUE, ...) {
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
  if (!include.w) {
    out <- predict.msPGOcc(object, X.0, ignore.RE, type)
  } else {
    if (tolower(type == 'occupancy')) {
      p.occ <- ncol(object$X)
      p.design <- p.occ
      if (object$psiRE & !ignore.RE) {
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
      if (q > 0) {
        lambda.samples <- array(object$lambda.samples, dim = c(n.post, N, q))
        w.samples <- object$w.samples
      }
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
      if (q > 0) {
        w.0.samples <- array(rnorm(n.post * q * J.str), dim = c(n.post, q, J.str))
        w.star.0.samples <- array(NA, dim = c(n.post, N, J.str))
        for (i in 1:n.post) {
          w.star.0.samples[i, , ] <- matrix(lambda.samples[i, , ], N, q) %*%
                                   matrix(w.0.samples[i, , ], q, J.str)
        }
      }

      out <- list()
      out$psi.0.samples <- array(NA, dim = c(n.post, N, nrow(X.fix)))
      out$z.0.samples <- array(NA, dim = c(n.post, N, nrow(X.fix)))
      # Make predictions
      for (i in 1:N) {
        for (j in 1:J.str) {
          if (q > 0) {
            out$psi.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
            				         t(beta.samples[, sp.indx == i]) + 
            				         w.star.0.samples[, i, j] + 
                                                   beta.star.sites.0.samples[, (j - 1) * N + i])
          } else {
            out$psi.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
            				         t(beta.samples[, sp.indx == i]) + 
                                                   beta.star.sites.0.samples[, (j - 1) * N + i])
          }
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
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.msPGOcc(object, X.0, ignore.RE, type)
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
    if (class(object) %in% c('stMsPGOcc', 'svcTMsPGOcc')) {
      cat("----------------------------------------\n");
      cat("\tSpatio-temporal Covariance: \n")
      cat("----------------------------------------\n");
    } else {
      cat("----------------------------------------\n");
      cat("\tSpatial Covariance\n");
      cat("----------------------------------------\n");
    }
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
			      ignore.RE = FALSE, type = 'occupancy', 
			      grid.index.0, ...) {

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
  # Check grid ------------------------------------------------------------
  if (missing(grid.index.0)) {
    grid.index.0 <- 1:nrow(X.0)
  }
  grid.index.0.c <- grid.index.0 - 1

  ptm <- proc.time()

  # Occurrence predictions ------------------------------------------------
  if (tolower(type == 'occupancy')) {
    n.post <- object$n.post * object$n.chains
    X <- object$X
    y <- object$y
    coords <- object$coords
    J <- nrow(X)
    J.w.0 <- nrow(coords.0)
    N <- dim(y)[1]
    q <- object$q
    p.occ <- ncol(X)
    theta.samples <- object$theta.samples
    beta.samples <- object$beta.samples
    lambda.samples <- object$lambda.samples
    w.samples <- object$w.samples
    n.neighbors <- object$n.neighbors
    cov.model.indx <- object$cov.model.indx
    std.by.sp <- object$std.by.sp
    species.sds <- object$species.sds
    species.means <- object$species.means
    sp.type <- object$type
    if (object$psiRE & !ignore.RE) {
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

    if (object$psiRE) {
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
      X.fix <- X.0
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0))
      p.occ.re <- 0
    }

    # Sub-sample previous
    theta.samples <- t(theta.samples)
    lambda.samples <- t(lambda.samples)
    beta.samples <- t(beta.samples)
    w.samples <- aperm(w.samples, c(2, 3, 1))
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

    # Logical indicating whether fitting a shared spatial model
    if (is(object, 'sfJSDM')) {
      if (object$shared.spatial) {
        shared.spatial <- TRUE
      } else {
        shared.spatial <- FALSE
      }
    } else {
      shared.spatial <- FALSE
    }

    # Get stuff for linking sampled sites to predicted sites
    sites.0.indx <- 0:(nrow(X.0) - 1)
    J.0 <- length(unique(sites.0.indx))
    sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
    sites.link <- rep(NA, J.0)
    sites.link[which(!is.na(match.indx))] <- coords.indx
    # For C
    sites.link <- sites.link - 1
    J.w <- nrow(coords)
    J.str <- nrow(X.0)

    # Get an individual set of covariates for each species to 
    # account for the potential of having different ranges
    # for each species.
    X.big <- array(NA, dim = c(J.str, ncol(X.fix), N))
    for (i in 1:N) {
      X.big[, , i] <- X.fix
      if (std.by.sp) {
        for (r in 1:ncol(X.fix)) {
          if (!is.na(species.sds[i, r])) {
            X.big[, r, i] <- (X.big[, r, i] - species.means[i, r]) / species.sds[i, r]
          }
        }
      }
    }

    if (sp.type == 'GP') {
      # Not currently implemented or accessed. 
    } else {
      # Get nearest neighbors
      # nn2 is a function from RANN.
      nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

      storage.mode(coords) <- "double"
      storage.mode(N) <- "integer"
      storage.mode(J) <- "integer"
      storage.mode(J.w) <- 'integer'
      storage.mode(J.w.0) <- 'integer'
      storage.mode(p.occ) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.big) <- "double"
      storage.mode(coords.0) <- "double"
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
      storage.mode(shared.spatial) <- 'integer'
      storage.mode(sites.link) <- "integer"
      storage.mode(sites.0.sampled) <- 'integer'
      storage.mode(grid.index.0.c) <- 'integer'

      out <- .Call("sfMsPGOccNNGPPredict", coords, J, N, q, p.occ, n.neighbors,
                   X.big, coords.0, J.str, nn.indx.0, beta.samples,
                   theta.samples, lambda.samples, w.samples, 
          	   beta.star.sites.0.samples, n.post,
                   cov.model.indx, n.omp.threads, verbose, n.report, shared.spatial, 
                   sites.link, sites.0.sampled, J.w.0, J.w, grid.index.0.c)

    }
    out$z.0.samples <- array(out$z.0.samples, dim = c(N, J.str, n.post))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$w.0.samples <- array(out$w.0.samples, dim = c(q, J.w.0, n.post))
    out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, J.str, n.post))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
  } # occurrence predictions
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.msPGOcc(object, X.0, ignore.RE, type)
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

# stPGOcc ----------------------------------------------------------------
print.stPGOcc <- function(x, ...) {
  print.spPGOcc(x)
}

summary.stPGOcc <- function(object,
			    quantiles = c(0.025, 0.5, 0.975), 
			    digits = max(3L, getOption("digits") - 3L), ...) {
  summary.spPGOcc(object, quantiles, digits)
}

fitted.stPGOcc <- function(object, ...) {
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
  if (!(class(object) %in% c('tPGOcc', 'stPGOcc', 'svcTPGOcc'))) {
    stop("error: object must be one of class tPGOcc or stPGOcc\n")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.years.max <- dim(y)[2]
  K.max <- dim(y)[3]
  J <- nrow(y)
  z.long.indx <- rep(1:(J * n.years.max), dim(y)[3])
  z.long.indx <- z.long.indx[!is.na(c(y))]
  y <- c(y)
  y <- y[!is.na(y)]
  z.samples <- object$z.samples
  z.samples <- aperm(z.samples, c(2, 3, 1))
  z.samples <- matrix(z.samples, nrow = J * n.years.max, ncol = n.post)
  alpha.samples <- object$alpha.samples
  tmp.samples <- matrix(0, n.post, length(y))
  if (object$pRE) {
    # Add 1 to get it to R indexing. 
    X.p.re <- object$X.p.re + 1
    for (i in 1:ncol(X.p.re)) {
      tmp.samples <- tmp.samples + object$alpha.star.samples[, X.p.re[, i]]
    }
    det.prob.samples <- logit.inv(X.p %*% t(alpha.samples) + t(tmp.samples))
  } else {
    det.prob.samples <- logit.inv(X.p %*% t(alpha.samples))
  }
  y.rep.samples <- t(apply(det.prob.samples * z.samples[z.long.indx, ], 
			   1, function(a) rbinom(n.post, 1, a)))

  tmp <- array(NA, dim = c(J * K.max * n.years.max, n.post))
  names.long <- which(!is.na(c(object$y)))
  tmp[names.long, ] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(J, n.years.max, K.max, n.post))
  y.rep.samples <- aperm(y.rep.samples, c(4, 1, 2, 3))
  tmp <- array(NA, dim = c(J * K.max * n.years.max, n.post))
  tmp[names.long, ] <- det.prob.samples
  det.prob.samples <- array(tmp, dim = c(J, n.years.max, K.max, n.post))
  det.prob.samples <- aperm(det.prob.samples, c(4, 1, 2, 3))
  out <- list()
  out$y.rep.samples <- y.rep.samples
  out$p.samples <- det.prob.samples
  return(out)
}

predict.stPGOcc <- function(object, X.0, coords.0, t.cols, n.omp.threads = 1,
                            verbose = TRUE, n.report = 100,
                            ignore.RE = FALSE, type = 'occupancy', 
                            forecast = FALSE, grid.index.0, ...) {

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

  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (length(dim(X.0)) != 3) {
    stop("error: X.0 must be an array with three dimensions corresponding to site, time, and covariate")
  }

  if (missing(t.cols) & forecast == FALSE) {
    stop("error: t.cols must be specified when forecast = FALSE\n")
  }
  if (missing(grid.index.0)) {
    grid.index.0 <- 1:nrow(X.0)
  }
  grid.index.0.c <- grid.index.0 - 1

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
    J <- dim(X)[1]
    J.w.0 <- nrow(coords.0)
    n.years.max <- dim(X.0)[2]
    p.occ <- dim(X)[3]
    theta.samples <- object$theta.samples
    beta.samples <- object$beta.samples
    w.samples <- object$w.samples
    n.neighbors <- object$n.neighbors
    cov.model.indx <- object$cov.model.indx
    sp.type <- object$type
    # Get AR1 random effect values for the corresponding years. 
    ar1 <- object$ar1
    if (forecast & ar1) {
      message("NOTE: forecasting in spOccupancy does not currently use AR(1) random effects\n")
    }
    if (ar1 & !forecast) {
      eta.samples <- object$eta.samples[, t.cols, drop = FALSE]
    } else {
      eta.samples <- matrix(0, n.post, n.years.max)
    }
    if (object$psiRE & !ignore.RE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }
    if (dim(X.0)[3] != p.occ + p.occ.re){
      stop(paste("error: the third dimension of X.0 must be ", p.occ + p.occ.re,"\n", sep = ''))
    }
    # Eliminate prediction sites that have already sampled been for now
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))

    if (object$psiRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get elements in design matrix with random effects
      x.re.names <- dimnames(object$X.re)[[3]]
      indx <- which(dimnames(X.0)[[3]] %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: dimnames(X.0)[[3]] must match variable names in data$occ.covs")
      }
      X.re <- X.0[, , indx, drop = FALSE]
      X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
		     ncol = dim(X.re)[3])
      X.fix <- X.0[, , -indx, drop = FALSE]
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
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
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.fix))
      p.occ.re <- 0
    }

    # Sub-sample previous
    if (ar1) {
      if (object$cov.model.indx == 2) { # matern == 2
        theta.samples <- t(theta.samples[, 1:3])
      } else {
        theta.samples <- t(theta.samples[, 1:2])
      }
    } else {
      theta.samples <- t(theta.samples)
    }
    beta.samples <- t(beta.samples)
    w.samples <- t(w.samples)
    eta.samples <- t(eta.samples)
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

    q <- nrow(X.fix) / n.years.max

    sites.0.indx <- 0:(nrow(X.0) - 1)
    J.0 <- length(unique(sites.0.indx))
    sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
    sites.link <- rep(NA, J.0)
    sites.link[which(!is.na(match.indx))] <- coords.indx
    # For C
    sites.link <- sites.link - 1
    J.w <- nrow(coords)
    
    # Check if sampled sites are included and make sure predicting across
    # all years.
    # if ((nrow(coords.0) != q) & n.years.max != dim(object$X)[2]) {
    #   stop("error: when predicting at sampled sites using stPGOcc, you must predict across all primary time periods")
    # }

    if (sp.type == 'GP') {
      stop("NNGP = FALSE is not currently supported for stPGOcc")
    } else {
      # Get nearest neighbors
      # nn2 is a function from RANN.
      nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

      storage.mode(coords) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(J.w) <- 'integer'
      storage.mode(J.w.0) <- 'integer'
      storage.mode(n.years.max) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.fix) <- "double"
      storage.mode(sites.link) <- "integer"
      storage.mode(sites.0.sampled) <- 'integer'
      storage.mode(coords.0) <- "double"
      storage.mode(grid.index.0.c) <- 'integer'
      storage.mode(q) <- "integer"
      storage.mode(beta.samples) <- "double"
      storage.mode(theta.samples) <- "double"
      storage.mode(w.samples) <- "double"
      storage.mode(beta.star.sites.0.samples) <- "double"
      storage.mode(eta.samples) <- "double"
      storage.mode(n.post) <- "integer"
      storage.mode(cov.model.indx) <- "integer"
      storage.mode(nn.indx.0) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      storage.mode(verbose) <- "integer"
      storage.mode(n.report) <- "integer"

      ptm <- proc.time()

      out <- .Call("stPGOccNNGPPredict", coords, J, n.years.max, p.occ, n.neighbors,
                   X.fix, coords.0, q, nn.indx.0, beta.samples,
                   theta.samples, w.samples, beta.star.sites.0.samples, eta.samples, 
		   sites.link, sites.0.sampled, n.post,
                   cov.model.indx, n.omp.threads, verbose, n.report, J.w.0, J.w,
                   grid.index.0.c)
    }
    out$z.0.samples <- array(out$z.0.samples, dim = c(q, n.years.max, n.post))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(q, n.years.max, n.post))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
    out$w.0.samples <- mcmc(t(out$w.0.samples))
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.tPGOcc(object, X.0, t.cols, ignore.RE, type)
  }
  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)
  class(out) <- "predict.stPGOcc"
  out
}
# tPGOcc ------------------------------------------------------------------
print.tPGOcc <- function(x, ...) {
  print.PGOcc(x)
}

summary.tPGOcc <- function(object, quantiles = c(0.025, 0.5, 0.975),
			   digits = max(3L, getOption("digits") - 3L), ...) {
  summary.PGOcc(object, quantiles, digits)
}

fitted.tPGOcc <- function(object, ...) {
  fitted.stPGOcc(object)
}

predict.tPGOcc <- function(object, X.0, t.cols, ignore.RE = FALSE,
                           type = 'occupancy', ...) {

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

  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (length(dim(X.0)) != 3) {
    stop("error: X.0 must be an array with three dimensions corresponding to site, time, and covariate")
  }
  if (missing(t.cols)) {
    stop("error: t.cols must be specified\n")
  }

  # Occurrence predictions ------------------------------------------------
  if (tolower(type) == 'occupancy') {
    n.post <- object$n.post * object$n.chains
    X <- object$X
    J <- dim(X)[1]
    n.years.max <- dim(X.0)[2]
    p.occ <- dim(X)[3]
    beta.samples <- object$beta.samples
    ar1 <- object$ar1
    if (object$psiRE & !ignore.RE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }
    if (dim(X.0)[3] != p.occ + p.occ.re){
      stop(paste("error: the third dimension of X.0 must be ", p.occ + p.occ.re,"\n", sep = ''))
    }

    if (object$psiRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get elements in design matrix with random effects
      x.re.names <- dimnames(object$X.re)[[3]]
      indx <- which(dimnames(X.0)[[3]] %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: dimnames(X.0)[[3]] must match variable names in data$occ.covs")
      }
      X.re <- X.0[, , indx, drop = FALSE]
      X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
		     ncol = dim(X.re)[3])
      X.fix <- X.0[, , -indx, drop = FALSE]
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
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
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.fix))
      p.occ.re <- 0
    }
    J.str <- nrow(X.fix) / n.years.max
    # Get AR1 REs if ar1 == TRUE
    if (ar1) {
      eta.samples <- object$eta.samples[, t.cols, drop = FALSE]
      t.indx <- rep(1:n.years.max, each = J.str)
    }

    out <- list()
    if (ar1) {
      out$psi.0.samples <- logit.inv(t(X.fix %*% t(beta.samples) +
            				t(beta.star.sites.0.samples) + 
					t(eta.samples[, t.indx])))
    } else {
      out$psi.0.samples <- logit.inv(t(X.fix %*% t(beta.samples) +
          				t(beta.star.sites.0.samples)))
    }
    out$z.0.samples <- t(matrix(rbinom(length(out$psi.0.samples), 1, c(out$psi.0.samples)),
    		                 nrow = n.post, ncol = nrow(X.fix)))
    out$psi.0.samples <- array(t(out$psi.0.samples), dim = c(J.str, n.years.max, n.post))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
    out$z.0.samples <- array(out$z.0.samples, dim = c(J.str, n.years.max, n.post))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE & !ignore.RE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (dim(X.0)[3] != p.design) {
      stop(paste("error: the third dimension of X.0 must be ", p.design, "\n", sep = ''))
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
      indx <- which(dimnames(X.0)[[3]] %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: dimnames(X.0)[[3]] must match variable names in data$det.covs")
      }
      X.re <- X.0[, , indx, drop = FALSE]
      X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
		     ncol = dim(X.re)[3])
      X.fix <- X.0[, , -indx, drop = FALSE]
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
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
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
      alpha.star.sites.0.samples <- matrix(0, n.post, nrow(X.fix))
      p.det.re <- 0
    }
    out$p.0.samples <- t(logit.inv(t(X.fix %*% t(alpha.samples) +
          				t(alpha.star.sites.0.samples))))
    out$p.0.samples <- array(out$p.0.samples, dim = c(nrow(X.0), ncol(X.0),
						      n.post))
    out$p.0.samples <- aperm(out$p.0.samples, c(3, 1, 2))

  }
  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)
  class(out) <- "predict.tPGOcc"
  out
}

# svcPGOcc ----------------------------------------------------------------
print.svcPGOcc <- function(x, ...) {
  print.spPGOcc(x)
}

fitted.svcPGOcc <- function(object, ...) {
  fitted.PGOcc(object)
}

residuals.svcPGOcc <- function(object, n.post.samples = 100, ...) {
  residuals.PGOcc(object, n.post.samples)  
}

summary.svcPGOcc <- function(object,
			    quantiles = c(0.025, 0.5, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {
  summary.spPGOcc(object, quantiles, digits)
}

predict.svcPGOcc <- function(object, X.0, coords.0, weights.0, n.omp.threads = 1,
			     verbose = TRUE, n.report = 100,
			     ignore.RE = FALSE, type = 'occupancy', 
			     grid.index.0, ...) {

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
  if (!(class(object) %in% c('svcPGOcc', 'svcPGBinom'))) {
    stop("error: requires an output object of class svcPGOcc or svcPGBinom\n")
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
  if (is(object, 'svcPGBinom')) {
    if (missing(weights.0)) {
      message('weights.0 not specified. Assuming weights = 1 for all prediction sites.')
      weights.0 <- rep(1, nrow(X.0))
    }
  } else {
    weights.0 <- rep(1, nrow(X.0))
  }
  if (missing(grid.index.0)) {
    grid.index.0 <- 1:nrow(X.0)
  }
  grid.index.0.c <- grid.index.0 - 1

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
    svc.cols <- object$svc.cols
    p.svc <- length(svc.cols)
    X.w.0 <- X.0[, svc.cols, drop = FALSE]
    X.w <- object$X.w
    coords <- object$coords
    J <- nrow(X)
    J.w.0 <- nrow(coords.0)
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
    # coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
    # X.0.new <- X.0[coords.0.indx, , drop = FALSE]
    # X.w.0.new <- X.w.0[coords.0.indx, , drop = FALSE]
    # weights.0.new <- weights.0[coords.0.indx]

    # if (length(coords.indx) == nrow(X.0)) {
    #   stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
    # }

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

    # Sub-sample previous
    theta.samples <- t(theta.samples)
    beta.samples <- t(beta.samples)
    # Order: site, svc within site, iteration within svc. 
    w.samples <- matrix(w.samples, n.post, J * p.svc)
    # Order: iteration, site within iteration, svc within site. 
    # Example: site 1, svc 1, iter 1, site 1, svc 2, iter 1, ..., site 2, svc 1, iter 1
    w.samples <- t(w.samples)
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

    J.str <- nrow(X.0)

    # Currently predict is only implemented for NNGP.
    # Get nearest neighbors
    # nn2 is a function from RANN.
    nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

    # Indicates whether a site has been sampled. 1 = sampled
    sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
    sites.link <- rep(NA, J.w.0)
    sites.link[which(!is.na(match.indx))] <- coords.indx
    # For C
    sites.link <- sites.link - 1
    # Number of spatial REs for model fit
    J.w <- nrow(coords)

    storage.mode(coords) <- "double"
    storage.mode(J) <- "integer"
    storage.mode(J.w) <- 'integer'
    storage.mode(J.w.0) <- 'integer'
    storage.mode(p.occ) <- "integer"
    storage.mode(p.svc) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(X.fix) <- "double"
    storage.mode(X.w.0) <- "double"
    storage.mode(coords.0) <- "double"
    storage.mode(weights.0) <- "double"
    storage.mode(grid.index.0.c) <- 'integer'
    storage.mode(J.str) <- "integer"
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
    storage.mode(sites.link) <- 'integer'
    storage.mode(sites.0.sampled) <- 'integer'

    ptm <- proc.time()

    out <- .Call("svcPGOccNNGPPredict", coords, J, p.occ, p.svc, n.neighbors,
                 X.fix, X.w.0, coords.0, weights.0, J.str, nn.indx.0, beta.samples,
                 theta.samples, w.samples, beta.star.sites.0.samples, n.post,
                 cov.model.indx, n.omp.threads, verbose, n.report, J.w.0, J.w, 
                 grid.index.0.c, sites.link, sites.0.sampled)

    if (is(object, 'svcPGOcc')) {
        out$z.0.samples <- mcmc(t(out$z.0.samples))
        out$psi.0.samples <- mcmc(t(out$psi.0.samples))
        out$w.0.samples <- array(out$w.0.samples, dim = c(p.svc, J.w.0, n.post))
        out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    }
    if (is(object, 'svcPGBinom')) {
        out$y.0.samples <- mcmc(t(out$z.0.samples))
        out$z.0.samples <- NULL
        out$psi.0.samples <- mcmc(t(out$psi.0.samples))
        out$w.0.samples <- array(out$w.0.samples, dim = c(p.svc, J.w.0, n.post))
        out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
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
  class(out) <- "predict.svcPGOcc"
  out
}

# svcPGBinom --------------------------------------------------------------
print.svcPGBinom <- function(x, ...) {
  print.spPGOcc(x)
}

fitted.svcPGBinom <- function(object, ...) {
  return(object$y.rep.samples)
}

summary.svcPGBinom <- function(object,
			       quantiles = c(0.025, 0.5, 0.975),
			       digits = max(3L, getOption("digits") - 3L), ...) {
  summary.spPGOcc(object, quantiles, digits)
}

predict.svcPGBinom <- function(object, X.0, coords.0, weights.0, n.omp.threads = 1,
			       verbose = TRUE, n.report = 100,
			       ignore.RE = FALSE, ...) {
  predict.svcPGOcc(object, X.0, coords.0, weights.0, n.omp.threads = n.omp.threads,
		   verbose, n.report, ignore.RE, type = 'occupancy')

}

# svcTPGOcc ---------------------------------------------------------------
print.svcTPGOcc <- function(x, ...) {
  print.spPGOcc(x)
}

summary.svcTPGOcc <- function(object,
			    quantiles = c(0.025, 0.5, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {
  summary.spPGOcc(object, quantiles, digits)
}

fitted.svcTPGOcc <- function(object, ...) {
  fitted.stPGOcc(object)
}

predict.svcTPGOcc <- function(object, X.0, coords.0, t.cols, weights.0, n.omp.threads = 1,
                              verbose = TRUE, n.report = 100,
                              ignore.RE = FALSE, type = 'occupancy', 
                              forecast = FALSE, grid.index.0, ...) {
  
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
  
  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }
  
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (length(dim(X.0)) != 3) {
    stop("error: X.0 must be an array with three dimensions corresponding to site, time, and covariate")
  }
  
  if (missing(t.cols) & forecast == FALSE) {
    stop("error: t.cols must be specified when forecast = FALSE")
  }
  
  if (is(object, 'svcTPGBinom')) {
    if (missing(weights.0)) {
      message('weights.0 not specified. Assuming weights = 1 for all prediction sites/times.')
      weights.0 <- matrix(1, nrow(X.0), ncol(X.0))
    }
  } else {
    weights.0 <- matrix(1, nrow(X.0), ncol(X.0))
  }
  if (missing(grid.index.0)) {
    grid.index.0 <- 1:nrow(X.0)
  }
  grid.index.0.c <- grid.index.0 - 1
  
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
    svc.cols <- object$svc.cols
    p.svc <- length(svc.cols)
    X.w.0 <- X.0[, , svc.cols, drop = FALSE]
    X.w <- object$X.w
    coords <- object$coords
    J <- nrow(X)
    J.w <- nrow(object$coords)
    J.w.0 <- nrow(coords.0)
    n.years.max <- dim(X.0)[2]
    p.occ <- dim(X)[3]
    theta.samples <- object$theta.samples
    beta.samples <- object$beta.samples
    w.samples <- object$w.samples
    n.neighbors <- object$n.neighbors
    cov.model.indx <- object$cov.model.indx
    sp.type <- object$type
    # Get AR1 random effect values for the corresponding years. 
    ar1 <- object$ar1
    if (forecast & ar1) {
      message("NOTE: forecasting in spOccupancy does not currently use AR(1) random effects\n")
    }
    if (ar1 & !forecast) {
      eta.samples <- object$eta.samples[, t.cols, drop = FALSE]
    } else {
      eta.samples <- matrix(0, n.post, n.years.max)
    }
    if (object$psiRE & !ignore.RE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }
    if (dim(X.0)[3] != p.occ + p.occ.re){
      stop(paste("error: the third dimension of X.0 must be ", p.occ + p.occ.re,"\n", sep = ''))
    }
    # Eliminate prediction sites that have already sampled been for now
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))
    
    if (object$psiRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get elements in design matrix with random effects
      x.re.names <- dimnames(object$X.re)[[3]]
      indx <- which(dimnames(X.0)[[3]] %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: dimnames(X.0)[[3]] must match variable names in data$occ.covs")
      }
      X.re <- X.0[, , indx, drop = FALSE]
      X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
                     ncol = dim(X.re)[3])
      X.fix <- X.0[, , -indx, drop = FALSE]
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
                      ncol = dim(X.fix)[3])
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
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
                      ncol = dim(X.fix)[3])
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.fix))
      p.occ.re <- 0
    }
    
    X.w.0 <- matrix(X.w.0, nrow = nrow(X.w.0) * ncol(X.w.0),
                    ncol = dim(X.w.0)[3])
    
    # Sub-sample previous
    if (ar1) {
      remove.indx <- (ncol(theta.samples) - 1):ncol(theta.samples)
      theta.samples <- t(theta.samples[, -c(remove.indx)])
    } else {
      theta.samples <- t(theta.samples)
    }
    beta.samples <- t(beta.samples)
    # Order: site, svc within site, iteration within svc.
    w.samples <- matrix(w.samples, n.post, J.w * p.svc)
    # Order: iteration, site within iteration, svc within site.
    # Example: site 1, svc 1, iter 1, site 1, svc 2, iter 1, ..., site 2, svc 1, iter 1
    w.samples <- t(w.samples)
    eta.samples <- t(eta.samples)
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)
    
    J.str <- nrow(X.fix) / n.years.max
    sites.0.indx <- 0:(nrow(X.0) - 1)
    J.0 <- length(unique(sites.0.indx))
    sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
    sites.link <- rep(NA, J.0)
    sites.link[which(!is.na(match.indx))] <- coords.indx
    # For C
    sites.link <- sites.link - 1
    J.w <- nrow(coords)
    
    # Check if sampled sites are included and make sure predicting across
    # all years.
    if ((sum(sites.0.sampled) > 0) & n.years.max != dim(object$X)[2]) {
      stop("error: when predicting at sampled sites using svcTPGOcc, you must predict across all primary time periods")
    }
    
    # Currently predict is only implemented for NNGP.
    # Get nearest neighbors
    # nn2 is a function from RANN.
    nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1
    
    storage.mode(coords) <- "double"
    storage.mode(J) <- "integer"
    storage.mode(J.w) <- 'integer'
    storage.mode(J.w.0) <- 'integer'
    storage.mode(n.years.max) <- "integer"
    storage.mode(p.occ) <- "integer"
    storage.mode(p.svc) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(X.fix) <- "double"
    storage.mode(X.w.0) <- "double"
    storage.mode(sites.link) <- "integer"
    storage.mode(sites.0.sampled) <- 'integer'
    storage.mode(grid.index.0.c) <- 'integer'
    storage.mode(coords.0) <- "double"
    storage.mode(weights.0) <- "double"
    storage.mode(J.str) <- "integer"
    storage.mode(beta.samples) <- "double"
    storage.mode(theta.samples) <- "double"
    storage.mode(w.samples) <- "double"
    storage.mode(beta.star.sites.0.samples) <- "double"
    storage.mode(eta.samples) <- "double"
    storage.mode(n.post) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(nn.indx.0) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    ptm <- proc.time()
    
    out <- .Call("svcTPGOccNNGPPredict", coords, J, n.years.max, p.occ, p.svc, n.neighbors,
                 X.fix, X.w.0, coords.0, weights.0, J.str, nn.indx.0, beta.samples,
                 theta.samples, w.samples, beta.star.sites.0.samples, eta.samples, 
                 sites.link, sites.0.sampled, n.post,
                 cov.model.indx, n.omp.threads, verbose, n.report, J.w.0, J.w, grid.index.0.c)
    
    if (class(object) %in% c('svcTPGOcc', 'svcTIntPGOcc')) {
      out$z.0.samples <- array(out$z.0.samples, dim = c(J.str, n.years.max, n.post))
      out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
      out$psi.0.samples <- array(out$psi.0.samples, dim = c(J.str, n.years.max, n.post))
      out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
      out$w.0.samples <- array(out$w.0.samples, dim = c(p.svc, J.w.0, n.post))
      out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    }
    if (is(object, 'svcTPGBinom')) {
      out$z.0.samples <- array(out$z.0.samples, dim = c(J.str, n.years.max, n.post))
      out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
      out$y.0.samples <- out$z.0.samples
      out$z.0.samples <- NULL
      out$psi.0.samples <- array(out$psi.0.samples, dim = c(J.str, n.years.max, n.post))
      out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))
      out$w.0.samples <- array(out$w.0.samples, dim = c(p.svc, J.w.0, n.post))
      out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
    }
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.tPGOcc(object, X.0, t.cols, ignore.RE, type)
  }
  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)
  class(out) <- "predict.svcTPGOcc"
  out
}

# svcTPGBinom -------------------------------------------------------------
print.svcTPGBinom <- function(x, ...) {
  print.spPGOcc(x)
}

summary.svcTPGBinom <- function(object,
			    quantiles = c(0.025, 0.5, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {
  summary.spPGOcc(object, quantiles, digits)
}

fitted.svcTPGBinom <- function(object, ...) {
  return(object$y.rep.samples)
}

predict.svcTPGBinom <- function(object, X.0, coords.0, t.cols, weights.0, n.omp.threads = 1,
			        verbose = TRUE, n.report = 100,
			        ignore.RE = FALSE, ...) {
  predict.svcTPGOcc(object, X.0, coords.0, t.cols, weights.0, n.omp.threads = n.omp.threads, 
		   verbose, n.report, ignore.RE, type = 'occupancy')

}

# intMsPGOcc --------------------------------------------------------------
print.intMsPGOcc <- function(x, ...) {
  print.msPGOcc(x)
}

summary.intMsPGOcc <- function(object,
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

  n.data <- length(object$y)
  p.det.long <- sapply(object$X.p, function(a) dim(a)[[2]])

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
    for (i in 1:n.data) {
      indx <- object$alpha.comm.indx == i
      cat("-----------------------------\n");
      cat(paste("\tData source ", i, "\n", sep = ''));
      cat("-----------------------------\n");
      cat("Detection Means (logit scale): \n")
      tmp.1 <- t(apply(object$alpha.comm.samples[, indx, drop = FALSE],
		       2, function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$alpha.comm.samples[, indx, drop = FALSE], 
		     2, function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$alpha.comm[indx], 
			round(object$ESS$alpha.comm[indx], 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')
      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
      cat("\n")
      cat("Detection Variances (logit scale): \n")
      tmp.1 <- t(apply(object$tau.sq.alpha.samples[, indx, drop = FALSE], 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$tau.sq.alpha.samples[, indx, drop = FALSE], 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$tau.sq.alpha[indx], round(object$ESS$tau.sq.alpha[indx], 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
      cat("\n")
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
    for (i in 1:n.data) {
      indx <- object$alpha.indx == i
      cat("-----------------------------\n");
      cat(paste("\tData source ", i, "\n", sep = ''));
      cat("-----------------------------\n");
      cat("Detection (logit scale): \n")
      tmp.1 <- t(apply(object$alpha.samples[, indx, drop = FALSE], 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$alpha.samples[, indx, drop = FALSE], 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$alpha[indx], round(object$ESS$alpha[indx], 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')
      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
      cat("\n")
    }
  }
}

predict.intMsPGOcc <- function(object, X.0, ignore.RE = FALSE, ...) {
  predict.msPGOcc(object, X.0, ignore.RE, type = 'occupancy', ...)
}

# postHocLM ---------------------------------------------------------------
print.postHocLM <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}
summary.postHocLM <- function(object,
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.samples <- object$n.samples
  n.chains <- object$n.chains
  n.post <- object$n.samples * object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  # Fixed Effects
  cat("----------------------------------------\n");
  cat("Fixed Effects: \n")
  cat("----------------------------------------\n");
  tmp.1 <- t(apply(object$beta.samples, 2,
		   function(x) c(mean(x), sd(x))))
  tmp.2 <- t(apply(object$tau.sq.samples, 2, function(x) c(mean(x), sd(x))))
  rownames(tmp.2) <- "Residual Variance"
  tmp.1 <- rbind(tmp.1, tmp.2)
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  tmp.2 <- t(apply(object$tau.sq.samples, 2,
		   function(x) quantile(x, prob = quantiles)))
  rownames(tmp.2) <- "Residual Variance"
  tmp <- rbind(tmp, tmp.2)
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  diags.2 <- matrix(c(object$rhat$tau.sq, round(object$ESS$tau.sq, 0)), ncol = 2)
  diags <- rbind(diags, diags.2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  if (object$RE) {
    cat("\n")
    cat("----------------------------------------\n");
    cat("Random Effects: \n")
    cat("----------------------------------------\n");
    tmp.1 <- t(apply(object$sigma.sq.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq, round(object$ESS$sigma.sq, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")

  cat("----------------------------------------\n");
  cat("Bayesian R2: \n")
  cat("----------------------------------------\n");
  tmp.1 <- t(apply(object$bayes.R2, 2,
        	   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  rownames(tmp.1) <- ""
  tmp <- t(apply(object$bayes.R2, 2,
        	 function(x) quantile(x, prob = quantiles)))
  rownames(tmp) <- ""
  print(noquote(round(cbind(tmp.1, tmp), digits)))
}

# svcMsPGOcc --------------------------------------------------------------
print.svcMsPGOcc <- function(x, ...) {
  print.sfMsPGOcc(x)
}

summary.svcMsPGOcc <- function(object,
			      level = 'both',
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {
  summary.sfMsPGOcc(object, level, quantiles, digits)
}

fitted.svcMsPGOcc <- function(object, ...) {
  fitted.msPGOcc(object)
}

predict.svcMsPGOcc <- function(object, X.0, coords.0, n.omp.threads = 1,
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
  if (!(class(object) %in% c('svcMsPGOcc'))) {
    stop("error: requires an output object of class svcMsPGOcc\n")
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
    if (missing(coords.0)) {
      stop("error: coords.0 must be specified\n")
    }
    if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
      stop("error: coords.0 must be a data.frame or matrix\n")
    }
    if (!ncol(coords.0) == 2){
      stop("error: coords.0 must have two columns\n")
    }
    n.post <- object$n.post * object$n.chains
    X <- object$X
    y <- object$y
    coords <- object$coords
    J <- nrow(X)
    N <- dim(y)[1]
    q <- object$q
    p.occ <- ncol(X)
    std.by.sp <- object$std.by.sp
    species.sds <- object$species.sds
    species.means <- object$species.means
    svc.cols <- object$svc.cols
    p.svc <- length(svc.cols)
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
    X.w.0 <- X.0[, svc.cols, drop = FALSE]
    X.w <- object$X.w

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

    if (object$psiRE) {
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
      X.fix <- X.0
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.0))
      p.occ.re <- 0
    }

    # Sub-sample previous
    theta.samples <- t(theta.samples)
    lambda.samples <- t(lambda.samples)
    beta.samples <- t(beta.samples)
    # Desired ordering: iteration, svc, site, factor
    w.samples <- aperm(w.samples, c(2, 3, 4, 1))
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

    J.str <- nrow(X.0)

    X.big <- array(NA, dim = c(J.str, ncol(X.fix), N))
    for (i in 1:N) {
      X.big[, , i] <- X.fix
      if (std.by.sp) {
        for (r in 1:ncol(X.fix)) {
          if (!is.na(species.sds[i, r])) {
            X.big[, r, i] <- (X.big[, r, i] - species.means[i, r]) / species.sds[i, r]
          }
        }
      }
    }
    X.w.big <- X.big[, svc.cols,  , drop = FALSE]

    if (sp.type == 'GP') {
      # Not currently implemented or accessed.
    } else {
      # Get nearest neighbors
      # nn2 is a function from RANN.
      nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1
      # Indicates whether a site has been sampled. 1 = sampled
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, J.str)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1

      storage.mode(coords) <- "double"
      storage.mode(N) <- "integer"
      storage.mode(J) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(p.svc) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.big) <- "double"
      storage.mode(X.w.big) <- "double"
      storage.mode(coords.0) <- "double"
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
      storage.mode(sites.link) <- 'integer'
      storage.mode(sites.0.sampled) <- 'integer'

      out <- .Call("svcMsPGOccNNGPPredict", coords, J, N, q, p.occ, p.svc, n.neighbors,
                   X.big, X.w.big, coords.0, J.str, nn.indx.0, beta.samples,
                   theta.samples, lambda.samples, w.samples,
          	   beta.star.sites.0.samples, n.post,
                   cov.model.indx, n.omp.threads, verbose, n.report, 
          	       sites.link, sites.0.sampled)

    }
    out$z.0.samples <- array(out$z.0.samples, dim = c(N, J.str, n.post))
    out$z.0.samples <- aperm(out$z.0.samples, c(3, 1, 2))
    out$w.0.samples <- array(out$w.0.samples, dim = c(q, J.str, p.svc, n.post))
    out$w.0.samples <- aperm(out$w.0.samples, c(4, 1, 2, 3))
    out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, J.str, n.post))
    out$psi.0.samples <- aperm(out$psi.0.samples, c(3, 1, 2))

  } # occurrence predictions
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.msPGOcc(object, X.0, ignore.RE, type)
  }

  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)

  class(out) <- "predict.svcMsPGOcc"

  out

}

# svcTMsPGOcc -------------------------------------------------------------
print.svcTMsPGOcc <- function(x, ...) {
  print.sfMsPGOcc(x)
}

summary.svcTMsPGOcc <- function(object,
			      level = 'both',
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {
  summary.sfMsPGOcc(object, level, quantiles, digits)
}

predict.svcTMsPGOcc <- function(object, X.0, coords.0,
				t.cols, n.omp.threads = 1,
			        verbose = TRUE, n.report = 100,
			        ignore.RE = FALSE, type = 'occupancy', 
				grid.index.0, ...) {

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
  if (!(class(object) %in% c('svcTMsPGOcc', 'stMsPGOcc'))) {
    stop("error: requires an output object of class svcTMsPGOcc, stMsPGOcc\n")
  }
  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (length(dim(X.0)) != 3) {
    stop("error: X.0 must be an array with three dimensions corresponding to site, time, and covariate.")
  }

  if (missing(t.cols)) {
    stop("error: t.cols must be specified\n")
  }

  if (missing(grid.index.0)) {
    grid.index.0 <- 1:nrow(X.0)
  }
  grid.index.0.c <- grid.index.0 - 1

  ptm <- proc.time()

  # Occurrence predictions ------------------------------------------------
  if (tolower(type == 'occupancy')) {
    if (missing(coords.0)) {
      stop("error: coords.0 must be specified\n")
    }
    if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
      stop("error: coords.0 must be a data.frame or matrix\n")
    }
    if (!ncol(coords.0) == 2){
      stop("error: coords.0 must have two columns\n")
    }
    n.post <- object$n.post * object$n.chains
    X <- object$X
    y <- object$y
    coords <- object$coords
    J <- nrow(X)
    J.w.0 <- nrow(coords.0)
    n.years.max <- dim(X.0)[2]
    N <- dim(y)[1]
    q <- object$q
    p.occ <- dim(X)[[3]]
    std.by.sp <- object$std.by.sp
    species.sds <- object$species.sds
    species.means <- object$species.means
    svc.cols <- object$svc.cols
    p.svc <- length(svc.cols)
    theta.samples <- object$theta.samples
    beta.samples <- object$beta.samples
    lambda.samples <- object$lambda.samples
    w.samples <- object$w.samples
    n.neighbors <- object$n.neighbors
    cov.model.indx <- object$cov.model.indx
    sp.type <- object$type
    # Get AR1 random effect values for the corresponding years.
    ar1 <- object$ar1
    if (ar1) {
      eta.samples <- object$eta.samples[, , t.cols, drop = FALSE]
    } else {
      eta.samples <- array(0, dim = c(n.post, N, n.years.max))
    }
    if (object$psiRE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }

    X.w.0 <- X.0[, , svc.cols, drop = FALSE]
    X.w <- object$X.w

    coords.0 <- as.matrix(coords.0)

    # Eliminate prediction sites that have already been sampled for now
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))

    if (object$psiRE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- dimnames(object$X.re)[[3]]
      x.names <- dimnames(object$X)[[3]]
      indx <- which(dimnames(X.0)[[3]] %in% x.re.names)
      X.re <- X.0[, , indx, drop = FALSE]
      X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
      	     ncol = dim(X.re)[3])
      remove.fix.indx <- indx[which(!(dimnames(X.0)[[3]][indx] %in% x.names))]
      if (length(remove.fix.indx)) {
        X.fix <- X.0[, , -remove.fix.indx, drop = FALSE]
      } else {
        X.fix <- X.0
      }
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
  		      ncol = dim(X.fix)[3])
      n.occ.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.occ.re)
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
      X.fix <- X.0
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.fix))
      p.occ.re <- 0
    }

    if (ar1) {
      if (object$cov.model.indx == 2) { # matern == 2
        theta.samples <- t(theta.samples[, 1:(q * p.svc * 2)])
      } else {
        theta.samples <- t(theta.samples[, 1:(q * p.svc)])
      }
    } else {
      theta.samples <- t(theta.samples)
    }
    beta.samples <- t(beta.samples)
    lambda.samples <- t(lambda.samples)
    # Desired ordering: iteration, svc, site, factor
    w.samples <- aperm(w.samples, c(2, 3, 4, 1))
    eta.samples <- aperm(eta.samples, c(2, 3, 1))
    beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

    J.str <- nrow(X.fix) / n.years.max

    X.big <- array(NA, dim = c(J.str, n.years.max, ncol(X.fix), N))
    for (i in 1:N) {
      X.big[, , , i] <- array(X.fix, dim = c(J.str, n.years.max, ncol(X.fix)))[, , , drop = FALSE]
      if (std.by.sp) {
        for (r in 1:ncol(X.fix)) {
          if (!is.na(species.sds[i, r])) {
            X.big[, , r, i] <- (X.big[, , r, i] - species.means[i, r]) / species.sds[i, r]
          }
        }
      }
    }
    X.big <- ifelse(is.na(X.big), 0, X.big)
    X.w.big <- X.big[, , svc.cols,  , drop = FALSE]

    # Get stuff for linking sampled sites to predicted sites
    sites.0.indx <- 0:(nrow(X.0) - 1)
    J.0 <- length(unique(sites.0.indx))
    sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
    sites.link <- rep(NA, J.0)
    sites.link[which(!is.na(match.indx))] <- coords.indx
    # For C
    sites.link <- sites.link - 1
    J.w <- nrow(coords)

    # Check if sampled sites are included and make sure predicting across
    # all years.
    if ((sum(sites.0.sampled) > 0) & n.years.max != dim(object$X)[2]) {
      stop("error: when predicting at sampled sites using svcTPGOcc, you must predict across all primary time periods")
    }

    if (sp.type == 'GP') {
      # Not currently implemented or accessed.
    } else {
      # Get nearest neighbors
      # nn2 is a function from RANN.
      nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

      storage.mode(coords) <- "double"
      storage.mode(N) <- "integer"
      storage.mode(J) <- "integer"
      storage.mode(J.w) <- 'integer'
      storage.mode(J.w.0) <- 'integer'
      storage.mode(n.years.max) <- "integer"
      storage.mode(p.occ) <- "integer"
      storage.mode(p.svc) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.big) <- "double"
      storage.mode(X.w.big) <- "double"
      storage.mode(coords.0) <- "double"
      storage.mode(sites.link) <- "integer"
      storage.mode(sites.0.sampled) <- 'integer'
      storage.mode(grid.index.0.c) <- 'integer'
      storage.mode(J.str) <- "integer"
      storage.mode(q) <- "integer"
      storage.mode(beta.samples) <- "double"
      storage.mode(theta.samples) <- "double"
      storage.mode(lambda.samples) <- "double"
      storage.mode(eta.samples) <- "double"
      storage.mode(beta.star.sites.0.samples) <- "double"
      storage.mode(w.samples) <- "double"
      storage.mode(n.post) <- "integer"
      storage.mode(cov.model.indx) <- "integer"
      storage.mode(nn.indx.0) <- "integer"
      storage.mode(n.omp.threads) <- "integer"
      storage.mode(verbose) <- "integer"
      storage.mode(n.report) <- "integer"

      out <- .Call("svcTMsPGOccNNGPPredict", coords, J, n.years.max, N, q, p.occ, p.svc, n.neighbors,
                   X.big, X.w.big, coords.0, J.str, nn.indx.0, beta.samples,
                   theta.samples, lambda.samples, w.samples,
          	   beta.star.sites.0.samples, eta.samples, 
		   sites.link, sites.0.sampled, n.post,
                   cov.model.indx, n.omp.threads, verbose, n.report, J.w.0, J.w, 
                   grid.index.0.c)

    }

      out$z.0.samples <- array(out$z.0.samples, dim = c(N, J.str, n.years.max, n.post))
      out$z.0.samples <- aperm(out$z.0.samples, c(4, 1, 2, 3))
      out$w.0.samples <- array(out$w.0.samples, dim = c(q, J.w.0, p.svc, n.post))
      out$w.0.samples <- aperm(out$w.0.samples, c(4, 1, 2, 3))
      out$psi.0.samples <- array(out$psi.0.samples, dim = c(N, J.str, n.years.max, n.post))
      out$psi.0.samples <- aperm(out$psi.0.samples, c(4, 1, 2, 3))
  } # occurrence predictions
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.tMsPGOcc(object, X.0, t.cols, ignore.RE = ignore.RE, type = type)
  }

  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)

  class(out) <- "predict.svcTMsPGOcc"

  out
}

fitted.svcTMsPGOcc <- function(object, ...) {
   fitted.tMsPGOcc(object)
}

# tMsPGOcc ----------------------------------------------------------------
print.tMsPGOcc <- function(x, ...) {
  print.msPGOcc(x)
}

summary.tMsPGOcc <- function(object, level = 'both', quantiles = c(0.025, 0.5, 0.975),
			     digits = max(3L, getOption("digits") - 3L), ...) {
  summary.msPGOcc(object, level, quantiles, digits)
}

fitted.tMsPGOcc <- function(object, ...) {
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
  if (!(class(object) %in% c('tMsPGOcc', 'stMsPGOcc', 'svcTMsPGOcc'))) {
    stop("error: object must be of class tMsPGOcc, stMsPGOcc, svcTMsPGOcc\n")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.years.max <- dim(y)[3]
  K.max <- dim(y)[4]
  J <- dim(y)[2]
  N <- dim(y)[1]
  z.long.indx <- rep(1:(J * n.years.max), K.max)
  z.long.indx <- z.long.indx[!is.na(c(y[1, , , ]))]
  z.samples <- object$z.samples
  alpha.samples <- object$alpha.samples
  n.obs <- nrow(X.p)
  det.prob.samples <- array(NA, dim = c(n.obs, N, n.post))
  sp.indx <- rep(1:N, ncol(X.p))
  y <- matrix(y, N, J * n.years.max * K.max)
  y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
  for (i in 1:N) {
    if (object$pRE) {
      sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
      # Add 1 to get it to R indexing.
      X.p.re <- object$X.p.re + 1
      tmp.samples <- matrix(0, n.post, n.obs)
      tmp.alpha.star <- object$alpha.star.samples[, sp.re.indx == i]
      for (j in 1:ncol(X.p.re)) {
        tmp.samples <- tmp.samples + tmp.alpha.star[, X.p.re[, j]]
      }
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + t(tmp.samples))
    } else {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }
  }

  out <- list()
  # Get detection probability
  # Need to be careful here that all arrays line up.
  det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
  tmp <- array(NA, dim = c(n.post, N, J * n.years.max * K.max))
  names.long <- which(!is.na(c(object$y[1, , , ])))
  tmp[, , names.long] <- det.prob.samples
  p.samples <- array(tmp, dim = c(n.post, N, J, n.years.max, K.max))
  out$p.samples <- p.samples
  # Get fitted values
  z.samples <- array(z.samples, dim = c(n.post, N, J * n.years.max))
  det.prob.samples <- det.prob.samples * z.samples[, , z.long.indx]
  y.rep.samples <- array(NA, dim = dim(det.prob.samples))
  for (i in 1:N) {
    y.rep.samples[, i, ] <- apply(det.prob.samples[, i, ], 2, function(a) rbinom(n.post, 1, a))
  }
  tmp <- array(NA, dim = c(n.post, N, J * n.years.max * K.max))
  names.long <- which(!is.na(c(object$y[1, , , ])))
  tmp[, , names.long] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(n.post, N, J, n.years.max, K.max))
  out$y.rep.samples <- y.rep.samples
  return(out)
}

predict.tMsPGOcc <- function(object, X.0, t.cols, ignore.RE = FALSE, 
			     type = 'occupancy', ...) {

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

  if (!(tolower(type) %in% c('occupancy', 'detection'))) {
    stop("error: prediction type must be either 'occupancy' or 'detection'")
  }

  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (length(dim(X.0)) != 3) {
    stop("error: X.0 must be an array with three dimensions corresponding to site, time, and covariate")
  }
  if (missing(t.cols)) {
    stop("error: t.cols must be specified\n")
  }

  # Occurrence predictions ------------------------------------------------
  if (tolower(type) == 'occupancy') {
    n.post <- object$n.post * object$n.chains
    X <- object$X
    J <- dim(X)[1]
    n.years.max <- dim(X.0)[2]
    p.occ <- dim(X)[3]
    beta.samples <- object$beta.samples
    ar1 <- object$ar1
    if (object$psiRE & !ignore.RE) {
      p.occ.re <- length(object$re.level.names)
    } else {
      p.occ.re <- 0
    }
    if (dim(X.0)[3] != p.occ + p.occ.re){
      stop(paste("error: the third dimension of X.0 must be ", p.occ + p.occ.re,"\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    sp.indx <- rep(1:N, p.occ)
    n.post <- object$n.post * object$n.chains
    beta.samples <- as.matrix(object$beta.samples)
    out <- list()
    out$psi.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0), n.years.max))
    out$z.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0), n.years.max))
    if (object$psiRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- dimnames(object$X.re)[[3]]
      indx <- which(dimnames(X.0)[[3]] %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: dimnames(X.0)[[3]] must match variable names in data$occ.covs")
      }
      X.re <- X.0[, , indx, drop = FALSE]
      X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
		     ncol = dim(X.re)[3])
      X.fix <- X.0[, , -indx, drop = FALSE]
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
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
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
      beta.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.fix))
      p.occ.re <- 0
    }
    J.str <- nrow(X.fix) / n.years.max
    site.indx <- rep(1:J.str, times = n.years.max)
    t.indx <- rep(1:n.years.max, each = J.str)
    # Get AR1 REs if ar1 == TRUE
    if (ar1) {
      eta.samples <- object$eta.samples[, , t.cols, drop = FALSE]
      # Make predictions
      for (i in 1:N) {
        for (j in 1:(J.str * n.years.max)) {
          out$psi.0.samples[, i, site.indx[j], t.indx[j]] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
            				     t(beta.samples[, sp.indx == i]) + 
                                                 beta.star.sites.0.samples[, (j - 1) * N + i] + 
	                                         c(eta.samples[, i, t.indx[j]]))
          out$z.0.samples[, i, site.indx[j], t.indx[j]] <- rbinom(n.post, 1, 
					      out$psi.0.samples[, i, site.indx[j], t.indx[j]])
        } # j
      } # i
    } else {
      # Make predictions
      for (i in 1:N) {
        for (j in 1:(J.str * n.years.max)) {
          out$psi.0.samples[, i, site.indx[j], t.indx[j]] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
            				     t(beta.samples[, sp.indx == i]) + 
                                                 beta.star.sites.0.samples[, (j - 1) * N + i])
          out$z.0.samples[, i, site.indx[j], t.indx[j]] <- rbinom(n.post, 1, out$psi.0.samples[, i, site.indx[j], t.indx[j]])
        } # j
      } # i
    }
  } # occurrence predictions
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE & !ignore.RE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (dim(X.0)[3] != p.design) {
      stop(paste("error: the third dimension of X.0 must be ", p.design, "\n", sep = ''))
    }
    # Composition sampling --------------------------------------------------
    N <- dim(object$y)[1]
    sp.indx <- rep(1:N, p.det)
    n.post <- object$n.post * object$n.chains
    alpha.samples <- as.matrix(object$alpha.samples)
    n.time.max <- dim(X.0)[[2]]
    out <- list()
    out$p.0.samples <- array(NA, dim = c(n.post, N, nrow(X.0), n.time.max))
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
      indx <- which(dimnames(X.0)[[3]] %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: dimnames(X.0)[[3]] must match variable names in data$det.covs")
      }
      X.re <- X.0[, , indx, drop = FALSE]
      X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
		     ncol = dim(X.re)[3])
      X.fix <- X.0[, , -indx, drop = FALSE]
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
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
      X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		      ncol = dim(X.fix)[3])
      alpha.star.sites.0.samples <- matrix(0, n.post, N * nrow(X.fix))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    site.indx <- rep(1:J.str, times = n.time.max)
    t.indx <- rep(1:n.time.max, each = J.str)
    # Make predictions
    for (i in 1:N) {
      for (j in 1:(J.str * n.time.max)) {
        out$p.0.samples[, i, site.indx[j], t.indx[j]] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				     t(alpha.samples[, sp.indx == i]) + 
                                               alpha.star.sites.0.samples[, (j - 1) * N + i])
      } # j
    } # i
  }
  out$call <- cl

  class(out) <- "predict.tMsPGOcc"
  out

}

# stMsPGOcc -------------------------------------------------------------
print.stMsPGOcc <- function(x, ...) {
  print.sfMsPGOcc(x)
}

summary.stMsPGOcc <- function(object,
			      level = 'both',
			      quantiles = c(0.025, 0.5, 0.975),
			      digits = max(3L, getOption("digits") - 3L), ...) {
  summary.sfMsPGOcc(object, level, quantiles, digits)
}

fitted.stMsPGOcc <- function(object, ...) {
   fitted.tMsPGOcc(object)
}

predict.stMsPGOcc <- function(object, X.0, coords.0, 
                              t.cols, n.omp.threads = 1,
                              verbose = TRUE, n.report = 100,
                              ignore.RE = FALSE, type = 'occupancy', 
			      grid.index.0, ...) {
  object$std.by.sp <- FALSE
  object$species.sds <- NA
  object$species.means <- NA
  object$svc.cols <- 1
  tmp <- array(NA, dim = c(object$n.post * object$n.chains, 
			   object$q, nrow(object$coords), 1))
  tmp[, , , 1] <- object$w.samples
  object$w.samples <- tmp
  object$X.w <- object$X[, , 1, drop = FALSE]
  out <- predict.svcTMsPGOcc(object, X.0, coords.0, 
			     t.cols, n.omp.threads, verbose, n.report,
			     ignore.RE, type, grid.index.0)
  out$w.0.samples <- out$w.0.samples[, , , 1]
  return(out)
}

# tIntPGOcc ---------------------------------------------------------------
predict.tIntPGOcc <- function(object, X.0, t.cols, ignore.RE = FALSE, 
                              type = 'occupancy', ...) {
  # Occupancy predictions -------------------------------------------------
  if (tolower(type == 'occupancy')) {	
    out <- predict.tPGOcc(object, X.0, t.cols, ignore.RE, type)
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type == 'detection')) {
  # TODO: this should be pretty easy. Just have an argument for which data
  #       source to predict with, and then send the model to predict.tPGOcc
    stop("detection prediction is not currently implemented.")
    out <- list()
  }
  class(out) <- "predict.tIntPGOcc"
  out
}

print.tIntPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.tIntPGOcc <- function(object, ...) {
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

  y <- object$y
  n.data <- length(y)
  sites <- object$sites
  seasons <- object$seasons
  X.p <- object$X.p
  y.rep.samples <- list()
  for (i in 1:n.data) {
    y.rep.samples[[i]] <- array(NA, dim = dim(object$p.samples[[i]]))
  }
  n.post <- object$n.post * object$n.chains

  for (q in 1:n.data) {
    for (j in 1:ncol(y.rep.samples[[q]])) {
      for (t in 1:dim(y.rep.samples[[q]])[3]) {
        for (k in 1:dim(y.rep.samples[[q]])[4]) {
          if (sum(!is.na(object$p.samples[[q]][, j, t, k])) != 0) {
            y.rep.samples[[q]][, j, t, k] <- rbinom(n.post, 1, 
                                                    object$z.samples[, sites[[q]][j], 
                                                                     seasons[[q]][t]] * 
                                                    object$p.samples[[q]][, j, t, k])
          } # if none na
        } # k
      } # t
    } # j
  } # q

  out <- list()
  out$y.rep.samples <- y.rep.samples
  out$p.samples <- object$p.samples
  return(out)
}

summary.tIntPGOcc <- function(object, quantiles = c(0.025, 0.5, 0.975),
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
  p.det.re.long <- sapply(object$X.p.re, function(a) dim(a)[[2]])

  # Occurrence ------------------------
  cat("----------------------------------------\n")
  cat("Occurrence\n")
  cat("----------------------------------------\n")
  cat("Fixed Effects (logit scale):\n")
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
    cat("Random Effect Variances (logit scale):\n")
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
  indx <- 1
  indx.re <- 1
  for (i in 1:n.data) {
    cat("----------------------------------------\n")
    cat(paste("Data source ", i, " Detection\n", sep = ""))
    cat("----------------------------------------\n")
    cat("Fixed Effects (logit scale):\n")
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
    if (object$pRELong[i]) {
      tmp.samples <- object$sigma.sq.p.samples[, indx.re:(indx.re+p.det.re.long[i] - 1), 
					       drop = FALSE]
      cat("Random Effect Variances (logit scale):\n")
      tmp.1 <- t(apply(tmp.samples, 2, function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(tmp.samples, 2, function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p[indx.re:(indx.re+p.det.re.long[i] - 1)], 
			round(object$ESS$sigma.sq.p[indx.re:(indx.re+p.det.re.long[i] - 1)],
			      0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
      cat("\n")
      indx.re <- indx.re + p.det.re.long[i]
    }
  }
  if (object$ar1 | class(object) %in% c('stIntPGOcc', 'svcTIntPGOcc')) {
    if (class(object) %in% c('tIntPGOcc')) {
      cat("----------------------------------------\n")
      cat("Occurrence AR(1) Temporal Covariance: \n")
      cat("----------------------------------------\n")
    } else {
      cat("----------------------------------------\n")
      cat("Occurrence Spatio-temporal Covariance: \n")
      cat("----------------------------------------\n")
    }
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

# stIntPGOcc --------------------------------------------------------------
print.stIntPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.stIntPGOcc <- function(object, quantiles = c(0.025, 0.5, 0.975),
                              digits = max(3L, getOption("digits") - 3L), ...) {
  summary.tIntPGOcc(object, quantiles, digits)  
}

fitted.stIntPGOcc <- function(object, ...) {
  fitted.tIntPGOcc(object)  
}

predict.stIntPGOcc <- function(object, X.0, coords.0, t.cols, 
                               n.omp.threads = 1, verbose = TRUE, n.report = 100, 
                               ignore.RE = FALSE,  type = 'occupancy', forecast = FALSE, ...) {
  # Occupancy predictions -------------------------------------------------
  if (tolower(type == 'occupancy')) {	
    out <- predict.stPGOcc(object, X.0, coords.0, t.cols, n.omp.threads, 
                           verbose, n.report, ignore.RE, type, forecast)
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type == 'detection')) {
    out <- predict.tIntPGOcc(object, X.0, t.cols, ignore.RE, type)
  }
  class(out) <- "predict.stIntPGOcc"
  out
}

# svcTIntPGOcc --------------------------------------------------------------
print.svcTIntPGOcc <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.svcTIntPGOcc <- function(object, quantiles = c(0.025, 0.5, 0.975),
                              digits = max(3L, getOption("digits") - 3L), ...) {
  summary.tIntPGOcc(object, quantiles, digits)  
}

fitted.svcTIntPGOcc <- function(object, ...) {
  fitted.tIntPGOcc(object)  
}

predict.svcTIntPGOcc <- function(object, X.0, coords.0, t.cols, 
                               n.omp.threads = 1, verbose = TRUE, n.report = 100, 
                               ignore.RE = FALSE,  type = 'occupancy', forecast = FALSE, ...) {
  # Occupancy predictions -------------------------------------------------
  if (tolower(type == 'occupancy')) {	
    out <- predict.svcTPGOcc(object = object, X.0 = X.0, coords.0 = coords.0, 
                             t.cols = t.cols, n.omp.threads = n.omp.threads, 
                             verbose = verbose, n.report = n.report, 
                             ignore.RE = ignore.RE, type = type, forecast = forecast)
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type == 'detection')) {
    out <- predict.tIntPGOcc(object, X.0, t.cols, ignore.RE, type)
  }
  class(out) <- "predict.svcTIntPGOcc"
  out
}
