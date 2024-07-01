updateMCMC <- function(object, n.batch, n.samples, n.burn = 0, n.thin, 
                       keep.orig = TRUE, verbose = TRUE, n.report = 100, 
                       save.fitted = TRUE, ...) {

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
  # TODO: temporary check until the function is implemented for all 
  #.      spOccupancy and spAbundance model types
  if (!class(object) %in% c('sfJSDM', 'msAbund', 'lfJSDM')) {
    stop("updateMCMC() is currently only implemented for sfJSDM, lfJSDM, msAbund model types.")
  }
  if (missing(n.batch)) {
    if (class(object) %in% c('spPGOcc', 'spMsPGOcc', 'spIntPGOcc', 
			     'sfMsPGOcc', 'sfJSDM', 'tPGOcc', 'stPGOcc', 
			     'svcPGOcc', 'svcPGBinom', 'svcMsPGBinom', 
			     'svcTPGOcc', 'svcTPGBinom', 'svcMsPGOcc', 
			     'svcTMsPGOcc', 'msAbund')) {
      n.batch <- object$update$n.batch
    }
  }
  if (missing(n.samples)) {
    if (class(object) %in% c('PGOcc', 'msPGOcc', 'intPGOcc', 
			     'lfMsPGOcc', 'lfJSDM')) {
      n.samples <- object$n.samples
    }
  }
  if (missing(n.thin)) {
    n.thin <- object$n.thin
  }

  # Updates for each specific function ------------------------------------
  n.chains <- object$n.chains
  n.post.samples <- object$n.post * n.chains
  n.post.one.chain <- object$n.post
  out.tmp <- list()
  seeds.new <- list()
  run.time.new <- 0 
  # msAbund ---------------------------------------------------------------
  if (is(object, 'msAbund')) {
    for (i in 1:n.chains) {
      # Set the random seed based on the previous set of the model
      assign(".Random.seed", object$update$final.seed[[i]], .GlobalEnv)
      n.sp <- nrow(object$y)
      p.abund <- ncol(object$beta.comm.samples)
      cov.model <- object$update$cov.model
      # Get initial values
      curr.inits <- n.post.one.chain * i
      inits <- list()
      # beta.comm, tau.sq.beta, beta, sigma.sq.mu, beta.star, kappa
      inits$beta.comm <- object$beta.comm.samples[curr.inits, ]
      inits$tau.sq.beta <- object$tau.sq.beta.samples[curr.inits, ]
      inits$beta <- matrix(object$beta.samples[curr.inits, ], n.sp, p.abund)
      if (object$muRE) {
        inits$sigma.sq.mu <- object$sigma.sq.mu.samples[curr.inits, ]
        inits$beta.star <- t(matrix(object$beta.star.samples[curr.inits, ], 
          			  ncol = n.sp))
      }
      if (object$dist == 'NB') {
        inits$kappa <- object$kappa.samples[curr.inits, ]
      }
      # Get tuning values
      tuning <- list()
      beta.tuning.indx <- 1:ncol(object$beta.samples)
      beta.star.start <- max(beta.tuning.indx) + 1
      beta.star.tuning.indx <- beta.star.start:(beta.star.start + ncol(object$beta.star.samples) - 1)
      if (object$dist == 'NB') {
        kappa.start <- max(beta.star.tuning.indx)
        kappa.tuning.indx <- kappa.start:(kappa.start + ncol(object$kappa.samples) - 1)
      }
      tuning$beta <- object$update$tuning[beta.tuning.indx, i]
      tuning$beta.star <- object$update$tuning[beta.star.tuning.indx, i]
      if (object$dist == 'NB') {
        tuning$kappa <- object$update$tuning[kappa.tuning.indx, i]
      }
      out.tmp[[i]] <- msAbund(formula = object$update$formula,
                              data = object$update$data,
                              inits = inits,
                              priors = object$update$priors,
                              tuning = tuning,
                              n.batch = n.batch,
                              family = object$dist,
                              batch.length = object$update$batch.length,
                              accept.rate = object$update$accept.rate,
                              n.omp.threads = object$update$n.omp.threads,
                              verbose = verbose,
                              n.report = n.report,
                              n.burn = n.burn,
                              n.thin = n.thin,
                              save.fitted = save.fitted,
                              n.chains = 1) # TODO: will make output look weird.
      run.time.new <- run.time.new + out.tmp[[i]]$run.time
      seeds.new[[i]] <- out.tmp[[i]]$update$final.seed[[1]]
    }
    # Put everything together
    beta.samples.new <- list()
    beta.comm.samples.new <- list()
    tau.sq.beta.samples.new <- list()
    if (object$dist == 'NB') {
      kappa.samples.new <- list()
    }
    sigma.sq.mu.samples.new <- list()
    beta.star.samples.new <- list()
    mu.samples.new <- list()
    y.rep.samples.new <- list()
    like.samples.new <- list()

    rhat.new <- list()
    ess.new <- list()
    n.samples.one.chain <- object$n.post
    for (i in 1:n.chains) {
      if (keep.orig) {
        beta.comm.samples.new[[i]] <- rbind(object$beta.comm.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.comm.samples)
        beta.samples.new[[i]] <- rbind(object$beta.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.samples)
        tau.sq.beta.samples.new[[i]] <- rbind(object$tau.sq.beta.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$tau.sq.beta.samples)
        if (object$dist == 'NB') {
          kappa.samples.new[[i]] <- rbind(object$kappa.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$kappa.samples)
        }
        if (object$muRE) {
          sigma.sq.mu.samples.new[[i]] <- rbind(object$sigma.sq.mu.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$sigma.sq.mu.samples)
          beta.star.samples.new[[i]] <- rbind(object$beta.star.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.star.samples)
        }
        mu.samples.new[[i]] <- abind(object$mu.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , , drop = FALSE], out.tmp[[i]]$mu.samples, along = 1)
        like.samples.new[[i]] <- abind(object$like.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , , drop = FALSE], out.tmp[[i]]$like.samples, along = 1)
        y.rep.samples.new[[i]] <- abind(object$y.rep.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , , drop = FALSE], out.tmp[[i]]$y.rep.samples, along = 1)
      } else {
        beta.comm.samples.new[[i]] <- out.tmp[[i]]$beta.comm.samples
        beta.samples.new[[i]] <- out.tmp[[i]]$beta.samples
        tau.sq.beta.samples.new[[i]] <- out.tmp[[i]]$tau.sq.beta.samples
        if (object$dist == 'NB') {
          kappa.samples.new[[i]] <- out.tmp[[i]]$kappa.samples
        }
	      if (object$muRE) {
          sigma.sq.mu.samples.new[[i]] <- out.tmp[[i]]$sigma.sq.mu.samples
          beta.star.samples.new[[i]] <- out.tmp[[i]]$beta.star.samples
	      }
        mu.samples.new[[i]] <- out.tmp[[i]]$mu.samples
        y.rep.samples.new[[i]] <- out.tmp[[i]]$y.rep.samples
        like.samples.new[[i]] <- out.tmp[[i]]$like.samples
      }
    }
    # Update Gelman-Rubin diagnostics. 
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable. 
      rhat.new$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(beta.comm.samples.new, function(a) 
        					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      rhat.new$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(tau.sq.beta.samples.new, function(a) 
      					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      rhat.new$beta <- as.vector(gelman.diag(mcmc.list(lapply(beta.samples.new, function(a) 
      					         mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      if (object$dist == 'NB') {
        rhat.new$kappa <- gelman.diag(mcmc.list(lapply(kappa.samples.new, function(a) 
        					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
      }
      if (object$muRE) {
        rhat.new$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(sigma.sq.mu.samples.new, function(a) 
        					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }
    } else {
      rhat.new$beta.comm <- rep(NA, p.abund)
      rhat.new$tau.sq.beta <- rep(NA, p.abund)
      rhat.new$beta <- rep(NA, p.abund * n.sp)
      if (object$dist == 'NB') {
        rhat.new$kappa <- rep(NA, n.sp)
      }
      if (object$muRE > 0) {
        rhat.new$sigma.sq.mu <- rep(NA, ncol(object$sigma.sq.mu.samples))
      }
    }
    object$rhat <- rhat.new

    object$beta.comm.samples <- mcmc(do.call(rbind, beta.comm.samples.new))
    object$tau.sq.beta.samples <- mcmc(do.call(rbind, tau.sq.beta.samples.new))
    object$beta.samples <- mcmc(do.call(rbind, beta.samples.new))
    if (object$dist == 'NB') {
      object$kappa.samples <- mcmc(do.call(rbind, kappa.samples.new))
    }
    if (object$muRE) {
      object$sigma.sq.mu.samples <- mcmc(do.call(rbind, sigma.sq.mu.samples.new))
      object$beta.star.samples <- mcmc(do.call(rbind, beta.star.samples.new))
    }
    object$mu.samples <- do.call(abind, list('...' = mu.samples.new, 
                                 along = 1))
    object$y.rep.samples <- do.call(abind, list('...' = y.rep.samples.new, 
                                    along = 1))
    object$like.samples <- do.call(abind, list('...' = like.samples.new, 
                                   along = 1))
    object$ESS <- list()
    # Calculate effective sample sizes
    object$ESS$beta.comm <- effectiveSize(object$beta.comm.samples)
    object$ESS$tau.sq.beta <- effectiveSize(object$tau.sq.beta.samples)
    object$ESS$beta <- effectiveSize(object$beta.samples)
    if (object$dist == 'NB') {
      object$ESS$kappa <- effectiveSize(object$kappa.samples)
    }
    if (object$muRE) {
      object$ESS$sigma.sq.mu <- effectiveSize(object$sigma.sq.mu.samples)
    }
    object$n.burn <- ifelse(keep.orig, object$n.burn + n.burn, object$n.samples + n.burn)
    object$n.samples <- object$n.samples + n.batch * object$update$batch.length
    n.post.new <- length(seq(from = n.burn + 1, 
                             to = n.batch * object$update$batch.length, 
                             by = as.integer(n.thin)))
    object$n.post <- ifelse(keep.orig, object$n.post + n.post.new, n.post.new)
    object$run.time <- object$run.time + run.time.new 
    object$update$tuning <- matrix(NA, nrow(out.tmp[[1]]$update$tuning), n.chains)
    for (i in 1:n.chains) {
      object$update$tuning[, i] <- out.tmp[[i]]$update$tuning
    }
    object$update$final.seed <- seeds.new
    object$update$n.batch <- n.batch + object$update$n.batch
  }
  # sfJSDM ----------------------------------------------------------------
  if (is(object, 'sfJSDM')) {
    for (i in 1:n.chains) {
      # Set the random seed based on the previous set of the model
      assign(".Random.seed", object$update$final.seed[[i]], .GlobalEnv)
      N <- nrow(object$y)
      p.occ <- ncol(object$beta.comm.samples)
      cov.model <- object$update$cov.model
      q <- object$q
      # Get initial values
      curr.inits <- n.post.one.chain * i
      inits <- list()
      # beta.comm, tau.sq.beta, beta, phi, lambda, nu, sigma.sq.psi, w
      inits$beta.comm <- object$beta.comm.samples[curr.inits, ]
      inits$tau.sq.beta <- object$tau.sq.beta.samples[curr.inits, ]
      inits$beta <- matrix(object$beta.samples[curr.inits, ], N, p.occ)
      if (cov.model != 'matern') {
        inits$phi <- object$theta.samples[curr.inits, ]
      } else {
        inits$phi <- object$theta.samples[curr.inits, 1:q]
        inits$nu <- object$theta.samples[curr.inits, (q+1):(q*2)]
      }
      inits$lambda <- matrix(object$lambda.samples[curr.inits, ], N, q)
      if (object$psiRE) {
        inits$sigma.sq.psi <- object$sigma.sq.psi.samples[curr.inits, ]
        inits$beta.star <- t(matrix(object$beta.star.samples[curr.inits, ], 
          			  ncol = N))
      }
      inits$w <- object$w.samples[curr.inits, , ]
      if (q == 1) {
        inits$w <- t(as.matrix(inits$w))
      }
      # Get tuning values
      tuning <- list()
      sigma.sq.indx <- 1
      phi.indx <- 2
      nu.indx <- 3
      tuning$phi <- object$update$tuning[((phi.indx - 1) * q + 1):(phi.indx * q), i]
      if (cov.model == 'matern') {
        tuning$nu <- object$update$tuning[((nu.indx - 1) * q + 1):(nu.indx * q), i]
      }
      out.tmp[[i]] <- sfJSDM(formula = object$update$formula,
                             data = object$update$data,
                             inits = inits,
                             priors = object$update$priors,
                             tuning = tuning,
                             cov.model = object$update$cov.model,
                             NNGP = ifelse(object$type == 'NNGP', TRUE, FALSE),
                             n.neighbors = object$n.neighbors,
                             search.type = object$update$search.type,
                             n.factors = object$q,
                             n.batch = n.batch,
                             batch.length = object$update$batch.length,
                             accept.rate = object$update$accept.rate,
                             n.omp.threads = object$update$n.omp.threads,
                             verbose = verbose,
                             n.report = n.report,
                             n.burn = n.burn,
                             n.thin = n.thin,
                             n.chains = 1) # TODO: will make output look weird.
      run.time.new <- run.time.new + out.tmp[[i]]$run.time
      seeds.new[[i]] <- out.tmp[[i]]$update$final.seed[[1]]
    }
    # Put everything together
    beta.samples.new <- list()
    beta.comm.samples.new <- list()
    tau.sq.beta.samples.new <- list()
    theta.samples.new <- list()
    lambda.samples.new <- list()
    sigma.sq.psi.samples.new <- list()
    beta.star.samples.new <- list()
    w.samples.new <- list()
    psi.samples.new <- list()
    z.samples.new <- list()
    like.samples.new <- list()

    rhat.new <- list()
    ess.new <- list()
    n.samples.one.chain <- object$n.post
    for (i in 1:n.chains) {
      if (keep.orig) {
        beta.comm.samples.new[[i]] <- rbind(object$beta.comm.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.comm.samples)
        beta.samples.new[[i]] <- rbind(object$beta.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.samples)
        tau.sq.beta.samples.new[[i]] <- rbind(object$tau.sq.beta.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$tau.sq.beta.samples)
        theta.samples.new[[i]] <- rbind(object$theta.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$theta.samples)
        lambda.samples.new[[i]] <- rbind(object$lambda.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$lambda.samples)
        w.samples.new[[i]] <- abind(object$w.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$w.samples, along = 1)
        if (object$psiRE) {
          sigma.sq.psi.samples.new[[i]] <- rbind(object$sigma.sq.psi.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$sigma.sq.psi.samples)
          beta.star.samples.new[[i]] <- rbind(object$beta.star.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.star.samples)
        }
        psi.samples.new[[i]] <- abind(object$psi.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$psi.samples, along = 1)
        like.samples.new[[i]] <- abind(object$like.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$like.samples, along = 1)
        z.samples.new[[i]] <- abind(object$z.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$z.samples, along = 1)
      } else {
        beta.comm.samples.new[[i]] <- out.tmp[[i]]$beta.comm.samples
        beta.samples.new[[i]] <- out.tmp[[i]]$beta.samples
        tau.sq.beta.samples.new[[i]] <- out.tmp[[i]]$tau.sq.beta.samples
        theta.samples.new[[i]] <- out.tmp[[i]]$theta.samples
        lambda.samples.new[[i]] <- out.tmp[[i]]$lambda.samples
        w.samples.new[[i]] <- out.tmp[[i]]$w.samples
	if (object$psiRE) {
          sigma.sq.psi.samples.new[[i]] <- out.tmp[[i]]$sigma.sq.psi.samples
          beta.star.samples.new[[i]] <- out.tmp[[i]]$beta.star.samples
	}
        psi.samples.new[[i]] <- out.tmp[[i]]$psi.samples
        z.samples.new[[i]] <- out.tmp[[i]]$z.samples
        like.samples.new[[i]] <- out.tmp[[i]]$like.samples
      }
    }
    # Update Gelman-Rubin diagnostics. 
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable. 
      rhat.new$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(beta.comm.samples.new, function(a) 
        					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      rhat.new$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(tau.sq.beta.samples.new, function(a) 
      					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      rhat.new$beta <- as.vector(gelman.diag(mcmc.list(lapply(beta.samples.new, function(a) 
      					         mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      rhat.new$theta <- gelman.diag(mcmc.list(lapply(theta.samples.new, function(a) 
      					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
      rhat.new$lambda.lower.tri <- as.vector(gelman.diag(mcmc.list(lapply(lambda.samples.new, function(a) 
        					       mcmc(a[, c(lower.tri(inits$lambda))]))), 
        					       autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      if (object$psiRE) {
        rhat.new$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(sigma.sq.psi.samples.new, function(a) 
        					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }
    } else {
      rhat.new$beta.comm <- rep(NA, p.occ)
      rhat.new$tau.sq.beta <- rep(NA, p.occ)
      rhat.new$beta <- rep(NA, p.occ * N)
      rhat.new$theta <- rep(NA, ifelse(cov.model == 'matern', 2 * q, q))
      if (object$psiRE > 0) {
        rhat.new$sigma.sq.psi <- rep(NA, ncol(object$sigma.sq.psi.samples))
      }
    }
    object$rhat <- rhat.new

    object$beta.comm.samples <- mcmc(do.call(rbind, beta.comm.samples.new))
    object$tau.sq.beta.samples <- mcmc(do.call(rbind, tau.sq.beta.samples.new))
    object$beta.samples <- mcmc(do.call(rbind, beta.samples.new))
    if (object$psiRE) {
      object$sigma.sq.psi.samples <- mcmc(do.call(rbind, sigma.sq.psi.samples.new))
      object$beta.star.samples <- mcmc(do.call(rbind, beta.star.samples.new))
    }
    object$lambda.samples <- mcmc(do.call(rbind, lambda.samples.new))
    object$theta.samples <- mcmc(do.call(rbind, theta.samples.new))
    object$w.samples <- do.call(abind, list('...' = w.samples.new, 
          				  along = 1))
    object$psi.samples <- do.call(abind, list('...' = psi.samples.new, 
          				  along = 1))
    object$z.samples <- do.call(abind, list('...' = z.samples.new, 
          				  along = 1))
    object$like.samples <- do.call(abind, list('...' = like.samples.new, 
          				  along = 1))
    object$ESS <- list()
    # Calculate effective sample sizes
    object$ESS$beta.comm <- effectiveSize(object$beta.comm.samples)
    object$ESS$tau.sq.beta <- effectiveSize(object$tau.sq.beta.samples)
    object$ESS$beta <- effectiveSize(object$beta.samples)
    object$ESS$theta <- effectiveSize(object$theta.samples)
    object$ESS$lambda <- effectiveSize(object$lambda.samples)
    if (object$psiRE) {
      object$ESS$sigma.sq.psi <- effectiveSize(object$sigma.sq.psi.samples)
    }
    object$n.burn <- ifelse(keep.orig, object$n.burn + n.burn, object$n.samples + n.burn)
    object$n.samples <- object$n.samples + n.batch * object$update$batch.length
    n.post.new <- length(seq(from = n.burn + 1, 
                             to = n.batch * object$update$batch.length, 
                             by = as.integer(n.thin)))
    object$n.post <- ifelse(keep.orig, object$n.post + n.post.new, n.post.new)
    # TODO: note the thinning rate may be different across models. Just ignoring
    #       this for now, but may want to update for summary. Note these values
    #       also might not be correct if keep.orig = FALSE, which is something
    #       you should look into. 
    object$run.time <- object$run.time + run.time.new 
    tmp.val <- ifelse(cov.model == 'matern', q * 3, q * 2)
    object$update$tuning <- matrix(NA, tmp.val, n.chains)
    for (i in 1:n.chains) {
      object$update$tuning[, i] <- out.tmp[[i]]$update$tuning
    }
    object$update$final.seed <- seeds.new
    object$update$n.batch <- n.batch + object$update$n.batch
  } # sfJSDM 
  # lfJSDM ----------------------------------------------------------------
  if (is(object, 'lfJSDM')) {
    for (i in 1:n.chains) {
      # Set the random seed based on the previous set of the model
      assign(".Random.seed", object$update$final.seed[[i]], .GlobalEnv)
      N <- nrow(object$y)
      p.occ <- ncol(object$beta.comm.samples)
      q <- object$q
      # Get initial values
      curr.inits <- n.post.one.chain * i
      inits <- list()
      # beta.comm, tau.sq.beta, beta, phi, lambda, nu, sigma.sq.psi, w
      inits$beta.comm <- object$beta.comm.samples[curr.inits, ]
      inits$tau.sq.beta <- object$tau.sq.beta.samples[curr.inits, ]
      inits$beta <- matrix(object$beta.samples[curr.inits, ], N, p.occ)
      if (q > 0) {
        inits$lambda <- matrix(object$lambda.samples[curr.inits, ], N, q)
      }
      if (object$psiRE) {
        inits$sigma.sq.psi <- object$sigma.sq.psi.samples[curr.inits, ]
        inits$beta.star <- t(matrix(object$beta.star.samples[curr.inits, ], 
          			  ncol = N))
      }
      if (q > 0) {
        inits$w <- object$w.samples[curr.inits, , ]
      }
      if (q == 1) {
        inits$w <- t(as.matrix(inits$w))
      }
      out.tmp[[i]] <- lfJSDM(formula = object$update$formula,
                             data = object$update$data,
                             inits = inits,
                             priors = object$update$priors,
                             n.factors = object$q,
                             n.samples = n.samples,
                             n.omp.threads = object$update$n.omp.threads,
                             verbose = verbose,
                             n.report = n.report,
                             n.burn = n.burn,
                             n.thin = n.thin,
                             n.chains = 1) # TODO: will make output look weird.
      run.time.new <- run.time.new + out.tmp[[i]]$run.time
      seeds.new[[i]] <- out.tmp[[i]]$update$final.seed[[1]]
    }
    # Put everything together
    beta.samples.new <- list()
    beta.comm.samples.new <- list()
    tau.sq.beta.samples.new <- list()
    theta.samples.new <- list()
    lambda.samples.new <- list()
    sigma.sq.psi.samples.new <- list()
    beta.star.samples.new <- list()
    w.samples.new <- list()
    psi.samples.new <- list()
    z.samples.new <- list()
    like.samples.new <- list()

    rhat.new <- list()
    ess.new <- list()
    n.samples.one.chain <- object$n.post
    for (i in 1:n.chains) {
      if (keep.orig) {
        beta.comm.samples.new[[i]] <- rbind(object$beta.comm.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.comm.samples)
        beta.samples.new[[i]] <- rbind(object$beta.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.samples)
        tau.sq.beta.samples.new[[i]] <- rbind(object$tau.sq.beta.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$tau.sq.beta.samples)
	if (q > 0) {
          lambda.samples.new[[i]] <- rbind(object$lambda.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$lambda.samples)
          w.samples.new[[i]] <- abind(object$w.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$w.samples, along = 1)
	}
        if (object$psiRE) {
          sigma.sq.psi.samples.new[[i]] <- rbind(object$sigma.sq.psi.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$sigma.sq.psi.samples)
          beta.star.samples.new[[i]] <- rbind(object$beta.star.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , drop = FALSE], out.tmp[[i]]$beta.star.samples)
        }
        psi.samples.new[[i]] <- abind(object$psi.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$psi.samples, along = 1)
        like.samples.new[[i]] <- abind(object$like.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$like.samples, along = 1)
        z.samples.new[[i]] <- abind(object$z.samples[((i - 1) * n.samples.one.chain + 1):(i * n.samples.one.chain), , , drop = FALSE], out.tmp[[i]]$z.samples, along = 1)
      } else {
        beta.comm.samples.new[[i]] <- out.tmp[[i]]$beta.comm.samples
        beta.samples.new[[i]] <- out.tmp[[i]]$beta.samples
        tau.sq.beta.samples.new[[i]] <- out.tmp[[i]]$tau.sq.beta.samples
	if (q > 0) {
          lambda.samples.new[[i]] <- out.tmp[[i]]$lambda.samples
          w.samples.new[[i]] <- out.tmp[[i]]$w.samples
	}
	if (object$psiRE) {
          sigma.sq.psi.samples.new[[i]] <- out.tmp[[i]]$sigma.sq.psi.samples
          beta.star.samples.new[[i]] <- out.tmp[[i]]$beta.star.samples
	}
        psi.samples.new[[i]] <- out.tmp[[i]]$psi.samples
        z.samples.new[[i]] <- out.tmp[[i]]$z.samples
        like.samples.new[[i]] <- out.tmp[[i]]$like.samples
      }
    }
    # Update Gelman-Rubin diagnostics. 
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable. 
      rhat.new$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(beta.comm.samples.new, function(a) 
        					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      rhat.new$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(tau.sq.beta.samples.new, function(a) 
      					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      rhat.new$beta <- as.vector(gelman.diag(mcmc.list(lapply(beta.samples.new, function(a) 
      					         mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      if (q > 0) {
        rhat.new$lambda.lower.tri <- as.vector(gelman.diag(mcmc.list(lapply(lambda.samples.new, function(a) 
          					       mcmc(a[, c(lower.tri(inits$lambda))]))), 
          					       autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }
      if (object$psiRE) {
        rhat.new$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(sigma.sq.psi.samples.new, function(a) 
        					      mcmc(a))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }
    } else {
      rhat.new$beta.comm <- rep(NA, p.occ)
      rhat.new$tau.sq.beta <- rep(NA, p.occ)
      rhat.new$beta <- rep(NA, p.occ * N)
      if (object$psiRE > 0) {
        rhat.new$sigma.sq.psi <- rep(NA, ncol(object$sigma.sq.psi.samples))
      }
    }
    object$rhat <- rhat.new

    object$beta.comm.samples <- mcmc(do.call(rbind, beta.comm.samples.new))
    object$tau.sq.beta.samples <- mcmc(do.call(rbind, tau.sq.beta.samples.new))
    object$beta.samples <- mcmc(do.call(rbind, beta.samples.new))
    if (object$psiRE) {
      object$sigma.sq.psi.samples <- mcmc(do.call(rbind, sigma.sq.psi.samples.new))
      object$beta.star.samples <- mcmc(do.call(rbind, beta.star.samples.new))
    }
    if (q > 0) {
      object$lambda.samples <- mcmc(do.call(rbind, lambda.samples.new))
      object$w.samples <- do.call(abind, list('...' = w.samples.new, along = 1))
    }
    object$psi.samples <- do.call(abind, list('...' = psi.samples.new, 
          				  along = 1))
    object$z.samples <- do.call(abind, list('...' = z.samples.new, 
          				  along = 1))
    object$like.samples <- do.call(abind, list('...' = like.samples.new, 
          				  along = 1))
    object$ESS <- list()
    # Calculate effective sample sizes
    object$ESS$beta.comm <- effectiveSize(object$beta.comm.samples)
    object$ESS$tau.sq.beta <- effectiveSize(object$tau.sq.beta.samples)
    object$ESS$beta <- effectiveSize(object$beta.samples)
    if (q > 0) {
      object$ESS$lambda <- effectiveSize(object$lambda.samples)
    }
    if (object$psiRE) {
      object$ESS$sigma.sq.psi <- effectiveSize(object$sigma.sq.psi.samples)
    }
    object$n.burn <- ifelse(keep.orig, object$n.burn + n.burn, object$n.samples + n.burn)
    object$n.samples <- object$n.samples + n.samples 
    n.post.new <- length(seq(from = n.burn + 1, 
			     to = n.samples,
                             by = as.integer(n.thin)))
    object$n.post <- ifelse(keep.orig, object$n.post + n.post.new, n.post.new)
    # TODO: note the thinning rate may be different across models. Just ignoring
    #       this for now, but may want to update for summary. Note these values
    #       also might not be correct if keep.orig = FALSE, which is something
    #       you should look into. 
    object$run.time <- object$run.time + run.time.new 
    object$update$final.seed <- seeds.new
    object$update$n.samples <- n.samples + object$update$n.samples
  } # lfJSDM 
  return(object)
}
