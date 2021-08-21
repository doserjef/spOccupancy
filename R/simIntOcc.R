simIntOcc <- function(n.data, J.x, J.y, J.obs, K.obs, beta, alpha, 
		      sigma.sq = 2, phi = 3/0.5, sp = FALSE) {

  # Subroutines -----------------------------------------------------------
  # MVN 
  rmvn <- function(n, mu=0, V = matrix(1)) {
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }

  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Form occupancy covariates (if any) ------------------------------------
  J <- J.x * J.y
  n.beta <- length(beta)
  X.psi <- matrix(1, nrow = J, ncol = n.beta) 
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      X.psi[, i] <- rnorm(J)
    } # i
  }

  # Get site ids for each of the data sets --------------------------------
  # Data sources can be obtained at multiple different sites. 
  sites <- list()
  for (i in 1:n.data) {
    sites[[i]] <- sort(sample(1:J, J.obs[i], replace = FALSE))
  }

  # Form detection covariates (if any) ------------------------------------
  X.p <- list()
  for (i in 1:n.data) {
    n.alpha.curr <- length(alpha[[i]])
    K.curr <- K.obs[[i]]
    J.curr <- J.obs[[i]]
    X.p[[i]] <- array(NA, dim = c(J.curr, max(K.curr), n.alpha.curr))
    X.p[[i]][, , 1] <- 1
    if (n.alpha.curr > 1) {
      for (q in 2:n.alpha.curr) {
        for (j in 1:J.curr) {
          X.p[[i]][j, 1:K.curr[j], q] <- rnorm(K.curr[j])
        } # j
      } # q
    }
  } # i

  # Simulate spatial random effect ----------------------------------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- expand.grid(s.x, s.y)
  # Distance matrix
  D <- as.matrix(dist(coords))
  # Exponential correlation function
  R <- exp(-phi * D)
  # Random spatial process
  w <- rmvn(1, rep(0, J), sigma.sq * R)

  # Latent Occupancy Process ----------------------------------------------
  if (sp) {
    psi <- logit.inv(X.psi %*% as.matrix(beta) + w)
  } else {
      psi <- logit.inv(X.psi %*% as.matrix(beta))
    }
  z <- rbinom(J, 1, psi)

  # Data Formation --------------------------------------------------------
  p <- list()
  y <- list()
  for (i in 1:n.data) {
    K.curr <- K.obs[[i]]
    J.curr <- J.obs[[i]]
    p[[i]] <- matrix(NA, nrow = J.curr, ncol = max(K.curr))
    y[[i]] <- matrix(NA, nrow = J.curr, ncol = max(K.curr))
    sites.curr <- sites[[i]]
    X.p.curr <- X.p[[i]]
    alpha.curr <- as.matrix(alpha[[i]])
    for (j in 1:J.curr) {
      p[[i]][j, 1:K.curr[j]] <- logit.inv(X.p.curr[j, 1:K.curr[j], ] %*% alpha.curr)
      y[[i]][j, 1:K.curr[j]] <- rbinom(K.curr[j], 1, p[[i]][j, 1:K.curr[j]] * z[sites.curr[j]])
    } # j
  } # i

  # Split up into observed and predicted ----------------------------------
  sites.obs <- sort(unique(unlist(sites))) 
  sites.pred <- (1:J)[!(1:J %in% sites.obs)]
  X.psi.obs <- X.psi[sites.obs, , drop = FALSE]
  X.psi.pred <- X.psi[sites.pred, , drop = FALSE]
  z.obs <- z[sites.obs]
  z.pred <- z[sites.pred]
  D.obs <- D[sites.obs, sites.obs, drop = FALSE]
  D.pred <- D[sites.pred, sites.pred, drop = FALSE]
  coords.obs <- coords[sites.obs,, drop = FALSE]
  coords.pred <- coords[sites.pred,, drop = FALSE]
  R.obs <- R[sites.obs, sites.obs, drop = FALSE]
  R.pred <- R[sites.pred, sites.pred, drop = FALSE]
  w.obs <- w[sites.obs, , drop = FALSE]
  w.pred <- w[sites.pred, , drop = FALSE]
  psi.obs <- psi[sites.obs, , drop = FALSE]
  psi.pred <- psi[sites.pred, , drop = FALSE]
  sites.vec <- unlist(sites)
  # Need to get index of each site value in only the sites.obs
  sites.new <- rep(0, length(sites.vec))
  for (i in 1:length(sites.vec)) {
    sites.new[i] <- which(sites.obs == sites.vec[i])
  }
  sites.return <- list()
  indx <- 1
  for (i in 1:n.data) {
    sites.return[[i]] <- sites.new[indx:(indx + J.obs[i] - 1)]
    indx <- indx + J.obs[i]
  }

  return(
    list(X.psi.obs = X.psi.obs, X.psi.pred = X.psi.pred, X.p = X.p, 
	 coords.obs = coords.obs, coords.pred = coords.pred, 
	 D.obs = D.obs, D.pred = D.pred, R.obs = R.obs, 
	 R.pred = R.pred,  w.obs = w.obs, w.pred = w.pred, 
	 psi.obs = psi.obs, psi.pred = psi.pred, z.obs = z.obs, 
	 z.pred = z.pred, p = p, y = y, sites = sites.return
	)
  )
}
