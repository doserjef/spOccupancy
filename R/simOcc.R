simOcc <- function(J.x, J.y, n.rep, beta, alpha, sigma.sq = 2, phi = 3/0.5,
                   sp = FALSE) {

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
  X <- matrix(1, nrow = J, ncol = n.beta)
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      X[, i] <- rnorm(J)
    } # i
  }

  # Form detection covariate (if any) -------------------------------------
  n.alpha <- length(alpha)
  X.p <- array(NA, dim = c(J, max(n.rep), n.alpha))
  X.p[, , 1] <- 1
  if (n.alpha > 1) {
    for (i in 2:n.alpha) {
      for (j in 1:J) {
        X.p[j, 1:n.rep[j], i] <- rnorm(n.rep[j])
      } # j
    } # i
  }

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
    psi <- logit.inv(X %*% as.matrix(beta) + w)
  } else {
    psi <- logit.inv(X %*% as.matrix(beta))
  }
  z <- rbinom(J, 1, psi)

  # Data Formation --------------------------------------------------------
  p <- matrix(NA, nrow = J, ncol = max(n.rep))
  y <- matrix(NA, nrow = J, ncol = max(n.rep))
  for (j in 1:J) {
    p[j, 1:n.rep[j]] <- logit.inv(X.p[j, 1:n.rep[j], ] %*% as.matrix(alpha))
    y[j, 1:n.rep[j]] <- rbinom(n.rep[j], 1, p[j, 1:n.rep[j]] * z[j])
  } # j

  return(
    list(X = X, X.p = X.p, coords = coords, D = D,
         w = w, psi = psi, z = z, p = p, y = y)
  )
}
