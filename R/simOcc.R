#' Simulate Single-Species Occupancy Data
#'
#' @param J.x Number of cells on x-axis
#' @param J.y Number of cells on y-axis
#' @param K Number of repeat visits
#' @param beta Occupancy covariate values
#' @param alpha Detection covariate values
#' @param sigma.sq Spatial variance parameter
#' @param phi Spatial range parameter
#' @param sp Logical indicating if simulate spatial model or not.
#'
#' @return stuff
#' @export
#'
#' @importFrom stats rnorm rbinom dist
#'
#' @examples J.x <- 10
#'J.y <- 10
#'K <- rep(4, J.x * J.y)
#'beta <- c(0.5, -0.15)
#'alpha <- c(0.7, 0.4)
#'phi <- 3 / .6
#'sigma.sq <- 2
#'dat <- simOcc(J.x = J.x, J.y = J.y, K = K, beta = beta, alpha = alpha,
#'              sigma.sq = sigma.sq, phi = phi, sp = TRUE)

simOcc <- function(J.x, J.y, K, beta, alpha, sigma.sq = 2, phi = 3/0.5,
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
  X.psi <- matrix(1, nrow = J, ncol = n.beta)
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      X.psi[, i] <- rnorm(J)
    } # i
  }

  # Form detection covariate (if any) -------------------------------------
  n.alpha <- length(alpha)
  X.p <- array(NA, dim = c(J, max(K), n.alpha))
  X.p[, , 1] <- 1
  if (n.alpha > 1) {
    for (i in 2:n.alpha) {
      for (j in 1:J) {
        X.p[j, 1:K[j], i] <- rnorm(K[j])
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
    psi <- logit.inv(X.psi %*% as.matrix(beta) + w)
  } else {
    psi <- logit.inv(X.psi %*% as.matrix(beta))
  }
  z <- rbinom(J, 1, psi)

  # Data Formation --------------------------------------------------------
  p <- matrix(NA, nrow = J, ncol = max(K))
  y <- matrix(NA, nrow = J, ncol = max(K))
  for (j in 1:J) {
    p[j, 1:K[j]] <- logit.inv(X.p[j, 1:K[j], ] %*% as.matrix(alpha))
    y[j, 1:K[j]] <- rbinom(K[j], 1, p[j, 1:K[j]] * z[j])
  } # j

  return(
    list(X.psi = X.psi, X.p = X.p, coords = coords, D = D, R = R,
         w = w, psi = psi, z = z, p = p, y = y, s.x = s.x, s.y = s.y)
  )
}
