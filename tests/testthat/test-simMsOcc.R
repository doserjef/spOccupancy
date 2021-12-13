# Test simMsOCc.R ---------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, size = J, replace = TRUE)
N <- 10
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.15)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.3)
# Detection
alpha.mean <- c(0.5, 0.2)
tau.sq.alpha <- c(0.2, 0.3)
p.det <- length(alpha.mean)
psi.RE <- list(levels = c(10), 
	       sigma.sq.psi = c(1.5))
p.RE <- list(levels = c(15), 
	     sigma.sq.p = 0.8)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
# Spatial parameters if desired
phi <- runif(N, 3/1, 3/.1)
sigma.sq <- runif(N, 0.3, 3)
sp <- TRUE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta,
		alpha = alpha, psi.RE = psi.RE, p.RE = p.RE, sp = TRUE,
		cov.model = 'exponential', phi = phi, sigma.sq = sigma.sq)

test_that("simulated object is a list", {
  expect_type(dat, "list")
})

test_that("X is correct", {
  expect_equal(dim(dat$X), c(J.x * J.y, length(beta.mean)))
})

test_that("X.p is correct", {
  expect_equal(dim(dat$X.p), c(J.x * J.y, max(n.rep), length(alpha.mean)))
})

test_that("y is correct", {
  expect_equal(dim(dat$y), c(N, J.x * J.y, max(n.rep)))
})

test_that("z is correct", {
  expect_equal(dim(dat$z), c(N, J.x * J.y))
})
