set.seed(400)
J.x <- 10
J.y <- 10
n.rep <- rep(4, J.x * J.y)
beta <- c(0.5, -0.15)
alpha <- c(0.7, 0.4)
phi <- 3 / .6
sigma.sq <- 2
psi.RE <- list(levels = 10, 
	       sigma.sq.psi = 1.2)
p.RE <- list(levels = 15, 
           sigma.sq.p = 0.8)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, cov.model = 'spherical', 
	      sigma.sq = sigma.sq, phi = phi)

test_that("simulated object is a list", {
  expect_type(dat, "list")
})

test_that("X is correct", {
  expect_equal(dim(dat$X), c(J.x * J.y, length(beta)))
})

test_that("X.p is correct", {
  expect_equal(dim(dat$X.p), c(J.x * J.y, max(n.rep), length(alpha)))
})

test_that("y is correct", {
  expect_equal(dim(dat$y), c(J.x * J.y, max(n.rep)))
})
