# Test simIntOcc.R  ----------------------------------------------------------

# Simulate Data -----------------------------------------------------------
J.x <- 15
J.y <- 15
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(1:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5, 1, -3)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(sample(1:4, 1), -1, 1)
}
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
sigma.sq <- 2
phi <- 3 / .5
sp <- TRUE

dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs,
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = TRUE,
		 cov.model = 'gaussian', sigma.sq = sigma.sq, phi = phi)

test_that("simulated object is a list", {
  expect_type(dat, "list")
})

test_that("X.obs and X.pred are correct", {
  expect_equal(nrow(dat$X.obs) + nrow(dat$X.pred), J.x * J.y)
  expect_equal(ncol(dat$X.obs), length(beta))
  expect_equal(ncol(dat$X.pred), length(beta))
})

test_that("X.p is correct", {
  expect_equal(length(dat$X.p), n.data)
  expect_equal(sapply(dat$X.p, function(a) dim(a)[1]), J.obs)
})

test_that("y is correct", {
  expect_equal(length(dat$y), n.data)
  expect_equal(sapply(dat$y, function(a) dim(a)[1]), J.obs)
  expect_equal(sapply(dat$y, function(a) dim(a)[2]), rep(max(unlist(n.rep)), n.data))
})
