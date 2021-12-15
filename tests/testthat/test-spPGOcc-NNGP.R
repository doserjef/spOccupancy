# Test spPGOcc.R  ---------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

# NNGP --------------------------------------------------------------------
set.seed(350)

J.x <- 12
J.y <- 12
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.5, -0.15)
p.occ <- length(beta)
alpha <- c(0.7, 0.4, -0.2)
p.det <- length(alpha)
phi <- 3 / .6
sigma.sq <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha, 
              sigma.sq = sigma.sq, phi = phi, sp = TRUE, cov.model = 'exponential')
y <- dat$y
X <- dat$X
X.p <- dat$X.p
coords <- as.matrix(dat$coords)

# Package all data into a list
occ.covs <- X[, -1, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs, 
		  coords = coords)

# Number of batches
n.batch <- 200
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length
# Priors 
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
		   alpha.normal = list(mean = 0, var = 2.72),
		   sigma.sq.ig = c(2, 2), 
		   phi.unif = c(3/1, 3/.1)) 
# Initial values
inits.list <- list(alpha = 0, beta = 0,
		   phi = 3 / .5, 
		   sigma.sq = 2,
		   w = rep(0, nrow(X)),
		   z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1) 

out <- spPGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov.1 + det.cov.2, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       cov.model = "exponential", 
	       tuning = tuning.list, 
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       search.type = 'cb', 
	       verbose = FALSE,
	       n.report = 10, 
	       n.burn = 4000, 
	       n.thin = 2,
	       n.chains = 1)

test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

test_that("samples are the right size", {
  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = out$n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_equal(dim(out$beta.samples), c(n.post.samples, length(beta)))
  expect_equal(dim(out$alpha.samples), c(n.post.samples, length(alpha)))
  expect_equal(dim(out$z.samples), c(n.post.samples, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, J))
})

test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
})

test_that("alphas are correct", {
  alpha.quants <- apply(out$alpha.samples, 2, quantile, c(0.025, 0.975)) 
  expect_true(sum(alpha.quants[1, ] < alpha) == length(alpha))
  expect_true(sum(alpha.quants[2, ] > alpha) == length(alpha))
})

test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = ~ occ.cov, 
	         det.formula = ~ det.cov.1 + det.cov.2, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)

  expect_s3_class(out, "spPGOcc")
})

test_that("all correlation functions work", {
  out <- spPGOcc(occ.formula = ~ occ.cov, 
	         det.formula = ~ det.cov.1 + det.cov.2, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "gaussian", 
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "spPGOcc")

  out <- spPGOcc(occ.formula = ~ occ.cov, 
	         det.formula = ~ det.cov.1 + det.cov.2, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "spherical", 
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "spPGOcc")

  out <- spPGOcc(occ.formula = ~ occ.cov, 
	         det.formula = ~ det.cov.1 + det.cov.2, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "matern", 
		 priors = list(nu.unif = c(0.5, 2)),
	         tuning = list(phi = 0.5, nu = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "spPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spPGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov.1 + det.cov.2, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       cov.model = "exponential", 
	       tuning = tuning.list, 
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       search.type = 'cb', 
	       n.report = 10, 
	       n.burn = 500, 
	       n.chains = 1))
})

test_that("cross-validation works", {
  
  out <- spPGOcc(occ.formula = ~ occ.cov, 
	         det.formula = ~ det.cov.1 + det.cov.2, 
	         data = data.list, 
	         inits = inits.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         priors = prior.list,
	         cov.model = "exponential", 
	         tuning = tuning.list, 
	         NNGP = TRUE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         verbose = FALSE,
	         n.burn = 400, 
	         n.thin = 2,
	         n.chains = 2, 
		 k.fold = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

test_that("random effects on detection work", {
set.seed(500)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(100), 
	     sigma.sq.p = c(2.2))
phi <- 3 / .7
sigma.sq <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
	      phi = phi, cov.model = 'exponential')
y <- dat$y
X <- dat$X
X.p <- dat$X.p
X.p.re <- dat$X.p.re
coords <- as.matrix(dat$coords)

det.covs <- list(det.factor = X.p.re[, , 1])
data.list <- list(y = y, coords = coords, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 100
n.samples <- batch.length * n.batch

out <- spPGOcc(occ.formula = ~ 1,
	       det.formula = ~ (1 | det.factor),
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length, 
	       n.batch = n.batch, 
	       priors = prior.list,
	       accept.rate = 0.43, 
	       cov.model = "exponential",
	       tuning = tuning.list, 
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE, 
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 400,
	       n.thin = 6, 
	       n.chains = 2,
	       k.fold = 2) 

  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = out$n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_s3_class(out, "spPGOcc")
  expect_equal(out$pRE, TRUE)
  expect_equal(dim(out$sigma.sq.p.samples), c(n.post.samples, length(p.RE$sigma.sq.p)))
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# For helper functions
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, -1)
p.occ <- length(beta)
alpha <- c(-0.5, 1)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
sigma.sq <- 2
phi <- 10
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, cov.model = 'matern', 
	      sigma.sq = 2, phi = 10, nu = 2)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

colnames(X) <- c('int', 'occ.cov')
det.covs <- list(det.cov = X.p[, , 2])
coords <- as.matrix(dat$coords)

data.list <- list(y = y, coords = coords, det.covs = det.covs, 
		  occ.covs = X)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72), 
		   nu.unif = c(0.3, 4))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 100

set.seed(400)
out <- spPGOcc(occ.formula = ~ occ.cov,
	       det.formula = ~ det.cov,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length, 
	       n.batch = n.batch, 
	       priors = prior.list,
	       accept.rate = 0.43, 
	       cov.model = "matern",
	       tuning = tuning.list, 
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE, 
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 400,
	       n.thin = 6, 
	       k.fold = 2, 
	       n.chains = 2,
	       k.fold.threads = 1)
n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = out$n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains

test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "array")
  expect_equal(dim(fitted.out), c(n.post.samples, J, max(n.rep)))
})

test_that("predict works for spPGOcc", {
  X.0 <- rbind(dat$X, matrix(c(1, rnorm(p.occ - 1)), nrow = 1))
  coords.0 <- rbind(dat$coords, matrix(c(0.538, 0.201), nrow = 1))
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, J + 1))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, J + 1))
})

test_that("posterior predictive checks work for spPGOcc", {
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

