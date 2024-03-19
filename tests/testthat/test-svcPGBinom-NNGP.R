# Test svcPGBinom.R  ---------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
# weights <- sample(2:4, J, replace = TRUE)
weights <- rep(2, J)
beta <- c(0.3)
p.occ <- length(beta)
psi.RE <- list()
svc.cols <- 1
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simBinom(J.x = J.x, J.y = J.y, weights = weights, beta = beta, 
	      psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
weights.0 <- weights[pred.indx]
weights.fit <- weights[-pred.indx]

data.list <- list(y = y, coords = coords, weights = weights.fit)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(a = 0.5, b = 2.5))
# Starting values
inits.list <- list(beta = 0, fix = TRUE)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ 1

out <- svcPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       svc.cols = svc.cols,
	       cov.model = "matern",
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcPGBinom", {
  expect_s3_class(out, "svcPGBinom")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- svcPGBinom(formula = formula, 
	       data = data.list, 
               n.thin = 13,
               n.batch = n.batch, 
               batch.length = batch.length, 
               accept.rate = 0.43, 
	       n.omp.threads = 1,
	       verbose = FALSE))
})
# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")
})

test_that("all correlation functions work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")

  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")


  out <- svcPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "matern", 
		 priors = list(nu.unif = list(0.5, 2)),
	         tuning = list(phi = 0.5, nu = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "svcPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcPGBinom(formula = formula, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(length(dim(fitted.out)), 2)
})

# Check predictions -------------------
test_that("predict works for svcPGBinom", {
  pred.out <- predict(out, X.0, coords.0, weights.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})

# Occurrence covariate only -----------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
# weights <- sample(2:4, J, replace = TRUE)
weights <- rep(2, J)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
psi.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- runif(p.svc, 3 /.9, 3 / .1)
sigma.sq <- runif(p.svc, 0.1, 2)
nu <- runif(p.svc, 1, 2)
dat <- simBinom(J.x = J.x, J.y = J.y, weights = weights, beta = beta, 
	      psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
weights.0 <- weights[pred.indx]
weights.fit <- weights[-pred.indx]

covs <- cbind(X)
colnames(covs) <- c('int', 'cov.1')
data.list <- list(y = y, coords = coords, weights = weights.fit, covs = covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(a = 0.5, b = 2.5))
# Starting values
inits.list <- list(beta = 0, fix = TRUE)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ cov.1

out <- svcPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       svc.cols = svc.cols,
	       cov.model = "matern",
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcPGBinom", {
  expect_s3_class(out, "svcPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")
})

test_that("all correlation functions work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")

  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")


  out <- svcPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "matern", 
		 priors = list(nu.unif = list(0.5, 2)),
	         tuning = list(phi = 0.5, nu = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "svcPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcPGBinom(formula = formula, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(length(dim(fitted.out)), 2)
})

# Check predictions -------------------
test_that("predict works for svcPGBinom", {
  pred.out <- predict(out, X.0, coords.0, weights.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$covs[3, ] <- NA
  expect_error(svcPGBinom(formula = formula, 
	               data = tmp.data, 
	               n.batch = 40, 
	               batch.length = batch.length, 
	               cov.model = "exponential", 
	               tuning = tuning.list, 
		       svc.cols = svc.cols,
	               NNGP = TRUE,
	               verbose = FALSE, 
	               n.neighbors = 5, 
	               search.type = 'cb', 
	               n.report = 10, 
	               n.burn = 500, 
	               n.chains = 1))
})

# Interactions ------------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
# weights <- sample(2:4, J, replace = TRUE)
weights <- rep(2, J)
beta <- c(0.3, 0.5, -0.5)
p.occ <- length(beta)
psi.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- runif(p.svc, 3 /.9, 3 / .1)
sigma.sq <- runif(p.svc, 0.1, 2)
nu <- runif(p.svc, 1, 2)
dat <- simBinom(J.x = J.x, J.y = J.y, weights = weights, beta = beta, 
	      psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
weights.0 <- weights[pred.indx]
weights.fit <- weights[-pred.indx]

covs <- cbind(X)
colnames(covs) <- c('int', 'cov.1', 'cov.2')
data.list <- list(y = y, coords = coords, weights = weights.fit, covs = covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(a = 0.5, b = 2.5))
# Starting values
inits.list <- list(beta = 0, fix = TRUE)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ cov.1 * cov.2

out <- svcPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       svc.cols = svc.cols,
	       cov.model = "matern",
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcPGBinom", {
  expect_s3_class(out, "svcPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")
})

test_that("all correlation functions work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")

  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")


  out <- svcPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "matern", 
		 priors = list(nu.unif = list(0.5, 2)),
	         tuning = list(phi = 0.5, nu = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "svcPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcPGBinom(formula = formula, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(length(dim(fitted.out)), 2)
})

# Check predictions -------------------
test_that("predict works for svcPGBinom", {
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0, coords.0, weights.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$covs[3, ] <- NA
  expect_error(svcPGBinom(formula = formula, 
	               data = tmp.data, 
	               n.batch = 40, 
	               batch.length = batch.length, 
	               cov.model = "exponential", 
	               tuning = tuning.list, 
		       svc.cols = svc.cols,
	               NNGP = TRUE,
	               verbose = FALSE, 
	               n.neighbors = 5, 
	               search.type = 'cb', 
	               n.report = 10, 
	               n.burn = 500, 
	               n.chains = 1))
})

# Random intercept on occurrence ------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
weights <- sample(1:4, J, replace = TRUE)
# weights <- rep(2, J)
beta <- c(0.3)
p.occ <- length(beta)
psi.RE <- list(levels = c(25), 
	       sigma.sq.psi = c(1.3))
svc.cols <- c(1)
p.svc <- length(svc.cols)
phi <- runif(p.svc, 3 /.9, 3 / .1)
sigma.sq <- runif(p.svc, 0.1, 2)
nu <- runif(p.svc, 1, 2)
dat <- simBinom(J.x = J.x, J.y = J.y, weights = weights, beta = beta, 
	      psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
weights.0 <- weights[pred.indx]
weights.fit <- weights[-pred.indx]

covs <- cbind(X, X.re)
colnames(covs) <- c('int', 'factor.1')
data.list <- list(y = y, coords = coords, weights = weights.fit, covs = covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(a = 0.5, b = 2.5))
# Starting values
inits.list <- list(beta = 0, fix = TRUE)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ (1 | factor.1)

out <- svcPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       svc.cols = svc.cols,
	       cov.model = "matern",
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcPGBinom", {
  expect_s3_class(out, "svcPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")
})

test_that("all correlation functions work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")

  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")


  out <- svcPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "matern", 
		 priors = list(nu.unif = list(0.5, 2)),
	         tuning = list(phi = 0.5, nu = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "svcPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcPGBinom(formula = formula, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(length(dim(fitted.out)), 2)
})

# Check predictions -------------------
test_that("predict works for svcPGBinom", {
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('(Intercept)', 'factor.1')
  pred.out <- predict(out, X.0, coords.0, weights.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$covs[3, ] <- NA
  expect_error(svcPGBinom(formula = formula, 
	               data = tmp.data, 
	               n.batch = 40, 
	               batch.length = batch.length, 
	               cov.model = "exponential", 
	               tuning = tuning.list, 
		       svc.cols = svc.cols,
	               NNGP = TRUE,
	               verbose = FALSE, 
	               n.neighbors = 5, 
	               search.type = 'cb', 
	               n.report = 10, 
	               n.burn = 500, 
	               n.chains = 1))
})

# Multiple random intercepts + covaraites ---------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
weights <- sample(1:4, J, replace = TRUE)
# weights <- rep(2, J)
beta <- c(0.3, 0.3, -0.5)
p.occ <- length(beta)
psi.RE <- list(levels = c(25, 10), 
	       sigma.sq.psi = c(1.3, 0.3))
svc.cols <- c(1, 3)
p.svc <- length(svc.cols)
phi <- runif(p.svc, 3 /.9, 3 / .1)
sigma.sq <- runif(p.svc, 0.1, 2)
nu <- runif(p.svc, 1, 2)
dat <- simBinom(J.x = J.x, J.y = J.y, weights = weights, beta = beta, 
	      psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
weights.0 <- weights[pred.indx]
weights.fit <- weights[-pred.indx]

covs <- cbind(X, X.re)
colnames(covs) <- c('int', 'cov.1', 'cov.2', 'factor.1', 'factor.2')
data.list <- list(y = y, coords = coords, weights = weights.fit, covs = covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(a = 0.5, b = 2.5))
# Starting values
inits.list <- list(beta = 0, fix = TRUE)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ cov.1 + cov.2 + (1 | factor.1) + (1 | factor.2)

out <- svcPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       svc.cols = svc.cols,
	       cov.model = "matern",
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcPGBinom", {
  expect_s3_class(out, "svcPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")
})

test_that("all correlation functions work", {
  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")

  out <- svcPGBinom(formula = formula, 
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
  expect_s3_class(out, "svcPGBinom")


  out <- svcPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "matern", 
		 priors = list(nu.unif = list(0.5, 2)),
	         tuning = list(phi = 0.5, nu = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "svcPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcPGBinom(formula = formula, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(length(dim(fitted.out)), 2)
})

# Check predictions -------------------
test_that("predict works for svcPGBinom", {
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('(Intercept)', 'cov.1', 'cov.2', 'factor.1', 'factor.2')
  pred.out <- predict(out, X.0, coords.0, weights.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$covs[3, ] <- NA
  expect_error(svcPGBinom(formula = formula, 
	               data = tmp.data, 
	               n.batch = 40, 
	               batch.length = batch.length, 
	               cov.model = "exponential", 
	               tuning = tuning.list, 
		       svc.cols = svc.cols,
	               NNGP = TRUE,
	               verbose = FALSE, 
	               n.neighbors = 5, 
	               search.type = 'cb', 
	               n.report = 10, 
	               n.burn = 500, 
	               n.chains = 1))
})

# Uniform sigma sq --------------------------------------------------------
test_that("svcPGBinom works with uniform prior on sigma.sq", {
  prior.list <- list(sigma.sq.unif = list(0, 5), 
                     nu.unif = list(0.1, 4))
  tuning.list <- list(phi = 0.5, nu = 0.6, sigma.sq = 0.7)
  out <- svcPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
		 svc.cols = svc.cols,
	         cov.model = "exponential", 
		 priors = prior.list,
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "svcPGBinom")
  out <- svcPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "matern", 
		 priors = prior.list,
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 500, 
	         n.chains = 1)
  expect_s3_class(out, "svcPGBinom")
})
