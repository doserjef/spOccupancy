# Test spPGOcc.R  ---------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(237)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

data.list <- list(y = y, coords = coords)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE), fix = TRUE)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ 1
det.formula <- ~ 1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
               n.thin = 13,
               n.batch = n.batch, 
               batch.length = batch.length, 
               accept.rate = 0.43, 
               NNGP = TRUE,
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
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Occurrence covariate only -----------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.9, 1.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, coords = coords, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ 1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spPGOcc(occ.formula = occ.formula, 
  #                      det.formula = det.formula, 
  #                      data = tmp.data, 
  #                      n.batch = 40, 
  #                      batch.length = batch.length, 
  #                      cov.model = "exponential", 
  #                      tuning = tuning.list, 
  #                      NNGP = TRUE,
  #                      verbose = FALSE, 
  #                      n.neighbors = 5, 
  #                      search.type = 'cb', 
  #                      n.report = 10, 
  #                      n.burn = 500, 
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- spPGOcc(occ.formula = occ.formula, 
  #                det.formula = det.formula, 
  #                data = tmp.data, 
  #                n.batch = 40, 
  #                batch.length = batch.length, 
  #                cov.model = "exponential", 
  #                tuning = tuning.list, 
  #                NNGP = TRUE,
  #                verbose = FALSE, 
  #                n.neighbors = 5, 
  #                search.type = 'cb', 
  #                n.report = 10, 
  #                n.burn = 500, 
  #                n.chains = 1)
  # expect_s3_class(out, "spPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Detection covariate only ------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5, 1.2, -0.4)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- X
colnames(occ.covs) <- c('int')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, 
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ 1
det.formula <- ~ det.cov.1 + det.cov.2

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs[3, ] <- NA
#   expect_error(spPGOcc(occ.formula = occ.formula, 
# 	               det.formula = det.formula, 
# 	               data = tmp.data, 
# 	               n.batch = 40, 
# 	               batch.length = batch.length, 
# 	               cov.model = "exponential", 
# 	               tuning = tuning.list, 
# 	               NNGP = TRUE,
# 	               verbose = FALSE, 
# 	               n.neighbors = 5, 
# 	               search.type = 'cb', 
# 	               n.report = 10, 
# 	               n.burn = 500, 
# 	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Covariates on both ------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.8)
p.occ <- length(beta)
alpha <- c(-0.5, 1.2, -0.4)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, 
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ occ.cov.1
det.formula <- ~ det.cov.1 + det.cov.2

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Interactions on both ---------------------------------------------------- 
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.8, 1.2)
p.occ <- length(beta)
alpha <- c(-0.5, 1.2, -0.4)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, 
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- ~ det.cov.1 * det.cov.2

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  X.p.0 <- cbind(X.p.0, X.p.0[, 2] * X.p.0[, 3])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Site covariate on detection ---------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.8)
p.occ <- length(beta)
alpha <- c(-0.5, 1.2, -0.4)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3], 
                 occ.cov.1 = X[, 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, 
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ 1 
det.formula <- ~ occ.cov.1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs[3, ] <- NA
#   expect_error(spPGOcc(occ.formula = occ.formula, 
# 	               det.formula = det.formula, 
# 	               data = tmp.data, 
# 	               n.batch = 40, 
# 	               batch.length = batch.length, 
# 	               cov.model = "exponential", 
# 	               tuning = tuning.list, 
# 	               NNGP = TRUE,
# 	               verbose = FALSE, 
# 	               n.neighbors = 5, 
# 	               search.type = 'cb', 
# 	               n.report = 10, 
# 	               n.burn = 500, 
# 	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[3]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  pred.out <- predict(out, X.0[, 1, drop = FALSE], coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(1, dat$X[, 2])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Random intercept on occurrence ------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(25), 
	       sigma.sq.psi = c(1.3))
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 10
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re) 
colnames(occ.covs) <- c('int', 'occ.factor.1')
data.list <- list(y = y, coords = coords, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ (1 | occ.factor.1)
det.formula <- ~ 1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spPGOcc(occ.formula = occ.formula, 
  #                      det.formula = det.formula, 
  #                      data = tmp.data, 
  #                      n.batch = 40, 
  #                      batch.length = batch.length, 
  #                      cov.model = "exponential", 
  #                      tuning = tuning.list, 
  #                      NNGP = TRUE,
  #                      verbose = FALSE, 
  #                      n.neighbors = 5, 
  #                      search.type = 'cb', 
  #                      n.report = 10, 
  #                      n.burn = 500, 
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- spPGOcc(occ.formula = occ.formula, 
  #                det.formula = det.formula, 
  #                data = tmp.data, 
  #                n.batch = 40, 
  #                batch.length = batch.length, 
  #                cov.model = "exponential", 
  #                tuning = tuning.list, 
  #                NNGP = TRUE,
  #                verbose = FALSE, 
  #                n.neighbors = 5, 
  #                search.type = 'cb', 
  #                n.report = 10, 
  #                n.burn = 500, 
  #                n.chains = 1)
  # expect_s3_class(out, "spPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- matrix(dat$X.p[, 1, ], ncol = 1)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Multiple random intercepts on occurrence --------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(10, 12), 
	       sigma.sq.psi = c(0.2, 3.5))
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re) 
colnames(occ.covs) <- c('int', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, coords = coords, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spPGOcc(occ.formula = occ.formula, 
  #                      det.formula = det.formula, 
  #                      data = tmp.data, 
  #                      n.batch = 40, 
  #                      batch.length = batch.length, 
  #                      cov.model = "exponential", 
  #                      tuning = tuning.list, 
  #                      NNGP = TRUE,
  #                      verbose = FALSE, 
  #                      n.neighbors = 5, 
  #                      search.type = 'cb', 
  #                      n.report = 10, 
  #                      n.burn = 500, 
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- spPGOcc(occ.formula = occ.formula, 
  #                det.formula = det.formula, 
  #                data = tmp.data, 
  #                n.batch = 40, 
  #                batch.length = batch.length, 
  #                cov.model = "exponential", 
  #                tuning = tuning.list, 
  #                NNGP = TRUE,
  #                verbose = FALSE, 
  #                n.neighbors = 5, 
  #                search.type = 'cb', 
  #                n.report = 10, 
  #                n.burn = 500, 
  #                n.chains = 1)
  # expect_s3_class(out, "spPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0, X.re.0)
  colnames(X.0.full) <- c('int', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- matrix(dat$X.p[, 1, ], ncol = 1)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Occurrence REs + covariates ---------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(20, 15), 
	       sigma.sq.psi = c(0.2, 3.5))
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re) 
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, coords = coords, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spPGOcc(occ.formula = occ.formula, 
  #                      det.formula = det.formula, 
  #                      data = tmp.data, 
  #                      n.batch = 40, 
  #                      batch.length = batch.length, 
  #                      cov.model = "exponential", 
  #                      tuning = tuning.list, 
  #                      NNGP = TRUE,
  #                      verbose = FALSE, 
  #                      n.neighbors = 5, 
  #                      search.type = 'cb', 
  #                      n.report = 10, 
  #                      n.burn = 500, 
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- spPGOcc(occ.formula = occ.formula, 
  #                det.formula = det.formula, 
  #                data = tmp.data, 
  #                n.batch = 40, 
  #                batch.length = batch.length, 
  #                cov.model = "exponential", 
  #                tuning = tuning.list, 
  #                NNGP = TRUE,
  #                verbose = FALSE, 
  #                n.neighbors = 5, 
  #                search.type = 'cb', 
  #                n.report = 10, 
  #                n.burn = 500, 
  #                n.chains = 1)
  # expect_s3_class(out, "spPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0, X.re.0)
  colnames(X.0.full) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- matrix(dat$X.p[, 1, ], ncol = 1)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Occurrence REs + covariates in everything -------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
alpha <- c(-0.5, 1.2, -0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(20, 15), 
	       sigma.sq.psi = c(0.2, 3.5))
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re) 
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ det.cov.1 + det.cov.2

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0, X.re.0)
  colnames(X.0.full) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Random intercept on detection -------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(50), 
	     sigma.sq.p = c(2.50))
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X) 
colnames(occ.covs) <- c('int')
det.covs <- list(det.factor.1 = X.p.re[, , 1])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ 1 
det.formula <- ~ (1 | det.factor.1)

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0)
  colnames(X.0.full) <- c('int')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1])
  colnames(X.p.0) <- c('intercept', 'det.factor.1')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Multiple random intercepts on detection ---------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(50, 25, 30), 
	     sigma.sq.p = c(2.50, 1.2, 0.3))
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X) 
colnames(occ.covs) <- c('int')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.factor.3 = X.p.re[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ 1 
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3)

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs[3, ] <- NA
#   expect_error(spPGOcc(occ.formula = occ.formula, 
# 	               det.formula = det.formula, 
# 	               data = tmp.data, 
# 	               n.batch = 40, 
# 	               batch.length = batch.length, 
# 	               cov.model = "exponential", 
# 	               tuning = tuning.list, 
# 	               NNGP = TRUE,
# 	               verbose = FALSE, 
# 	               n.neighbors = 5, 
# 	               search.type = 'cb', 
# 	               n.report = 10, 
# 	               n.burn = 500, 
# 	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0)
  colnames(X.0.full) <- c('int')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.factor.1', 'det.factor.2', 'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Multiple random intercepts with covariates on detection -----------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5, 0.3)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(50, 25, 30), 
	     sigma.sq.p = c(2.50, 1.2, 0.3))
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X) 
colnames(occ.covs) <- c('int')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.factor.3 = X.p.re[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ 1 
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3) + det.cov.1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs[3, ] <- NA
#   expect_error(spPGOcc(occ.formula = occ.formula, 
# 	               det.formula = det.formula, 
# 	               data = tmp.data, 
# 	               n.batch = 40, 
# 	               batch.length = batch.length, 
# 	               cov.model = "exponential", 
# 	               tuning = tuning.list, 
# 	               NNGP = TRUE,
# 	               verbose = FALSE, 
# 	               n.neighbors = 5, 
# 	               search.type = 'cb', 
# 	               n.report = 10, 
# 	               n.burn = 500, 
# 	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0)
  colnames(X.0.full) <- c('int')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2', 'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Multiple random intercepts with covariates on both ----------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.4)
p.occ <- length(beta)
alpha <- c(-0.5, 0.3)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(50, 25, 30), 
	     sigma.sq.p = c(2.50, 1.2, 0.3))
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X) 
colnames(occ.covs) <- c('int', 'occ.cov.1')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.factor.3 = X.p.re[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ occ.cov.1
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3) + det.cov.1

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0)
  colnames(X.0.full) <- c('int', 'occ.cov.1')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2', 'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Random intercepts on both -----------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(20), 
               sigma.sq.psi = c(2.5))
p.RE <- list(levels = c(50, 25, 30), 
	     sigma.sq.p = c(2.50, 1.2, 0.3))
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re) 
colnames(occ.covs) <- c('int', 'occ.factor.1')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.factor.3 = X.p.re[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ (1 | occ.factor.1) 
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3)

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0, X.re.0)
  colnames(X.0.full) <- c('int', 'occ.factor.1')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.factor.1', 'det.factor.2', 'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Random intercepts on both with covariates -------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
alpha <- c(-0.5, -1.2, 0.9)
p.det <- length(alpha)
psi.RE <- list(levels = c(20), 
               sigma.sq.psi = c(2.5))
p.RE <- list(levels = c(50, 25, 30), 
	     sigma.sq.p = c(2.50, 1.2, 0.3))
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 42
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 42
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re) 
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3],
		 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.factor.3 = X.p.re[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ (1 | occ.factor.1) + occ.cov.1 
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3) + det.cov.1 + det.cov.2

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
                              tuning = tuning.list, 
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 400, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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

  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_output(spPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
	               det.formula = det.formula, 
	               data = tmp.data, 
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula, 
                       det.formula = det.formula, 
                       data = tmp.data, 
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = tmp.data, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  X.0.full <- cbind(X.0, X.re.0)
  colnames(X.0.full) <- c('int', 'occ.cov.1', 'occ.factor.1')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.cov.2', 
		       'det.factor.1', 'det.factor.2', 'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
})

# Fixed sigma sq ----------------------------------------------------
test_that("spPGOcc works with fixed sigma.sq", {
  prior.list <- list(sigma.sq.ig = 'fixed')
  out <- spPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.batch = 40,
	         batch.length = batch.length,
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
  expect_s3_class(out, "spPGOcc")
  expect_equal(length(unique(out$theta.samples[, 1])), 1)
})

# Uniform sigma sq --------------------------------------------------------
test_that("spPGOcc works with uniform prior on sigma.sq", {
  prior.list <- list(sigma.sq.unif = c(0, 5), 
                     nu.unif = c(0.1, 4))
  tuning.list <- list(phi = 0.5, nu = 0.6, sigma.sq = 0.7)
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
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
  expect_s3_class(out, "spPGOcc")
  out <- spPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_s3_class(out, "spPGOcc")
})

# Third dimension of y != max(n.rep) --------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
n.rep.max <- 7
beta <- c(0.3, 0.8)
p.occ <- length(beta)
alpha <- c(-0.5, 1.2, -0.4)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
	      phi = phi, cov.model = 'matern', nu = nu, n.rep.max = n.rep.max)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1')
det.covs <- list(det.cov.1 = X.p[, , 2],
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs,
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 40
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ occ.cov.1
det.formula <- ~ det.cov.1 + det.cov.2

out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
	       n.burn = 500,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class spPGOcc", {
  expect_s3_class(out, "spPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
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
  out <- spPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
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

  out <- spPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
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

  out <- spPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
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
  expect_output(spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula,
	               det.formula = det.formula,
	               data = tmp.data,
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
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spPGOcc(occ.formula = occ.formula,
                       det.formula = det.formula,
                       data = tmp.data,
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
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- spPGOcc(occ.formula = occ.formula,
                 det.formula = det.formula,
                 data = tmp.data,
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

# Check waicOcc -----------------------
test_that("waicOCC works for spPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spPGOcc", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for spPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.rep.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.rep.max))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.rep.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.rep.max))
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 2)
})

