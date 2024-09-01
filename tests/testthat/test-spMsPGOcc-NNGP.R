# Test spMsPGOcc.R  -------------------------------------------------------
# NNGP --------------------------------------------------------------------

skip_on_cran()

# Intercept Only ----------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
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
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      fix = TRUE,
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ 1

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  out.k.fold <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
		             k.fold.only = TRUE,
                 k.fold.threads = 1)
  expect_equal(length(out.k.fold$k.fold.deviance), N)
  expect_type(out.k.fold$k.fold.deviance, "double")
  expect_equal(sum(out.k.fold$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
               n.thin = 13,
               n.batch = n.batch, 
               batch.length = batch.length, 
               accept.rate = 0.43, 
	       n.omp.threads = 1,
	       verbose = FALSE))
})
# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Occurrence coviarate only -----------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2.3)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
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
data.list <- list(y = y, coords = coords, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1
det.formula <- ~ 1

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
	               n.burn = 100,
	               n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
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
  #                      n.burn = 100,
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[, 1, 1] <- NA
  # out <- spMsPGOcc(occ.formula = occ.formula,
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
  #                n.burn = 100,
  #                n.chains = 1)
  # expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Detection covariate only ------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0, -0.5)
tau.sq.alpha <- c(1, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
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
det.covs <- list(det.cov.1 = X.p[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ 1 
det.formula <- ~ det.cov.1

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs[3, ] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
  #                        det.formula = det.formula,
  #                        data = tmp.data,
  #                        n.batch = 40,
  #                        batch.length = batch.length,
  #                        cov.model = "exponential",
  #                        tuning = tuning.list,
  #                        NNGP = TRUE,
  #                        verbose = FALSE,
  #                        n.neighbors = 5,
  #                        search.type = 'cb',
  #                        n.report = 10,
  #                        n.burn = 100,
  #                        n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Covariates on both ------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5, 1.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 2.3)
# Detection
alpha.mean <- c(0, -0.5)
tau.sq.alpha <- c(1, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
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
det.covs <- list(det.cov.1 = X.p[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ det.cov.1

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Interactions on both ----------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5, 1.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 2.3)
# Detection
alpha.mean <- c(0, -0.5, 1.2)
tau.sq.alpha <- c(1, 2.3, 1.5)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
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
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0,
		      beta.comm = 0,
		      beta = 0,
		      alpha = 0,
		      tau.sq.beta = 1,
		      tau.sq.alpha = 1,
		      z = apply(y, c(1, 2), max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- ~ det.cov.1 * det.cov.2

out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
		 n.thin = 2,
		 n.chains = 2,
                 k.fold = 2,
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   cov.model = "exponential",
                   tuning = tuning.list,
                   NNGP = TRUE,
		   verbose = FALSE,
                   n.neighbors = 5,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   inits = inits.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   priors = prior.list,
                   cov.model = "gaussian",
                   tuning = list(phi = 0.3),
                   n.omp.threads = 1,
                   verbose = FALSE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   search.type = 'cb',
                   n.report = 10,
                   n.burn = 100,
                   n.thin = 1,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   inits = inits.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   priors = prior.list,
                   cov.model = "spherical",
                   tuning = list(phi = 0.3),
                   n.omp.threads = 1,
                   verbose = FALSE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   search.type = 'cb',
                   n.report = 10,
                   n.burn = 100,
                   n.thin = 1,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   inits = inits.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   priors = list(nu.unif = list(0.4, 3)),
                   cov.model = "matern",
                   tuning = list(phi = 0.3, nu = 0.2),
                   n.omp.threads = 1,
                   verbose = FALSE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   search.type = 'cb',
                   n.report = 10,
                   n.burn = 100,
                   n.thin = 1,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula,
                 det.formula = det.formula,
                 data = data.list,
                 inits = inits.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 accept.rate = 0.43,
                 priors = prior.list,
                 cov.model = "exponential",
                 tuning = tuning.list,
                 n.omp.threads = 1,
                 verbose = TRUE,
                 NNGP = TRUE,
                 n.neighbors = 5,
                 search.type = 'cb',
                 n.report = 10,
                 n.burn = 100,
                 n.thin = 1,
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.cov.1:occ.cov.2')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  X.p.0 <- cbind(X.p.0, X.p.0[, 2] * X.p.0[, 3])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Site covariate on detection ---------------------------------------------
set.seed(400)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5, 1.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 2.3)
# Detection
alpha.mean <- c(0, -0.5)
tau.sq.alpha <- c(1, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
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
                 occ.cov.1 = X[, 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ occ.cov.1

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs[3, ] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
  #                        det.formula = det.formula,
  #                        data = tmp.data,
  #                        n.batch = 40,
  #                        batch.length = batch.length,
  #                        cov.model = "exponential",
  #                        tuning = tuning.list,
  #                        NNGP = TRUE,
  #                        verbose = FALSE,
  #                        n.neighbors = 5,
  #                        search.type = 'cb',
  #                        n.report = 10,
  #                        n.burn = 100,
  #                        n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[2]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- X.0[, 1, drop = FALSE]
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(1, data.list$det.covs$occ.cov.1)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, nrow(X.p.0)))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Random intercept on occurrence -------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(45), 
               sigma.sq.psi = c(1.3))
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
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
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1)
det.formula <- ~ 1 

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
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
  #                      n.burn = 100,
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[, 1, 1] <- NA
  # out <- spMsPGOcc(occ.formula = occ.formula,
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
  #                n.burn = 100,
  #                n.chains = 1)
  # expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Multiple random intercepts on occurrence --------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(45, 20), 
               sigma.sq.psi = c(1.3, 3.4))
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, coords = coords, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1 

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
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
  #                      n.burn = 100,
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[, 1, 1] <- NA
  # out <- spMsPGOcc(occ.formula = occ.formula,
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
  #                n.burn = 100,
  #                n.chains = 1)
  # expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Occurrence REs + covariates ---------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.2, -1.9)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 2.3)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(20, 10), 
               sigma.sq.psi = c(1.3, 3.4))
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, coords = coords, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1 

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
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
  #                      n.burn = 100,
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$y[, 1, 1] <- NA
  # out <- spMsPGOcc(occ.formula = occ.formula,
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
  #                n.burn = 100,
  #                n.chains = 1)
  # expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Occurrence REs + covariates in everything -------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.2, -1.9)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 2.3)
# Detection
alpha.mean <- c(0, 1.2)
tau.sq.alpha <- c(1, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(20, 10), 
               sigma.sq.psi = c(1.3, 3.4))
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.factor.1', 'occ.factor.2')
det.covs <- list(det.cov.1 = X.p[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ det.cov.1

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  # data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Random intercept on detection -------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(50), 
	     sigma.sq.p = c(2.50))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
det.covs <- list(det.factor.1 = X.p.re[, , 1])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ 1 
det.formula <- ~ (1 | det.factor.1) 

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs[3, ] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
  #                        det.formula = det.formula,
  #                        data = tmp.data,
  #                        n.batch = 40,
  #                        batch.length = batch.length,
  #                        cov.model = "exponential",
  #                        tuning = tuning.list,
  #                        NNGP = TRUE,
  #                        verbose = FALSE,
  #                        n.neighbors = 5,
  #                        search.type = 'cb',
  #                        n.report = 10,
  #                        n.burn = 100,
  #                        n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0)
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, 1], dat$X.p.re[, 1, 1])
  colnames(X.p.0) <- c('intercept', 'det.factor.1')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Multiple random intercepts on detection ---------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(50, 10), 
	     sigma.sq.p = c(2.50, 1.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ 1 
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2)

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs[3, ] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
  #                        det.formula = det.formula,
  #                        data = tmp.data,
  #                        n.batch = 40,
  #                        batch.length = batch.length,
  #                        cov.model = "exponential",
  #                        tuning = tuning.list,
  #                        NNGP = TRUE,
  #                        verbose = FALSE,
  #                        n.neighbors = 5,
  #                        search.type = 'cb',
  #                        n.report = 10,
  #                        n.burn = 100,
  #                        n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0)
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, 1], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Detection random effects with covariates --------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0, 0.5)
tau.sq.alpha <- c(1, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(50, 10), 
	     sigma.sq.p = c(2.50, 1.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ 1 
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs[3, ] <- NA
  # expect_error(spMsPGOcc(occ.formula = occ.formula,
  #                        det.formula = det.formula,
  #                        data = tmp.data,
  #                        n.batch = 40,
  #                        batch.length = batch.length,
  #                        cov.model = "exponential",
  #                        tuning = tuning.list,
  #                        NNGP = TRUE,
  #                        verbose = FALSE,
  #                        n.neighbors = 5,
  #                        search.type = 'cb',
  #                        n.report = 10,
  #                        n.burn = 100,
  #                        n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0)
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Detection random effects with covariates on all -------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.5, 0.3)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.5, 3.3)
# Detection
alpha.mean <- c(0, 0.5)
tau.sq.alpha <- c(1, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(50, 10), 
	     sigma.sq.p = c(2.50, 1.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
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

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})


# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Random intercepts on both -----------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(20), 
               sigma.sq.psi = c(2.5))
p.RE <- list(levels = c(50, 10), 
	     sigma.sq.p = c(2.50, 1.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1)
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2)

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
	       sort(unique(c(X.re))))
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Random intercepts on both with covariates -------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.2, -0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2.4, 0.3)
# Detection
alpha.mean <- c(0, 0.5)
tau.sq.alpha <- c(1, 3.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(20), 
               sigma.sq.psi = c(2.5))
p.RE <- list(levels = c(50, 10), 
	     sigma.sq.p = c(2.50, 1.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
		phi = phi, nu = nu, cov.model = 'matern')

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.factor.1')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2 + (1 | occ.factor.1)
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

out <- spMsPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
		 n.thin = 2,
		 n.chains = 2, 
                 k.fold = 2, 
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spMsPGOcc(occ.formula = occ.formula,
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
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
	       sort(unique(c(X.re))))
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   NNGP = TRUE,
		   verbose = FALSE, 
                   n.neighbors = 5, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "gaussian", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "spherical", 
                   tuning = list(phi = 0.3), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula, 
                   det.formula = det.formula, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = list(nu.unif = list(0.4, 3)), 
                   cov.model = "matern", 
                   tuning = list(phi = 0.3, nu = 0.2), 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 100, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 100, 
                 n.thin = 1, 
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.factor.1')
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})

# Model with fixed sigma.sq -----------------------------------------------
test_that("spMsPGOcc works with fixed sigma.sq", {
  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   inits = inits.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   priors = list(sigma.sq.ig = "fixed"),
                   cov.model = "exponential",
                   tuning = list(phi = 0.3, nu = 0.2),
                   n.omp.threads = 1,
                   verbose = FALSE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   search.type = 'cb',
                   n.report = 10,
                   n.burn = 100,
                   n.thin = 1,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
  expect_equal(length(unique(out$theta.samples[, 1])), 1)
})

# Uniform sigma sq --------------------------------------------------------
test_that("spMsPGOcc works with uniform prior on sigma.sq", {
  prior.list <- list(sigma.sq.unif = list(a = 0, b = 5),
                     nu.unif = list(a = 0.1, b = 4))
  tuning.list <- list(phi = 0.5, nu = 0.6, sigma.sq = 0.7)
  out <- spMsPGOcc(occ.formula = occ.formula,
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
	         n.burn = 100,
	         n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
  out <- spMsPGOcc(occ.formula = occ.formula,
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
	         n.burn = 100,
	         n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Third dimension of y != max(n.rep)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
n.rep.max <- 7
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5, 1.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 2.3)
# Detection
alpha.mean <- c(0, -0.5)
tau.sq.alpha <- c(1, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
alpha.true <- alpha
phi <- rep(3 / .7, N)
sigma.sq <- rep(2, N)
nu <- rep(2, N)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
		phi = phi, nu = nu, cov.model = 'matern', n.rep.max = n.rep.max)

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
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
det.covs <- list(det.cov.1 = X.p[, , 2])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha.comm = 0,
		      beta.comm = 0,
		      beta = 0,
		      alpha = 0,
		      tau.sq.beta = 1,
		      tau.sq.alpha = 1,
		      z = apply(y, c(1, 2), max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 10
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ det.cov.1

out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
		 n.thin = 2,
		 n.chains = 2,
                 k.fold = 2,
                 k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                         n.burn = 100,
                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- 1
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(spMsPGOcc(occ.formula = occ.formula,
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
                       n.burn = 100,
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[, 1, 1] <- NA
  out <- spMsPGOcc(occ.formula = occ.formula,
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
                 n.burn = 100,
                 n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

# Check default values ----------------
test_that("default priors, inits, burn, thin work", {
  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   cov.model = "exponential",
                   tuning = tuning.list,
                   NNGP = TRUE,
		   verbose = FALSE,
                   n.neighbors = 5,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   inits = inits.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   priors = prior.list,
                   cov.model = "gaussian",
                   tuning = list(phi = 0.3),
                   n.omp.threads = 1,
                   verbose = FALSE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   search.type = 'cb',
                   n.report = 10,
                   n.burn = 100,
                   n.thin = 1,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   inits = inits.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   priors = prior.list,
                   cov.model = "spherical",
                   tuning = list(phi = 0.3),
                   n.omp.threads = 1,
                   verbose = FALSE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   search.type = 'cb',
                   n.report = 10,
                   n.burn = 100,
                   n.thin = 1,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  out <- spMsPGOcc(occ.formula = occ.formula,
                   det.formula = det.formula,
                   data = data.list,
                   inits = inits.list,
                   n.batch = n.batch,
                   batch.length = batch.length,
                   accept.rate = 0.43,
                   priors = list(nu.unif = list(0.4, 3)),
                   cov.model = "matern",
                   tuning = list(phi = 0.3, nu = 0.2),
                   n.omp.threads = 1,
                   verbose = FALSE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   search.type = 'cb',
                   n.report = 10,
                   n.burn = 100,
                   n.thin = 1,
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spMsPGOcc(occ.formula = occ.formula,
                 det.formula = det.formula,
                 data = data.list,
                 inits = inits.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 accept.rate = 0.43,
                 priors = prior.list,
                 cov.model = "exponential",
                 tuning = tuning.list,
                 n.omp.threads = 1,
                 verbose = TRUE,
                 NNGP = TRUE,
                 n.neighbors = 5,
                 search.type = 'cb',
                 n.report = 10,
                 n.burn = 100,
                 n.thin = 1,
		 n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for spMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  J.fit <- nrow(X)
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.rep.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.rep.max))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J.fit))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J.fit))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.rep.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.rep.max))
})

