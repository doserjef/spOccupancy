# Test lfJSDM.R  -------------------------------------------------------

skip_on_cran()

# Intercept only ----------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, factor.model = TRUE, n.factors = 3)
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

y <- apply(y, c(1, 2), max, na.rm = TRUE)
data.list <- list(y = y, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta.comm = 0, 
		   beta = 0, 
		   fix = TRUE,
	           tau.sq.beta = 1) 

n.samples <- 1000
n.report <- 100
formula <- ~ 1

out <- lfJSDM(formula = formula, 
	      data = data.list, 
	      inits = inits.list, 
	      n.samples = n.samples, 
	      n.factors = 3,
	      priors = prior.list, 
              n.omp.threads = 1, 
	      verbose = FALSE, 
	      n.report = n.report, 
	      n.burn = 400,
	      n.thin = 2, 
	      n.chains = 2,
	      k.fold = 2,
              k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class lfJSDM", {
  expect_s3_class(out, "lfJSDM")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  out.k.fold <- lfJSDM(formula = formula, 
	      data = data.list, 
	      inits = inits.list, 
	      n.samples = n.samples, 
	      n.factors = 3,
	      priors = prior.list, 
              n.omp.threads = 1, 
	      verbose = FALSE, 
	      n.report = n.report, 
	      n.burn = 400,
	      n.thin = 2, 
	      n.chains = 2,
	      k.fold = 2,
              k.fold.threads = 1, k.fold.only = TRUE)
  expect_equal(length(out.k.fold$k.fold.deviance), N)
  expect_type(out.k.fold$k.fold.deviance, "double")
  expect_equal(sum(out.k.fold$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfJSDM(formula = formula, 
	        data = data.list, 
		n.factors = 3,
	        n.samples = n.samples,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "lfJSDM")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- lfJSDM(formula = ~ 1, 
	       data = data.list, 
               n.thin = 13,
               n.factors = 3,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
})

test_that("verbose prints to the screen", {
  expect_output(lfJSDM(formula = formula, 
	               data = data.list, 
	               n.samples = 100,
	               n.omp.threads = 1,
		       n.factors = 3,
	               verbose = TRUE,
	               n.report = 100, 
	               n.burn = 1, 
	               n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfJSDM", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfJSDM", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$z.samples), "array")
  expect_equal(class(fitted.out$psi.samples), "array")
  expect_equal(dim(fitted.out$z.samples), dim(fitted.out$psi.samples))
})

test_that("predict works for lfJSDM", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})

test_that("posterior predictive checks work for lfJSDM", {
  expect_error(ppcOcc(out, 'chi-square', 2))
})

# Covariates --------------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.2, -0.4)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.3, 2.5)
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, factor.model = TRUE, n.factors = 3)
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

y <- apply(y, c(1, 2), max, na.rm = TRUE)
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, coords = coords, covs = occ.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta.comm = 0, 
		   beta = 0, 
	           tau.sq.beta = 1) 

n.samples <- 1000
n.report <- 100
formula <- ~ occ.cov.1 + occ.cov.2

out <- lfJSDM(formula = formula, 
	      data = data.list, 
	      inits = inits.list, 
	      n.samples = n.samples, 
	      n.factors = 3,
	      priors = prior.list, 
              n.omp.threads = 1, 
	      verbose = FALSE, 
	      n.report = n.report, 
	      n.burn = 400,
	      n.thin = 2, 
	      n.chains = 2,
	      k.fold = 2,
              k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class lfJSDM", {
  expect_s3_class(out, "lfJSDM")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfJSDM(formula = formula, 
	        data = data.list, 
		n.factors = 3,
	        n.samples = n.samples,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "lfJSDM")
})

test_that("verbose prints to the screen", {
  expect_output(lfJSDM(formula = formula, 
	               data = data.list, 
	               n.samples = 100,
	               n.omp.threads = 1,
		       n.factors = 3,
	               verbose = TRUE,
	               n.report = 100, 
	               n.burn = 1, 
	               n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfJSDM", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfJSDM", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$z.samples), "array")
  expect_equal(class(fitted.out$psi.samples), "array")
  expect_equal(dim(fitted.out$z.samples), dim(fitted.out$psi.samples))
})

test_that("predict works for lfJSDM", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2')
  pred.out <- predict(out, X.0, coords.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})

test_that("posterior predictive checks work for lfJSDM", {
  expect_error(ppcOcc(out, 'chi-square', 2))
})

# Random intercept --------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
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
	       sigma.sq.psi = c(0.5))
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, factor.model = TRUE, n.factors = 3)
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

y <- apply(y, c(1, 2), max, na.rm = TRUE)
occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1')
data.list <- list(y = y, coords = coords, covs = occ.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta.comm = 0, 
		   beta = 0, 
	           tau.sq.beta = 1) 

n.samples <- 1000
n.report <- 100
formula <- ~ (1 | occ.factor.1)

out <- lfJSDM(formula = formula, 
	      data = data.list, 
	      inits = inits.list, 
	      n.samples = n.samples, 
	      n.factors = 3,
	      priors = prior.list, 
              n.omp.threads = 1, 
	      verbose = FALSE, 
	      n.report = n.report, 
	      n.burn = 400,
	      n.thin = 2, 
	      n.chains = 2,
	      k.fold = 2,
              k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class lfJSDM", {
  expect_s3_class(out, "lfJSDM")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs <- as.data.frame(data.list$covs)
  data.list$covs$occ.factor.1 <- factor(data.list$covs$occ.factor.1)
  expect_error(out <- lfJSDM(formula = formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
			    n.factors = 3,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  data.list$covs$occ.factor.1 <- as.character(factor(data.list$covs$occ.factor.1))
  expect_error(out <- lfJSDM(formula = formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
			    n.factors = 3,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfJSDM(formula = formula, 
	        data = data.list, 
		n.factors = 3,
	        n.samples = n.samples,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "lfJSDM")
})

test_that("verbose prints to the screen", {
  expect_output(lfJSDM(formula = formula, 
	               data = data.list, 
	               n.samples = 100,
	               n.omp.threads = 1,
		       n.factors = 3,
	               verbose = TRUE,
	               n.report = 100, 
	               n.burn = 1, 
	               n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfJSDM", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfJSDM", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$z.samples), "array")
  expect_equal(class(fitted.out$psi.samples), "array")
  expect_equal(dim(fitted.out$z.samples), dim(fitted.out$psi.samples))
})

test_that("predict works for lfJSDM", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1')
  pred.out <- predict(out, X.0, coords.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})

test_that("posterior predictive checks work for lfJSDM", {
  expect_error(ppcOcc(out, 'chi-square', 2))
})

# Multiple random intercepts ----------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
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
psi.RE <- list(levels = c(20, 30), 
	       sigma.sq.psi = c(0.5, 1.5))
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, factor.model = TRUE, n.factors = 3)
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

y <- apply(y, c(1, 2), max, na.rm = TRUE)
occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, coords = coords, covs = occ.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta.comm = 0, 
		   beta = 0, 
	           tau.sq.beta = 1) 

n.samples <- 1000
n.report <- 100
formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2)

out <- lfJSDM(formula = formula, 
	      data = data.list, 
	      inits = inits.list, 
	      n.samples = n.samples, 
	      n.factors = 3,
	      priors = prior.list, 
              n.omp.threads = 1, 
	      verbose = FALSE, 
	      n.report = n.report, 
	      n.burn = 400,
	      n.thin = 2, 
	      n.chains = 2,
	      k.fold = 2,
              k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class lfJSDM", {
  expect_s3_class(out, "lfJSDM")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs <- as.data.frame(data.list$covs)
  data.list$covs$occ.factor.1 <- factor(data.list$covs$occ.factor.1)
  expect_error(out <- lfJSDM(formula = formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
			    n.factors = 3,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  data.list$covs$occ.factor.1 <- as.character(factor(data.list$covs$occ.factor.1))
  expect_error(out <- lfJSDM(formula = formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
			    n.factors = 3,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfJSDM(formula = formula, 
	        data = data.list, 
		n.factors = 3,
	        n.samples = n.samples,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "lfJSDM")
})

test_that("verbose prints to the screen", {
  expect_output(lfJSDM(formula = formula, 
	               data = data.list, 
	               n.samples = 100,
	               n.omp.threads = 1,
		       n.factors = 3,
	               verbose = TRUE,
	               n.report = 100, 
	               n.burn = 1, 
	               n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfJSDM", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfJSDM", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$z.samples), "array")
  expect_equal(class(fitted.out$psi.samples), "array")
  expect_equal(dim(fitted.out$z.samples), dim(fitted.out$psi.samples))
})

test_that("predict works for lfJSDM", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})

test_that("posterior predictive checks work for lfJSDM", {
  expect_error(ppcOcc(out, 'chi-square', 2))
})

# Random intercepts + covariates ------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -1)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2.3)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(20, 30), 
	       sigma.sq.psi = c(0.5, 1.5))
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, factor.model = TRUE, n.factors = 3)
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

y <- apply(y, c(1, 2), max, na.rm = TRUE)
occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, coords = coords, covs = occ.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta.comm = 0, 
		   beta = 0, 
	           tau.sq.beta = 1) 

n.samples <- 1000
n.report <- 100
formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)

out <- lfJSDM(formula = formula, 
	      data = data.list, 
	      inits = inits.list, 
	      n.samples = n.samples, 
	      n.factors = 3,
	      priors = prior.list, 
              n.omp.threads = 1, 
	      verbose = FALSE, 
	      n.report = n.report, 
	      n.burn = 400,
	      n.thin = 2, 
	      n.chains = 1,
	      k.fold = 2,
              k.fold.threads = 1)

# To make sure it worked --------------
test_that("out is of class lfJSDM", {
  expect_s3_class(out, "lfJSDM")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs <- as.data.frame(data.list$covs)
  data.list$covs$occ.factor.1 <- factor(data.list$covs$occ.factor.1)
  expect_error(out <- lfJSDM(formula = formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
			    n.factors = 3,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  data.list$covs$occ.factor.1 <- as.character(factor(data.list$covs$occ.factor.1))
  expect_error(out <- lfJSDM(formula = formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
			    n.factors = 3,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfJSDM(formula = formula, 
	        data = data.list, 
		n.factors = 3,
	        n.samples = n.samples,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "lfJSDM")
})

test_that("verbose prints to the screen", {
  expect_output(lfJSDM(formula = formula, 
	               data = data.list, 
	               n.samples = 100,
	               n.omp.threads = 1,
		       n.factors = 3,
	               verbose = TRUE,
	               n.report = 100, 
	               n.burn = 1, 
	               n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfJSDM", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfJSDM", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$z.samples), "array")
  expect_equal(class(fitted.out$psi.samples), "array")
  expect_equal(dim(fitted.out$z.samples), dim(fitted.out$psi.samples))
})

test_that("predict works for lfJSDM", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})

test_that("posterior predictive checks work for lfJSDM", {
  expect_error(ppcOcc(out, 'chi-square', 2))
})

