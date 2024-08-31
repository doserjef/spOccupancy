# Test lfMsPGOcc.R  -------------------------------------------------------

skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(20)
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

data.list <- list(y = y, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      fix = TRUE,
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ 1

out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  out.k.fold <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
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
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	           det.formula = det.formula, 
	           data = data.list, 
		   n.factors = 3,
	           n.samples = n.samples,
	           n.omp.threads = 1,
	           verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
               n.thin = 13,
               n.factors = 3, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
	                n.omp.threads = 1,
                  parallel.chains = TRUE,
			n.factors = 3,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
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


# Occurrence covariate only -----------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5, 0.3)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.3, 2.8)
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, factor.model = TRUE, n.factor = 3)
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
data.list <- list(y = y, occ.covs = occ.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ 1

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3, 
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
                  parallel.chains = TRUE,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
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

# Detection covariate only -----------------------------------------------
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
alpha.mean <- c(0, 0.5, 1.2)
tau.sq.alpha <- c(1, 2, 3)
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1 
det.formula <- ~ det.cov.1 + det.cov.2

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
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
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2.3)
# Detection
alpha.mean <- c(0, 0.5, 1.2)
tau.sq.alpha <- c(1, 2, 3)
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1
det.formula <- ~ det.cov.1 + det.cov.2

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
           parallel.chains = TRUE,
	         n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  colnames(X.0) <- c('int', 'occ.cov.1')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.5, -0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2.3, 1.2)
# Detection
alpha.mean <- c(0, 0.5, 1.2)
tau.sq.alpha <- c(1, 2, 3)
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- ~ det.cov.1 * det.cov.2

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
	         n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.cov.1:occ.cov.2')
  pred.out <- predict(out, X.0, coords.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, J))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  X.p.0 <- cbind(X.p.0, X.p.0[, 2] * X.p.0[, 3])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2.3)
# Detection
alpha.mean <- c(0, 0.5, 1.2)
tau.sq.alpha <- c(1, 2, 3)
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3], 
                 occ.cov.1 = X[, 2])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1
det.formula <- ~ occ.cov.1 

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       inits = inits.list, 
	       n.samples = n.samples, 
         parallel.chains = TRUE,
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  colnames(X.0) <- c('int', 'occ.cov.1')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

# Random intercept on occurrence ------------------------------------------
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
# det.covs <- list(det.cov.1 = X.p[, , 2], 
# 		 det.cov.2 = X.p[, , 3]) 
data.list <- list(y = y, occ.covs = occ.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1) 
det.formula <- ~ 1 

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
           parallel.chains = TRUE,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('int', 'occ.factor.1')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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
psi.RE <- list(levels = c(45, 15), 
               sigma.sq.psi = c(1.3, 0.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
# det.covs <- list(det.cov.1 = X.p[, , 2], 
# 		 det.cov.2 = X.p[, , 3]) 
data.list <- list(y = y, occ.covs = occ.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1 

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
                  parallel.chains = TRUE,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(45, 15), 
               sigma.sq.psi = c(1.3, 0.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
# det.covs <- list(det.cov.1 = X.p[, , 2], 
# 		 det.cov.2 = X.p[, , 3]) 
data.list <- list(y = y, occ.covs = occ.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1 

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.str <- cbind(X.0, X.re.0)
  colnames(X.str) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.str, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

# Occurrence REs + covariates in all --------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2)
# Detection
alpha.mean <- c(0, 1)
tau.sq.alpha <- c(1, 2.5)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(45, 15), 
               sigma.sq.psi = c(1.3, 0.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
det.covs <- list(det.cov.1 = X.p[, , 2]) 
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ det.cov.1 

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  # data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

# Random intercepts on detection ------------------------------------------
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
alpha.mean <- c(0, 1)
tau.sq.alpha <- c(1, 2.5)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(30), 
             sigma.sq.p = c(2))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1) 

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0)
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0, 1)
tau.sq.alpha <- c(1, 2.5)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(30, 45), 
             sigma.sq.p = c(2, 1.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2)

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

# Multiple random intercepts with covariate -------------------------------
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
alpha.mean <- c(0, 1)
tau.sq.alpha <- c(1, 2.5)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(30, 45), 
             sigma.sq.p = c(2, 1.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
                 det.factor.2 = X.p.re[, , 2], 
                 det.cov.1 = X.p[, , 2]) 
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  colnames(X.0) <- c('int')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

# Multiple random intercepts with covariate on both -----------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5)
# Detection
alpha.mean <- c(0, 1)
tau.sq.alpha <- c(1, 2.5)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list(levels = c(30, 45), 
             sigma.sq.p = c(2, 1.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
colnames(occ.covs) <- c('int', 'occ.cov.1')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.cov.1 = X.p[, , 2]) 
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0)
  colnames(X.0) <- c('int', 'occ.cov.1')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

# Multiple random intercepts on both --------------------------------------
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
p.RE <- list(levels = c(30, 45), 
             sigma.sq.p = c(2, 1.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1)
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2)

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.factor.1')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.2)
# Detection
alpha.mean <- c(0, -0.5)
tau.sq.alpha <- c(1, 1.5)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = c(20, 15), 
	       sigma.sq.psi = c(0.5, 2.4))
p.RE <- list(levels = c(30, 45, 30), 
             sigma.sq.p = c(2, 1.5, 0.5))
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
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
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.cov.1 = X.p[, , 2], 
                 det.factor.3 = X.p.re[, , 3]) 
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3) + det.cov.1

out <- lfMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.samples = n.samples,
		 n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula, 
	                det.formula = det.formula, 
	                data = data.list, 
	                n.samples = 100,
			n.factors = 3, 
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, coords.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, nrow(X.0)))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2', 
                       'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J))
})

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

# Third dimension of y != max(n.rep) --------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
n.rep.max <- max(n.rep)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 2.3)
# Detection
alpha.mean <- c(0, 0.5, 1.2)
tau.sq.alpha <- c(1, 2, 3)
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
	        psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, n.rep.max = n.rep.max)
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
det.covs <- list(det.cov.1 = X.p[, , 2],
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0,
		      beta.comm = 0,
		      beta = 0,
		      alpha = 0,
		      tau.sq.beta = 1,
		      tau.sq.alpha = 1,
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1
det.formula <- ~ det.cov.1 + det.cov.2

out <- lfMsPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
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
test_that("out is of class lfMsPGOcc", {
  expect_s3_class(out, "lfMsPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
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

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- lfMsPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.samples = n.samples,
	         n.factors = 3,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "lfMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(lfMsPGOcc(occ.formula = occ.formula,
	                det.formula = det.formula,
	                data = data.list,
	                n.samples = 100,
			n.factors = 3,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100,
	                n.burn = 1,
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for lfMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

test_that("fitted works for lfMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
  expect_equal(class(fitted.out$y.rep.samples), "array")
  expect_equal(class(fitted.out$p.samples), "array")
  expect_equal(dim(fitted.out$y.rep.samples), dim(fitted.out$p.samples))
})

test_that("predict works for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  colnames(X.0) <- c('int', 'occ.cov.1')
  pred.out <- predict(out, X.0, coords.0)
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

test_that("posterior predictive checks work for lfMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  J.fit <- nrow(X)
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

