# Test spMsPGOcc.R  -------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()
set.seed(101)
J.x <- 7
J.y <- 7
J <- J.x * J.y
n.rep <- sample(2:4, size = J, replace = TRUE)
N <- 5
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.15)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.3)
# Detection
alpha.mean <- c(0.5, 0.2, -.2)
tau.sq.alpha <- c(0.2, 0.3, 0.8)
p.det <- length(alpha.mean)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
phi <- runif(N, 3/1, 3/.4)
sigma.sq <- runif(N, 0.3, 3)
sp <- TRUE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                phi = phi, sigma.sq = sigma.sq, sp = TRUE, cov.model = 'exponential')

# Number of batches
n.batch <- 40
# Batch length
batch.length <- 25
n.samples <- n.batch * batch.length

y <- dat$y
X <- dat$X
X.p <- dat$X.p
coords <- as.matrix(dat$coords)

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2], 
                 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                   alpha.comm.normal = list(mean = 0, var = 2.72), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3/1, b = 3/.1), 
                   sigma.sq.ig = list(a = 2, b = 2)) 
# Initial values
inits.list <- list(alpha.comm = 0, 
                   beta.comm = 0, 
                   beta = 0, 
                   alpha = 0,
                   tau.sq.beta = 1, 
                   tau.sq.alpha = 1, 
                   phi = 3 / .5, 
                   sigma.sq = 2,
                   w = matrix(0, nrow = N, ncol = nrow(X)),
                   z = apply(y, c(1, 2), max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1) 

out <- spMsPGOcc(occ.formula = ~ occ.cov, 
                 det.formula = ~ det.cov.1 + det.cov.2, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = FALSE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 500, 
                 n.thin = 1, 
		 n.chains = 1)
n.post.samples <- length(seq(from = out$n.burn + 1, 
			     to = n.samples, 
			     by = as.integer(out$n.thin))) * out$n.chains

test_that("out is of class spMsPGOcc", {
  expect_s3_class(out, "spMsPGOcc")
})

test_that("samples are the right size", {
  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = out$n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_equal(dim(out$beta.comm.samples), c(n.post.samples, length(beta.mean)))
  expect_equal(dim(out$alpha.comm.samples), c(n.post.samples, length(alpha.mean)))
  expect_equal(dim(out$beta.samples), c(n.post.samples, length(beta)))
  expect_equal(dim(out$alpha.samples), c(n.post.samples, length(alpha)))
  expect_equal(dim(out$z.samples), c(n.post.samples, N, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, N, J))
})

test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
})

test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

test_that("default priors and inits work", {

  set.seed(1010)
  out <- spMsPGOcc(occ.formula = ~ occ.cov, 
                   det.formula = ~ det.cov.1 + det.cov.2, 
                   data = data.list,
                   #inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   #priors = prior.list, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   n.omp.threads = 1, 
                   verbose = TRUE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 500, 
                   n.thin = 1, 
		   n.chains = 2)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("all correlation functions work", {
  set.seed(555)
  out <- spMsPGOcc(occ.formula = ~ occ.cov, 
                   det.formula = ~ det.cov.1 + det.cov.2, 
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
                   n.burn = 500, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")

  set.seed(557)
  out <- spMsPGOcc(occ.formula = ~ occ.cov, 
                   det.formula = ~ det.cov.1 + det.cov.2, 
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
                   n.burn = 500, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
  set.seed(556)
  out <- spMsPGOcc(occ.formula = ~ occ.cov, 
                   det.formula = ~ det.cov.1 + det.cov.2, 
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
                   n.burn = 500, 
                   n.thin = 1, 
		   n.chains = 1)
  expect_s3_class(out, "spMsPGOcc")
})

test_that("verbose prints to the screen", {

  set.seed(1010)		  
  expect_output(spMsPGOcc(occ.formula = ~ occ.cov, 
                 det.formula = ~ det.cov.1 + det.cov.2, 
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
                 n.burn = 500, 
                 n.thin = 1, 
		 n.chains = 2))
})

test_that("cross-validation works", {
  set.seed(117) 
  out <- spMsPGOcc(occ.formula = ~ occ.cov, 
                   det.formula = ~ det.cov.1 + det.cov.2, 
                   data = data.list,
                   inits = inits.list, 
                   n.batch = n.batch, 
                   batch.length = batch.length, 
                   accept.rate = 0.43, 
                   priors = prior.list, 
                   cov.model = "exponential", 
                   tuning = tuning.list, 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   NNGP = TRUE, 
                   n.neighbors = 5, 
                   search.type = 'cb', 
                   n.report = 10, 
                   n.burn = 500, 
                   n.thin = 1, 
		   n.chains = 2, 
		   k.fold = 2)

  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

test_that("random effects on detection work", {
  set.seed(500)
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
  p.RE <- list(levels = c(45), 
  	     sigma.sq.p = c(1.2))
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
  
  dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
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
  # Tuning
  tuning.list <- list(phi = 1, nu = 1)
  
  batch.length <- 25
  n.batch <- 40
  n.report <- 100


  out <- spMsPGOcc(occ.formula = ~ 1,
	         det.formula = ~ (1 | det.factor),
	         data = data.list,
	         inits = inits.list,
		 batch.length = batch.length,
		 n.batch = n.batch,
	         priors = prior.list,
		 accept.rate = 0.43,
		 cov.model = 'exponential',
		 tuning = tuning.list,
	         n.omp.threads = 1,
	         verbose = FALSE,
		 NNGP = TRUE, 
		 n.neighbors = 5,
	         n.report = n.report,
	         n.burn = 400,
	         n.thin = 6, 
		 n.chains = 2,
	         k.fold = 2, 
	         k.fold.threads = 2)

  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_s3_class(out, "spMsPGOcc")
  expect_equal(out$pRE, TRUE)
  expect_equal(dim(out$sigma.sq.p.samples), c(n.post.samples, length(p.RE$sigma.sq.p)))
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})


# For helper functions ----------------------------------------------------
test_that("waicOCC works for spMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for spMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "array")
  expect_equal(dim(fitted.out), c(n.post.samples, N, J, max(n.rep)))
})

test_that("predict works for spMsPGOcc", {
  X.0 <- rbind(dat$X, matrix(c(1, rnorm(p.occ - 1)), nrow = 1))
  coords.0 <- rbind(dat$coords, matrix(c(0.538, 0.201), nrow = 1))
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, J + 1))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, J + 1))
})

test_that("posterior predictive checks work for spMsPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, max(n.rep)))
})
