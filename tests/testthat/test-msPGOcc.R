# Test msPGOcc.R  --------------------------------------------------------

skip_on_cran()

set.seed(400)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.3)
# Detection
alpha.mean <- c(0.5, 0.2, -0.1)
tau.sq.alpha <- c(0.2, 0.3, 1)
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
		sp = FALSE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3]
		 )
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  det.covs = det.covs)

# Occupancy initial values
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Initial values
inits.list <- list(alpha.comm = 0, 
		   beta.comm = 0, 
		   beta = 0, 
		   alpha = 0,
		   tau.sq.beta = 1, 
		   tau.sq.alpha = 1, 
		   z = apply(y, c(1, 2), max, na.rm = TRUE))

n.samples <- 5000
n.burn <- 3000
n.thin <- 2

out <- msPGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov.1 + det.cov.2, 
	       data = data.list, 
	       inits = inits.list, 
	       n.samples = n.samples, 
	       priors = prior.list, 
               n.omp.threads = 1, 
	       verbose = FALSE, 
	       n.report = 1000, 
	       n.burn = n.burn, 
	       n.chains = 2,
	       n.thin = n.thin)

test_that("out is of class msPGOcc", {
  expect_s3_class(out, "msPGOcc")
})

test_that("samples are the right size", {
  n.post.samples <- length(seq(from = n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(n.thin))) * out$n.chains
  expect_equal(dim(out$beta.comm.samples), c(n.post.samples, length(beta.mean)))
  expect_equal(dim(out$alpha.comm.samples), c(n.post.samples, length(alpha.mean)))
  expect_equal(dim(out$beta.samples), c(n.post.samples, length(beta)))
  expect_equal(dim(out$alpha.samples), c(n.post.samples, length(alpha)))
  expect_equal(dim(out$z.samples), c(n.post.samples, N, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, N, J))
})

test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

test_that("betas are correct", {
  beta.comm.quants <- apply(out$beta.comm.samples, 2, quantile, c(0.025, 0.975)) 
  expect_true(sum(beta.comm.quants[1, ] < beta.mean) == length(beta.mean))
  expect_true(sum(beta.comm.quants[2, ] > beta.mean) == length(beta.mean))
})

test_that("alphas are correct", {
  alpha.comm.quants <- apply(out$alpha.comm.samples, 2, quantile, c(0.025, 0.975)) 
  expect_true(sum(alpha.comm.quants[1, ] < alpha.mean) == length(alpha.mean))
  expect_true(sum(alpha.comm.quants[2, ] > alpha.mean) == length(alpha.mean))
})

test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

test_that("default priors and inits work", {

  out <- msPGOcc(occ.formula = ~ occ.cov, 
	         det.formula = ~ det.cov.1 + det.cov.2, 
	         data = data.list, 
	         n.samples = n.samples,
	         n.omp.threads = 1,
	         verbose = FALSE,
	         n.burn = n.burn, 
	         n.thin = n.thin)
  expect_s3_class(out, "msPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(msPGOcc(occ.formula = ~ occ.cov, 
	                det.formula = ~ det.cov.1 + det.cov.2, 
	                data = data.list, 
	                n.samples = 100,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100, 
	                n.burn = 1, 
	                n.thin = 1))
})

test_that("cross-validation works", {
  
  out <- msPGOcc(occ.formula = ~ occ.cov, 
	         det.formula = ~ det.cov.1 + det.cov.2, 
	         data = data.list, 
	         inits = inits.list, 
	         n.samples = n.samples, 
	         priors = prior.list, 
                 n.omp.threads = 1, 
	         verbose = FALSE, 
	         n.report = 1000, 
	         n.burn = n.burn, 
	         n.thin = n.thin, 
		 k.fold = 2)

  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

test_that("random effects on occurrence work", {
  set.seed(500)
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
  y <- dat$y
  X <- dat$X
  X.re <- dat$X.re
  X.p <- dat$X.p
  
  colnames(X.re) <- c("occ.factor")
  data.list <- list(y = y, occ.covs = X.re)
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
  n.report <- 200

  out <- msPGOcc(occ.formula = ~ (1 | occ.factor),
	         det.formula = ~ 1,
	         data = data.list,
	         inits = inits.list,
	         n.samples = n.samples,
	         priors = prior.list,
	         n.omp.threads = 1,
	         verbose = FALSE,
	         n.report = n.report,
	         n.burn = 400,
	         n.thin = 6, 
		 n.chains = 2,
	         k.fold = 2) 
  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_s3_class(out, "msPGOcc")
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$pRE, FALSE)
  expect_equal(dim(out$sigma.sq.psi.samples), c(n.post.samples, length(psi.RE$sigma.sq.psi)))
})

test_that("random effects on detection work", {
  #set.seed(500)
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
  p.RE <- list(levels = c(45), 
  	     sigma.sq.p = c(1.3))
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
  y <- dat$y
  X <- dat$X
  X.p.re <- dat$X.p.re
  X.p <- dat$X.p
  
  det.covs <- list(det.factor = X.p.re[, , 1])
  data.list <- list(y = y, det.covs = det.covs)
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
  n.report <- 200

  out <- msPGOcc(occ.formula = ~ 1,
	       det.formula = ~ (1 | det.factor),
	       data = data.list,
	       inits = inits.list,
	       n.samples = n.samples,
	       priors = prior.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = n.report,
	       n.burn = 400,
	       n.thin = 6, 
	       n.chains = 2,
	       k.fold = 2) 

  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_s3_class(out, "msPGOcc")
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
  expect_equal(dim(out$sigma.sq.p.samples), c(n.post.samples, length(p.RE$sigma.sq.p)))
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

test_that("random effects on occurrence and detection work", {
  set.seed(500)
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
  alpha.mean <- c(0, -0.3)
  tau.sq.alpha <- c(1, 1.6)
  p.det <- length(alpha.mean)
  # Random effects
  psi.RE <- list(levels = c(30), 
  	       sigma.sq.psi = c(2.1))
  p.RE <- list(levels = c(45), 
  	     sigma.sq.p = c(1.3))
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
  y <- dat$y
  X <- dat$X
  X.re <- dat$X.re
  X.p.re <- dat$X.p.re
  X.p <- dat$X.p
  
  occ.covs <- cbind(X[, 2, drop = FALSE], X.re)
  colnames(occ.covs) <- c("occ.cov", "occ.factor")
  det.covs <- list(det.cov = X.p[, , 2], 
  		 det.factor = X.p.re[, , 1])
  data.list <- list(y = y, det.covs = det.covs, occ.covs = occ.covs)
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
  n.report <- 200

  out <- msPGOcc(occ.formula = ~ occ.cov + (1 | occ.factor),
	         det.formula = ~ det.cov + (1 | det.factor),
	         data = data.list,
	         inits = inits.list,
	         n.samples = n.samples,
	         priors = prior.list,
	         n.omp.threads = 1,
	         verbose = FALSE,
	         n.report = n.report,
	         n.burn = 400,
	         n.thin = 6, 
	         k.fold = 2) 
  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_s3_class(out, "msPGOcc")
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
  expect_equal(dim(out$sigma.sq.p.samples), c(n.post.samples, length(p.RE$sigma.sq.p)))
  expect_equal(dim(out$sigma.sq.psi.samples), c(n.post.samples, length(psi.RE$sigma.sq.psi)))
  expect_equal(length(out$k.fold.deviance), N)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# For helper functions ----------------------------------------------------
#set.seed(500)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 8
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.4)
# Detection
alpha.mean <- c(0, 0.5)
tau.sq.alpha <- c(1, 1.4)
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
y <- dat$y
X <- dat$X
X.p <- dat$X.p

colnames(X) <- c('int', 'occ.cov')
det.covs <- list(det.cov = X.p[, , 2])
data.list <- list(y = y, occ.covs = X, det.covs = det.covs)
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
n.report <- 200

out <- msPGOcc(occ.formula = ~ occ.cov,
	       det.formula = ~ det.cov,
	       data = data.list,
	       inits = inits.list,
	       n.samples = n.samples,
	       priors = prior.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = n.report,
	       n.burn = 400,
	       n.thin = 6, 
	       k.fold = 2) 
n.post.samples <- length(seq(from = out$n.burn + 1, 
			     to = n.samples, 
			     by = as.integer(out$n.thin))) * out$n.chains
test_that("waicOCC works for msPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for msPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "array")
  expect_equal(dim(fitted.out), c(n.post.samples, N, J, max(n.rep)))
})

test_that("predict works for msPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, N, J))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, N, J))
})

test_that("posterior predictive checks work for msPGOcc", {
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
