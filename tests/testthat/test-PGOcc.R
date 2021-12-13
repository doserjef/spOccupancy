# Test PGOcc.R  ----------------------------------------------------------

skip_on_cran()

set.seed(207)

J.x <- 25
J.y <- 25
J <- J.x * J.y
n.rep <- sample(4, J, replace = TRUE)
beta <- c(0.5, -0.15)
p.occ <- length(beta)
alpha <- c(0.7, 0.4)
p.det <- length(alpha)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      sp = FALSE)
occ.covs <- dat$X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov = dat$X.p[, , 2])
# Data bundle
data.list <- list(y = dat$y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Initial values
inits.list <- list(alpha = 0, beta = 0,
		      z = apply(data.list$y, 1, max, na.rm = TRUE))

n.samples <- 5000
n.thin <- 1
n.burn <- 4000
n.report <- 1000

out <- PGOcc(occ.formula = ~ occ.cov, 
	     det.formula = ~ det.cov, 
	     data = data.list, 
	     inits = inits.list,
	     n.samples = n.samples,
	     priors = prior.list,
	     n.omp.threads = 1,
	     verbose = FALSE,
	     n.report = n.report, 
	     n.burn = n.burn, 
	     n.thin = n.thin)

test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

test_that("samples are the right size", {
  n.post.samples <- length(seq(from = n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(n.thin)))
  expect_equal(dim(out$beta.samples), c(n.post.samples, length(beta)))
  expect_equal(dim(out$alpha.samples), c(n.post.samples, length(alpha)))
  expect_equal(dim(out$z.samples), c(n.post.samples, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, J))
})

test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

test_that("betas are correct", {
  beta.quants <- apply(out$beta.samples, 2, quantile, c(0.025, 0.975)) 
  expect_true(sum(beta.quants[1, ] < beta) == length(beta))
  expect_true(sum(beta.quants[2, ] > beta) == length(beta))
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

  out <- PGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = n.report, 
	       n.burn = n.burn, 
	       n.thin = n.thin)
  expect_s3_class(out, "PGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

test_that("cross-validation works", {
  
  out <- PGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov, 
	       data = data.list, 
	       inits = inits.list,
	       n.samples = n.samples,
	       priors = prior.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = n.report, 
	       n.burn = n.burn, 
	       n.thin = n.thin, 
	       k.fold = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

test_that("PGOccREOcc works", {
  J.x <- 15
  J.y <- 15
  J <- J.x * J.y
  n.rep <- sample(2:4, J, replace = TRUE)
  beta <- c(0.3)
  p.occ <- length(beta)
  alpha <- c(-0.5)
  p.det <- length(alpha)
  psi.RE <- list(levels = c(45), 
  	       sigma.sq.psi = c(1.3))
  p.RE <- list()
  dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
  	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
  y <- dat$y
  X <- dat$X
  X.re <- dat$X.re
  X.p <- dat$X.p
  
  X <- cbind(X, X.re)
  colnames(X) <- c("int", "occ.factor")
  occ.covs <- as.data.frame(X)
  data.list <- list(y = y, occ.covs = occ.covs)
  # Priors
  prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
  		   alpha.normal = list(mean = 0, var = 2.72))
  # Starting values
  inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
  
  n.samples <- 2000
  n.report <- 1
  
  set.seed(400)
  out <- PGOcc(occ.formula = ~ (1 | occ.factor),
	     det.formula = ~ 1,
	     data = data.list,
	     inits = inits.list,
	     n.samples = n.samples,
	     priors = prior.list,
	     n.omp.threads = 1,
	     verbose = FALSE,
	     n.report = n.report,
	     n.burn = 1000,
	     n.thin = 2, 
	     k.fold = 2)

  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin)))
  expect_s3_class(out, "PGOcc")
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$pRE, FALSE)
  expect_equal(dim(out$sigma.sq.psi.samples), c(n.post.samples, length(psi.RE$sigma.sq.psi)))
})

test_that("random effects on detection work", {
  J.x <- 15
  J.y <- 15
  J <- J.x * J.y
  n.rep <- sample(2:4, J, replace = TRUE)
  beta <- c(0.3)
  p.occ <- length(beta)
  alpha <- c(-0.5)
  p.det <- length(alpha)
  psi.RE <- list()
  p.RE <- list(levels = c(100), 
  	     sigma.sq.p = c(2.2))
  dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
  	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
  y <- dat$y
  X <- dat$X
  X.p.re <- dat$X.p.re
  X.p <- dat$X.p
  
  det.covs <- list(det.factor = X.p.re[, , 1])
  data.list <- list(y = y, det.covs = det.covs)
  # Priors
  prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
  		   alpha.normal = list(mean = 0, var = 2.72))
  # Starting values
  inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
  
  n.samples <- 2000
  n.report <- 1
  
  set.seed(400)
  out <- PGOcc(occ.formula = ~ 1,
	       det.formula = ~ (1 | det.factor),
	       data = data.list,
	       inits = inits.list,
	       n.samples = n.samples,
	       priors = prior.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = n.report,
	       n.burn = 1000,
	       n.thin = 2, 
	       k.fold = 2) 
  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin)))
  expect_s3_class(out, "PGOcc")
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
  expect_equal(dim(out$sigma.sq.p.samples), c(n.post.samples, length(p.RE$sigma.sq.p)))
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

test_that("random effects on occurrence and detection work", {
  J.x <- 15
  J.y <- 15
  J <- J.x * J.y
  n.rep <- sample(2:4, J, replace = TRUE)
  beta <- c(0.3, 0.5)
  p.occ <- length(beta)
  alpha <- c(-0.5, 0.3)
  p.det <- length(alpha)
  psi.RE <- list(levels = c(35), 
  	       sigma.sq.psi = c(0.4))
  p.RE <- list(levels = c(100), 
  	     sigma.sq.p = c(2.3))
  dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
  	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
  y <- dat$y
  X <- dat$X
  X.re <- dat$X.re
  X.p.re <- dat$X.p.re
  X.p <- dat$X.p
  
  X <- cbind(X, X.re)
  colnames(X) <- c("int", "occ.cov", "occ.factor")
  occ.covs <- as.data.frame(X)
  det.covs <- list(det.cov = X.p[, , 2], 
  		 det.factor = X.p.re[, , 1], 
  		 occ.cov = X[, 2])
  data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
  # Priors
  prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
  		   alpha.normal = list(mean = 0, var = 2.72))
  # Starting values
  inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))
  
  n.samples <- 2000
  n.report <- 1
  
  set.seed(400)
  out <- PGOcc(occ.formula = ~ occ.cov + (1 | occ.factor),
	     det.formula = ~ det.cov + (1 | det.factor),
	     data = data.list,
	     inits = inits.list,
	     n.samples = n.samples,
	     priors = prior.list,
	     n.omp.threads = 1,
	     verbose = FALSE,
	     n.report = n.report,
	     n.burn = 1000,
	     n.thin = 2, 
	     k.fold = 2) 

  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin)))
  expect_s3_class(out, "PGOcc")
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
  expect_equal(dim(out$sigma.sq.p.samples), c(n.post.samples, length(p.RE$sigma.sq.p)))
  expect_equal(dim(out$sigma.sq.psi.samples), c(n.post.samples, length(psi.RE$sigma.sq.psi)))
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# For helper functions
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(4, J, replace = TRUE)
beta <- c(0.5, -0.15)
p.occ <- length(beta)
alpha <- c(0.7, 0.4)
p.det <- length(alpha)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      sp = FALSE)
occ.covs <- dat$X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov = dat$X.p[, , 2])
# Data bundle
data.list <- list(y = dat$y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Initial values
inits.list <- list(alpha = 0, beta = 0,
		      z = apply(data.list$y, 1, max, na.rm = TRUE))

n.samples <- 5000
n.thin <- 1
n.burn <- 4000
n.report <- 1000

out <- PGOcc(occ.formula = ~ occ.cov, 
           det.formula = ~ det.cov, 
           data = data.list, 
           inits = inits.list,
           n.samples = n.samples,
           priors = prior.list,
           n.omp.threads = 1,
           verbose = FALSE,
           n.report = n.report, 
           n.burn = n.burn, 
           n.thin = n.thin)
n.post.samples <- length(seq(from = n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(n.thin)))

test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "array")
  expect_equal(dim(fitted.out), c(n.post.samples, J, max(n.rep)))
})

test_that("predict works for PGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, J))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, J))
})

test_that("posterior predictive checks work for PGOcc", {
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
