# Test intPGOcc.R  ----------------------------------------------------------

skip_on_cran()

set.seed(1006)

# Simulate Data -----------------------------------------------------------
J.x <- 25
J.y <- 25
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
beta <- c(0.5, 1)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(2, -1, 1)
}
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)

# Simulate occupancy data. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = FALSE)

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2]) 
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2]) 
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , 2]) 
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2]) 
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  det.covs = det.covs, 
		  sites = sites)

J <- length(dat$z.obs)
# Initial values
inits.list <- list(alpha = list(0, 0, 0, 0), 
		      beta = 0, 
		      z = rep(1, J))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
		   alpha.normal = list(mean = list(0, 0, 0, 0), 
			               var = list(2.72, 2.72, 2.72, 2.72)))
n.samples <- 10000
out <- intPGOcc(occ.formula = ~ occ.cov, 
                det.formula = list(f.1 = ~ det.cov.1.1, 
		                   f.2 = ~ det.cov.2.1, 
		                   f.3 = ~ det.cov.3.1, 
		                   f.4 = ~ det.cov.4.1), 
		data = data.list,
		inits = inits.list,
		n.samples = n.samples, 
		priors = prior.list, 
		n.omp.threads = 1, 
		verbose = FALSE, 
		n.report = 1000, 
		n.burn = 5000, 
		n.thin = 5, 
		n.chains = 1)
n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains

test_that("out is of class intPGOcc", {
  expect_s3_class(out, "intPGOcc")
})

test_that("samples are the right size", {
  n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains
  expect_equal(dim(out$beta.samples), c(n.post.samples, p.occ))
  expect_equal(dim(out$alpha.samples), c(n.post.samples, p.det))
  expect_equal(dim(out$z.samples), c(n.post.samples, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, J))
})

test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

test_that("default priors and inits work", {

  out <- intPGOcc(occ.formula = ~ occ.cov, 
                  det.formula = list(f.1 = ~ det.cov.1.1, 
		                     f.2 = ~ det.cov.2.1, 
		                     f.3 = ~ det.cov.3.1, 
		                     f.4 = ~ det.cov.4.1), 
		  data = data.list,
		  #inits = inits.list,
		  n.samples = 100, 
		  #priors = prior.list, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.report = 1000, 
		  n.burn = 1, 
		  n.thin = 1, 
		  n.chains = 1)
  expect_s3_class(out, "intPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(out <- intPGOcc(occ.formula = ~ occ.cov, 
                                det.formula = list(f.1 = ~ det.cov.1.1, 
		                                   f.2 = ~ det.cov.2.1, 
		                                   f.3 = ~ det.cov.3.1, 
		                                   f.4 = ~ det.cov.4.1), 
		                data = data.list,
		                inits = inits.list,
		                n.samples = 100, 
		                priors = prior.list, 
		                n.omp.threads = 1, 
		                verbose = TRUE, 
		                n.report = 25, 
		                n.burn = 1, 
		                n.thin = 1, 
		                n.chains = 1))
})

test_that("cross-validation works", {
  
  out <- intPGOcc(occ.formula = ~ occ.cov, 
                  det.formula = list(f.1 = ~ det.cov.1.1, 
		                     f.2 = ~ det.cov.2.1, 
		                     f.3 = ~ det.cov.3.1, 
		                     f.4 = ~ det.cov.4.1), 
		  data = data.list,
		  #inits = inits.list,
		  n.samples = 100, 
		  #priors = prior.list, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.report = 1000, 
		  n.burn = 1, 
		  n.thin = 1, 
		  n.chains = 1, 
		  k.fold = 2)
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)

  out <- intPGOcc(occ.formula = ~ occ.cov, 
                  det.formula = list(f.1 = ~ det.cov.1.1, 
		                     f.2 = ~ det.cov.2.1, 
		                     f.3 = ~ det.cov.3.1, 
		                     f.4 = ~ det.cov.4.1), 
		  data = data.list,
		  inits = inits.list,
		  n.samples = 100, 
		  priors = prior.list, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.report = 1000, 
		  n.burn = 1, 
		  n.thin = 1, 
		  n.chains = 2, 
		  k.fold = 2, 
		  k.fold.data = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

test_that("waicOCC works for intPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
})

test_that("fitted works for intPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(sapply(fitted.out, function(a) dim(apply(a, c(2, 3), mean))), 
	       sapply(dat$y, dim))
})

test_that("predict works for intPGOcc", {
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for PGOcc", {
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(sapply(ppc.out$fit.y, length), rep(n.post.samples, n.data))
  expect_equal(sapply(ppc.out$fit.y.rep, length), rep(n.post.samples, n.data))
  
  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(sapply(ppc.out$fit.y, length), rep(n.post.samples, n.data))
  expect_equal(sapply(ppc.out$fit.y.rep, length), rep(n.post.samples, n.data))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(sapply(ppc.out$fit.y, length), rep(n.post.samples, n.data))
  expect_equal(sapply(ppc.out$fit.y.rep, length), rep(n.post.samples, n.data))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(sapply(ppc.out$fit.y, length), rep(n.post.samples, n.data))
  expect_equal(sapply(ppc.out$fit.y.rep, length), rep(n.post.samples, n.data))
})
