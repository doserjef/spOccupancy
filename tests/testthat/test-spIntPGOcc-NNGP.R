# Test spIntPGOcc.R  ------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

set.seed(1006)

# Simulate Data -----------------------------------------------------------
J.x <- 8
J.y <- 8
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
beta <- c(0.5, 0.5)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(2, 0, 1)
alpha[[2]] <- runif(3, 0, 1)
alpha[[3]] <- runif(2, -1, 1)
alpha[[4]] <- runif(4, -1, 1)
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
sigma.sq <- 2
phi <- 3 / .5
sp <- TRUE

# Simulate occupancy data from multiple data sources. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = sp, 
		 sigma.sq = sigma.sq, phi = phi, cov.model = 'exponential')

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites
X.0 <- dat$X.pred
psi.0 <- dat$psi.pred
coords <- as.matrix(dat$coords.obs)
coords.0 <- as.matrix(dat$coords.pred)

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2], 
		      det.cov.2.2 = X.p[[2]][, , 3])
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , 2])
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2], 
		      det.cov.4.2 = X.p[[4]][, , 3], 
		      det.cov.4.3 = X.p[[4]][, , 4])
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  det.covs = det.covs, 
		  sites = sites, 
		  coords = coords)

J <- length(dat$z.obs)

# Initial values
inits.list <- list(alpha = list(0, 0, 0, 0), 
		      beta = 0, 
		      phi = 3 / .5, 
		      sigma.sq = 2, 
		      w = rep(0, J), 
		      z = rep(1, J))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
		   alpha.normal = list(mean = list(0, 0, 0, 0), 
				       var = list(2.72, 2.72, 2.72, 2.72)),
		   phi.unif = c(3/1, 3/.1), 
		   sigma.sq.ig = c(2, 2))
# Tuning
tuning.list <- list(phi = 0.3) 

# Number of batches
n.batch <- 40
# Batch length
batch.length <- 25
n.samples <- n.batch * batch.length

out <- spIntPGOcc(occ.formula = ~ occ.cov, 
		  det.formula = list(f.1 = ~ det.cov.1.1, 
				     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
				     f.3 = ~ det.cov.3.1, 
				     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
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
		  n.thin = 1)

n.post.samples <- length(seq(from = out$n.burn + 1, 
			       to = n.samples, 
			       by = as.integer(out$n.thin))) * out$n.chains

test_that("out is of class spIntPGOcc", {
  expect_s3_class(out, "spIntPGOcc")
})

test_that("samples are the right size", {
  expect_equal(dim(out$beta.samples), c(n.post.samples, p.occ))
  expect_equal(dim(out$alpha.samples), c(n.post.samples, p.det))
  expect_equal(dim(out$z.samples), c(n.post.samples, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, J))
})

test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

test_that("default priors and inits work", {

  out <- spIntPGOcc(occ.formula = ~ occ.cov, 
		  det.formula = list(f.1 = ~ det.cov.1.1, 
				     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
				     f.3 = ~ det.cov.3.1, 
				     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
		  data = data.list,  
		  #inits = inits.list, 
		  n.batch = n.batch, 
		  batch.length = batch.length, 
		  accept.rate = 0.43, 
		  #priors = prior.list, 
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
		  n.chains = 2)  
  expect_s3_class(out, "spIntPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(spIntPGOcc(occ.formula = ~ occ.cov, 
		  det.formula = list(f.1 = ~ det.cov.1.1, 
				     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
				     f.3 = ~ det.cov.3.1, 
				     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
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
		  n.chains = 1))
})

test_that("cross-validation works", {
  
  out <- spIntPGOcc(occ.formula = ~ occ.cov, 
		  det.formula = list(f.1 = ~ det.cov.1.1, 
				     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
				     f.3 = ~ det.cov.3.1, 
				     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
		  data = data.list,  
		  #inits = inits.list, 
		  n.batch = n.batch, 
		  batch.length = batch.length, 
		  accept.rate = 0.43, 
		  #priors = prior.list, 
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
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)

  # out <- spIntPGOcc(occ.formula = ~ occ.cov, 
  #       	  det.formula = list(f.1 = ~ det.cov.1.1, 
  #       			     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
  #       			     f.3 = ~ det.cov.3.1, 
  #       			     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
  #       	  data = data.list,  
  #       	  #inits = inits.list, 
  #       	  n.batch = n.batch, 
  #       	  batch.length = batch.length, 
  #       	  accept.rate = 0.43, 
  #       	  #priors = prior.list, 
  #       	  cov.model = "exponential", 
  #       	  tuning = tuning.list, 
  #       	  n.omp.threads = 1, 
  #       	  verbose = FALSE, 
  #       	  NNGP = TRUE, 
  #       	  n.neighbors = 5, 
  #       	  search.type = 'cb', 
  #       	  n.report = 10, 
  #       	  n.burn = 500, 
  #       	  n.thin = 1, 
  #       	  n.chains = 2, 
  #       	  k.fold = 2, 
  #       	  k.fold.data = 1)  
  # expect_equal(length(out$k.fold.deviance), 1)
  # expect_type(out$k.fold.deviance, "double")
  # expect_equal(sum(out$k.fold.deviance < 0), 0)
})

test_that("waicOCC works for spIntPGOcc", {
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
  coords.0 <- dat$coords.pred
  pred.out <- predict(out, X.0, coords.0)
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
