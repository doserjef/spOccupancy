# Test spIntPGOcc.R  ------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(1010)
J.x <- 10
J.y <- 10
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
beta <- c(0.5)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(1, 0, 1)
alpha[[2]] <- runif(1, 0, 1)
alpha[[3]] <- runif(1, -1, 1)
alpha[[4]] <- runif(1, -1, 1)
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
occ.covs <- X
colnames(occ.covs) <- c('int')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int.1 = X.p[[1]][, , 1])
det.covs[[2]] <- list(int.2 = X.p[[2]][, , 1]) 
det.covs[[3]] <- list(int.3 = X.p[[3]][, , 1])
det.covs[[4]] <- list(int.4 = X.p[[4]][, , 1]) 
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
		   fix = TRUE,
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
occ.formula <- ~ 1
det.formula <- list(f.1 = ~ 1, f.2 = ~ 1, f.3 = ~ 1, f.4 = ~ 1)

out <- spIntPGOcc(occ.formula = occ.formula, 
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
		  verbose = FALSE, 
		  NNGP = TRUE, 
		  n.neighbors = 5, 
		  search.type = 'cb', 
		  n.report = 10, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2,
                  k.fold = 2,
                  k.fold.threads = 1)

# Test to make sure it worked --------- 
test_that("out is of class spIntPGOcc", {
  expect_s3_class(out, "spIntPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.burn = 500, 
		    n.thin = 2, 
                    k.fold = 2,
                    k.fold.data = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check output data is correct --------
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors and inits work", {

  out <- spIntPGOcc(occ.formula = occ.formula, 
		    det.formula = det.formula,
		    data = data.list,  
		    n.batch = n.batch, 
		    batch.length = batch.length, 
		    accept.rate = 0.43, 
		    cov.model = "exponential", 
		    tuning = tuning.list, 
		    n.omp.threads = 1, 
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.chains = 1)  
  expect_s3_class(out, "spIntPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(spIntPGOcc(occ.formula = occ.formula,
			   det.formula = det.formula, 
		           data = data.list,  
		           n.batch = n.batch, 
		           batch.length = batch.length, 
		           accept.rate = 0.43, 
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

# Check all correlation functions -----
test_that("all correlation functions work", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spIntPGOcc(occ.formula = occ.formula, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spIntPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for spIntPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  coords.0 <- dat$coords.pred
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Occurrence covariate only -----------------------------------------------
J.x <- 10
J.y <- 10
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
beta <- c(0.5, 1.2, -0.5)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(1, 0, 1)
alpha[[2]] <- runif(1, 0, 1)
alpha[[3]] <- runif(1, -1, 1)
alpha[[4]] <- runif(1, -1, 1)
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
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int.1 = X.p[[1]][, , 1])
det.covs[[2]] <- list(int.2 = X.p[[2]][, , 1]) 
det.covs[[3]] <- list(int.3 = X.p[[3]][, , 1])
det.covs[[4]] <- list(int.4 = X.p[[4]][, , 1]) 
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
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- list(f.1 = ~ 1, f.2 = ~ 1, f.3 = ~ 1, f.4 = ~ 1)

out <- spIntPGOcc(occ.formula = occ.formula, 
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
		  verbose = FALSE, 
		  NNGP = TRUE, 
		  n.neighbors = 5, 
		  search.type = 'cb', 
		  n.report = 10, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2,
                  k.fold = 2,
                  k.fold.threads = 1)

# Test to make sure it worked --------- 
test_that("out is of class spIntPGOcc", {
  expect_s3_class(out, "spIntPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.burn = 500, 
		    n.thin = 2, 
                    k.fold = 2,
                    k.fold.data = 1)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check output data is correct --------
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  # tmp.data$det.covs[[1]][[1]][1] <- NA
  # expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  # tmp.data$y[[1]][1, 1] <- NA
  # out <- spIntPGOcc(occ.formula = occ.formula,
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
  # expect_s3_class(out, "spIntPGOcc")
})

# Check default priors ----------------
test_that("default priors and inits work", {

  out <- spIntPGOcc(occ.formula = occ.formula, 
		    det.formula = det.formula,
		    data = data.list,  
		    n.batch = n.batch, 
		    batch.length = batch.length, 
		    accept.rate = 0.43, 
		    cov.model = "exponential", 
		    tuning = tuning.list, 
		    n.omp.threads = 1, 
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.chains = 1)  
  expect_s3_class(out, "spIntPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(spIntPGOcc(occ.formula = occ.formula,
			   det.formula = det.formula, 
		           data = data.list,  
		           n.batch = n.batch, 
		           batch.length = batch.length, 
		           accept.rate = 0.43, 
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

# Check all correlation functions -----
test_that("all correlation functions work", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spIntPGOcc(occ.formula = occ.formula, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spIntPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for spIntPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  coords.0 <- dat$coords.pred
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Detection covariate only ------------------------------------------------
J.x <- 10
J.y <- 10
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(2:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(2, 0, 1)
alpha[[2]] <- runif(1, 0, 1)
alpha[[3]] <- runif(3, -1, 1)
alpha[[4]] <- runif(2, -1, 1)
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
occ.covs <- X
colnames(occ.covs) <- c('int')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int.1 = X.p[[1]][, , 1], 
                      det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(int.2 = X.p[[2]][, , 1]) 
det.covs[[3]] <- list(int.3 = X.p[[3]][, , 1], 
                      det.cov.3.1 = X.p[[3]][, , 2], 
                      det.cov.3.2 = X.p[[3]][, , 3])
det.covs[[4]] <- list(int.4 = X.p[[4]][, , 1], 
                      det.cov.4.1 = X.p[[4]][, , 2]) 
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
occ.formula <- ~ 1 
det.formula <- list(f.1 = ~ det.cov.1.1, 
		    f.2 = ~ 1, 
		    f.3 = ~ det.cov.3.1 + det.cov.3.2, 
		    f.4 = ~ det.cov.4.1)

out <- spIntPGOcc(occ.formula = occ.formula, 
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
		  verbose = FALSE, 
		  NNGP = TRUE, 
		  n.neighbors = 5, 
		  search.type = 'cb', 
		  n.report = 10, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2,
                  k.fold = 2,
                  k.fold.threads = 1)

# Test to make sure it worked --------- 
test_that("out is of class spIntPGOcc", {
  expect_s3_class(out, "spIntPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.burn = 500, 
		    n.thin = 2, 
                    k.fold = 2,
                    k.fold.data = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check output data is correct --------
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs[3, ] <- NA
  # expect_error(spIntPGOcc(occ.formula = occ.formula,
  #                         det.formula = det.formula,
  #                         data = tmp.data,
  #                         n.batch = 40,
  #                         batch.length = batch.length,
  #                         cov.model = "exponential",
  #                         tuning = tuning.list,
  #                         NNGP = TRUE,
  #                         verbose = FALSE,
  #                         n.neighbors = 5,
  #                         search.type = 'cb',
  #                         n.report = 10,
  #                         n.burn = 500,
  #                         n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[2]][1] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  tmp.data$y[[1]][1, 1] <- NA
  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")
})

# Check default priors ----------------
test_that("default priors and inits work", {

  out <- spIntPGOcc(occ.formula = occ.formula, 
		    det.formula = det.formula,
		    data = data.list,  
		    n.batch = n.batch, 
		    batch.length = batch.length, 
		    accept.rate = 0.43, 
		    cov.model = "exponential", 
		    tuning = tuning.list, 
		    n.omp.threads = 1, 
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.chains = 1)  
  expect_s3_class(out, "spIntPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(spIntPGOcc(occ.formula = occ.formula,
			   det.formula = det.formula, 
		           data = data.list,  
		           n.batch = n.batch, 
		           batch.length = batch.length, 
		           accept.rate = 0.43, 
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

# Check all correlation functions -----
test_that("all correlation functions work", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spIntPGOcc(occ.formula = occ.formula, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spIntPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for spIntPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  coords.0 <- dat$coords.pred
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Covariates on both ------------------------------------------------------
J.x <- 10
J.y <- 10
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(2:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5, 1.2, -0.3)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(2, 0, 1)
alpha[[2]] <- runif(1, 0, 1)
alpha[[3]] <- runif(3, -1, 1)
alpha[[4]] <- runif(2, -1, 1)
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
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int.1 = X.p[[1]][, , 1], 
                      det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(int.2 = X.p[[2]][, , 1]) 
det.covs[[3]] <- list(int.3 = X.p[[3]][, , 1], 
                      det.cov.3.1 = X.p[[3]][, , 2], 
                      det.cov.3.2 = X.p[[3]][, , 3])
det.covs[[4]] <- list(int.4 = X.p[[4]][, , 1], 
                      det.cov.4.1 = X.p[[4]][, , 2]) 
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
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1, 
		    f.2 = ~ 1, 
		    f.3 = ~ det.cov.3.1 + det.cov.3.2, 
		    f.4 = ~ det.cov.4.1)

out <- spIntPGOcc(occ.formula = occ.formula, 
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
		  verbose = FALSE, 
		  NNGP = TRUE, 
		  n.neighbors = 5, 
		  search.type = 'cb', 
		  n.report = 10, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2,
                  k.fold = 2,
                  k.fold.threads = 1)

# Test to make sure it worked --------- 
test_that("out is of class spIntPGOcc", {
  expect_s3_class(out, "spIntPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.burn = 500, 
		    n.thin = 2, 
                    k.fold = 2,
                    k.fold.data = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check output data is correct --------
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[2]][1] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  tmp.data$y[[1]][1, 1] <- NA
  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")
})

# Check default priors ----------------
test_that("default priors and inits work", {

  out <- spIntPGOcc(occ.formula = occ.formula, 
		    det.formula = det.formula,
		    data = data.list,  
		    n.batch = n.batch, 
		    batch.length = batch.length, 
		    accept.rate = 0.43, 
		    cov.model = "exponential", 
		    tuning = tuning.list, 
		    n.omp.threads = 1, 
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.chains = 1)  
  expect_s3_class(out, "spIntPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(spIntPGOcc(occ.formula = occ.formula,
			   det.formula = det.formula, 
		           data = data.list,  
		           n.batch = n.batch, 
		           batch.length = batch.length, 
		           accept.rate = 0.43, 
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

# Check all correlation functions -----
test_that("all correlation functions work", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spIntPGOcc(occ.formula = occ.formula, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spIntPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for spIntPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  coords.0 <- dat$coords.pred
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Interactions on both ----------------------------------------------------
J.x <- 10
J.y <- 10
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(2:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5, 1.2, -0.3)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(2, 0, 1)
alpha[[2]] <- runif(1, 0, 1)
alpha[[3]] <- runif(3, -1, 1)
alpha[[4]] <- runif(2, -1, 1)
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
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int.1 = X.p[[1]][, , 1],
                      det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(int.2 = X.p[[2]][, , 1])
det.covs[[3]] <- list(int.3 = X.p[[3]][, , 1],
                      det.cov.3.1 = X.p[[3]][, , 2],
                      det.cov.3.2 = X.p[[3]][, , 3])
det.covs[[4]] <- list(int.4 = X.p[[4]][, , 1],
                      det.cov.4.1 = X.p[[4]][, , 2])
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
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1,
		    f.2 = ~ 1,
		    f.3 = ~ det.cov.3.1 * det.cov.3.2,
		    f.4 = ~ det.cov.4.1)

out <- spIntPGOcc(occ.formula = occ.formula,
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
		  verbose = FALSE,
		  NNGP = TRUE,
		  n.neighbors = 5,
		  search.type = 'cb',
		  n.report = 10,
		  n.burn = 500,
		  n.thin = 2,
		  n.chains = 2,
                  k.fold = 2,
                  k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class spIntPGOcc", {
  expect_s3_class(out, "spIntPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- spIntPGOcc(occ.formula = occ.formula,
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
		    verbose = FALSE,
		    NNGP = TRUE,
		    n.neighbors = 5,
		    search.type = 'cb',
		    n.report = 10,
		    n.burn = 500,
		    n.thin = 2,
                    k.fold = 2,
                    k.fold.data = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check output data is correct --------
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[2]][1] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  tmp.data$y[[1]][1, 1] <- NA
  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")
})

# Check default priors ----------------
test_that("default priors and inits work", {

  out <- spIntPGOcc(occ.formula = occ.formula,
		    det.formula = det.formula,
		    data = data.list,
		    n.batch = n.batch,
		    batch.length = batch.length,
		    accept.rate = 0.43,
		    cov.model = "exponential",
		    tuning = tuning.list,
		    n.omp.threads = 1,
		    verbose = FALSE,
		    NNGP = TRUE,
		    n.neighbors = 5,
		    search.type = 'cb',
		    n.report = 10,
		    n.chains = 1)
  expect_s3_class(out, "spIntPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(spIntPGOcc(occ.formula = occ.formula,
			   det.formula = det.formula,
		           data = data.list,
		           n.batch = n.batch,
		           batch.length = batch.length,
		           accept.rate = 0.43,
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

# Check all correlation functions -----
test_that("all correlation functions work", {
  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spIntPGOcc(occ.formula = occ.formula,
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

# Check waicOcc -----------------------
test_that("waicOCC works for spIntPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for spIntPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  coords.0 <- dat$coords.pred
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Model with fixed sigma.sq -----------------------------------------------
test_that("spIntPGOcc works with fixed sigma.sq", {
  out <- spIntPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.batch = 40,
	         batch.length = batch.length,
	         cov.model = "spherical",
		 priors = list(sigma.sq.ig = "fixed"),
	         tuning = list(phi = 0.5),
	         NNGP = TRUE,
		 verbose = FALSE,
	         n.neighbors = 5,
	         search.type = 'cb',
	         n.report = 10,
	         n.burn = 500,
	         n.chains = 1)
  expect_s3_class(out, "spIntPGOcc")
  expect_equal(length(unique(out$theta.samples[, 1])), 1)
})

# Uniform sigma sq --------------------------------------------------------
test_that("spIntPGOcc works with uniform prior on sigma.sq", {
  prior.list <- list(sigma.sq.unif = c(0, 5),
                     nu.unif = c(0.1, 4))
  tuning.list <- list(phi = 0.5, nu = 0.6, sigma.sq = 0.7)
  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")
  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")
})

# Different max(n.rep) ----------------------------------------------------
J.x <- 10
J.y <- 10
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(2:4, size = J.obs[i], replace = TRUE)
}
n.rep.max <- rep(5, n.data)
# Occupancy covariates
beta <- c(0.5, 1.2, -0.3)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(2, 0, 1)
alpha[[2]] <- runif(1, 0, 1)
alpha[[3]] <- runif(3, -1, 1)
alpha[[4]] <- runif(2, -1, 1)
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
sigma.sq <- 2
phi <- 3 / .5
sp <- TRUE

# Simulate occupancy data from multiple data sources. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = sp, 
		 sigma.sq = sigma.sq, phi = phi, cov.model = 'exponential', n.rep.max = n.rep.max)

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites
X.0 <- dat$X.pred
psi.0 <- dat$psi.pred
coords <- as.matrix(dat$coords.obs)
coords.0 <- as.matrix(dat$coords.pred)

# Package all data into a list
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int.1 = X.p[[1]][, , 1], 
                      det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(int.2 = X.p[[2]][, , 1]) 
det.covs[[3]] <- list(int.3 = X.p[[3]][, , 1], 
                      det.cov.3.1 = X.p[[3]][, , 2], 
                      det.cov.3.2 = X.p[[3]][, , 3])
det.covs[[4]] <- list(int.4 = X.p[[4]][, , 1], 
                      det.cov.4.1 = X.p[[4]][, , 2]) 
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
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1, 
		    f.2 = ~ 1, 
		    f.3 = ~ det.cov.3.1 + det.cov.3.2, 
		    f.4 = ~ det.cov.4.1)

out <- spIntPGOcc(occ.formula = occ.formula, 
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
		  verbose = FALSE, 
		  NNGP = TRUE, 
		  n.neighbors = 5, 
		  search.type = 'cb', 
		  n.report = 10, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2,
                  k.fold = 2,
                  k.fold.threads = 1)

# Test to make sure it worked --------- 
test_that("out is of class spIntPGOcc", {
  expect_s3_class(out, "spIntPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.burn = 500, 
		    n.thin = 2, 
                    k.fold = 2,
                    k.fold.data = 2)
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check output data is correct --------
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[2]][1] <- NA
  expect_error(spIntPGOcc(occ.formula = occ.formula,
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
  tmp.data$y[[1]][1, 1] <- NA
  out <- spIntPGOcc(occ.formula = occ.formula,
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
  expect_s3_class(out, "spIntPGOcc")
})

# Check default priors ----------------
test_that("default priors and inits work", {

  out <- spIntPGOcc(occ.formula = occ.formula, 
		    det.formula = det.formula,
		    data = data.list,  
		    n.batch = n.batch, 
		    batch.length = batch.length, 
		    accept.rate = 0.43, 
		    cov.model = "exponential", 
		    tuning = tuning.list, 
		    n.omp.threads = 1, 
		    verbose = FALSE, 
		    NNGP = TRUE, 
		    n.neighbors = 5, 
		    search.type = 'cb', 
		    n.report = 10, 
		    n.chains = 1)  
  expect_s3_class(out, "spIntPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(spIntPGOcc(occ.formula = occ.formula,
			   det.formula = det.formula, 
		           data = data.list,  
		           n.batch = n.batch, 
		           batch.length = batch.length, 
		           accept.rate = 0.43, 
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

# Check all correlation functions -----
test_that("all correlation functions work", {
  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")

  out <- spIntPGOcc(occ.formula = occ.formula, 
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
  expect_s3_class(out, "spIntPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spIntPGOcc(occ.formula = occ.formula, 
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

# Check waicOcc -----------------------
test_that("waicOCC works for spIntPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for spIntPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  coords.0 <- dat$coords.pred
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for spIntPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

