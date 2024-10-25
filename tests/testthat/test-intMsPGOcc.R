# test-intMsPGOcc.R -------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
# Simulate data -----------------------------------------------------------
set.seed(883)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
# for (i in 1:n.data) {
#   n.rep[[i]] <- sample(1:4, size = J.obs[i], replace = TRUE)
# }
n.rep[[1]] <- sample(1, size = J.obs[1], replace = TRUE)
n.rep[[2]] <- sample(2, size = J.obs[2], replace = TRUE)
n.rep[[3]] <- sample(3, size = J.obs[3], replace = TRUE)

N <- c(7, 4, 5)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
p.det.long <- c(1, 1, 1)

for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
p.det <- sum(p.det.long)

# Random effects
psi.RE <- list()
# Not currently supported
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE

dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
		   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	           psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model)

J <- nrow(dat$coords.obs)
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
X.re <- dat$X.re.obs
X.p.re <- dat$X.p.re
sites <- dat$sites
species <- dat$species

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
#colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int = X.p[[1]][, , 1])
det.covs[[2]] <- list(int = X.p[[2]][, , 1])
det.covs[[3]] <- list(int = X.p[[3]][, , 1])

data.list <- list(y = y,
		  occ.covs = occ.covs,
		  det.covs = det.covs,
                  sites = sites,
                  species = species)

prior.list <- list(beta.comm.normal = list(mean = 0,
				      var = 2.73),
		   alpha.comm.normal = list(mean = list(0, 0, 0),
			                    var = list(2.73, 2.73, 2.73)),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = list(0.1, 0.1, 0.1),
					  b = list(0.1, 0.1, 0.1)))
inits.list <- list(alpha.comm = list(0, 0, 0),
		   beta.comm = 0,
                   tau.sq.beta = 1,
		   tau.sq.alpha = list(1, 1, 1),
                   alpha = list(a = matrix(rnorm(p.det.long[1] * N[1]), N[1], p.det.long[1]),
				b = matrix(rnorm(p.det.long[2] * N[2]), N[2], p.det.long[2]),
				c = matrix(rnorm(p.det.long[3] * N[3]), N[3], p.det.long[3])))
n.samples <- 1000
occ.formula <- ~ 1
det.formula <- list(f.1 = ~ 1, 
		    f.2 = ~ 1, 
		    f.3 = ~ 1)
out <- intMsPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
		  inits = inits.list,
		  priors = prior.list,
	          data = data.list,
	          n.samples = n.samples,
                  n.omp.threads = 1,
	          verbose = FALSE,
	          n.report = 100,
	          n.burn = 400,
	          n.thin = 3,
	          n.chains = 2)

# To make sure it worked --------------
test_that("out is of class intMsPGOcc", {
  expect_s3_class(out, "intMsPGOcc")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- intMsPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.samples = n.samples,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "intMsPGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- intMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
               n.thin = 13,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
})

test_that("verbose prints to the screen", {
  expect_output(intMsPGOcc(occ.formula = occ.formula,
	                det.formula = det.formula,
	                data = data.list,
	                n.samples = 100,
                  parallel.chains = FALSE,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100,
	                n.burn = 1,
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for intMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check predictions -------------------
test_that("predict works for intMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
})

# Occurrence covariate only -----------------------------------------------
# Simulate data -----------------------------------------------------------
# set.seed(883)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
# for (i in 1:n.data) {
#   n.rep[[i]] <- sample(1:4, size = J.obs[i], replace = TRUE)
# }
n.rep[[1]] <- sample(1, size = J.obs[1], replace = TRUE)
n.rep[[2]] <- sample(2, size = J.obs[2], replace = TRUE)
n.rep[[3]] <- sample(3, size = J.obs[3], replace = TRUE)

N <- c(7, 4, 5)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.3)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4, 0.5)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
p.det.long <- c(1, 1, 1)

for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
p.det <- sum(p.det.long)

# Random effects
psi.RE <- list()
# Not currently supported
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE

dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
		   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	           psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model)

J <- nrow(dat$coords.obs)
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
X.re <- dat$X.re.obs
X.p.re <- dat$X.p.re
sites <- dat$sites
species <- dat$species

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1')
#colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int = X.p[[1]][, , 1])
det.covs[[2]] <- list(int = X.p[[2]][, , 1])
det.covs[[3]] <- list(int = X.p[[3]][, , 1])

data.list <- list(y = y,
		  occ.covs = occ.covs,
		  det.covs = det.covs,
                  sites = sites,
                  species = species)

prior.list <- list(beta.comm.normal = list(mean = 0,
				      var = 2.73),
		   alpha.comm.normal = list(mean = list(0, 0, 0),
			                    var = list(2.73, 2.73, 2.73)),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = list(0.1, 0.1, 0.1),
					  b = list(0.1, 0.1, 0.1)))
inits.list <- list(alpha.comm = list(0, 0, 0),
		   beta.comm = 0,
                   tau.sq.beta = 1,
		   tau.sq.alpha = list(1, 1, 1),
                   alpha = list(a = matrix(rnorm(p.det.long[1] * N[1]), N[1], p.det.long[1]),
				b = matrix(rnorm(p.det.long[2] * N[2]), N[2], p.det.long[2]),
				c = matrix(rnorm(p.det.long[3] * N[3]), N[3], p.det.long[3])))
n.samples <- 1000
occ.formula <- ~ occ.cov.1
det.formula <- list(f.1 = ~ 1, 
		    f.2 = ~ 1, 
		    f.3 = ~ 1)
out <- intMsPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
		 inits = inits.list,
		 priors = prior.list,
	         data = data.list,
	         n.samples = n.samples,
                 n.omp.threads = 1,
	         verbose = FALSE,
	         n.report = 100,
	         n.burn = 400,
	         n.thin = 3,
	         n.chains = 2)

# To make sure it worked --------------
test_that("out is of class intMsPGOcc", {
  expect_s3_class(out, "intMsPGOcc")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- intMsPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.samples = n.samples,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "intMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(intMsPGOcc(occ.formula = occ.formula,
	                det.formula = det.formula,
	                data = data.list,
	                n.samples = 100,
                  parallel.chains = FALSE,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100,
	                n.burn = 1,
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for intMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check predictions -------------------
test_that("predict works for intMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
})

# Detection covariate only ------------------------------------------------
# Simulate data -----------------------------------------------------------
# set.seed(883)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
n.rep[[1]] <- sample(1, size = J.obs[1], replace = TRUE)
n.rep[[2]] <- sample(2, size = J.obs[2], replace = TRUE)
n.rep[[3]] <- sample(3, size = J.obs[3], replace = TRUE)

N <- c(7, 4, 5)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
p.det.long <- c(3, 2, 2)

for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
p.det <- sum(p.det.long)

# Random effects
psi.RE <- list()
# Not currently supported
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE

dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
		   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	           psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model)

J <- nrow(dat$coords.obs)
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
X.re <- dat$X.re.obs
X.p.re <- dat$X.p.re
sites <- dat$sites
species <- dat$species

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
#colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int = X.p[[1]][, , 1], 
                      det.cov.1.1 = X.p[[1]][, , 2],
                      det.cov.2.1 = X.p[[1]][, , 3])
det.covs[[2]] <- list(int = X.p[[2]][, , 1], 
                      det.cov.1.2 = X.p[[2]][, , 2])
det.covs[[3]] <- list(int = X.p[[3]][, , 1], 
                      det.cov.1.3 = X.p[[3]][, , 2])

data.list <- list(y = y,
		  occ.covs = occ.covs,
		  det.covs = det.covs,
                  sites = sites,
                  species = species)

prior.list <- list(beta.comm.normal = list(mean = 0,
				      var = 2.73),
		   alpha.comm.normal = list(mean = list(0, 0, 0),
			                    var = list(2.73, 2.73, 2.73)),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = list(0.1, 0.1, 0.1),
					  b = list(0.1, 0.1, 0.1)))
inits.list <- list(alpha.comm = list(0, 0, 0),
		   beta.comm = 0,
                   tau.sq.beta = 1,
		   tau.sq.alpha = list(1, 1, 1),
                   alpha = list(a = matrix(rnorm(p.det.long[1] * N[1]), N[1], p.det.long[1]),
				b = matrix(rnorm(p.det.long[2] * N[2]), N[2], p.det.long[2]),
				c = matrix(rnorm(p.det.long[3] * N[3]), N[3], p.det.long[3])))
n.samples <- 1000
occ.formula <- ~ 1 
det.formula <- list(f.1 = ~ det.cov.1.1 + det.cov.2.1, 
		    f.2 = ~ det.cov.1.2, 
		    f.3 = ~ det.cov.1.3)
out <- intMsPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
		 inits = inits.list,
		 priors = prior.list,
	         data = data.list,
	         n.samples = n.samples,
                 n.omp.threads = 1,
	         verbose = FALSE,
	         n.report = 100,
	         n.burn = 400,
	         n.thin = 3,
	         n.chains = 2)

# To make sure it worked --------------
test_that("out is of class intMsPGOcc", {
  expect_s3_class(out, "intMsPGOcc")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- intMsPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.samples = n.samples,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "intMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(intMsPGOcc(occ.formula = occ.formula,
	                det.formula = det.formula,
	                data = data.list,
	                n.samples = 100,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100,
	                n.burn = 1,
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for intMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check predictions -------------------
test_that("predict works for intMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
})

# Covariates on both ------------------------------------------------------
# Simulate data -----------------------------------------------------------
# set.seed(883)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
n.rep[[1]] <- sample(1, size = J.obs[1], replace = TRUE)
n.rep[[2]] <- sample(2, size = J.obs[2], replace = TRUE)
n.rep[[3]] <- sample(3, size = J.obs[3], replace = TRUE)

N <- c(7, 4, 5)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.3, 0.4)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4, 0.5, 1.0)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
p.det.long <- c(3, 2, 2)

for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
p.det <- sum(p.det.long)

# Random effects
psi.RE <- list()
# Not currently supported
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE

dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
		   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	           psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model)

J <- nrow(dat$coords.obs)
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
X.re <- dat$X.re.obs
X.p.re <- dat$X.p.re
sites <- dat$sites
species <- dat$species

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int = X.p[[1]][, , 1], 
                      det.cov.1.1 = X.p[[1]][, , 2],
                      det.cov.2.1 = X.p[[1]][, , 3])
det.covs[[2]] <- list(int = X.p[[2]][, , 1], 
                      det.cov.1.2 = X.p[[2]][, , 2])
det.covs[[3]] <- list(int = X.p[[3]][, , 1], 
                      det.cov.1.3 = X.p[[3]][, , 2])

data.list <- list(y = y,
		  occ.covs = occ.covs,
		  det.covs = det.covs,
                  sites = sites,
                  species = species)

prior.list <- list(beta.comm.normal = list(mean = 0,
				      var = 2.73),
		   alpha.comm.normal = list(mean = list(0, 0, 0),
			                    var = list(2.73, 2.73, 2.73)),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = list(0.1, 0.1, 0.1),
					  b = list(0.1, 0.1, 0.1)))
inits.list <- list(alpha.comm = list(0, 0, 0),
		   beta.comm = 0,
                   tau.sq.beta = 1,
		   tau.sq.alpha = list(1, 1, 1),
                   alpha = list(a = matrix(rnorm(p.det.long[1] * N[1]), N[1], p.det.long[1]),
				b = matrix(rnorm(p.det.long[2] * N[2]), N[2], p.det.long[2]),
				c = matrix(rnorm(p.det.long[3] * N[3]), N[3], p.det.long[3])))
n.samples <- 1000
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1 + det.cov.2.1, 
		    f.2 = ~ det.cov.1.2, 
		    f.3 = ~ det.cov.1.3)
out <- intMsPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
		 inits = inits.list,
		 priors = prior.list,
	         data = data.list,
	         n.samples = n.samples,
                 n.omp.threads = 1,
	         verbose = FALSE,
	         n.report = 100,
	         n.burn = 400,
	         n.thin = 3,
	         n.chains = 2)

# To make sure it worked --------------
test_that("out is of class intMsPGOcc", {
  expect_s3_class(out, "intMsPGOcc")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- intMsPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.samples = n.samples,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "intMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(intMsPGOcc(occ.formula = occ.formula,
	                det.formula = det.formula,
	                data = data.list,
	                n.samples = 100,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100,
	                n.burn = 1,
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for intMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check predictions -------------------
test_that("predict works for intMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
})

# Interactions on both ----------------------------------------------------
# Simulate data -----------------------------------------------------------
# set.seed(883)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
n.rep[[1]] <- sample(1, size = J.obs[1], replace = TRUE)
n.rep[[2]] <- sample(2, size = J.obs[2], replace = TRUE)
n.rep[[3]] <- sample(3, size = J.obs[3], replace = TRUE)

N <- c(7, 4, 5)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.3, 0.4)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4, 0.5, 1.0)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
p.det.long <- c(3, 2, 2)

for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
p.det <- sum(p.det.long)

# Random effects
psi.RE <- list()
# Not currently supported
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE

dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
		   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	           psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model)

J <- nrow(dat$coords.obs)
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
X.re <- dat$X.re.obs
X.p.re <- dat$X.p.re
sites <- dat$sites
species <- dat$species

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(int = X.p[[1]][, , 1], 
                      det.cov.1.1 = X.p[[1]][, , 2],
                      det.cov.2.1 = X.p[[1]][, , 3])
det.covs[[2]] <- list(int = X.p[[2]][, , 1], 
                      det.cov.1.2 = X.p[[2]][, , 2])
det.covs[[3]] <- list(int = X.p[[3]][, , 1], 
                      det.cov.1.3 = X.p[[3]][, , 2])

data.list <- list(y = y,
		  occ.covs = occ.covs,
		  det.covs = det.covs,
                  sites = sites,
                  species = species)

prior.list <- list(beta.comm.normal = list(mean = 0,
				      var = 2.73),
		   alpha.comm.normal = list(mean = list(0, 0, 0),
			                    var = list(2.73, 2.73, 2.73)),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = list(0.1, 0.1, 0.1),
					  b = list(0.1, 0.1, 0.1)))
inits.list <- list(alpha.comm = list(0, 0, 0),
		   beta.comm = 0,
                   tau.sq.beta = 1,
		   tau.sq.alpha = list(1, 1, 1),
                   alpha = list(a = 0,
				b = matrix(rnorm(p.det.long[2] * N[2]), N[2], p.det.long[2]),
				c = matrix(rnorm(p.det.long[3] * N[3]), N[3], p.det.long[3])))
n.samples <- 1000
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1 * det.cov.2.1, 
		    f.2 = ~ det.cov.1.2, 
		    f.3 = ~ det.cov.1.3)
out <- intMsPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
		 inits = inits.list,
		 priors = prior.list,
	         data = data.list,
	         n.samples = n.samples,
                 n.omp.threads = 1,
	         verbose = FALSE,
	         n.report = 100,
	         n.burn = 400,
	         n.thin = 3,
	         n.chains = 2)

# To make sure it worked --------------
test_that("out is of class intMsPGOcc", {
  expect_s3_class(out, "intMsPGOcc")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- intMsPGOcc(occ.formula = occ.formula,
	         det.formula = det.formula,
	         data = data.list,
	         n.samples = n.samples,
	         n.omp.threads = 1,
	         verbose = FALSE)
  expect_s3_class(out, "intMsPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(intMsPGOcc(occ.formula = occ.formula,
	                det.formula = det.formula,
	                data = data.list,
	                n.samples = 100,
	                n.omp.threads = 1,
	                verbose = TRUE,
	                n.report = 100,
	                n.burn = 1,
	                n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for intMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check predictions -------------------
test_that("predict works for intMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- cbind(dat$X.pred, dat$X.pred[, 3] * dat$X.pred[, 2])
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, max(N), nrow(dat$X.pred)))
})
