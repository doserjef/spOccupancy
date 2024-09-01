# Test svcTPGBinom.R  ---------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(150)
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
weights <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  weights[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
svc.cols <- 1
psi.RE <- list()
phi <- rep(3 / .7, length(svc.cols))
sigma.sq <- rep(2, length(svc.cols))
nu <- rep(2, length(svc.cols))
dat <- simTBinom(J.x = J.x, J.y = J.y, n.time = n.time, weights = weights,
	         beta = beta, sp.only = sp.only, trend = trend,
	         psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq, 
                 phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , , drop = FALSE]
# X.p.re.0 <- dat$X.p.re[pred.indx, , , , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx, ]
weights.0 <- weights[pred.indx, ]
weights.fit <- weights[-pred.indx, ]

covs <- list(int = X[, , 1]) 

data.list <- list(y = y, coords = coords, covs = covs, weights = weights.fit)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(a = 0.5, b = 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0)
# Tuning
tuning.list <- list(phi = 1, nu = 1, rho = 1)

batch.length <- 25
n.batch <- 10
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ 1 

out <- svcTPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       svc.cols = svc.cols,
	       cov.model = "matern",
	       tuning = tuning.list,
	       ar1 = TRUE,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 100,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcTPGBinom", {
  expect_s3_class(out, "svcTPGBinom")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- svcTPGBinom(formula =formula, 
	       data = data.list, 
               n.thin = 13,
               n.batch = n.batch, 
               batch.length = batch.length, 
               accept.rate = 0.43, 
	       n.omp.threads = 1,
	       verbose = FALSE))
})
# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
		 svc.cols = svc.cols,
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("all correlation functions work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "gaussian", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")

  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "spherical", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcTPGBinom(formula = formula, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       svc.cols = svc.cols,
	       cov.model = "exponential", 
	       tuning = tuning.list, 
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       search.type = 'cb', 
	       n.report = 10, 
	       n.burn = 100, 
	       n.chains = 1))
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcTPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcTPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(fitted.out, out$y.rep.samples)
})

# Check predictions -------------------
test_that("predict works for svcTPGBinom", {
  # X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  # dimnames(X.0.full)[[3]] <- c('(Intercept)', 'trend', 'occ.factor.1')
  X.0.full <- X.0
  pred.out <- predict(out, X.0.full, coords.0, weights.0, verbose = FALSE, t.cols = 1:n.time.max, 
                      weights = weights.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
})

# Occurrence covariate only -----------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
weights <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  weights[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.8)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
psi.RE <- list()
svc.cols <- c(1, 2)
phi <- rep(3 / .7, length(svc.cols))
sigma.sq <- rep(2, length(svc.cols))
nu <- rep(2, length(svc.cols))
dat <- simTBinom(J.x = J.x, J.y = J.y, n.time = n.time, weights = weights,
	       beta = beta, sp.only = sp.only, trend = trend,
	       psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq, 
               phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx, ]
weights.0 <- weights[pred.indx, ]
weights.fit <- weights[-pred.indx, ]


covs <- list(int = X[, , 1], 
             trend = X[, , 2]) 

data.list <- list(y = y, coords = coords, covs = covs, weights = weights.fit)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0)
# Tuning
tuning.list <- list(phi = 1, nu = 1, rho = 1)

batch.length <- 25
n.batch <- 10
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ trend

out <- svcTPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       svc.cols = svc.cols,
	       cov.model = "matern",
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 100,
	       n.chains = 1,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcTPGBinom", {
  expect_s3_class(out, "svcTPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
		 svc.cols = svc.cols,
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("all correlation functions work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
		 svc.cols = svc.cols,
	         cov.model = "gaussian", 
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")

  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
		 svc.cols = svc.cols,
	         cov.model = "spherical", 
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcTPGBinom(formula = formula, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       svc.cols = svc.cols,
	       cov.model = "exponential", 
	       tuning = tuning.list, 
	       ar1 = TRUE,
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       search.type = 'cb', 
	       n.report = 10, 
	       n.burn = 100, 
	       n.chains = 1))
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$covs$trend[3, ] <- NA
  expect_error(svcTPGBinom(formula = formula, 
	               data = tmp.data, 
	               n.batch = 40, 
	               batch.length = batch.length, 
		       svc.cols = svc.cols,
	               cov.model = "exponential", 
	               tuning = tuning.list, 
	               NNGP = TRUE,
	               verbose = FALSE, 
	               n.neighbors = 5, 
	               search.type = 'cb', 
	               n.report = 10, 
	               n.burn = 100, 
	               n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcTPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcTPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(fitted.out, out$y.rep.samples)
})

# Check predictions -------------------
test_that("predict works for svcTPGBinom", {
  X.0.full <- X.0
  # X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  dimnames(X.0.full)[[3]] <- c('(Intercept)', 'trend')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE, t.cols = 1:n.time.max, 
                      weights = weights.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
})

# Random intercept on occurrence ------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
weights <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  weights[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
psi.RE <- list(levels = c(10),
               sigma.sq.psi = c(0.8))
svc.cols = c(1)
phi <- rep(3 / .7, length(svc.cols))
sigma.sq <- rep(2, length(svc.cols))
nu <- rep(2, length(svc.cols))
dat <- simTBinom(J.x = J.x, J.y = J.y, n.time = n.time, weights = weights,
	       beta = beta, sp.only = sp.only, trend = trend,
	       psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq, 
               phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx, ]
weights.0 <- weights[pred.indx, ]
weights.fit <- weights[-pred.indx, ]

covs <- list(int = X[, , 1], 
                 occ.factor.1 = X.re[, , 1]) 

data.list <- list(y = y, coords = coords, covs = covs, weights = weights.fit)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0)
# Tuning
tuning.list <- list(phi = 1, nu = 1, rho = 1)

batch.length <- 25
n.batch <- 10
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ (1 | occ.factor.1)
det.formula <- ~ 1 

out <- svcTPGBinom(formula = formula,
	           data = data.list,
	           inits = inits.list,
	           batch.length = batch.length,
	           n.batch = n.batch,
	           priors = prior.list,
	           accept.rate = 0.43,
	           cov.model = "matern",
	           svc.cols = svc.cols,
	           tuning = tuning.list,
	           n.omp.threads = 1,
	           verbose = FALSE,
	           NNGP = TRUE,
	           n.neighbors = 10,
	           n.report = n.report,
	           n.burn = 100,
	           n.chains = 2,
	           n.thin = 2,
	           k.fold = 2,
	           k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcTPGBinom", {
  expect_s3_class(out, "svcTPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$occ.factor.1 <- factor(data.list$covs$occ.factor.1)
  expect_error(out <- svcTPGBinom(formula = formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
			      svc.cols = svc.cols,
                              tuning = tuning.list, 
			      ar1 = TRUE,
                              n.omp.threads = 1, 
                              verbose = FALSE, 
                              NNGP = TRUE, 
                              n.neighbors = 5, 
                              search.type = 'cb', 
                              n.report = 10, 
                              n.burn = 100, 
		              n.thin = 2,
		              n.chains = 1))
  data.list$covs$occ.factor.1 <- as.character(factor(data.list$covs$occ.factor.1))
  expect_error(out <- svcTPGBinom(formula = formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
			      svc.cols = svc.cols,
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
#  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
#	       sort(unique(c(X.p.re))))
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
		 svc.cols = svc.cols,
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("all correlation functions work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "gaussian", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")

  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "spherical", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcTPGBinom(formula = formula, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       cov.model = "exponential", 
	       svc.cols = svc.cols,
	       tuning = tuning.list, 
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       search.type = 'cb', 
	       n.report = 10, 
	       n.burn = 100, 
	       n.chains = 1))
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcTPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcTPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(fitted.out, out$y.rep.samples)
})

# Check predictions -------------------
test_that("predict works for svcTPGBinom", {
  X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  dimnames(X.0.full)[[3]] <- c('(Intercept)', 'occ.factor.1')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE, t.cols = 1:n.time.max, 
                      weights = weights.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
})

# Multiple random intercepts on occurrence --------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
weights <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  weights[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
psi.RE <- list(levels = c(10, 15),
               sigma.sq.psi = c(0.8, 0.3))
svc.cols <- c(1)
phi <- rep(3 / .7, length(svc.cols))
sigma.sq <- rep(2, length(svc.cols))
nu <- rep(2, length(svc.cols))
dat <- simTBinom(J.x = J.x, J.y = J.y, n.time = n.time, weights = weights,
	       beta = beta, sp.only = sp.only, trend = trend,
	       psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq, 
               phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx, ]
weights.0 <- weights[pred.indx, ]
weights.fit <- weights[-pred.indx, ]

covs <- list(int = X[, , 1], 
                 occ.factor.1 = X.re[, , 1], 
                 occ.factor.2 = X.re[, , 2]) 

data.list <- list(y = y, coords = coords, covs = covs, weights = weights.fit)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0)
# Tuning
tuning.list <- list(phi = 1, nu = 1, rho = 1)

batch.length <- 25
n.batch <- 10
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2)

out <- svcTPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       cov.model = "matern",
	       svc.cols = svc.cols,
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 100,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcTPGBinom", {
  expect_s3_class(out, "svcTPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$occ.factor.1 <- factor(data.list$covs$occ.factor.1)
  expect_error(out <- svcTPGBinom(formula = formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
			      svc.cols = svc.cols,
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
  data.list$covs$occ.factor.1 <- as.character(factor(data.list$covs$occ.factor.1))
  expect_error(out <- svcTPGBinom(formula = formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
			      svc.cols = svc.cols,
                              tuning = tuning.list, 
			      ar1 = TRUE,
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

test_that("default priors and inits work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
		 svc.cols = svc.cols,
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("all correlation functions work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "gaussian", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")

  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "spherical", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcTPGBinom(formula = formula, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       cov.model = "exponential", 
	       svc.cols = svc.cols,
	       tuning = tuning.list, 
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       search.type = 'cb', 
	       n.report = 10, 
	       n.burn = 100, 
	       n.chains = 1))
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcTPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcTPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(fitted.out, out$y.rep.samples)
})

# Check predictions -------------------
test_that("predict works for svcTPGBinom", {
  X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  dimnames(X.0.full)[[3]] <- c('(Intercept)', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE, t.cols = 1:n.time.max, 
                      weights = weights.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
})

# Occurrence REs + covariates ---------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
weights <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  weights[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.5)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
psi.RE <- list(levels = c(10, 15),
               sigma.sq.psi = c(0.8, 0.3))
svc.cols <- c(1, 2)
phi <- rep(3 / .7, length(svc.cols))
sigma.sq <- rep(2, length(svc.cols))
nu <- rep(2, length(svc.cols))
dat <- simTBinom(J.x = J.x, J.y = J.y, n.time = n.time, weights = weights,
	       beta = beta, sp.only = sp.only, trend = trend,
	       psi.RE = psi.RE, sp = TRUE, sigma.sq = sigma.sq, 
               phi = phi, cov.model = 'matern', nu = nu, svc.cols = svc.cols)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx, ]
weights.0 <- weights[pred.indx, ]
weights.fit <- weights[-pred.indx, ]

covs <- list(int = X[, , 1], 
		 trend = X[, , 2],
                 occ.factor.1 = X.re[, , 1], 
                 occ.factor.2 = X.re[, , 2]) 

data.list <- list(y = y, coords = coords, covs = covs, weights = weights.fit)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(alpha = 0, beta = 0)
# Tuning
tuning.list <- list(phi = 1, nu = 1, rho = 1)

batch.length <- 25
n.batch <- 10
n.report <- 10
n.samples <- batch.length * n.batch
formula <- ~ trend + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1 

out <- svcTPGBinom(formula = formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       cov.model = "matern",
	       svc.cols = svc.cols,
	       tuning = tuning.list,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       NNGP = TRUE,
	       n.neighbors = 10,
	       n.report = n.report,
	       n.burn = 100,
	       n.chains = 2,
	       n.thin = 2,
	       k.fold = 2,
	       k.fold.threads = 2)

# Test to make sure it worked ---------
test_that("out is of class svcTPGBinom", {
  expect_s3_class(out, "svcTPGBinom")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$occ.factor.1 <- factor(data.list$covs$occ.factor.1)
  expect_error(out <- svcTPGBinom(formula = formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
			      svc.cols = svc.cols,
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
  data.list$covs$occ.factor.1 <- as.character(factor(data.list$covs$occ.factor.1))
  expect_error(out <- svcTPGBinom(formula = formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
			      svc.cols = svc.cols,
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

test_that("default priors and inits work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
		 svc.cols = svc.cols,
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("all correlation functions work", {
  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "gaussian", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")

  out <- svcTPGBinom(formula = formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "spherical", 
		 svc.cols = svc.cols,
	         tuning = list(phi = 0.5), 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "svcTPGBinom")
})

test_that("verbose prints to the screen", {
  expect_output(svcTPGBinom(formula = formula, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       cov.model = "exponential", 
	       svc.cols = svc.cols,
	       tuning = tuning.list, 
	       ar1 = TRUE,
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       search.type = 'cb', 
	       n.report = 10, 
	       n.burn = 100, 
	       n.chains = 1))
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$covs$trend[3, ] <- NA
  expect_error(svcTPGBinom(formula = formula, 
                       data = tmp.data, 
                       n.batch = 40, 
                       batch.length = batch.length, 
                       cov.model = "exponential", 
		       svc.cols = svc.cols,
                       tuning = tuning.list, 
                       NNGP = TRUE,
                       verbose = FALSE, 
                       n.neighbors = 5, 
                       search.type = 'cb', 
                       n.report = 10, 
                       n.burn = 100, 
                       n.chains = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcTPGBinom", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcTPGBinom", {
  fitted.out <- fitted(out)
  expect_equal(fitted.out, out$y.rep.samples)
})

# Check predictions -------------------
test_that("predict works for svcTPGBinom", {
  X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  dimnames(X.0.full)[[3]] <- c('(Intercept)', 'trend', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE, t.cols = 1:n.time.max, 
                      weights = weights.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
})
