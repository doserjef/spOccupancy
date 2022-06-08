# Test spPGOcc.R  ---------------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(150)
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
psi.RE <- list()
alpha <- c(-1)
p.det <- length(alpha)
p.RE <- list()
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep,
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
	       psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
               phi = phi, cov.model = 'matern', nu = nu)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
# X.re <- dat$X.re[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
# X.re.0 <- dat$X.re[pred.indx, , , drop = FALSE]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , , drop = FALSE]
# X.p.re <- dat$X.p.re[-pred.indx, , , , drop = FALSE]
X.p.0 <- dat$X.p[pred.indx, , , , drop = FALSE]
# X.p.re.0 <- dat$X.p.re[pred.indx, , , , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx]

occ.covs <- list(int = X[, , 1]) 
det.covs <- list(int = X.p[, , , 1]) 

data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 10
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ 1 
det.formula <- ~ 1 

out <- spTPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       cov.model = "matern",
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
test_that("out is of class spTPGOcc", {
  expect_s3_class(out, "spTPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE error ----------------------
# test_that("random effect gives error when non-numeric", {
#   data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
#   data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
#   expect_error(out <- spTPGOcc(occ.formula = occ.formula,
# 	                      det.formula = det.formula,
#                               data = data.list,
#                               inits = inits.list, 
#                               n.batch = n.batch, 
#                               batch.length = batch.length, 
#                               accept.rate = 0.43, 
#                               priors = prior.list, 
#                               cov.model = "matern", 
#                               tuning = tuning.list, 
#                               n.omp.threads = 1, 
#                               verbose = FALSE, 
#                               NNGP = TRUE, 
#                               n.neighbors = 5, 
#                               search.type = 'cb', 
#                               n.report = 10, 
#                               n.burn = 100, 
# 		              n.thin = 2,
# 		              n.chains = 1))
#   data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
#   data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
#   expect_error(out <- spTPGOcc(occ.formula = occ.formula,
# 	                      det.formula = det.formula,
#                               data = data.list,
#                               inits = inits.list, 
#                               n.batch = n.batch, 
#                               batch.length = batch.length, 
#                               accept.rate = 0.43, 
#                               priors = prior.list, 
#                               cov.model = "matern", 
#                               tuning = tuning.list, 
#                               n.omp.threads = 1, 
#                               verbose = FALSE, 
#                               NNGP = TRUE, 
#                               n.neighbors = 5, 
#                               search.type = 'cb', 
#                               n.report = 10, 
#                               n.burn = 100, 
# 		              n.thin = 2,
# 		              n.chains = 1))
# })

# Check RE levels ---------------------
# test_that("random effect levels are correct", {
#   expect_equal(sort(unique(unlist(out$p.re.level.names))), 
# 	       sort(unique(c(X.p.re))))
#   expect_equal(sort(unique(unlist(out$re.level.names))), 
# 	       sort(unique(c(X.re))))
# })

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spTPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "spTPGOcc")
})

test_that("all correlation functions work", {
  out <- spTPGOcc(occ.formula = occ.formula, 
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
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "spTPGOcc")

  out <- spTPGOcc(occ.formula = occ.formula, 
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
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "spTPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spTPGOcc(occ.formula = occ.formula, 
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
	       n.burn = 100, 
	       n.chains = 1))
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs$trend[3, ] <- NA
  # expect_error(spTPGOcc(occ.formula = occ.formula, 
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
  #                      n.burn = 100, 
  #                      n.chains = 1))
  # tmp.data <- data.list
  # tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  # expect_error(spTPGOcc(occ.formula = occ.formula, 
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
  #                      n.burn = 100, 
  #                      n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- spTPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
                 n.chains = 1)
  expect_s3_class(out, "spTPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for spTPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spTPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spTPGOcc", {
  # X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  # dimnames(X.0.full)[[3]] <- c('(Intercept)', 'trend', 'occ.factor.1')
  X.0.full <- X.0
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
})
test_that("detection prediction works", {
  # X.p.0 <- abind::abind(dat$X.p[, , 1, ], dat$X.p.re[, , 1, ], along = 3)
  # dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  X.p.0 <- array(dat$X.p[, , 1, ], dim = c(J, n.time.max, p.det))
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
  J.fit <- nrow(X)
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Random intercepts on both with covariates -------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.8)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
psi.RE <- list(levels = c(10),
               sigma.sq.psi = c(0.8))
alpha <- c(-1, 0.5)
p.det <- length(alpha)
p.RE <- list(levels = c(20, 15),
             sigma.sq.p = c(1, 0.5))
phi <- 3 / .7
sigma.sq <- 2
nu <- 2
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep,
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
	       psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
               phi = phi, cov.model = 'matern', nu = nu)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , , drop = FALSE]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , , drop = FALSE]
X.p.re <- dat$X.p.re[-pred.indx, , , , drop = FALSE]
X.p.0 <- dat$X.p[pred.indx, , , , drop = FALSE]
X.p.re.0 <- dat$X.p.re[pred.indx, , , , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx]

occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.factor.1 = X.re[, , 1]) 
det.covs <- list(int = X.p[, , , 1], 
		 det.cov.1 = X.p[, , , 2],
                 det.factor.1 = X.p.re[, , , 1], 
                 det.factor.2 = X.p.re[, , , 2]) 

data.list <- list(y = y, coords = coords, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72),
		   nu.unif = c(0.5, 2.5))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)
# Tuning
tuning.list <- list(phi = 1, nu = 1)

batch.length <- 25
n.batch <- 10
n.report <- 10
n.samples <- batch.length * n.batch
occ.formula <- ~ trend + (1 | occ.factor.1)
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2) 

out <- spTPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = data.list,
	       inits = inits.list,
	       batch.length = batch.length,
	       n.batch = n.batch,
	       priors = prior.list,
	       accept.rate = 0.43,
	       cov.model = "matern",
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
test_that("out is of class spTPGOcc", {
  expect_s3_class(out, "spTPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- spTPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
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
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- spTPGOcc(occ.formula = occ.formula,
	                      det.formula = det.formula,
                              data = data.list,
                              inits = inits.list, 
                              n.batch = n.batch, 
                              batch.length = batch.length, 
                              accept.rate = 0.43, 
                              priors = prior.list, 
                              cov.model = "matern", 
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
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

test_that("default priors and inits work", {
  out <- spTPGOcc(occ.formula = occ.formula, 
	         det.formula = det.formula, 
	         data = data.list, 
	         n.batch = 40, 
	         batch.length = batch.length, 
	         cov.model = "exponential", 
	         tuning = tuning.list, 
	         NNGP = TRUE,
		 verbose = FALSE, 
	         n.neighbors = 5, 
	         search.type = 'cb', 
	         n.report = 10, 
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "spTPGOcc")
})

test_that("all correlation functions work", {
  out <- spTPGOcc(occ.formula = occ.formula, 
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
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "spTPGOcc")

  out <- spTPGOcc(occ.formula = occ.formula, 
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
	         n.burn = 100, 
	         n.chains = 1)
  expect_s3_class(out, "spTPGOcc")
})

test_that("verbose prints to the screen", {
  expect_output(spTPGOcc(occ.formula = occ.formula, 
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
  tmp.data$occ.covs$trend[3, ] <- NA
  expect_error(spTPGOcc(occ.formula = occ.formula, 
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
	               n.burn = 100, 
	               n.chains = 1))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(spTPGOcc(occ.formula = occ.formula, 
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
                       n.burn = 100, 
                       n.chains = 1))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- spTPGOcc(occ.formula = occ.formula, 
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
                 n.burn = 100, 
                 n.chains = 1)
  expect_s3_class(out, "spTPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for spTPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spTPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for spTPGOcc", {
  X.0.full <- abind::abind(X.0, X.re.0, along = 3)
  dimnames(X.0.full)[[3]] <- c('(Intercept)', 'trend', 'occ.factor.1')
  pred.out <- predict(out, X.0.full, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, nrow(X.0.full), n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind::abind(dat$X.p[, , 1, ], dat$X.p.re[, , 1, ], along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
  J.fit <- nrow(X)
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

