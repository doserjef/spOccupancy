# Test tPGOcc.R -----------------------------------------------------------

# skip_on_cran()

# Intercept Only ----------------------------------------------------------
set.seed(100)
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1)
p.det <- length(alpha)
p.RE <- list()
ar1 <- TRUE
rho <- 0.9
sigma.sq.t <- 2.4
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE, ar1 = TRUE, rho = rho, sigma.sq.t = sigma.sq.t)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1])
det.covs <- list(int = X.p[, , , 1])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

out <- tPGOcc(occ.formula = ~ 1,
	      det.formula = ~ 1,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      ar1 = TRUE,
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = ~ 1, 
	        det.formula = ~ 1, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = ~ 1, 
	       det.formula = ~ 1, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- array(1, dim = c(J.str, 1, p.det))
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, 1))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Occurrence covariate only -----------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.5)
p.occ <- length(beta)
trend <- TRUE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2])
det.covs <- list(int = X.p[, , , 1])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 100

occ.formula <- ~ trend
det.formula <- ~ 1
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- array(1, dim = c(J.str, 1, p.det))
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, 1))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Detection covariate only -----------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1, 0.5)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 2])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ 1 
det.formula <- ~ det.cov.1
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       ar1 = TRUE,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- dat$X.p[, , 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Covariates on both ------------------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.6, 0.5)
p.occ <- length(beta)
trend <- TRUE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1, 0.5, 0.3, -0.4)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.cov.1 = X[, , 3])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 2], 
                 det.cov.2 = X.p[, , , 3], 
                 det.cov.3 = X.p[, , , 4])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ trend + occ.cov.1
det.formula <- ~ det.cov.1 + det.cov.1 + det.cov.2 + det.cov.3
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$trend[3, ] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- dat$X.p[, , 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Intercations in both ----------------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.6, 0.5)
p.occ <- length(beta)
trend <- TRUE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1, 0.5, 0.3, -0.4)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.cov.1 = X[, , 3])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 2], 
                 det.cov.2 = X.p[, , , 3], 
                 det.cov.3 = X.p[, , , 4])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ trend * occ.cov.1
det.formula <- ~ det.cov.1 + det.cov.1 + det.cov.2 * det.cov.3
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      ar1 = TRUE,
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$trend[3, ] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       ar1 = TRUE,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- out$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- abind::abind(dat$X.p[, , 1, ], dat$X.p[, , 1, 3] * dat$X.p[, , 1, 4], along = 3)
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Site/year covariate on detection ----------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.6, 0.5)
p.occ <- length(beta)
trend <- TRUE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1, 0.5, 0.3, -0.4)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.cov.1 = X[, , 3])
det.covs <- list(occ.cov.1 = X[, , 3])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ trend + occ.cov.1
det.formula <- ~ occ.cov.1 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       ar1 = TRUE,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$trend[3, ] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$occ.cov.1[1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- dat$X.p[, , 1, 1:2]
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Site covariate on occurrence/detection ----------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.6)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 2
psi.RE <- list()
alpha <- c(-1, 0.5, 0.3, -0.4)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 occ.cov.1 = X[, , 2]) 
det.covs <- list(occ.cov.1 = X[, , 2])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ occ.cov.1
det.formula <- ~ occ.cov.1 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$occ.cov.1[1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$occ.cov.1[1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- dat$X.p[, , 1, 1:2]
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Random intercept on occurrence ------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.6, 0.5)
p.occ <- length(beta)
trend <- TRUE 
sp.only <- 0
psi.RE <- list(levels = c(10), 
               sigma.sq.psi = c(1))
alpha <- c(-1, 0.5, 0.3, -0.4)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.cov.1 = X[, , 3], 
                 occ.factor.1 = X.re[, , 1])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 2], 
                 det.cov.2 = X.p[, , , 3], 
                 det.cov.3 = X.p[, , , 4])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ trend + occ.cov.1 + (1 | occ.factor.1)
det.formula <- ~ det.cov.1 + det.cov.1 + det.cov.2 + det.cov.3
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))

})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$occ.cov.1[1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- abind(dat$X, dat$X.re, along = 3)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'trend', 'occ.cov.1', 'occ.factor.1')
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- dat$X.p[, , 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Multiple random intercepts on occurrence --------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.6, 0.5)
p.occ <- length(beta)
trend <- TRUE 
sp.only <- 0
psi.RE <- list(levels = c(10, 20), 
               sigma.sq.psi = c(1, 0.5))
alpha <- c(-1, 0.5, 0.3, -0.4)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.cov.1 = X[, , 3], 
                 occ.factor.1 = X.re[, , 1], 
                 occ.factor.2 = X.re[, , 2])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 2], 
                 det.cov.2 = X.p[, , , 3], 
                 det.cov.3 = X.p[, , , 4])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ trend + occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ det.cov.1 + det.cov.1 + det.cov.2 + det.cov.3
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))

})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$occ.cov.1[1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- abind(dat$X, dat$X.re, along = 3)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'trend', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- dat$X.p[, , 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Occurrence REs only -----------------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 0
psi.RE <- list(levels = c(10, 20), 
               sigma.sq.psi = c(1, 0.5))
alpha <- c(-1)
p.det <- length(alpha)
p.RE <- list()
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 occ.factor.1 = X.re[, , 1], 
                 occ.factor.2 = X.re[, , 2])
det.covs <- list(int = X.p[, , , 1]) 

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))

})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- abind(dat$X, dat$X.re, along = 3)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'occ.factor.1', 'occ.factor.2')
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- array(dat$X.p[, , 1, ], dim = c(J, n.time.max, 1))
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Random intercept on detection -------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1)
p.det <- length(alpha)
p.RE <- list(levels = c(20), 
             sigma.sq.p = c(1))
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
occ.covs <- list(int = X[, , 1]) 
det.covs <- list(int = X.p[, , , 1], 
                 det.factor.1 = X.p.re[, , , 1]) 

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ 1 
det.formula <- ~ (1 | det.factor.1) 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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
  expect_equal(out$psiRE, FALSE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind(array(dat$X.p[, , 1, ], dim = c(J, n.time.max, 1)), 
		 array(dat$X.p.re[, , 1, ], dim = c(J, n.time.max, 1)), along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.factor.1')  
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Multiple random intercepts on detection ---------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1)
p.det <- length(alpha)
p.RE <- list(levels = c(20, 15), 
             sigma.sq.p = c(1, 0.5))
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
occ.covs <- list(int = X[, , 1]) 
det.covs <- list(int = X.p[, , , 1], 
                 det.factor.1 = X.p.re[, , , 1], 
                 det.factor.2 = X.p.re[, , , 2]) 

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ 1 
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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
  expect_equal(out$psiRE, FALSE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       ar1 = TRUE,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind::abind(array(dat$X.p[, , 1, ], dim = c(J, n.time.max, 1)), 
		 array(dat$X.p.re[, , 1, ], dim = c(J, n.time.max, 2)), along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.factor.1', 'det.factor.2')  
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Multiple random intercepts with covariates ------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1, 0.5)
p.det <- length(alpha)
p.RE <- list(levels = c(20, 15), 
             sigma.sq.p = c(1, 0.5))
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
occ.covs <- list(int = X[, , 1]) 
det.covs <- list(int = X.p[, , , 1], 
		 det.cov.1 = X.p[, , , 2],
                 det.factor.1 = X.p.re[, , , 1], 
                 det.factor.2 = X.p.re[, , , 2]) 

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ 1 
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2) 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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
  expect_equal(out$psiRE, FALSE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
		ar1 = TRUE,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs$trend[3, ] <- NA
#   expect_error(tPGOcc(occ.formula = occ.formula,
# 	       det.formula = det.formula,
# 	       data = tmp.data,
#	       n.batch = n.batch, 
#	       batch.length = batch.length,
#	       tuning = list(rho = 1),
# 	       n.omp.threads = 1,
# 	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind::abind(dat$X.p[, , 1, ], dat$X.p.re[, , 1, ], along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')  
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Multiple random intercepts with covariates on both ----------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4, 0.8)
p.occ <- length(beta)
trend <- TRUE 
sp.only <- 0
psi.RE <- list()
alpha <- c(-1, 0.5)
p.det <- length(alpha)
p.RE <- list(levels = c(20, 15), 
             sigma.sq.p = c(1, 0.5))
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2]) 
det.covs <- list(int = X.p[, , , 1], 
		 det.cov.1 = X.p[, , , 2],
                 det.factor.1 = X.p.re[, , , 1], 
                 det.factor.2 = X.p.re[, , , 2]) 

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ trend 
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2) 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      ar1 = TRUE,
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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
  expect_equal(out$psiRE, FALSE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       ar1 = TRUE,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$trend[3, ] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind::abind(dat$X.p[, , 1, ], dat$X.p.re[, , 1, ], along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')  
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Random intercepts on both -----------------------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
n.time.max <- max(n.time)
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
beta <- c(0.4)
p.occ <- length(beta)
trend <- FALSE 
sp.only <- 0
psi.RE <- list(levels = c(10), 
               sigma.sq.psi = c(0.8))
alpha <- c(-1)
p.det <- length(alpha)
p.RE <- list(levels = c(20, 15), 
             sigma.sq.p = c(1, 0.5))
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
occ.covs <- list(int = X[, , 1], 
                 occ.factor.1 = X.re[, , 1]) 
det.covs <- list(int = X.p[, , , 1], 
                 det.factor.1 = X.p.re[, , , 1], 
                 det.factor.2 = X.p.re[, , , 2]) 

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ (1 | occ.factor.1)
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      ar1 = TRUE,
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  # data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  # data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       ar1 = TRUE,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, ignore.RE = TRUE, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind::abind(array(dat$X.p[, , 1, ], dim = c(J, n.time.max, 1)), 
		 array(dat$X.p.re[, , 1, ], dim = c(J, n.time.max, 2)), along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.factor.1', 'det.factor.2')  
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Random intercepts and covariates on both --------------------------------
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.time <- sample(2:10, J, replace = TRUE)
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
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
	       beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
	       psi.RE = psi.RE, p.RE = p.RE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.factor.1 = X.re[, , 1]) 
det.covs <- list(int = X.p[, , , 1], 
		 det.cov.1 = X.p[, , , 2],
                 det.factor.1 = X.p.re[, , , 1], 
                 det.factor.2 = X.p.re[, , , 2]) 

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha = 0, beta = 0, z = z.init)

n.batch <- 40
batch.length <- 25
n.samples <- n.batch * batch.length 
n.report <- 100

occ.formula <- ~ trend + (1 | occ.factor.1)
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2) 
out <- tPGOcc(occ.formula = occ.formula,
	      det.formula = det.formula,
	      data = data.list,
	      inits = inits.list,
	      n.batch = n.batch, 
	      batch.length = batch.length,
	      tuning = list(rho = 1),
	      priors = prior.list,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = n.report,
	      n.burn = 400,
	      n.thin = 6, 
	      n.chains = 2,
	      k.fold = 2, 
	      k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class tPGOcc", {
  expect_s3_class(out, "tPGOcc")
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))),
	       sort(unique(c(X.p.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
			    ar1 = TRUE,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- tPGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.batch = n.batch, 
	                    batch.length = batch.length,
	                    tuning = list(rho = 1),
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tPGOcc(occ.formula = occ.formula, 
	        det.formula = det.formula, 
	        data = data.list, 
	        n.batch = n.batch, 
	        batch.length = batch.length,
	        tuning = list(rho = 1),
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs$trend[3, ] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$det.covs$det.cov.1[1, 1, 1] <- NA
  expect_error(tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1, 1] <- NA
  out <- tPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.batch = n.batch, 
	       batch.length = batch.length,
	       tuning = list(rho = 1),
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "tPGOcc")
})

# Check waicOcc -----------------------
test_that("waicOCC works for tPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for tPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tPGOcc", {
  X.0 <- abind::abind(dat$X, dat$X.re, along = 3)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'trend', 'occ.factor.1')
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind::abind(dat$X.p[, , 1, ], dat$X.p.re[, , 1, ], along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')  
  pred.out <- predict(out, X.p.0, type = 'detection', t.cols = 1:n.time.max)
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tPGOcc", {
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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, n.time.max))
  
  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.time.max, max(n.rep, na.rm = TRUE)))
})
