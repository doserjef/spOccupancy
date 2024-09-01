# Test tMsPGOcc.R ---------------------------------------------------------

skip_on_cran()

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
N <- 6
beta.mean <- c(0.4)
tau.sq.beta <- c(1)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list()
alpha.mean <- c(-1)
tau.sq.alpha <- c(0.5)
p.det <- length(alpha.mean)
p.RE <- list()
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1])
det.covs <- list(int = X.p[, , , 1])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                   rho.unif = list(a = -1, b = 1), 
                   sigma.sq.t.ig = list(a = 2, b = 1)) 
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ 1
det.formula <- ~ 1

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- tMsPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
               n.thin = 13,
               n.batch = n.batch, 
               batch.length = batch.length, 
               accept.rate = 0.43, 
	       n.omp.threads = 1,
	       verbose = FALSE))
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
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
          parallel.chains = TRUE,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- array(1, dim = c(J.str, 1, p.det))
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str, 1))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Occurrence covariate only -----------------------------------------------
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
N <- 6
beta.mean <- c(0.4, 0.5, -0.2)
tau.sq.beta <- c(1, 0.5, 1)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list()
alpha.mean <- c(-1)
tau.sq.alpha <- c(0.5)
p.det <- length(alpha.mean)
p.RE <- list()
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 occ.cov.1 = X[, , 2],
                 occ.cov.2 = X[, , 3])
det.covs <- list(int = X.p[, , , 1])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                   rho.unif = list(a = -1, b = 1), 
                   sigma.sq.t.ig = list(a = 2, b = 1)) 
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ 1

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
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
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- array(1, dim = c(J.str, 1, p.det))
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str, 1))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Detection covariate only ------------------------------------------------
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
N <- 6
beta.mean <- c(0.4)
tau.sq.beta <- c(1)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list()
alpha.mean <- c(-1, 0.2, 0.4)
tau.sq.alpha <- c(0.5, 1, 1.2)
p.det <- length(alpha.mean)
p.RE <- list()
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 1], 
                 det.cov.2 = X.p[, , , 2])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                   rho.unif = list(a = -1, b = 1), 
                   sigma.sq.t.ig = list(a = 2, b = 1)) 
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ 1
det.formula <- ~ det.cov.1 + det.cov.2

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
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
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, , 1, ]
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Covariates on both ------------------------------------------------------
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
N <- 6
beta.mean <- c(0.4, 0.3)
tau.sq.beta <- c(1, 0.5)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list()
alpha.mean <- c(-1, 0.2, 0.4)
tau.sq.alpha <- c(0.5, 1, 1.2)
p.det <- length(alpha.mean)
p.RE <- list()
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 occ.cov.1 = X[, , 2])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 1], 
                 det.cov.2 = X.p[, , , 2])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                   rho.unif = list(a = -1, b = 1), 
                   sigma.sq.t.ig = list(a = 2, b = 1)) 
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ occ.cov.1
det.formula <- ~ det.cov.1 + det.cov.2

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
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
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, , 1, ]
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Interactions on both ----------------------------------------------------
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
N <- 6
beta.mean <- c(0.4, 0.3, -0.3)
tau.sq.beta <- c(1, 0.5, 1.2)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list()
alpha.mean <- c(-1, 0.2, 0.4)
tau.sq.alpha <- c(0.5, 1, 1.2)
p.det <- length(alpha.mean)
p.RE <- list()
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 occ.cov.1 = X[, , 2], 
                 occ.cov.2 = X[, , 3])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X.p[, , , 1], 
                 det.cov.2 = X.p[, , , 2])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                   rho.unif = list(a = -1, b = 1), 
                   sigma.sq.t.ig = list(a = 2, b = 1)) 
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- ~ det.cov.1 * det.cov.2

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
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
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- dat$X
  X.0 <- abind(X.0, X.0[, , 2] * X.0[, , 3], along = 3)
  dimnames(X.0)[[3]] <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.cov.1:occ.cov.2')
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, , 1, ]
  X.p.0 <- abind(X.p.0, X.p.0[, , 2] * X.p.0[, , 3], along = 3)
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Site covariate on detection ---------------------------------------------
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
N <- 6
beta.mean <- c(0.4, 0.5, -0.2)
tau.sq.beta <- c(1, 0.5, 1)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list()
alpha.mean <- c(-1)
tau.sq.alpha <- c(0.5)
p.det <- length(alpha.mean)
p.RE <- list()
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 occ.cov.1 = X[, , 2],
                 occ.cov.2 = X[, , 3])
det.covs <- list(int = X.p[, , , 1], 
                 det.cov.1 = X[, , 2])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                   rho.unif = list(a = -1, b = 1), 
                   sigma.sq.t.ig = list(a = 2, b = 1)) 
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ 1
det.formula <- ~ det.cov.1

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
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
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- dat$X[, , 1, drop = FALSE]
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- abind(dat$X[, , 1], dat$X[, , 2], along = 3)
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Random intercept on occurrence -------------------------------------------
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
N <- 6
beta.mean <- c(0.4)
tau.sq.beta <- c(1)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list(levels = c(45), 
               sigma.sq.psi = c(1.3))
alpha.mean <- c(-1)
tau.sq.alpha <- c(0.5)
p.det <- length(alpha.mean)
p.RE <- list()
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
occ.covs <- list(int = X[, , 1], 
                 occ.factor.1 = X.re[, , 1])
det.covs <- list(int = X.p[, , , 1])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   rho.unif = list(a = -1, b = 1),
                   sigma.sq.t.ig = list(a = 2, b = 1))
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ (1 | occ.factor.1)
det.formula <- ~ 1

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- tMsPGOcc(occ.formula = occ.formula,
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
  expect_error(out <- tMsPGOcc(occ.formula = occ.formula,
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

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- abind(dat$X, dat$X.re, along = 3)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'occ.factor.1')
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- array(1, dim = c(J.str, 1, p.det))
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J.str, 1))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})

# Multiple random intercepts on occurrence --------------------------------
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
N <- 6
beta.mean <- c(0.4, 0.3)
tau.sq.beta <- c(1, 0.5)
p.occ <- length(beta.mean)
trend <- FALSE
sp.only <- 0
psi.RE <- list(levels = c(20),
               sigma.sq.psi = c(2.5))
p.RE <- list(levels = c(50, 10),
	     sigma.sq.p = c(2.50, 1.5))
alpha.mean <- c(-1, 0.2, 0.4)
tau.sq.alpha <- c(0.5, 1, 1.2)
p.det <- length(alpha.mean)
sp <- FALSE
factor.model <- FALSE
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.1, 5)
rho <- runif(N, 0, 1)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
X.re <- dat$X.re
X.p.re <- dat$X.p.re
occ.covs <- list(int = X[, , 1],
                 occ.cov.1 = X[, , 2], 
                 occ.factor.1 = X.re[, , 1])
det.covs <- list(int = X.p[, , , 1],
                 det.cov.1 = X.p[, , , 1],
                 det.cov.2 = X.p[, , , 2], 
                 det.factor.1 = X.p.re[, , , 1],
                 det.factor.2 = X.p.re[, , , 2])

data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   rho.unif = list(a = -1, b = 1),
                   sigma.sq.t.ig = list(a = 2, b = 1))
# Starting values
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   sigma.sq.p = 0.5,
		   z = z.init)

n.batch <- 10
batch.length <- 25
n.samples <- n.batch * batch.length
n.report <- 10
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1)
det.formula <- ~ det.cov.1 + det.cov.2 + (1 | det.factor.1) + (1 | det.factor.2)

out <- tMsPGOcc(occ.formula = occ.formula,
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
	        n.burn = 100,
	        n.thin = 2,
	        n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class tMsPGOcc", {
  expect_s3_class(out, "tMsPGOcc")
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
  expect_error(out <- tMsPGOcc(occ.formula = occ.formula,
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
	                    n.burn = 0,
	                    n.thin = 1,
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- tMsPGOcc(occ.formula = occ.formula,
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
	                    n.burn = 0,
	                    n.thin = 1,
	                    n.chains = 1))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- tMsPGOcc(occ.formula = occ.formula,
	        det.formula = det.formula,
	        data = data.list,
	        n.batch = n.batch,
		ar1 = FALSE,
	        batch.length = batch.length,
	        n.omp.threads = 1,
	        verbose = FALSE)
  expect_s3_class(out, "tMsPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(tMsPGOcc(occ.formula = occ.formula,
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
test_that("waicOCC works for tMsPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})
test_that("waicOCC works for multiple species", {
  # as.vector gets rid of names
  waic.out <- waicOcc(out, by.sp = TRUE)
  expect_equal(nrow(waic.out), N)
})

# Check fitted ------------------------
test_that("fitted works for tMsPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for tMsPGOcc", {
  X.0 <- abind(dat$X, dat$X.re, along = 3)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'occ.cov.1', 'occ.factor.1')
  pred.out <- predict(out, X.0, t.cols = 1:n.time.max)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})
test_that("detection prediction works", {
  X.p.0 <- abind(dat$X.p[, , 1, ], dat$X.p.re[, , 1, ], along = 3)
  dimnames(X.p.0)[[3]] <- c('(Intercept)', 'det.cov.1', 'det.cov.2', 
			    'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, t.cols = 1:n.time.max, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, N, J, n.time.max))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for tMsPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))

  ppc.out <- ppcOcc(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, J, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, J, n.time.max))

  ppc.out <- ppcOcc(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.rep), c(n.post.samples, N, n.time.max))
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, N, n.time.max, max(n.rep, na.rm = TRUE)))
})


