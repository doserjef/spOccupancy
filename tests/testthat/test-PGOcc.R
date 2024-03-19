# Test PGOcc.R  ----------------------------------------------------------

skip_on_cran()

# Intercept Only ----------------------------------------------------------
set.seed(100)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

data.list <- list(y = y)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE), 
                   fix = TRUE)

n.samples <- 1000
n.report <- 100

out <- PGOcc(occ.formula = ~ 1,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
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
  out <- PGOcc(occ.formula = ~ 1, 
	       det.formula = ~ 1, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- PGOcc(occ.formula = ~ 1, 
	       det.formula = ~ 1, 
	       data = data.list, 
               n.thin = 13,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
})
# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = ~ 1, 
	       det.formula = ~ 1, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Occurrence covariate only -----------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

colnames(X) <- c('int', 'occ.cov.1')

data.list <- list(y = y, occ.covs = X)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100

out <- PGOcc(occ.formula = ~ occ.cov.1,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
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
  out <- PGOcc(occ.formula = ~ occ.cov.1, 
	       det.formula = ~ 1, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = ~ occ.cov.1, 
	       det.formula = ~ 1, 
	       data = tmp.data, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = ~ occ.cov.1, 
	       det.formula = ~ 1, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Detection covariate only ------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5, 1)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

data.list <- list(y = y, det.covs = list(det.cov.1 = X.p[, , 2]))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ det.cov.1

out <- PGOcc(occ.formula = ~ 1,
	     det.formula = ~ det.cov.1,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
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
  out <- PGOcc(occ.formula = ~ 1, 
	       det.formula = ~ det.cov.1, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs$det.cov.1[1, 1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = tmp.data, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = tmp.data, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc") 
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = ~ 1, 
	       det.formula = ~ det.cov.1, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Covariates on both ------------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.2, -1.0)
p.occ <- length(beta)
alpha <- c(-0.5, 1)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

colnames(X) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, det.covs = list(det.cov.1 = X.p[, , 2]), 
                  occ.covs = X)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ det.cov.1

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs$det.cov.1[1, 1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Interactions in both ----------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.2, -1.0)
p.occ <- length(beta)
alpha <- c(-0.5, 1, 0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

colnames(X) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, det.covs = list(det.cov.1 = X.p[, , 2], 
					 det.cov.2 = X.p[, , 3]), 
                  occ.covs = X)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- ~ det.cov.1 * det.cov.2

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$det.covs$det.cov.1[1, 1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  X.0 <- dat$X
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  X.p.0 <- cbind(X.p.0, X.p.0[, 2] * X.p.0[, 3])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Site covariate on detection ---------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.2, -1.0)
p.occ <- length(beta)
alpha <- c(-0.5, 1)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

colnames(X) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, det.covs = list(det.cov.1 = X[, 2]), 
                  occ.covs = X)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ det.cov.1

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs$det.cov.1[1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Random intercept on occurrence ------------------------------------------
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
X.re <- dat$X.re + 15
X.p <- dat$X.p

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1')
data.list <- list(y = y, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1)
det.formula <- ~ 1

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)		 
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$det.covs$det.cov.1[1] <- NA
  # expect_error(PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE)
  # expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- matrix(dat$X.p[, 1, ], ncol = 1)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})


# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Multiple random intercepts on occurrence --------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(45, 15), 
	       sigma.sq.psi = c(1.3, 0.5))
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p <- dat$X.p

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)		 
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$det.covs$det.cov.1[1] <- NA
  # expect_error(PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE)
  # expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- matrix(dat$X.p[, 1, ], ncol = 1)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Occurrence REs + covariates ---------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.7)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(45, 15), 
	       sigma.sq.psi = c(1.3, 0.5))
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p <- dat$X.p

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, occ.covs = occ.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ 1

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)		 
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$det.covs$det.cov.1[1] <- NA
  # expect_error(PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE)
  # expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- matrix(dat$X.p[, 1, ], ncol = 1)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Occurrence REs + covariates in all --------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.7)
p.occ <- length(beta)
alpha <- c(-0.5, 0.8, 0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(45, 15), 
	       sigma.sq.psi = c(1.3, 0.5))
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p <- dat$X.p

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
det.formula <- ~ det.cov.1 + det.cov.2

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)		 
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs$det.cov.1[1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Random intercept on detection ------------------------------------------
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
	     sigma.sq.p = c(2.50))
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p.re <- dat$X.p.re + 12
X.p <- dat$X.p

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
det.covs <- list(det.factor.1 = X.p.re[, , 1])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1)
data <- data.list
inits <- inits.list
priors <- prior.list
n.omp.threads <- 1
verbose <- FALSE
n.burn <- 400
n.thin <- 6
n.chains <- 2
k.fold <- 2
k.fold.threads <- 1
k.fold.only <- FALSE
k.fold.seed <- 100

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs[3, ] <- NA
#   expect_error(PGOcc(occ.formula = occ.formula,
# 	       det.formula = det.formula,
# 	       data = tmp.data,
# 	       n.samples = n.samples,
# 	       n.omp.threads = 1,
# 	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1])
  colnames(X.p.0) <- c('intercept', 'det.factor.1')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Multiple random intercepts on detection ---------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(100, 52), 
	     sigma.sq.p = c(4, 0.3))
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p.re <- dat$X.p.re + 12
X.p <- dat$X.p

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2)

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)
# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs[3, ] <- NA
#   expect_error(PGOcc(occ.formula = occ.formula,
# 	       det.formula = det.formula,
# 	       data = tmp.data,
# 	       n.samples = n.samples,
# 	       n.omp.threads = 1,
# 	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Multiple random intercepts with covariate -------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5, 1.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(100, 52), 
	     sigma.sq.p = c(4, 0.3))
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p.re <- dat$X.p.re + 12
X.p <- dat$X.p

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.cov.1 = X.p[, , 2])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)
# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
#   tmp.data <- data.list
#   tmp.data$occ.covs[3, ] <- NA
#   expect_error(PGOcc(occ.formula = occ.formula,
# 	       det.formula = det.formula,
# 	       data = tmp.data,
# 	       n.samples = n.samples,
# 	       n.omp.threads = 1,
# 	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Multiple random intercepts with covariate on both -----------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, -1)
p.occ <- length(beta)
alpha <- c(-0.5, 1.5)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list(levels = c(100, 52), 
	     sigma.sq.p = c(4, 0.3))
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p.re <- dat$X.p.re + 12
X.p <- dat$X.p

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.cov.1 = X.p[, , 2])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)
# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
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
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6, 
	                    n.chains = 1))
  # data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:2])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Random intercepts on both -----------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
alpha <- c(-0.5)
p.det <- length(alpha)
psi.RE <- list(levels = c(45, 30, 15), 
	       sigma.sq.psi = c(1.3, 2.5, 0.5))
p.RE <- list(levels = c(100, 52, 15), 
	     sigma.sq.p = c(4, 0.3, 1.5))
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p.re <- dat$X.p.re + 12
X.p <- dat$X.p

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1', 'occ.factor.2', 'occ.factor.3')
det.covs <- list(det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.factor.3 = X.p.re[, , 3])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ (1 | occ.factor.1) + (1 | occ.factor.2) + (1 | occ.factor.3)
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3)

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)
# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.factor.1', 'det.factor.2', 'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Random intercepts on both with covariates -------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 1.2)
p.occ <- length(beta)
alpha <- c(-0.5, 0.4)
p.det <- length(alpha)
psi.RE <- list(levels = c(45, 30, 15), 
	       sigma.sq.psi = c(1.3, 2.5, 0.5))
p.RE <- list(levels = c(100, 52, 15), 
	     sigma.sq.p = c(4, 0.3, 1.5))
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.p.re <- dat$X.p.re + 12
X.p <- dat$X.p

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2', 'occ.factor.3')
det.covs <- list(det.cov.1 = X.p.re[, , 2], 
		 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2], 
                 det.factor.3 = X.p.re[, , 3])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2) + (1 | occ.factor.3)
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2) + (1 | det.factor.3)

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)
# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, TRUE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.p.re))))
  expect_equal(sort(unique(unlist(out$re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
	                    priors = prior.list,
	                    n.omp.threads = 1,
	                    verbose = FALSE,
	                    n.report = n.report,
	                    n.burn = 400,
	                    n.thin = 6,
	                    n.chains = 1))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- PGOcc(occ.formula = occ.formula,
	                    det.formula = det.formula,
	                    data = data.list,
	                    inits = inits.list,
	                    n.samples = n.samples,
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
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1 
  tmp.data$det.covs[[1]][1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1:3])
  colnames(X.p.0) <- c('intercept', 'det.cov.1', 'det.factor.1', 'det.factor.2', 'det.factor.3')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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

# Dimension of y is not equal to max(n.rep) -------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
n.rep.max <- 7
beta <- c(0.3, 0.2, -1.0)
p.occ <- length(beta)
alpha <- c(-0.5, 1)
p.det <- length(alpha)
psi.RE <- list()
p.RE <- list()
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, n.rep.max = n.rep.max)
y <- dat$y
X <- dat$X
X.p <- dat$X.p

colnames(X) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, det.covs = list(det.cov.1 = X.p[, , 2]),
                  occ.covs = X)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- ~ det.cov.1

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2,
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
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
  out <- PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = data.list,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$det.covs$det.cov.1[1, 1] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  out <- PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = data.list,
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report,
	       n.burn = 1,
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcOcc(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.rep.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.rep.max))

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
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.rep.max))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.rep.max))
})

# Random site-level intercept on detection --------------------------------
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
X.re <- dat$X.re + 15
X.p <- dat$X.p

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int')
data.list <- list(y = y, occ.covs = occ.covs,
                  det.covs = list(det.factor.1 = X.re[, 1]))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   alpha.normal = list(mean = 0, var = 2.72))
# Starting values
inits.list <- list(alpha = 0, beta = 0, z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 1000
n.report <- 100
occ.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1)

out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
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
	     k.fold = 2, 
	     k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class PGOcc", {
  expect_s3_class(out, "PGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), 1)
  expect_type(out$k.fold.deviance, "double")
  expect_gt(out$k.fold.deviance, 0)
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$psiRE, FALSE)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$p.re.level.names))), 
	       sort(unique(c(X.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE)
  expect_s3_class(out, "PGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(PGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = tmp.data,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$det.covs$det.cov.1[1] <- NA
  # expect_error(PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$y[1, 1] <- NA
  # out <- PGOcc(occ.formula = occ.formula,
  #              det.formula = det.formula,
  #              data = tmp.data,
  #              n.samples = n.samples,
  #              n.omp.threads = 1,
  #              verbose = FALSE)
  # expect_s3_class(out, "PGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {

  expect_output(PGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
	       n.samples = 100,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = n.report, 
	       n.burn = 1, 
	       n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for PGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))  
})

# Check fitted ------------------------
test_that("fitted works for PGOcc", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for PGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.re[, 1])
  colnames(X.p.0) <- c('(Intercept)', 'det.factor.1')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for PGOcc", {
  n.post.samples <- out$n.post * out$n.chains
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
