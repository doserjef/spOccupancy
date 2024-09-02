# Test intPGOcc.R  ----------------------------------------------------------

skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(111)
J.x <- 15
J.y <- 15
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
for (i in 1:n.data) {
  alpha[[i]] <- runif(1, -1, 1)
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
occ.covs <- X
colnames(occ.covs) <- c('int')
det.covs <- list()
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  sites = sites)

J <- length(dat$z.obs)
# Initial values
inits.list <- list(alpha = list(0, 0, 0, 0), 
		   beta = 0, 
		   fix = TRUE,
		   z = rep(1, J))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
		   alpha.normal = list(mean = list(0, 0, 0, 0), 
			               var = list(2.72, 2.72, 2.72, 2.72)))
n.samples <- 1000
occ.formula <- ~ 1
det.formula <- list(f.1 = ~ 1, f.2 = ~ 1, f.3 = ~ 1, f.4 = ~ 1)
out <- intPGOcc(occ.formula = occ.formula, 
                det.formula = det.formula, 
		data = data.list,
		inits = inits.list,
		n.samples = n.samples, 
		priors = prior.list, 
		n.omp.threads = 1, 
		verbose = FALSE, 
		n.report = 100, 
		n.burn = 500, 
		n.thin = 2, 
		n.chains = 2, 
                k.fold = 2, 
                k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class intPGOcc", {
  expect_s3_class(out, "intPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- intPGOcc(occ.formula = occ.formula, 
                  det.formula = det.formula, 
		  data = data.list,
		  inits = inits.list,
		  n.samples = n.samples, 
		  priors = prior.list, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.report = 100, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check default priors ---------------- 
test_that("default priors, inits, burn, thin, work", {

  out <- intPGOcc(occ.formula = occ.formula, 
		  det.formula = det.formula,
		  data = data.list,
		  n.samples = 100, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.chains = 1)
  expect_s3_class(out, "intPGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- intPGOcc(occ.formula = occ.formula, 
	       det.formula = det.formula, 
	       data = data.list, 
               n.thin = 13,
	       n.samples = n.samples,
	       n.omp.threads = 1,
	       verbose = FALSE))
})

test_that("verbose prints to the screen", {

  expect_output(out <- intPGOcc(occ.formula = occ.formula, 
                                det.formula = det.formula,  
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

# Check waicOcc -----------------------
test_that("waicOCC works for intPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for intPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for intPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for intPGOcc", {
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
J.x <- 15
J.y <- 15
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
beta <- c(0.5, 1.2, -0.8)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(1, -1, 1)
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
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
data.list <- list(y = y, 
		  occ.covs = occ.covs,
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
n.samples <- 1000
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- list(f.1 = ~ 1, f.2 = ~ 1, f.3 = ~ 1, f.4 = ~ 1)
out <- intPGOcc(occ.formula = occ.formula, 
                det.formula = det.formula, 
		data = data.list,
		inits = inits.list,
		n.samples = n.samples, 
		priors = prior.list, 
		n.omp.threads = 1, 
		verbose = FALSE, 
		n.report = 100, 
		n.burn = 500, 
		n.thin = 2, 
		n.chains = 2, 
                k.fold = 2, 
                k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class intPGOcc", {
  expect_s3_class(out, "intPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- intPGOcc(occ.formula = occ.formula, 
                  det.formula = det.formula, 
		  data = data.list,
		  inits = inits.list,
		  n.samples = n.samples, 
		  priors = prior.list, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.report = 100, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 1, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
                        det.formula = det.formula,
                        data = tmp.data,
                        n.samples = n.samples,
                        n.omp.threads = 1,
                        verbose = FALSE))
  # tmp.data <- data.list
  # tmp.data$det.covs[[1]][1] <- NA
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

# Check default priors ---------------- 
test_that("default priors, inits, burn, thin, work", {

  out <- intPGOcc(occ.formula = occ.formula, 
		  det.formula = det.formula,
		  data = data.list,
		  n.samples = 100, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.chains = 1)
  expect_s3_class(out, "intPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(out <- intPGOcc(occ.formula = occ.formula, 
                                det.formula = det.formula,  
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

# Check waicOcc -----------------------
test_that("waicOCC works for intPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for intPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for intPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for intPGOcc", {
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

# Detection covariate only -----------------------------------------------
J.x <- 15
J.y <- 15
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
occ.covs <- X
colnames(occ.covs) <- c('int')
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
n.samples <- 1000
occ.formula <- ~ 1 
det.formula <- list(f.1 = ~ det.cov.1.1, 
		    f.2 = ~ det.cov.2.1, 
		    f.3 = ~ det.cov.3.1, 
		    f.4 = ~ det.cov.4.1)
out <- intPGOcc(occ.formula = occ.formula, 
                det.formula = det.formula, 
		data = data.list,
		inits = inits.list,
		n.samples = n.samples, 
		priors = prior.list, 
		n.omp.threads = 1, 
		verbose = FALSE, 
		n.report = 100, 
		n.burn = 500, 
		n.thin = 2, 
		n.chains = 1, 
                k.fold = 2, 
                k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class intPGOcc", {
  expect_s3_class(out, "intPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- intPGOcc(occ.formula = occ.formula, 
                  det.formula = det.formula, 
		  data = data.list,
		  inits = inits.list,
		  n.samples = n.samples, 
		  priors = prior.list, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.report = 100, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  # tmp.data <- data.list
  # tmp.data$occ.covs[3, ] <- NA
  # expect_error(intPGOcc(occ.formula = occ.formula,
  #                       det.formula = det.formula,
  #                       data = tmp.data,
  #                       n.samples = n.samples,
  #                       n.omp.threads = 1,
  #                       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[1]][1] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- NA
  out <- intPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
                  data = tmp.data,
                  n.samples = n.samples,
                  n.omp.threads = 1,
                  verbose = FALSE)
  expect_s3_class(out, "intPGOcc")
})

# Check default priors ---------------- 
test_that("default priors, inits, burn, thin, work", {

  out <- intPGOcc(occ.formula = occ.formula, 
		  det.formula = det.formula,
		  data = data.list,
		  n.samples = 100, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.chains = 1)
  expect_s3_class(out, "intPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(out <- intPGOcc(occ.formula = occ.formula, 
                                det.formula = det.formula,  
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

# Check waicOcc -----------------------
test_that("waicOCC works for intPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for intPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for intPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for intPGOcc", {
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
J.x <- 15
J.y <- 15
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
beta <- c(0.5, 1.2, -0.8)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- c(1, -0.5)
alpha[[2]] <- c(0, -0.2, 0.5)
alpha[[3]] <- c(0.25)
alpha[[4]] <- c(-1, 1, 0.3)
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
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2], 
                      det.cov.2.2 = X.p[[2]][, , 3])
det.covs[[3]] <- list(int = X.p[[3]][, , 1])
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2], 
                      det.cov.4.2 = X.p[[4]][, , 3])
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
n.samples <- 1000
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1, 
		    f.2 = ~ det.cov.2.1 + det.cov.2.2, 
		    f.3 = ~ 1, 
		    f.4 = ~ det.cov.4.1 + det.cov.4.2)
out <- intPGOcc(occ.formula = occ.formula, 
                det.formula = det.formula, 
		data = data.list,
		inits = inits.list,
		n.samples = n.samples, 
		priors = prior.list, 
		n.omp.threads = 1, 
		verbose = FALSE, 
		n.report = 100, 
		n.burn = 500, 
		n.thin = 2, 
		n.chains = 2, 
                k.fold = 2, 
                k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class intPGOcc", {
  expect_s3_class(out, "intPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- intPGOcc(occ.formula = occ.formula, 
                  det.formula = det.formula, 
		  data = data.list,
		  inits = inits.list,
		  n.samples = n.samples, 
		  priors = prior.list, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.report = 100, 
		  n.burn = 500, 
		  n.thin = 2, 
		  n.chains = 2, 
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
                        det.formula = det.formula,
                        data = tmp.data,
                        n.samples = n.samples,
                        n.omp.threads = 1,
                        verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[1]][1] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- NA
  out <- intPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
                  data = tmp.data,
                  n.samples = n.samples,
                  n.omp.threads = 1,
                  verbose = FALSE)
  expect_s3_class(out, "intPGOcc")
})

# Check default priors ---------------- 
test_that("default priors, inits, burn, thin, work", {

  out <- intPGOcc(occ.formula = occ.formula, 
		  det.formula = det.formula,
		  data = data.list,
		  n.samples = 100, 
		  n.omp.threads = 1, 
		  verbose = FALSE, 
		  n.chains = 1)
  expect_s3_class(out, "intPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(out <- intPGOcc(occ.formula = occ.formula, 
                                det.formula = det.formula,  
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

# Check waicOcc -----------------------
test_that("waicOCC works for intPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for intPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for intPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for intPGOcc", {
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
J.x <- 15
J.y <- 15
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
beta <- c(0.5, 1.2, -0.8)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- c(1, -0.5)
alpha[[2]] <- c(0, -0.2, 0.5)
alpha[[3]] <- c(0.25)
alpha[[4]] <- c(-1, 1, 0.3)
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
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2],
                      det.cov.2.2 = X.p[[2]][, , 3])
det.covs[[3]] <- list(int = X.p[[3]][, , 1])
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2],
                      det.cov.4.2 = X.p[[4]][, , 3])
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
n.samples <- 1000
occ.formula <- ~ occ.cov.1 * occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1,
		    f.2 = ~ det.cov.2.1 * det.cov.2.2,
		    f.3 = ~ 1,
		    f.4 = ~ det.cov.4.1 + det.cov.4.2)
out <- intPGOcc(occ.formula = occ.formula,
                det.formula = det.formula,
		data = data.list,
		inits = inits.list,
		n.samples = n.samples,
		priors = prior.list,
		n.omp.threads = 1,
		verbose = FALSE,
		n.report = 100,
		n.burn = 500,
		n.thin = 2,
		n.chains = 2,
                k.fold = 2,
                k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class intPGOcc", {
  expect_s3_class(out, "intPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- intPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
		  data = data.list,
		  inits = inits.list,
		  n.samples = n.samples,
		  priors = prior.list,
		  n.omp.threads = 1,
		  verbose = FALSE,
		  n.report = 100,
		  n.burn = 500,
		  n.thin = 2,
		  n.chains = 2,
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
                        det.formula = det.formula,
                        data = tmp.data,
                        n.samples = n.samples,
                        n.omp.threads = 1,
                        verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[1]][1] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- NA
  out <- intPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
                  data = tmp.data,
                  n.samples = n.samples,
                  n.omp.threads = 1,
                  verbose = FALSE)
  expect_s3_class(out, "intPGOcc")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, work", {

  out <- intPGOcc(occ.formula = occ.formula,
		  det.formula = det.formula,
		  data = data.list,
		  n.samples = 100,
		  n.omp.threads = 1,
		  verbose = FALSE,
		  n.chains = 1)
  expect_s3_class(out, "intPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(out <- intPGOcc(occ.formula = occ.formula,
                                det.formula = det.formula,
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

# Check waicOcc -----------------------
test_that("waicOCC works for intPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for intPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for intPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for intPGOcc", {
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
J.x <- 15
J.y <- 15
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
n.rep.max <- c(6, 4, 4, 5)
# Occupancy covariates
beta <- c(0.5, 1.2, -0.8)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- c(1, -0.5)
alpha[[2]] <- c(0, -0.2, 0.5)
alpha[[3]] <- c(0.25)
alpha[[4]] <- c(-1, 1, 0.3)
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)

# Simulate occupancy data.
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs,
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = FALSE, n.rep.max = n.rep.max)

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites

# Package all data into a list
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2],
                      det.cov.2.2 = X.p[[2]][, , 3])
det.covs[[3]] <- list(int = X.p[[3]][, , 1])
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2],
                      det.cov.4.2 = X.p[[4]][, , 3])
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
n.samples <- 1000
occ.formula <- ~ occ.cov.1 + occ.cov.2
det.formula <- list(f.1 = ~ det.cov.1.1,
		    f.2 = ~ det.cov.2.1 + det.cov.2.2,
		    f.3 = ~ 1,
		    f.4 = ~ det.cov.4.1 + det.cov.4.2)
out <- intPGOcc(occ.formula = occ.formula,
                det.formula = det.formula,
		data = data.list,
		inits = inits.list,
		n.samples = n.samples,
		priors = prior.list,
		n.omp.threads = 1,
		verbose = FALSE,
		n.report = 100,
		n.burn = 500,
		n.thin = 2,
		n.chains = 2,
                k.fold = 2,
                k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class intPGOcc", {
  expect_s3_class(out, "intPGOcc")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$k.fold.deviance), n.data)
  expect_type(out$k.fold.deviance, "double")
  expect_equal(sum(out$k.fold.deviance < 0), 0)
})

# Check individual data set cv --------
test_that("individual data set cv works", {
  out <- intPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
		  data = data.list,
		  inits = inits.list,
		  n.samples = n.samples,
		  priors = prior.list,
		  n.omp.threads = 1,
		  verbose = FALSE,
		  n.report = 100,
		  n.burn = 500,
		  n.thin = 2,
		  n.chains = 2,
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

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
                        det.formula = det.formula,
                        data = tmp.data,
                        n.samples = n.samples,
                        n.omp.threads = 1,
                        verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- 1
  tmp.data$det.covs[[1]][[1]][1] <- NA
  expect_error(intPGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = tmp.data,
               n.samples = n.samples,
               n.omp.threads = 1,
               verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[[1]][1, 1] <- NA
  out <- intPGOcc(occ.formula = occ.formula,
                  det.formula = det.formula,
                  data = tmp.data,
                  n.samples = n.samples,
                  n.omp.threads = 1,
                  verbose = FALSE)
  expect_s3_class(out, "intPGOcc")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, work", {

  out <- intPGOcc(occ.formula = occ.formula,
		  det.formula = det.formula,
		  data = data.list,
		  n.samples = 100,
		  n.omp.threads = 1,
		  verbose = FALSE,
		  n.chains = 1)
  expect_s3_class(out, "intPGOcc")
})

test_that("verbose prints to the screen", {

  expect_output(out <- intPGOcc(occ.formula = occ.formula,
                                det.formula = det.formula,
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

# Check waicOcc -----------------------
test_that("waicOCC works for intPGOcc", {
  waic.out <- waicOcc(out)
  expect_equal(ncol(waic.out), 3)
  expect_equal(nrow(waic.out), n.data)
  expect_equal(waic.out[, 3], -2 * (waic.out[, 1] - waic.out[, 2]))
})

test_that("fitted works for intPGOcc", {
  fitted.out <- fitted(out)
  expect_equal(class(fitted.out), "list")
  expect_equal(length(fitted.out), 2)
  expect_equal(length(fitted.out[[1]]), n.data)
})

test_that("predict works for intPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0 <- dat$X.pred
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})

test_that("posterior predictive checks work for intPGOcc", {
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

