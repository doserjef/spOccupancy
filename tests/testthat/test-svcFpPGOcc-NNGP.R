# Test svcFpPGOcc.R  -----------------------------------------------------
# NNGP --------------------------------------------------------------------
skip_on_cran()

# Intercept Only ----------------------------------------------------------
set.seed(100)
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- 1
p.svc <- length(svc.cols)
phi <- 3 / .6
sigma.sq <- 2
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

data.list <- list(y = y, coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0),
                   fix = TRUE)
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)

n.batch <- 40
batch.length <- 25
n.report <- 100
occ.formula <- ~ 1
tp.det.formula <- ~ 1
fp.det.formula <- ~ 1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
  expect_equal(out$n.val, n.val)
  expect_equal(out$d, d)
})

# Check output dimensions -------------
test_that("posterior sample dimensions are correct", {
  n.post.samples <- out$n.post * out$n.chains
  expect_equal(n.post.samples, 200)
  expect_equal(dim(out$beta.samples), c(n.post.samples, p.occ))
  expect_equal(dim(out$gamma.samples), c(n.post.samples, p.tp.det))
  expect_equal(dim(out$delta.samples), c(n.post.samples, p.fp.det))
  expect_equal(dim(out$z.samples), c(n.post.samples, J.fit))
  expect_equal(dim(out$psi.samples), c(n.post.samples, J.fit))
  expect_equal(dim(out$like.samples), c(n.post.samples, J.fit))
  expect_equal(length(out$rhat$beta), p.occ)
  expect_equal(length(out$rhat$gamma), p.tp.det)
  expect_equal(length(out$rhat$delta), p.fp.det)
})

# Check spatial output ----------------
test_that("spatial parameter dimensions are correct", {
  n.post.samples <- out$n.post * out$n.chains
  expect_equal(dim(out$theta.samples), c(n.post.samples, 2 * p.svc))
  expect_equal(colnames(out$theta.samples),
               c('sigma.sq-(Intercept)', 'phi-(Intercept)'))
  expect_equal(dim(out$w.samples), c(n.post.samples, p.svc, J.fit))
  expect_equal(length(out$rhat$theta), 2 * p.svc)
  expect_equal(length(out$ESS$theta), 2 * p.svc)
  expect_equal(out$type, "NNGP")
  expect_equal(out$n.neighbors, 10)
  expect_equal(out$svc.cols, svc.cols)
  expect_true(all(out$theta.samples[, 'phi-(Intercept)'] >= 3 / 1 &
                  out$theta.samples[, 'phi-(Intercept)'] <= 3 / 0.05))
})
test_that("coordinates are returned in the original order", {
  expect_equal(unname(out$coords), unname(coords))
  expect_equal(unname(out$X.w), unname(X[, svc.cols, drop = FALSE]))
})

# Check fixed initial values ----------
test_that("validated true positives force z = 1", {
  z.fixed <- which(apply(d, 1, sum, na.rm = TRUE) > 0)
  expect_true(all(out$z.samples[, z.fixed] == 1))
})

# Check equivalence with spFpPGOcc ----
test_that("a spatially-varying intercept reproduces spFpPGOcc", {
  # With every initial value supplied, both functions use the random number
  # stream only in the sampler, so the chains should be identical.
  sp.inits <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
  set.seed(1001)
  out.svc <- svcFpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        inits = sp.inits,
                        priors = prior.list,
                        tuning = tuning.list,
                        svc.cols = 1,
                        n.neighbors = 10,
                        n.batch = n.batch,
                        batch.length = batch.length,
                        n.omp.threads = 1,
                        verbose = FALSE,
                        n.burn = 400,
                        n.thin = 6)
  set.seed(1001)
  out.sp <- spFpPGOcc(occ.formula = occ.formula,
                      tp.det.formula = tp.det.formula,
                      fp.det.formula = fp.det.formula,
                      data = data.list,
                      inits = sp.inits,
                      priors = list(beta.normal = list(mean = 0, var = 2.72),
                                    gamma.normal = list(mean = 0, var = 10),
                                    delta.normal = list(mean = 0, var = 10),
                                    phi.unif = c(3 / 1, 3 / 0.05),
                                    sigma.sq.ig = c(2, 1)),
                      tuning = tuning.list,
                      n.neighbors = 10,
                      n.batch = n.batch,
                      batch.length = batch.length,
                      n.omp.threads = 1,
                      verbose = FALSE,
                      n.burn = 400,
                      n.thin = 6)
  expect_equal(as.matrix(out.svc$beta.samples), as.matrix(out.sp$beta.samples))
  expect_equal(as.matrix(out.svc$gamma.samples), as.matrix(out.sp$gamma.samples))
  expect_equal(as.matrix(out.svc$delta.samples), as.matrix(out.sp$delta.samples))
  expect_equal(unname(as.matrix(out.svc$theta.samples)),
               unname(as.matrix(out.sp$theta.samples)))
  expect_equal(unname(out.svc$w.samples[, 1, ]), unname(as.matrix(out.sp$w.samples)))
  expect_equal(as.matrix(out.svc$z.samples), as.matrix(out.sp$z.samples))
  expect_equal(as.matrix(out.svc$psi.samples), as.matrix(out.sp$psi.samples))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- svcFpPGOcc(occ.formula = occ.formula,
                                 tp.det.formula = tp.det.formula,
                                 fp.det.formula = fp.det.formula,
                                 data = data.list,
                                 n.thin = 13,
                                 n.batch = n.batch,
                                 batch.length = batch.length,
                                 svc.cols = svc.cols,
                                 n.omp.threads = 1,
                                 verbose = FALSE))
})

# Check correlation functions ---------
test_that("all correlation functions work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    tuning = list(gamma = 0.3, delta = 0.3, phi = 0.5),
                    svc.cols = svc.cols,
                    cov.model = "gaussian",
                    n.neighbors = 5,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")

  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    tuning = list(gamma = 0.3, delta = 0.3, phi = 0.5),
                    svc.cols = svc.cols,
                    cov.model = "spherical",
                    n.neighbors = 5,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")

  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    priors = list(nu.unif = list(0.5, 2)),
                    tuning = list(gamma = 0.3, delta = 0.3, phi = 0.5, nu = 0.5),
                    svc.cols = svc.cols,
                    cov.model = "matern",
                    n.neighbors = 5,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE,
                    n.chains = 2)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(dim(out$theta.samples), c(out$n.post * out$n.chains, 3 * p.svc))
  expect_equal(colnames(out$theta.samples),
               c('sigma.sq-(Intercept)', 'phi-(Intercept)', 'nu-(Intercept)'))
  expect_equal(length(out$rhat$theta), 3 * p.svc)
  expect_true(all(out$theta.samples[, 'nu-(Intercept)'] >= 0.5 &
                  out$theta.samples[, 'nu-(Intercept)'] <= 2))
})

# Check search type -------------------
test_that("brute force neighbor search works", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    n.neighbors = 5,
                    search.type = 'brute',
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(out$n.neighbors, 5)
})

# Fixed sigma sq ----------------------
test_that("svcFpPGOcc works with fixed sigma.sq", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    priors = list(sigma.sq.ig = 'fixed'),
                    inits = list(sigma.sq = 1.5),
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(length(unique(out$theta.samples[, 1])), 1)
  expect_equal(unique(as.vector(out$theta.samples[, 1])), 1.5)
  # The fixed value must be supplied
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = data.list,
                          priors = list(sigma.sq.ig = 'fixed'),
                          tuning = tuning.list,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE),
               "sigma.sq is specified as fixed")
})

# Uniform sigma sq --------------------
test_that("svcFpPGOcc works with uniform prior on sigma.sq", {
  prior.list <- list(sigma.sq.unif = list(0, 5),
                     nu.unif = list(0.1, 4))
  tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 0.5, nu = 0.6, sigma.sq = 0.7)
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    priors = prior.list,
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    cov.model = "exponential",
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_true(all(out$theta.samples[, 'sigma.sq-(Intercept)'] >= 0 &
                  out$theta.samples[, 'sigma.sq-(Intercept)'] <= 5))
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    priors = prior.list,
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    cov.model = "matern",
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  # sigma.sq tuning is required with a uniform prior
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = data.list,
                          priors = list(sigma.sq.unif = list(0, 5)),
                          tuning = list(gamma = 0.3, delta = 0.3, phi = 0.5),
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE),
               "sigma.sq must be specified in tuning")
})

# Fixed phi ---------------------------
test_that("svcFpPGOcc works with fixed phi and nu", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    priors = list(phi.unif = 'fixed'),
                    inits = list(phi = 4),
                    tuning = list(gamma = 0.3, delta = 0.3),
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(unique(as.vector(out$theta.samples[, 'phi-(Intercept)'])), 4)
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    priors = list(phi.unif = 'fixed', nu.unif = 'fixed'),
                    inits = list(phi = 4, nu = 1.5),
                    tuning = list(gamma = 0.3, delta = 0.3),
                    svc.cols = svc.cols,
                    cov.model = 'matern',
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(unique(as.vector(out$theta.samples[, 'phi-(Intercept)'])), 4)
  expect_equal(unique(as.vector(out$theta.samples[, 'nu-(Intercept)'])), 1.5)
  # The fixed value must be supplied
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = data.list,
                          priors = list(phi.unif = 'fixed'),
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE),
               "phi is specified as fixed")
})

# Check spatial errors ----------------
test_that("spatial error handling works", {
  fit.sp <- function(data = data.list, ...) {
    svcFpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data,
               svc.cols = svc.cols,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE, ...)
  }
  # coords must be supplied, a matrix, and unique
  tmp.data <- data.list
  tmp.data$coords <- NULL
  expect_error(fit.sp(tmp.data), "coords must be specified")
  tmp.data <- data.list
  tmp.data$coords <- c(coords)
  expect_error(fit.sp(tmp.data), "coords must be a matrix or data frame")
  tmp.data <- data.list
  tmp.data$coords[2, ] <- tmp.data$coords[1, ]
  expect_error(fit.sp(tmp.data), "not all unique")
  tmp.data <- data.list
  tmp.data$coords <- tmp.data$coords[-1, ]
  expect_error(fit.sp(tmp.data), "grid.index must be specified")
  # Covariance model and search type
  expect_error(fit.sp(cov.model = 'linear'), "is not a valid option")
  expect_error(fit.sp(search.type = 'kd'), "is not a valid option")
  # Matern requires a prior on nu
  expect_error(fit.sp(cov.model = 'matern'), "nu.unif must be specified")
  # Spatial priors are lists with elements of length 1 or p.svc
  expect_error(fit.sp(priors = list(phi.unif = c(3, 60))),
               "phi.unif must be a list of length 2")
  expect_error(fit.sp(priors = list(phi.unif = list(c(1, 2), 60))),
               "phi.unif\\[\\[1\\]\\] must be a vector of length 1")
  expect_error(fit.sp(priors = list(sigma.sq.ig = list(2, 1, 1))),
               "sigma.sq.ig must be a list of length 2")
  expect_error(fit.sp(priors = list(sigma.sq.ig = list(2, 1), sigma.sq.unif = list(0, 5))),
               "cannot specify both")
  # Spatial initial values
  expect_error(fit.sp(inits = list(phi = c(1, 2))), "initial values for phi must be of length")
  expect_error(fit.sp(inits = list(sigma.sq = c(1, 2))),
               "initial values for sigma.sq must be of length")
  expect_error(fit.sp(inits = list(w = rep(0, J.fit))), "initial values for w")
  expect_error(fit.sp(inits = list(w = matrix(0, p.svc, J.fit - 1))), "initial values for w")
  # Spatial tuning
  expect_error(fit.sp(tuning = list(gamma = 0.3, delta = 0.3)),
               "phi must be specified in tuning")
  expect_error(fit.sp(tuning = list(gamma = 0.3, delta = 0.3, phi = c(1, 1))),
               "phi tuning must be")
  expect_error(fit.sp(cov.model = 'matern', priors = list(nu.unif = list(0.5, 2)),
                      tuning = list(gamma = 0.3, delta = 0.3, phi = 1)),
               "nu must be specified in tuning")
  # Unused arguments give a warning
  expect_warning(fit.sp(k.fold = 2), "is not an argument")
  # Initial values for w are a p.svc x J matrix
  out <- fit.sp(inits = list(w = matrix(0, p.svc, J.fit)))
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
  expect_output(summary(out), "Spatial Covariance")
  expect_output(print(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$mu.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$omega.samples), c(n.post.samples, dim(y)))
  # Missing cells in y stay missing
  expect_true(all(is.na(fitted.out$y.rep.samples[1, , ][is.na(y)])))
  # Replicate validated true positives cannot exceed validated signals
  expect_true(all(sweep(fitted.out$d.rep.samples, c(2, 3), n.val, '<=') |
                  is.na(fitted.out$d.rep.samples)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_s3_class(pred.out, "predict.svcFpPGOcc")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$w.0.samples), c(n.post.samples, p.svc, nrow(X.0)))
  expect_true(all(pred.out$psi.0.samples >= 0 & pred.out$psi.0.samples <= 1))
})
test_that("predict at observed locations returns the fitted w", {
  pred.out <- predict(out, X, coords, verbose = FALSE)
  expect_equal(pred.out$w.0.samples, out$w.samples, check.attributes = FALSE)
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
  expect_true(all(pred.out$mu.0.samples > 0))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
  expect_true(all(pred.out$omega.0.samples > 0))
})
test_that("predict errors on bad input", {
  expect_error(predict(out, X.0, coords.0, type = 'detection'))
  expect_error(predict(out, X.0, verbose = FALSE), "coords.0 must be specified")
  expect_error(predict(out, X.0, coords.0[, 1, drop = FALSE], verbose = FALSE))
  expect_error(predict(out, matrix(1, nrow(X.0), 3), coords.0, verbose = FALSE))
})

# Check SVC samples -------------------
test_that("getSVCSamples works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  svc.samples <- getSVCSamples(out)
  expect_type(svc.samples, 'list')
  expect_equal(length(svc.samples), p.svc)
  expect_equal(names(svc.samples), colnames(out$beta.samples)[svc.cols])
  expect_equal(dim(svc.samples[[1]]), c(n.post.samples, J.fit))
  expect_equal(as.matrix(svc.samples[[1]]),
               as.matrix(out$beta.samples)[, 1] + out$w.samples[, 1, ],
               check.attributes = FALSE)
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  svc.samples <- getSVCSamples(out, pred.out)
  expect_equal(dim(svc.samples[[1]]), c(n.post.samples, nrow(X.0)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(length(ppc.out$fit.d.rep), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.rep.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(length(ppc.out$fit.d.rep), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.rep.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(length(ppc.out$fit.d.rep), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.rep.group.quants), c(5, max(n.rep)))
  }
  expect_output(summary(ppc.out))
  expect_error(ppcOcc(out, 'chi-square', 3))
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
  expect_equal(dim(out.resids$occ.resids), c(10, J.fit))
  expect_equal(dim(out.resids$det.resids), c(10, dim(y)))
  expect_equal(dim(out.resids$val.resids), c(10, dim(y)))
  for (type in c('pearson', 'response')) {
    out.resids <- residuals(out, n.post.samples = 10, type = type)
    expect_equal(dim(out.resids$det.resids), c(10, dim(y)))
  }
  expect_error(residuals(out, n.post.samples = 10, type = 'deviance'))
  expect_error(residuals(out, n.post.samples = out$n.post * out$n.chains + 1))
})

# Check plot --------------------------
test_that("plot works for svcFpPGOcc", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'beta', density = FALSE))
  expect_no_error(plot(out, 'gamma', density = FALSE))
  expect_no_error(plot(out, 'delta', density = FALSE))
  expect_no_error(plot(out, 'theta', density = FALSE))
  expect_error(plot(out, 'beta.comm'))
})

# Occurrence covariate only -----------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

colnames(X) <- c('int', 'occ.cov.1')

data.list <- list(y = y, occ.covs = X, coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ 1
fp.det.formula <- ~ 1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check covariates are in original order
test_that("design matrices are returned in the original order", {
  expect_equal(unname(out$X), unname(X))
  expect_equal(unname(out$X.w), unname(X[, svc.cols, drop = FALSE]))
  expect_equal(unname(out$coords), unname(coords))
})

# Check spatial output ----------------
test_that("spatial parameter dimensions are correct", {
  n.post.samples <- out$n.post * out$n.chains
  expect_equal(dim(out$theta.samples), c(n.post.samples, 2 * p.svc))
  expect_equal(colnames(out$theta.samples),
               c('sigma.sq-(Intercept)', 'sigma.sq-occ.cov.1',
                 'phi-(Intercept)', 'phi-occ.cov.1'))
  expect_equal(dim(out$w.samples), c(n.post.samples, p.svc, J.fit))
  expect_equal(length(out$rhat$theta), 2 * p.svc)
  expect_equal(length(out$ESS$theta), 2 * p.svc)
  expect_equal(out$svc.cols, svc.cols)
})

# Check svc.cols ----------------------
test_that("svc.cols can be given as column names", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    tuning = tuning.list,
                    svc.cols = c('(Intercept)', 'occ.cov.1'),
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(out$svc.cols, c(1, 2))
  # A single SVC that is not the intercept
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    tuning = tuning.list,
                    svc.cols = 'occ.cov.1',
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(out$svc.cols, 2)
  expect_equal(colnames(out$theta.samples), c('sigma.sq-occ.cov.1', 'phi-occ.cov.1'))
  expect_equal(dim(out$w.samples), c(out$n.post * out$n.chains, 1, J.fit))
})
test_that("svc.cols error handling works", {
  fit.svc <- function(svc.cols) {
    svcFpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               svc.cols = svc.cols,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE)
  }
  expect_error(fit.svc(c(1, 3)), "not in design matrix columns")
  expect_error(fit.svc(c('(Intercept)', 'occ.cov.2')), "not in occurrence covariates")
  expect_error(fit.svc(list(1)), "must be a numeric or character vector")
  expect_error(fit.svc(c(1, 1)), "must not contain duplicate values")
})

# Check SVC-specific spatial inputs ---
test_that("spatial priors, inits, and tuning can differ across SVCs", {
  fit.sp <- function(...) {
    svcFpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               svc.cols = svc.cols,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE, ...)
  }
  out <- fit.sp(priors = list(phi.unif = list(a = c(3 / 1, 3 / 0.8), b = c(3 / 0.05, 3 / 0.1)),
                              sigma.sq.ig = list(a = c(2, 3), b = c(1, 0.5))),
                inits = list(phi = c(3, 6), sigma.sq = c(1, 0.5)),
                tuning = list(gamma = 0.3, delta = 0.3, phi = c(1, 0.5)))
  expect_s3_class(out, "svcFpPGOcc")
  expect_true(all(out$theta.samples[, 'phi-(Intercept)'] >= 3 / 1 &
                  out$theta.samples[, 'phi-(Intercept)'] <= 3 / 0.05))
  expect_true(all(out$theta.samples[, 'phi-occ.cov.1'] >= 3 / 0.8 &
                  out$theta.samples[, 'phi-occ.cov.1'] <= 3 / 0.1))
  # Fixed spatial parameters, with a different value for each SVC
  out <- fit.sp(priors = list(sigma.sq.ig = 'fixed', phi.unif = 'fixed'),
                inits = list(sigma.sq = c(1.5, 0.5), phi = c(4, 8)),
                tuning = list(gamma = 0.3, delta = 0.3))
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(unname(apply(out$theta.samples, 2, function(a) length(unique(a)))),
               rep(1, 2 * p.svc))
  expect_equal(unname(out$theta.samples[1, ]), c(1.5, 0.5, 4, 8))
  # Lengths must be 1 or p.svc
  expect_error(fit.sp(priors = list(phi.unif = list(a = c(1, 2, 3), b = 60))),
               "phi.unif\\[\\[1\\]\\] must be a vector of length 2")
  expect_error(fit.sp(priors = list(sigma.sq.ig = list(a = 2, b = c(1, 1, 1)))),
               "sigma.sq.ig\\[\\[2\\]\\] must be a vector of length 2")
  expect_error(fit.sp(inits = list(phi = c(1, 2, 3))),
               "initial values for phi must be of length 2 or 1")
  expect_error(fit.sp(inits = list(sigma.sq = c(1, 2, 3))),
               "initial values for sigma.sq must be of length 2 or 1")
  expect_error(fit.sp(inits = list(w = matrix(0, 1, J.fit))),
               "initial values for w must be a matrix with dimensions 2")
  expect_error(fit.sp(tuning = list(gamma = 0.3, delta = 0.3, phi = c(1, 1, 1))),
               "phi tuning must be either a single value or a vector of length 2")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = tmp.data,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, ] <- NA
  tmp.data$n.val[1, ] <- NA
  tmp.data$d[1, ] <- NA
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = tmp.data,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE),
               "all missing detection histories")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$w.0.samples), c(n.post.samples, p.svc, nrow(X.0)))
})
test_that("predict at observed locations returns the fitted w", {
  pred.out <- predict(out, X, coords, verbose = FALSE)
  expect_equal(pred.out$w.0.samples, out$w.samples, check.attributes = FALSE)
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check SVC samples -------------------
test_that("getSVCSamples works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  svc.samples <- getSVCSamples(out)
  expect_type(svc.samples, 'list')
  expect_equal(length(svc.samples), p.svc)
  expect_equal(names(svc.samples), c('(Intercept)', 'occ.cov.1'))
  for (i in 1:p.svc) {
    expect_equal(dim(svc.samples[[i]]), c(n.post.samples, J.fit))
    expect_equal(as.matrix(svc.samples[[i]]),
                 as.matrix(out$beta.samples)[, svc.cols[i]] + out$w.samples[, i, ],
                 check.attributes = FALSE)
  }
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  svc.samples <- getSVCSamples(out, pred.out)
  expect_equal(length(svc.samples), p.svc)
  expect_equal(dim(svc.samples[[2]]), c(n.post.samples, nrow(X.0)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for svcFpPGOcc with multiple SVCs", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'beta', density = FALSE))
  expect_no_error(plot(out, 'theta', density = FALSE))
})

# True positive covariate only --------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- 1
p.svc <- length(svc.cols)
phi <- 3 / .6
sigma.sq <- 2
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

data.list <- list(y = y, tp.det.covs = list(tp.cov.1 = X.mu[, , 2]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ 1
tp.det.formula <- ~ tp.cov.1
fp.det.formula <- ~ 1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check observation-level covariates --
test_that("detection design matrices are returned in the original order", {
  n.obs <- sum(!is.na(y))
  expect_equal(nrow(out$X.mu), n.obs)
  expect_equal(unname(out$X.mu[, 2]), unname(c(X.mu[, , 2])[!is.na(c(y))]))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  # Missing covariate where y is observed
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$n.val[1, 1] <- 0
  tmp.data$d[1, 1] <- 0
  tmp.data$tp.det.covs$tp.cov.1[1, 1] <- NA
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = tmp.data,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE))
  # Missing y without matching missing n.val and d
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = tmp.data,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE))
  # Missing y, n.val, and d in the same cell
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = tmp.data,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(nrow(out$X.mu), sum(!is.na(tmp.data$y)))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X.mu[pred.indx, 1, ]
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.mu.0)))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# False positive covariate only -------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- 1
p.svc <- length(svc.cols)
phi <- 3 / .6
sigma.sq <- 2
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

data.list <- list(y = y, fp.det.covs = list(fp.cov.1 = X.omega[, , 2]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ 1
tp.det.formula <- ~ 1
fp.det.formula <- ~ fp.cov.1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check observation-level covariates --
test_that("detection design matrices are returned in the original order", {
  n.obs <- sum(!is.na(y))
  expect_equal(nrow(out$X.omega), n.obs)
  expect_equal(unname(out$X.omega[, 2]), unname(c(X.omega[, , 2])[!is.na(c(y))]))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$n.val[1, 1] <- 0
  tmp.data$d[1, 1] <- 0
  tmp.data$fp.det.covs$fp.cov.1[1, 1] <- NA
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = tmp.data,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = tmp.data,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(nrow(out$X.omega), sum(!is.na(tmp.data$y)))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X.omega[pred.indx, 1, ]
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.omega.0)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Covariates on all -------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6, -0.3)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

colnames(X) <- c('int', 'occ.cov.1')
colnames(X.0) <- c('int', 'occ.cov.1')
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.cov.2 = X.omega[, , 3]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ tp.cov.1
fp.det.formula <- ~ fp.cov.1 + fp.cov.2

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check parameter names ---------------
test_that("parameter names are correct", {
  expect_equal(colnames(out$beta.samples), c('(Intercept)', 'occ.cov.1'))
  expect_equal(colnames(out$gamma.samples), c('(Intercept)', 'tp.cov.1'))
  expect_equal(colnames(out$delta.samples), c('(Intercept)', 'fp.cov.1', 'fp.cov.2'))
  expect_equal(colnames(out$theta.samples),
               c('sigma.sq-(Intercept)', 'sigma.sq-occ.cov.1',
                 'phi-(Intercept)', 'phi-occ.cov.1'))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = tmp.data,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$n.val[1, 1] <- 0
  tmp.data$d[1, 1] <- 0
  tmp.data$fp.det.covs$fp.cov.2[1, 1] <- NA
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = tmp.data,
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          n.omp.threads = 1,
                          verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = tmp.data,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check validation data ---------------
test_that("validation data error handling works", {
  fit.fp <- function(data, ...) {
    svcFpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data,
               svc.cols = svc.cols,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE, ...)
  }
  # n.val and d must be supplied
  tmp.data <- data.list
  tmp.data$n.val <- NULL
  expect_error(fit.fp(tmp.data), "n.val must be specified")
  tmp.data <- data.list
  tmp.data$d <- NULL
  expect_error(fit.fp(tmp.data), "d must be specified")
  # Dimensions must match y
  tmp.data <- data.list
  tmp.data$n.val <- tmp.data$n.val[, -1]
  expect_error(fit.fp(tmp.data), "dimension mismatch")
  tmp.data <- data.list
  tmp.data$d <- tmp.data$d[-1, ]
  expect_error(fit.fp(tmp.data), "dimension mismatch")
  # Missingness pattern must match y
  tmp.data <- data.list
  obs.indx <- which(!is.na(tmp.data$y))[1]
  tmp.data$d[obs.indx] <- NA
  expect_error(fit.fp(tmp.data), "mismatch between missing values")
  # Counts must be non-negative integers
  tmp.data <- data.list
  tmp.data$y[obs.indx] <- 1.5
  expect_error(fit.fp(tmp.data), "non-negative integer")
  tmp.data <- data.list
  tmp.data$n.val[obs.indx] <- -1
  expect_error(fit.fp(tmp.data), "non-negative integer")
  # d <= n.val <= y
  tmp.data <- data.list
  tmp.data$n.val[obs.indx] <- tmp.data$y[obs.indx] + 1
  expect_error(fit.fp(tmp.data), "n.val cannot be larger than y")
  tmp.data <- data.list
  tmp.data$d[obs.indx] <- tmp.data$n.val[obs.indx] + 1
  expect_error(fit.fp(tmp.data), "d cannot be larger than n.val")
  # z inits must be 1 where a validated true positive exists
  z.inits <- rep(0, J.fit)
  expect_error(fit.fp(data.list, inits = list(z = z.inits)),
               "initial values for latent occurrence")
  # Missing formulas and data
  expect_error(svcFpPGOcc(occ.formula = occ.formula, fp.det.formula = fp.det.formula,
                          data = data.list, svc.cols = svc.cols, n.batch = n.batch,
                          batch.length = batch.length, verbose = FALSE),
               "tp.det.formula must be specified")
  expect_error(svcFpPGOcc(occ.formula = occ.formula, tp.det.formula = tp.det.formula,
                          data = data.list, svc.cols = svc.cols, n.batch = n.batch,
                          batch.length = batch.length, verbose = FALSE),
               "fp.det.formula must be specified")
  tmp.data <- data.list
  tmp.data$fp.det.covs <- NULL
  expect_error(fit.fp(tmp.data), "fp.det.covs must be specified")
})

# Check priors, inits, and tuning -----
test_that("prior, initial value, and tuning error handling works", {
  fit.fp <- function(...) {
    svcFpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               svc.cols = svc.cols,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE, ...)
  }
  expect_error(fit.fp(priors = list(gamma.normal = list(0, rep(10, 3)))))
  expect_error(fit.fp(priors = list(delta.normal = list(rep(0, 2), 10))))
  expect_error(fit.fp(priors = list(delta.normal = 10)))
  expect_error(fit.fp(inits = list(gamma = c(0, 0, 0))))
  expect_error(fit.fp(inits = list(delta = c(0, 0))))
  expect_error(fit.fp(inits = list(z = rep(2, J.fit))))
  expect_error(fit.fp(tuning = list(gamma = 0.3, phi = 1)), "delta must be specified")
  expect_error(fit.fp(tuning = list(delta = 0.3, phi = 1)), "gamma must be specified")
  expect_error(fit.fp(tuning = list(gamma = c(0.3, 0.3, 0.3), delta = 0.3, phi = 1)),
               "gamma tuning must be")
  # Vector-valued priors, inits, and tuning work
  out <- fit.fp(priors = list(gamma.normal = list(rep(0, p.tp.det), rep(10, p.tp.det)),
                              delta.normal = list(rep(0, p.fp.det), rep(10, p.fp.det))),
                inits = list(gamma = rep(0, p.tp.det), delta = rep(0, p.fp.det)),
                tuning = list(gamma = rep(0.3, p.tp.det), delta = rep(0.3, p.fp.det),
                              phi = c(1, 1)))
  expect_s3_class(out, "svcFpPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$w.0.samples), c(n.post.samples, p.svc, nrow(X.0)))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X.mu[pred.indx, 1, ]
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.mu.0)))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X.omega[pred.indx, 1, ]
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.omega.0)))
  expect_error(predict(out, dat$X.omega[pred.indx, 1, 1:2], type = 'false positive'))
})

# Check SVC samples -------------------
test_that("getSVCSamples works for svcFpPGOcc", {
  svc.samples <- getSVCSamples(out)
  expect_equal(length(svc.samples), p.svc)
  expect_equal(dim(svc.samples[[2]]), c(out$n.post * out$n.chains, J.fit))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Interactions in all -----------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5, -0.4)
p.occ <- length(beta)
gamma <- c(1.5, -0.4, 0.2)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6, -0.3)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

colnames(X) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2],
                                     tp.cov.2 = X.mu[, , 3]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.cov.2 = X.omega[, , 3]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ occ.cov.1 * occ.cov.2
tp.det.formula <- ~ tp.cov.1 * tp.cov.2
fp.det.formula <- ~ fp.cov.1 * fp.cov.2
# SVC on the intercept and the interaction term
svc.cols.fit <- c('(Intercept)', 'occ.cov.1:occ.cov.2')

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols.fit,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check parameter names ---------------
test_that("parameter names are correct", {
  expect_equal(colnames(out$beta.samples),
               c('(Intercept)', 'occ.cov.1', 'occ.cov.2', 'occ.cov.1:occ.cov.2'))
  expect_equal(colnames(out$gamma.samples),
               c('(Intercept)', 'tp.cov.1', 'tp.cov.2', 'tp.cov.1:tp.cov.2'))
  expect_equal(colnames(out$delta.samples),
               c('(Intercept)', 'fp.cov.1', 'fp.cov.2', 'fp.cov.1:fp.cov.2'))
  expect_equal(out$svc.cols, c(1, 4))
  expect_equal(colnames(out$theta.samples),
               c('sigma.sq-(Intercept)', 'sigma.sq-occ.cov.1:occ.cov.2',
                 'phi-(Intercept)', 'phi-occ.cov.1:occ.cov.2'))
  expect_equal(unname(out$X.w[, 2]), unname(X[, 2] * X[, 3]))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols.fit,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols.fit,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  X.0.int <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0.int, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$w.0.samples), c(n.post.samples, p.svc, nrow(X.0)))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X.mu[pred.indx, 1, ]
  X.mu.0 <- cbind(X.mu.0, X.mu.0[, 2] * X.mu.0[, 3])
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.mu.0)))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X.omega[pred.indx, 1, ]
  X.omega.0 <- cbind(X.omega.0, X.omega.0[, 2] * X.omega.0[, 3])
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.omega.0)))
})

# Check SVC samples -------------------
test_that("getSVCSamples works for svcFpPGOcc", {
  svc.samples <- getSVCSamples(out)
  expect_equal(names(svc.samples), c('(Intercept)', 'occ.cov.1:occ.cov.2'))
  expect_equal(as.matrix(svc.samples[[2]]),
               as.matrix(out$beta.samples)[, 4] + out$w.samples[, 2, ],
               check.attributes = FALSE)
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Site covariate on detection ---------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- 2
delta <- c(0.5)
p.fp.det <- 2
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

colnames(X) <- c('int', 'occ.cov.1')
# Site-level detection covariates supplied as length J vectors
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(site.cov = X[, 2]),
                  fp.det.covs = list(site.cov = X[, 2]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ site.cov
fp.det.formula <- ~ site.cov

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check site covariate expansion ------
test_that("site-level detection covariates are expanded to all surveys", {
  n.obs <- sum(!is.na(y))
  expect_equal(nrow(out$X.mu), n.obs)
  expect_equal(nrow(out$X.omega), n.obs)
  site.long <- rep(1:J.fit, ncol(y))[!is.na(c(y))]
  expect_equal(unname(out$X.mu[, 2]), unname(X[site.long, 2]))
  expect_equal(unname(out$X.omega[, 2]), unname(X[site.long, 2]))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})
test_that("true positive prediction works", {
  pred.out <- predict(out, X.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})
test_that("false positive prediction works", {
  pred.out <- predict(out, X.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.0)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Random intercept on occurrence ------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list(levels = c(20),
               sigma.sq.psi = c(1.3))
mu.RE <- list()
omega.RE <- list()
svc.cols <- 1
p.svc <- length(svc.cols)
phi <- 3 / .6
sigma.sq <- 2
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE] + 15
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE] + 15
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1')
occ.covs.0 <- cbind(X.0, X.re.0)
colnames(occ.covs.0) <- c('int', 'occ.factor.1')
data.list <- list(y = y, occ.covs = occ.covs, coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.psi.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.psi = 1,
                   sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ (1 | occ.factor.1)
tp.det.formula <- ~ 1
fp.det.formula <- ~ 1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
  expect_equal(dim(out$sigma.sq.psi.samples), c(out$n.post * out$n.chains, 1))
  expect_equal(ncol(out$beta.star.samples), length(unique(c(X.re))))
  expect_equal(length(out$rhat$sigma.sq.psi), 1)
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
  expect_error(out <- svcFpPGOcc(occ.formula = occ.formula,
                                 tp.det.formula = tp.det.formula,
                                 fp.det.formula = fp.det.formula,
                                 data = data.list,
                                 svc.cols = svc.cols,
                                 n.batch = n.batch,
                                 batch.length = batch.length,
                                 n.omp.threads = 1,
                                 verbose = FALSE))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- svcFpPGOcc(occ.formula = occ.formula,
                                 tp.det.formula = tp.det.formula,
                                 fp.det.formula = fp.det.formula,
                                 data = data.list,
                                 svc.cols = svc.cols,
                                 n.batch = n.batch,
                                 batch.length = batch.length,
                                 n.omp.threads = 1,
                                 verbose = FALSE))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, occ.covs.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
  # Ignoring the random effects
  pred.out <- predict(out, occ.covs.0[, 1, drop = FALSE], coords.0,
                      ignore.RE = TRUE, verbose = FALSE)
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  # New random effect levels
  X.0.new <- occ.covs.0
  X.0.new[, 'occ.factor.1'] <- max(X.re) + 1
  pred.out <- predict(out, X.0.new, coords.0, verbose = FALSE)
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  # Column names must match
  X.0.bad <- occ.covs.0
  colnames(X.0.bad) <- c('int', 'wrong.name')
  expect_error(predict(out, X.0.bad, coords.0, verbose = FALSE))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for random effect variances", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'sigma.sq.psi', density = FALSE))
})

# Multiple random intercepts on occurrence with covariate -----------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list(levels = c(20, 10),
               sigma.sq.psi = c(1.3, 0.5))
mu.RE <- list()
omega.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
occ.covs.0 <- cbind(X.0, X.re.0)
colnames(occ.covs.0) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, occ.covs = occ.covs, coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.psi.ig = list(a = c(0.1, 0.1), b = c(0.1, 0.1)),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.psi = c(1, 1),
                   sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
tp.det.formula <- ~ 1
fp.det.formula <- ~ 1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
  expect_equal(dim(out$sigma.sq.psi.samples), c(out$n.post * out$n.chains, 2))
  expect_equal(colnames(out$sigma.sq.psi.samples), c('occ.factor.1', 'occ.factor.2'))
  expect_equal(ncol(out$beta.star.samples),
               length(unique(X.re[, 1])) + length(unique(X.re[, 2])))
  expect_equal(length(out$rhat$sigma.sq.psi), 2)
})

# Check RE prior and inits errors -----
test_that("random effect prior and initial value errors work", {
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = data.list,
                          priors = list(sigma.sq.psi.ig = list(a = rep(0.1, 3), b = 0.1)),
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          verbose = FALSE))
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = data.list,
                          inits = list(sigma.sq.psi = c(1, 1, 1)),
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          verbose = FALSE))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, occ.covs.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$w.0.samples), c(n.post.samples, p.svc, nrow(X.0)))
  # Ignoring the random effects
  pred.out <- predict(out, X.0, coords.0, ignore.RE = TRUE, verbose = FALSE)
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Random intercept on true positive detection -----------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list(levels = c(12),
              sigma.sq.mu = c(0.75))
omega.RE <- list()
svc.cols <- 1
p.svc <- length(svc.cols)
phi <- 3 / .6
sigma.sq <- 2
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.mu.re <- dat$X.mu.re[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

data.list <- list(y = y,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2],
                                     tp.factor.1 = X.mu.re[, , 1]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.mu = 0.5,
                   sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, gamma.star = 0.3, phi = 1)
occ.formula <- ~ 1
tp.det.formula <- ~ tp.cov.1 + (1 | tp.factor.1)
fp.det.formula <- ~ 1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, TRUE)
  expect_equal(out$omegaRE, FALSE)
  expect_equal(dim(out$sigma.sq.mu.samples), c(out$n.post * out$n.chains, 1))
  expect_equal(ncol(out$gamma.star.samples), length(unique(c(na.omit(c(X.mu.re))))))
  expect_equal(length(out$rhat$sigma.sq.mu), 1)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$mu.re.level.names))),
               sort(unique(c(X.mu.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$tp.det.covs$tp.factor.1 <- factor(data.list$tp.det.covs$tp.factor.1)
  expect_error(out <- svcFpPGOcc(occ.formula = occ.formula,
                                 tp.det.formula = tp.det.formula,
                                 fp.det.formula = fp.det.formula,
                                 data = data.list,
                                 svc.cols = svc.cols,
                                 n.batch = n.batch,
                                 batch.length = batch.length,
                                 n.omp.threads = 1,
                                 verbose = FALSE))
  data.list$tp.det.covs$tp.factor.1 <- as.character(factor(data.list$tp.det.covs$tp.factor.1))
  expect_error(out <- svcFpPGOcc(occ.formula = occ.formula,
                                 tp.det.formula = tp.det.formula,
                                 fp.det.formula = fp.det.formula,
                                 data = data.list,
                                 svc.cols = svc.cols,
                                 n.batch = n.batch,
                                 batch.length = batch.length,
                                 n.omp.threads = 1,
                                 verbose = FALSE))
})

# Check RE tuning ---------------------
test_that("gamma.star tuning is required with true positive random effects", {
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = data.list,
                          tuning = list(gamma = 0.3, delta = 0.3, phi = 1),
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          verbose = FALSE),
               "gamma.star must be specified")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = tmp.data,
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(nrow(out$X.mu.re), sum(!is.na(tmp.data$y)))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})
test_that("true positive prediction works", {
  X.mu.0 <- cbind(dat$X.mu[pred.indx, 1, ], dat$X.mu.re[pred.indx, 1, 1])
  colnames(X.mu.0) <- c('(Intercept)', 'tp.cov.1', 'tp.factor.1')
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.mu.0)))
  # Ignoring the random effects
  pred.out <- predict(out, X.mu.0[, 1:2], type = 'true positive', ignore.RE = TRUE)
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.mu.0)))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for random effect variances", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'sigma.sq.mu', density = FALSE))
})

# Random intercept on false positive detection ----------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list(levels = c(10),
                 sigma.sq.omega = c(0.5))
svc.cols <- 1
p.svc <- length(svc.cols)
phi <- 3 / .6
sigma.sq <- 2
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]
X.omega.re <- dat$X.omega.re[-pred.indx, , , drop = FALSE]

data.list <- list(y = y,
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.factor.1 = X.omega.re[, , 1]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.omega.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.omega = 0.5,
                   sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, delta.star = 0.3, phi = 1)
occ.formula <- ~ 1
tp.det.formula <- ~ 1
fp.det.formula <- ~ fp.cov.1 + (1 | fp.factor.1)

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, TRUE)
  expect_equal(dim(out$sigma.sq.omega.samples), c(out$n.post * out$n.chains, 1))
  expect_equal(ncol(out$delta.star.samples), length(unique(c(na.omit(c(X.omega.re))))))
  expect_equal(length(out$rhat$sigma.sq.omega), 1)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$omega.re.level.names))),
               sort(unique(c(X.omega.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$fp.det.covs$fp.factor.1 <- factor(data.list$fp.det.covs$fp.factor.1)
  expect_error(out <- svcFpPGOcc(occ.formula = occ.formula,
                                 tp.det.formula = tp.det.formula,
                                 fp.det.formula = fp.det.formula,
                                 data = data.list,
                                 svc.cols = svc.cols,
                                 n.batch = n.batch,
                                 batch.length = batch.length,
                                 n.omp.threads = 1,
                                 verbose = FALSE))
  data.list$fp.det.covs$fp.factor.1 <- as.character(factor(data.list$fp.det.covs$fp.factor.1))
  expect_error(out <- svcFpPGOcc(occ.formula = occ.formula,
                                 tp.det.formula = tp.det.formula,
                                 fp.det.formula = fp.det.formula,
                                 data = data.list,
                                 svc.cols = svc.cols,
                                 n.batch = n.batch,
                                 batch.length = batch.length,
                                 n.omp.threads = 1,
                                 verbose = FALSE))
})

# Check RE tuning ---------------------
test_that("delta.star tuning is required with false positive random effects", {
  expect_error(svcFpPGOcc(occ.formula = occ.formula,
                          tp.det.formula = tp.det.formula,
                          fp.det.formula = fp.det.formula,
                          data = data.list,
                          tuning = list(gamma = 0.3, delta = 0.3, phi = 1),
                          svc.cols = svc.cols,
                          n.batch = n.batch,
                          batch.length = batch.length,
                          verbose = FALSE),
               "delta.star must be specified")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = tmp.data,
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
  expect_equal(nrow(out$X.omega.re), sum(!is.na(tmp.data$y)))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  X.omega.0 <- cbind(dat$X.omega[pred.indx, 1, ], dat$X.omega.re[pred.indx, 1, 1])
  colnames(X.omega.0) <- c('(Intercept)', 'fp.cov.1', 'fp.factor.1')
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.omega.0)))
  # Ignoring the random effects
  pred.out <- predict(out, X.omega.0[, 1:2], type = 'false positive', ignore.RE = TRUE)
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.omega.0)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for random effect variances", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'sigma.sq.omega', density = FALSE))
})

# Random intercepts on all with covariates --------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6, -0.3)
p.fp.det <- length(delta)
psi.RE <- list(levels = c(20, 15),
               sigma.sq.psi = c(1.3, 0.2))
mu.RE <- list(levels = c(12),
              sigma.sq.mu = c(0.75))
omega.RE <- list(levels = c(10),
                 sigma.sq.omega = c(0.5))
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.mu.re <- dat$X.mu.re[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]
X.omega.re <- dat$X.omega.re[-pred.indx, , , drop = FALSE]

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
occ.covs.0 <- cbind(X.0, X.re.0)
colnames(occ.covs.0) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, occ.covs = occ.covs,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2],
                                     tp.factor.1 = X.mu.re[, , 1]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.cov.2 = X.omega[, , 3],
                                     fp.factor.1 = X.omega.re[, , 1]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.psi.ig = list(a = 0.1, b = 0.1),
                   sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
                   sigma.sq.omega.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.psi = 1,
                   sigma.sq.mu = 0.5, sigma.sq.omega = 0.5,
                   sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, gamma.star = 0.3, delta.star = 0.3,
                    phi = 1)
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
tp.det.formula <- ~ tp.cov.1 + (1 | tp.factor.1)
fp.det.formula <- ~ fp.cov.1 + fp.cov.2 + (1 | fp.factor.1)

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$muRE, TRUE)
  expect_equal(out$omegaRE, TRUE)
  n.post.samples <- out$n.post * out$n.chains
  expect_equal(dim(out$sigma.sq.psi.samples), c(n.post.samples, 2))
  expect_equal(dim(out$sigma.sq.mu.samples), c(n.post.samples, 1))
  expect_equal(dim(out$sigma.sq.omega.samples), c(n.post.samples, 1))
  expect_equal(dim(out$theta.samples), c(n.post.samples, 2 * p.svc))
  expect_equal(dim(out$w.samples), c(n.post.samples, p.svc, J.fit))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
               sort(unique(c(X.re))))
  expect_equal(sort(unique(unlist(out$mu.re.level.names))),
               sort(unique(c(X.mu.re))))
  expect_equal(sort(unique(unlist(out$omega.re.level.names))),
               sort(unique(c(X.omega.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check fixed initial values ----------
test_that("fixed initial values across chains work", {
  tmp.inits <- inits.list
  tmp.inits$fix <- TRUE
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    inits = tmp.inits,
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.chains = 2,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = tmp.data,
                    tuning = tuning.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, occ.covs.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$w.0.samples), c(n.post.samples, p.svc, nrow(X.0)))
})
test_that("true positive prediction works", {
  X.mu.0 <- cbind(dat$X.mu[pred.indx, 1, ], dat$X.mu.re[pred.indx, 1, 1])
  colnames(X.mu.0) <- c('(Intercept)', 'tp.cov.1', 'tp.factor.1')
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.mu.0)))
})
test_that("false positive prediction works", {
  X.omega.0 <- cbind(dat$X.omega[pred.indx, 1, ], dat$X.omega.re[pred.indx, 1, 1])
  colnames(X.omega.0) <- c('(Intercept)', 'fp.cov.1', 'fp.cov.2', 'fp.factor.1')
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.omega.0)))
})

# Check SVC samples -------------------
test_that("getSVCSamples works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  svc.samples <- getSVCSamples(out)
  expect_equal(length(svc.samples), p.svc)
  expect_equal(dim(svc.samples[[1]]), c(n.post.samples, J.fit))
  pred.out <- predict(out, occ.covs.0, coords.0, verbose = FALSE)
  svc.samples <- getSVCSamples(out, pred.out)
  expect_equal(dim(svc.samples[[2]]), c(n.post.samples, nrow(X.0)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for all parameters", {
  pdf(NULL)
  on.exit(dev.off())
  for (param in c('beta', 'gamma', 'delta', 'theta', 'sigma.sq.psi',
                  'sigma.sq.mu', 'sigma.sq.omega')) {
    expect_no_error(plot(out, param, density = FALSE))
  }
})

# Dimension of y is not equal to max(n.rep) -------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
n.rep.max <- 7
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, n.rep.max = n.rep.max,
                beta = beta, gamma = gamma, delta = delta, psi.RE = psi.RE,
                mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', prop.val = 0.2)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]
J.fit <- nrow(y)
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]

colnames(X) <- c('int', 'occ.cov.1')
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2]),
                  coords = coords, n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ tp.cov.1
fp.det.formula <- ~ fp.cov.1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
  expect_equal(ncol(out$y), n.rep.max)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- svcFpPGOcc(occ.formula = occ.formula,
                    tp.det.formula = tp.det.formula,
                    fp.det.formula = fp.det.formula,
                    data = data.list,
                    svc.cols = svc.cols,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = 1,
                    verbose = FALSE)
  expect_s3_class(out, "svcFpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(svcFpPGOcc(occ.formula = occ.formula,
                           tp.det.formula = tp.det.formula,
                           fp.det.formula = fp.det.formula,
                           data = data.list,
                           svc.cols = svc.cols,
                           n.batch = 4,
                           batch.length = 25,
                           n.omp.threads = 1,
                           verbose = TRUE,
                           n.report = n.report,
                           n.burn = 1,
                           n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, J.fit, n.rep.max))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, J.fit, n.rep.max))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, nrow(X.0)))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, nrow(X.0)))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X.mu[pred.indx, 1, ]
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, nrow(X.mu.0)))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X.omega[pred.indx, 1, ]
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, nrow(X.omega.0)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  # Replicate groups are the columns of y with any observed surveys. simFPOcc
  # spreads surveys across all n.rep.max columns, so this is not max(n.rep).
  n.rep.obs <- sum(colSums(!is.na(y)) > 0)
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit, n.rep.max))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit, n.rep.max))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J.fit))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J.fit))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.rep.obs))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, n.rep.obs))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
  expect_equal(dim(out.resids$det.resids), c(10, J.fit, n.rep.max))
})

# Spatial random effect on a coarser grid (grid.index) --------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
phi <- c(3 / .6, 3 / .4)
sigma.sq <- c(2, 0.8)
# Pairs of neighboring sites share a spatial random effect
grid <- rep(1:(J / 2), each = 2)
J.w <- length(unique(grid))
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, svc.cols = svc.cols, sigma.sq = sigma.sq, phi = phi,
                cov.model = 'exponential', grid = grid, prop.val = 0.2)
y <- dat$y
n.val <- dat$n.val
d <- dat$d
X <- dat$X
coords <- as.matrix(dat$coords)
X.mu <- dat$X.mu
X.omega <- dat$X.omega

colnames(X) <- c('int', 'occ.cov.1')
data.list <- list(y = y, occ.covs = X, coords = coords, grid.index = grid,
                  n.val = n.val, d = d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.05),
                   sigma.sq.ig = list(a = 2, b = 1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq = 1, phi = 3 / 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, phi = 1)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ 1
fp.det.formula <- ~ 1

out <- svcFpPGOcc(occ.formula = occ.formula,
                  tp.det.formula = tp.det.formula,
                  fp.det.formula = fp.det.formula,
                  data = data.list,
                  inits = inits.list,
                  priors = prior.list,
                  tuning = tuning.list,
                  svc.cols = svc.cols,
                  cov.model = 'exponential',
                  n.neighbors = 10,
                  search.type = 'cb',
                  n.batch = n.batch,
                  batch.length = batch.length,
                  accept.rate = 0.43,
                  n.omp.threads = 1,
                  verbose = FALSE,
                  n.report = n.report,
                  n.burn = 400,
                  n.thin = 6,
                  n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class svcFpPGOcc", {
  expect_s3_class(out, "svcFpPGOcc")
})

# Check output dimensions -------------
test_that("spatial random effects are at the grid level", {
  n.post.samples <- out$n.post * out$n.chains
  expect_equal(dim(out$w.samples), c(n.post.samples, p.svc, J.w))
  expect_equal(dim(out$z.samples), c(n.post.samples, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, J))
  expect_equal(unname(out$coords), unname(coords))
  # X.w is at the site level, not the grid level
  expect_equal(unname(out$X.w), unname(X[, svc.cols]))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
  expect_equal(out$n.val, n.val)
  expect_equal(out$d, d)
  expect_equal(unname(out$X), unname(X))
})

# Check grid.index errors -------------
test_that("grid.index error handling works", {
  fit.sp <- function(data) {
    svcFpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data,
               svc.cols = svc.cols,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE)
  }
  tmp.data <- data.list
  tmp.data$grid.index <- NULL
  expect_error(fit.sp(tmp.data), "grid.index must be specified")
  tmp.data <- data.list
  tmp.data$grid.index <- grid[-1]
  expect_error(fit.sp(tmp.data), "length\\(data\\$grid.index\\)")
  tmp.data <- data.list
  tmp.data$grid.index[tmp.data$grid.index == 1] <- 2
  # The level check compares vectors of different lengths, so R also warns
  expect_error(suppressWarnings(fit.sp(tmp.data)), "must contain each integer")
  tmp.data <- data.list
  tmp.data$grid.index <- as.character(grid)
  expect_error(fit.sp(tmp.data), "must be a numeric vector")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check waicOcc -----------------------
test_that("waicOCC works for svcFpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for svcFpPGOcc with grid.index.0", {
  n.post.samples <- out$n.post * out$n.chains
  pred.out <- predict(out, X, coords, grid.index.0 = grid, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(n.post.samples, J))
  expect_equal(dim(pred.out$z.0.samples), c(n.post.samples, J))
  expect_equal(dim(pred.out$w.0.samples), c(n.post.samples, p.svc, J.w))
})

# Check SVC samples -------------------
test_that("getSVCSamples works with grid.index", {
  svc.samples <- getSVCSamples(out)
  expect_equal(length(svc.samples), p.svc)
  expect_equal(dim(svc.samples[[1]]), c(out$n.post * out$n.chains, J.w))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for svcFpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
  expect_equal(dim(out.resids$occ.resids), c(10, J))
})
