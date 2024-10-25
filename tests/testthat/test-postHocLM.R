# Test postHocLM.R --------------------------------------------------------

skip_on_cran()


# Fixed effects only ------------------------------------------------------
set.seed(100)
N <- 100
beta <- c(0, 0.5, 1.2)
tau.sq <- 2 
p <- length(beta)
X <- matrix(1, nrow = N, ncol = p)
if (p > 1) {
  for (i in 2:p) {
    X[, i] <- rnorm(N)
  } # i
}
mu <- X %*% as.matrix(beta)
y <- rnorm(N, mu, sqrt(tau.sq))
y.small <- y
# Replicate y n.samples times 
n.samples <- 1000
y <- matrix(y, n.samples, N, byrow = TRUE)
y <- y + rnorm(length(y), 0, 0.5)

covs <- cbind(X)
colnames(covs) <- c('int', 'cov.1', 'cov.2')
data.list <- list(y = y, 
		  covs = covs)
inits <- list(beta = 0, tau.sq = 1, sigma.sq = 1)
priors <- list(beta.normal = list(mean = 0, var = 10000), 
	       tau.sq.ig = c(0.001, 0.001), 
	       sigma.sq.ig = list(a = 0.1, b = 0.1))
out <- postHocLM(formula = ~ cov.1 + cov.2, 
		 inits = inits, 
		 data = data.list, 
		 priors = priors, 
		 verbose = FALSE, 
		 n.samples = 1000,
		 n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class postHocLM", {
  expect_s3_class(out, "postHocLM")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$RE, FALSE)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- postHocLM(formula = ~ cov.1 + cov.2, 
		 data = data.list, 
		 verbose = FALSE, 
		 n.samples = 10000,
		 n.chains = 1)
  expect_s3_class(out, "postHocLM")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
expect_output(postHocLM(formula = ~ cov.1 + cov.2, 
		 inits = inits, 
		 data = data.list, 
		 priors = priors, 
		 verbose = TRUE, 
		 n.samples = 1000,
		 n.chains = 1))
})

# Fixed and random effects ------------------------------------------------
set.seed(100)
N <- 100
beta <- c(0, 0.5, 1.2)
tau.sq <- 2 
RE <- list(levels = c(20, 10),
           sigma.sq = c(2, 0.5))
p <- length(beta)
X <- matrix(1, nrow = N, ncol = p)
if (p > 1) {
  for (i in 2:p) {
    X[, i] <- rnorm(N)
  } # i
}
if (length(RE) > 0) {
  p.re <- length(RE$levels)
  sigma.sq <- rep(NA, p.re)
  n.re.long <- RE$levels
  n.re <- sum(n.re.long)
  beta.star.indx <- rep(1:p.re, n.re.long)
  beta.star <- rep(0, n.re)
  X.re <- matrix(NA, N, p.re)
  for (i in 1:p.re) {
    X.re[, i] <- sample(1:RE$levels[i], N, replace = TRUE)
    beta.star[which(beta.star.indx == i)] <- rnorm(RE$levels[i], 0, sqrt(RE$sigma.sq[i]))
  }
  if (p.re > 1) {
    for (j in 2:p.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
    }
  }
  beta.star.sites <- apply(X.re, 1, function(a) sum(beta.star[a]))
} else {
  X.re <- NA
  beta.star <- NA
  beta.star.sites <- rep(0, N)
}
mu <- X %*% as.matrix(beta) + beta.star.sites
y <- rnorm(N, mu, sqrt(tau.sq))
y.small <- y
# Replicate y n.samples times 
n.samples <- 1000
y <- matrix(y, n.samples, N, byrow = TRUE)
y <- y + rnorm(length(y), 0, 0.5)

covs <- cbind(X, X.re)
colnames(covs) <- c('int', 'cov.1', 'cov.2', 'factor.1', 'factor.2')
data.list <- list(y = y, 
		  covs = covs)
inits <- list(beta = 0, tau.sq = 1, sigma.sq = 1)
priors <- list(beta.normal = list(mean = 0, var = 10000), 
	       tau.sq.ig = c(0.001, 0.001), 
	       sigma.sq.ig = list(a = 0.1, b = 0.1))
out <- postHocLM(formula = ~ cov.1 + cov.2 + (1 | factor.1) + (1 | factor.2), 
		 inits = inits, 
		 data = data.list, 
		 priors = priors, 
		 verbose = FALSE, 
		 n.samples = 10000,
		 n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class postHocLM", {
  expect_s3_class(out, "postHocLM")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$RE, TRUE)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- postHocLM(formula = ~ cov.1 + cov.2 + (1 | factor.1) + (1 | factor.2), 
		 data = data.list, 
		 verbose = FALSE, 
		 n.samples = 10000,
		 n.chains = 1)
  expect_s3_class(out, "postHocLM")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
expect_output(postHocLM(formula = ~ cov.1 + cov.2 + (1 | factor.1) + (1 | factor.2), 
		 inits = inits, 
		 data = data.list, 
		 priors = priors, 
		 verbose = TRUE, 
		 n.samples = 1000,
		 n.chains = 1))
})
