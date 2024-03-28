\name{postHocLM}
\alias{postHocLM}
\title{Function for Fitting Linear Mixed Models with Previous Model Estimates}

\usage{
postHocLM(formula, data, inits, priors, verbose = FALSE, 
          n.report = 100, n.samples, n.chains = 1, ...)
}

\description{
  Function for fitting a linear (mixed) model as a second-stage model where
  the response variable itself comes from a previous model fit and has 
  uncertainty associated with it. The response variable is assumed to be 
  a set of estimates from a previous model fit, where each value in the 
  response variable has a posterior MCMC sample of estimates. This function
  is useful for doing "posthoc" analyses of model estimates (e.g., 
  exploring how species traits relate to species-specific parameter estimates
  from a multi-species occupancy model). Such analyses are sometimes referred
  to as "two-stage" analyses.
}

\arguments{
\item{formula}{a symbolic description of the model to be fit
  for the model using R's model syntax. Only right-hand side of 
  formula is specified. See example below. Random intercepts
  are allowed using lme4 syntax (Bates et al. 2015).}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y} and \code{covs}. \code{y}
  is a matrix or data frame with first dimension equal to the number of 
  posterior samples of each value in the response variable and the 
  second dimension is equal to the number of values in the response
  variable. For example, if the response is species-specific covariate 
  effect estimates from a multi-species occupancy model, the rows correspond
  to the posterior MCMC samples and the columns correspond to species. \code{covs}
  is a matrix or data frame containing the independent variables used in the model. 
  Note the number of rows of \code{covs} should be equal to the number of columns
  in \code{y}.} 
  
\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{beta}, \code{tau.sq}, and \code{sigma.sq}.
  The value portion of each tag is the 
  parameter's initial value. \code{sigma.sq} is 
  only relevant when including random effects in the model. See \code{priors} 
  description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal}, \code{tau.sq.ig}, and \code{sigma.sq.ig}. 
  Regression coefficients (\code{beta}) are assumed to follow a normal distribution. 
  The hyperparameters of the normal distribution are passed as a list of 
  length two with the first and second elements corresponding to the mean 
  and variance of the normal distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not specified, prior means 
  are set to 0 and prior variances set to 100. \code{tau.sq} is the residual variance, 
  and is assumed to follow an inverse-Gamma distribution. 
  The hyperparameters of the inverse-Gamma distribution
  are passed as a vector of length two with first and second elements corresponding
  to the shape and scale parameters, respectively. \code{sigma.sq} are 
  the variances of any random intercepts included in the model, which similarly to 
  \code{tau.sq} follow an inverse-Gamma distribution. The hyperparameters of the 
  inverse-Gamma distribution are passed as a list of length two with first and second
  elements corresponding to the shape and scale parameters, respectively, which are 
  each specified as vectors of length equal to the number of random intercepts or 
  of length one if priors are the same for all random effect variances.}

\item{verbose}{if \code{TRUE}, messages about data preparation, 
  model specification, and progress of the sampler are printed to the screen. 
  Otherwise, no messages are printed.}

\item{n.report}{the interval to report MCMC progress.}

\item{n.samples}{the number of posterior samples to collect in each chain. Note that
  by default, the same number of MCMC samples fit in the first stage model is 
  assumed to be fit for the second stage model. If \code{n.samples} is specified, 
  it must be a multiple of the number of samples fit in the first stage, otherwise
  an error will be reported.}

\item{n.chains}{the number of chains to run in sequence.}

\item{...}{currently no additional arguments}
}

\references{

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  An object of class \code{postHocLM} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the regression coefficients.}

  \item{tau.sq.samples}{a \code{coda} object of posterior samples
    for the residual variances.}

  \item{y.hat.samples}{a \code{coda} object of posterior samples 
    of fitted values.}

  \item{sigma.sq.samples}{a \code{coda} object of posterior samples
    for the random effect variances if any random intercepts were included
    in the model.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the random effects. Only included if random intercepts
    are specified in \code{formula}.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  \item{bayes.R2}{a \code{coda} object of posterior samples 
    of the Bayesian R-squared as a measure of model fit. Note that when 
    random intercepts are included in the model, this is the conditional
    Bayesian R-squared, not the marginal Bayesian R-squared.}

  The return object will include additional objects used for 
  subsequent summarization. 

}

\examples{
# Simulate Data -----------------------------------------------------------
set.seed(100)
N <- 100
beta <- c(0, 0.5, 1.2)
tau.sq <- 1 
p <- length(beta)
X <- matrix(1, nrow = N, ncol = p)
if (p > 1) {
  for (i in 2:p) {
    X[, i] <- rnorm(N)
  } # i
}
mu <- X[, 1] * beta[1] + X[, 2] * beta[2] + X[, 3] * beta[3]
y <- rnorm(N, mu, sqrt(tau.sq))
# Replicate y n.samples times and add a small amount of noise that corresponds
# to uncertainty from a first stage model.
n.samples <- 1000
y <- matrix(y, n.samples, N, byrow = TRUE)
y <- y + rnorm(length(y), 0, 0.25)

# Package data for use with postHocLM -------------------------------------
colnames(X) <- c('int', 'cov.1', 'cov.2')
data.list <- list(y = y, covs = X)
data <- data.list
inits <- list(beta = 0, tau.sq = 1)
priors <- list(beta.normal = list(mean = 0, var = 10000),
               tau.sq.ig = c(0.001, 0.001))

# Run the model -----------------------------------------------------------
out <- postHocLM(formula = ~ cov.1 + cov.2, 
                 inits = inits, 
                 data = data.list, 
                 priors = priors, 
                 verbose = FALSE, 
                 n.chains = 1)
summary(out)
}
