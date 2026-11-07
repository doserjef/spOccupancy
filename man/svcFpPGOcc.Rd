\name{svcFpPGOcc}
\alias{svcFpPGOcc}
\title{Function for Fitting Single-Species Spatially Varying Coefficient False Positive Occupancy Models}

\usage{
svcFpPGOcc(occ.formula, tp.det.formula, fp.det.formula, data, inits, priors, tuning, 
           svc.cols = 1, cov.model = 'exponential', n.neighbors = 15, search.type = 'cb', 
           n.batch, batch.length, accept.rate = 0.43, n.omp.threads = 1, verbose = TRUE, 
           n.report = 100, n.burn = round(.10 * n.batch * batch.length), 
           n.thin = 1, n.chains = 1, ...)
}

\description{
  Function for fitting single-species spatially varying coefficient false positive occupancy models using Polya-Gamma 
  latent variables. The function assumes data arise from some automated algorithm that 
  processed remote sensing data (e.g., acoustic recording units, camera traps) to generate 
  a number of observations for each site x replicate combination (i.e., count data). 
  A human observer is assumed to verify a subset of the observations to help separately 
  identify the false positive and true positive detection rates. Spatial random effects
  are estimated using a Nearest Neighbor Gaussian Process. 
}

\arguments{
  \item{occ.formula}{a symbolic description of the model to be fit
    for the occurrence portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random intercepts
    are allowed using lme4 syntax (Bates et al. 2015).}
  
  \item{tp.det.formula}{a symbolic description of the model to be fit
    for the true positive rate detection portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random intercepts
    are allowed using lme4 syntax (Bates et al. 2015).}
  
  \item{fp.det.formula}{a symbolic description of the model to be fit
    for the false positive rate detection portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random intercepts
    are allowed using lme4 syntax (Bates et al. 2015).}
  
  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{occ.covs}, \code{tp.det.covs}, \code{fp.det.covs}, 
    \code{n.val}, \code{d}, \code{coords}, and \code{grid.index}. \code{y} is a matrix or data frame with 
    first dimension equal to the number of sites (\eqn{J}{J}) and second dimension equal to the maximum number of 
    replicates at a given site. \code{y} contains the number of automated detections
    of the target species at each site x replicate combination. \code{occ.covs} is a matrix or data frame
    containing the variables used in the occurrence portion of the model, with
    \eqn{J}{J} rows for each column (variable). \code{tp.det.covs} is a list of
    variables included in the true positive rate detection portion of the model. Each list element
    is a different detection covariate, which can be site-level or 
    observational-level. Site-level covariates are specified as a vector 
    of length \eqn{J}{J} while observation-level covariates are specified 
    as a matrix or data frame with the number of rows equal to \eqn{J}{J} and 
    number of columns equal to the maximum number of replicates at a given site. 
    \code{fp.det.covs} takes the same form as \code{tp.det.covs} but contains covarites for the 
    false positive rate detection portion of the model. \code{n.val} is a two-dimensional 
    matrix with rows corresponding to sites and columns corresponding to replicates that contains
    the number of observations for each site x replicate cell that were validated by a human 
    observer. \code{d} is a two-dimensional matrix with rows corresponding to sites and columns
    corresponding to replicates that contains the number of the \code{n.val} observations that 
    were validated to be true during the manual validation process. \code{coords} is a matrix of the 
    observation coordinates used to estimate the spatial random effect for each site. \code{coords} has two columns for the 
    easting and northing coordinate, respectively. Typically, each site in the data 
    set will have it's own coordinate, such that \code{coords} is a \eqn{J \times 2}{J x 2} 
    matrix and \code{grid.index} should not be specified. If you desire to estimate 
    SVCs at some larger spatial level, e.g., if points fall within 
    grid cells and you want to estimate a spatial random effect for each grid cell 
    instead of each point, \code{coords} can be specified as the coordinate for each 
    grid cell. In such a case, \code{grid.index} is an indexing vector of length J, where each 
    value of \code{grid.index} indicates the corresponding row in \code{coords} that the given 
    site corresponds to. Note that \code{spOccupancy} assumes coordinates are specified 
    in a projected coordinate system. }
  
  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{z}, \code{beta}, \code{gamma}, \code{sigma.sq.psi}, 
    \code{sigma.sq.mu}, \code{sigma.sq.omega}, \code{sigma.sq}, \code{phi}, and \code{nu}. 
    The value portion of each tag is the parameter's initial value. \code{nu} is only specified if \code{cov.model = "matern"}. 
    \code{sigma.sq.psi}, \code{sigma.sq.mu}, and \code{sigma.sq.omega} are 
    only relevant when including random effects in the occurrence, true positive rate, and 
    false positive rate portions of the occupancy model, respectively. See \code{priors} 
    description for definition of each parameter name.
    Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}
  
  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal}, \code{gamma.normal}, \code{delta.normal}, 
    \code{sigma.sq.psi.ig}, \code{sigma.sq.mu.ig}, \code{sigma.sq.omega.ig}, 
    \code{sigma.sq.ig}, \code{sigma.sq.unif}, \code{nu.unif}, and \code{phi.unif}. 
    Occupancy (\code{beta}), true positive rate (\code{omega}), and false positive rate (\code{delta}) 
    regression coefficients are assumed to follow a normal distribution. 
    The hyperparameters of the normal distribution are passed as a list of 
    length two with the first and second elements corresponding to the mean 
    and variance of the normal distribution, which are each specified as vectors of 
    length equal to the number of coefficients to be estimated or of length one
    if priors are the same for all coefficients. If not specified, prior means 
    are set to 0 and prior variances set to 2.72 for the occupancy coefficients and 10 for 
    the true positive and false positive detection rate coefficients. \code{sigma.sq.psi}, \code{sigma.sq.mu}, and  
    \code{sigma.sq.omega} are the random effect variances for any occurrence, true positive rate, 
    or false positive rate random effects, respectively, and are assumed to follow an inverse
    Gamma distribution. The hyperparameters of the inverse-Gamma distribution
    are passed as a list of length two with first and second elements corresponding
    to the shape and scale parameters, respectively, which are each specified as
    vectors of length equal to the number of random intercepts or of length one 
    if priors are the same for all random effect variances. The spatial variance parameter, \code{sigma.sq}, is assumed to follow an 
    inverse-Gamma distribution or a uniform distribution (default is inverse-Gamma). 
    \code{sigma.sq} can also be fixed at its initial value by setting the prior value to 
    \code{"fixed"}. The spatial decay \code{phi} and 
    smoothness \code{nu} parameters are assumed to follow Uniform 
    distributions. The hyperparameters of the inverse-Gamma for \code{sigma.sq}
    are passed as a list with two elements corresponding to the 
    shape and scale parameters, respetively, with each element comprised of 
    a vector equal to the number of spatially-varying coefficients to be estimated
    or of length one if priors are the same for all coefficients.
    The hyperparameters of any uniform priors are also passed as a list of 
    length two with the first and second elements corresponding to 
    the lower and upper support, respectively, which can be passed as
    a vector equal to the total number of spatially-varying coefficients to 
    be estimated or of length one if priors are the same for all coefficients.}
  
  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{gamma}, \code{delta}, \code{sigma.sq.mu}, 
    \code{sigma.sq.omega}, \code{gamma.star} (the true-positive rate random intercept values), 
    \code{delta.star} (the false-positive rate random intercept values), \code{phi}, \code{nu},
    and \code{sigma.sq} (if using a uniform prior for sigma.sq).
    The value portion of each tag defines the initial variance of the Adaptive sampler. Note 
    that the introduction of count data through \code{y} and \code{d} eliminates the use of 
    Polya-Gamma data augmentation for the detection model parameter updates, and so all detection 
    model parameters are updated with Metropolis (or Metropolis-Hastings) updates using an 
    Adaptive sampler. See Roberts and Rosenthal (2009) for details.}

  \item{svc.cols}{a vector indicating the variables whose effects will be 
    estimated as spatially-varying coefficients. \code{svc.cols} can be an 
    integer vector with values indicating the order of covariates specified 
    in the model formula (with 1 being the intercept if specified), or it can 
    be specified as a character vector with names corresponding to variable 
    names in \code{occ.covs} (for the intercept, use \code{'(Intercept)'}). \code{svc.cols}
    default argument of 1 results in a spatial occupancy model analogous to 
    \code{spFpPGOcc} (assuming an intercept is included in the model).} 

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}

  \item{n.neighbors}{number of neighbors used in the NNGP. Datta et al. (2016) showed that 15 neighbors is usually 
    sufficient, but that as few as 5 neighbors can be adequate for certain data
    sets, which can lead to even greater decreases in run time. We recommend
    starting with 15 neighbors (the default) and if additional gains in computation
    time are desired, subsequently compare the results with a smaller number
    of neighbors using WAIC or k-fold cross-validation.}
  
  \item{search.type}{a quoted keyword that specifies the type of nearest
    neighbor search algorithm. Supported method key words are: \code{"cb"} and
    \code{"brute"}. The \code{"cb"} should generally be much
    faster. If locations do not have identical coordinate values on
    the axis used for the nearest neighbor ordering then \code{"cb"} 
    and \code{"brute"} should produce identical neighbor sets. 
    However, if there are identical coordinate values on the axis used 
    for nearest neighbor ordering, then \code{"cb"} and \code{"brute"} 
    might produce different, but equally valid, neighbor sets, 
    e.g., if data are on a grid. }

  \item{n.batch}{the number of MCMC batches in each chain to run for the Adaptive MCMC
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch in each chain to run for the Adaptive
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is
    0.43. See Roberts and Rosenthal (2009) for details.}
  
  \item{n.omp.threads}{a positive integer indicating the number of threads
    to use for SMP parallel processing within-chains. This will have no impact
    on model run time for non-spatial models. The package must be compiled for 
    OpenMP support. For most Intel-based machines, we recommend setting 
    \code{n.omp.threads} up to the number of hypterthreaded cores. Note, 
    \code{n.omp.threads} > 1 might not work on some systems. Currently only 
    relevant for spatial models.}
  
  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}
  
  \item{n.report}{the interval to report MCMC progress.}
  
  \item{n.burn}{the number of samples out of the total \code{n.samples} to 
    discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}
  
  \item{n.chains}{the number of chains to run.}
  
  \item{...}{currently no additional arguments}
}

\note{
 Some of the underlying code used for generating random numbers from the
 Polya-Gamma distribution is taken from the \pkg{pgdraw} package written
 by Daniel F. Schmidt and Enes Makalic. Their code implements Algorithm
 6 in PhD thesis of Jesse Bennett Windle (2013) \url{https://repositories.lib.utexas.edu/handle/2152/21842}.
}

\references{

  Chambert, T., Waddle, J. H., Miller, D. A., Walls, S. C., & Nichols, J. D. (2018). 
  A new framework for analysing automated acoustic species detection data: 
  Occupancy estimation and optimization of recordings post‐processing. 
  Methods in Ecology and Evolution, 9(3), 560-570.

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
    Association}, \doi{10.1080/01621459.2015.1044091}.

  Doser, J. W., Finley, A. O., Saunders, S. P., Kery, M., Weed, A. S., &
  Zipkin, E. F. (2025). Modeling complex species-environment
  relationships through spatially-varying coefficient occupancy models.
  Journal of Agricultural, Biological and Environmental Statistics.
  \doi{10.1007/s13253-023-00595-6}.

  Doser, J. W., Kery, M., Saunders, S. P., Finley, A. O., Bateman, B. L.,
  Grand, J., Reault, S., Weed, A. S., & Zipkin, E. F. (2024). Guidelines
  for the use of spatially varying coefficients in species distribution
  models. Global Ecology and Biogeography, 33(4), e13814.
  \doi{10.1111/geb.13814}.
  
  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
    Statistics}, \doi{10.1080/10618600.2018.1537924}.
  
  MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, 
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy 
  Rates When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. 
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC. 
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.

}

\author{
  Jeffrey W. Doser \email{jwdoser@ncsu.edu}
}

\value{
  An object of class \code{svcFpPGOcc} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the occupancy regression coefficients.}

  \item{gamma.samples}{a \code{coda} object of posterior samples
    for the true positive rate regression coefficients.}
  
  \item{delta.samples}{a \code{coda} object of posterior samples
    for the false positive rate regression coefficients.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for covariance parameters.}

  \item{w.samples}{a three-dimensional array of posterior samples for the 
        latent spatial random effects for all spatially-varying coefficients. 
        Dimensions correspond to MCMC sample, coefficient, and sites.}

  \item{z.samples}{a \code{coda} object of posterior samples 
    for the latent occupancy values}

  \item{psi.samples}{a \code{coda} object of posterior samples
    for the latent occupancy probability values}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occupancy portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the true positive rate portion 
    of the model. Only included if random intercepts are specified in 
    \code{tp.det.formula}.}

  \item{sigma.sq.omega.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the false positive rate portion 
    of the model. Only included if random intercepts are specified in 
    \code{fp.det.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the occurrence random effects. Only included if random intercepts
    are specified in \code{occ.formula}.}

  \item{gamma.star.samples}{a \code{coda} object of posterior samples
    for the true positive rate random effects. Only included if random intercepts 
    are specified in \code{tp.det.formula}.}
  
  \item{delta.star.samples}{a \code{coda} object of posterior samples
    for the false positive rate random effects. Only included if random intercepts 
    are specified in \code{fp.det.formula}.}

  \item{like.samples}{a \code{coda} object of posterior samples
    for the likelihood value associated with each site. Used for calculating
    WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 
}

\examples{
set.seed(8372)
J.x <- 5
J.y <- 5
J <- J.x * J.y
n.rep <- sample(3, J, replace = TRUE)
n.rep.max <- max(n.rep)
beta <- c(0, 1.3)
p.occ <- length(beta)
# Covariates on true positive detection rate 
gamma <- c(1, -0.5)
p.tp <- length(gamma)
# False positive rate
delta <- c(2, 0.3, -0.8)
p.fp <- length(delta)

# Proportion of samples manually validated
prop.val <- 0.05

# Spatial information
svc.cols <- c(1, 2)
phi <- c(3/.6, 3 / 0.8)
sigma.sq <- c(2, 0.8)
sp <- TRUE
cov.model <- "exponential"

# Simulate the data
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, sp = TRUE, prop.val = prop.val, family = "Poisson", 
                phi = phi, sigma.sq = sigma.sq, cov.model = cov.model, svc.cols = svc.cols)

# Extract necessary data from the data list
y <- dat$y
# Occupancy covariates
X <- dat$X
# True positive rate covariates
X.mu <- dat$X.mu
# False positive rate covariates
X.omega <- dat$X.omega
# Human validation data
n.val <- dat$n.val
# True positive validated data
d <- dat$d
# Coordinates
coords <- dat$coords

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1')
# True positive covariates
tp.det.covs <- list(tp.cov.1 = X.mu[, , 2])
# False positive covariates
fp.det.covs <- list(fp.cov.1 = X.omega[, , 2], fp.cov.2 = X.omega[, , 3])
data.list <- list(y = y,
                  occ.covs = occ.covs,
                  tp.det.covs = tp.det.covs,
                  fp.det.covs = fp.det.covs,
                  coords = dat$coords,
                  n.val = n.val,
                  d = d)

# Priors
prior.list <- list(beta.normal = list(0, 2.72),
                   gamma.normal = list(mean = rep(0, p.tp),
                                       var = rep(10, p.tp)),
                   delta.normal = list(mean = rep(0, p.fp),
                                       var = rep(10, p.fp)), 
                   sigma.sq.ig = list(2, 1), 
                   phi.unif = list(3 / 1, 3 / .1))

# Initial values
inits.list <- list(gamma = 0, beta = 0, delta = 0, phi = 3 / 0.5, 
                   sigma.sq = 1)

# Tuning
tuning.list <- list(delta = 0.5, gamma = 0.5, phi = 0.5)

# Run a very short test run of the model. The full model would need
# to run longer for convergence. 
out <- svcFpPGOcc(occ.formula = ~ occ.cov.1,
               tp.det.formula = ~ tp.cov.1,
               fp.det.formula = ~ fp.cov.1 + fp.cov.2,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               n.batch = 5,
               tuning = tuning.list,
               n.neighbors = 15,
               svc.cols = c(1, 2),
               cov.model = "exponential", 
               batch.length = 25,
               accept.rate = 0.43,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = 1,
               n.chains = 1,
               n.burn = 0,
               n.thin = 1)

summary(out)

}
