\name{simBinom}
\alias{simBinom}
\title{Simulate Single-Species Binomial Data}

\description{
  The function \code{simBinom} simulates single-species binomial data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the model. Non-spatial random intercepts can also be included in the model.
}

\usage{
simBinom(J.x, J.y, weights, beta, psi.RE = list(), 
         sp = FALSE, svc.cols = 1, cov.model, sigma.sq, phi, nu, 
         x.positive = FALSE, ...)
}

\arguments{
\item{J.x}{a single numeric value indicating the number of sites to simulate data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}.}

\item{J.y}{a single numeric value indicating the number of sites to simulate data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}.}

\item{weights}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of Bernoulli trials at each of the \eqn{J}{J} sites.}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the model.}

\item{psi.RE}{a list used to specify the non-spatial random intercepts included in the model. The list must have two tags: \code{levels} and \code{sigma.sq.psi}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.psi} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the model.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process. By default set to \code{FALSE}.}

\item{svc.cols}{a vector indicating the variables whose effects will be
  estimated as spatially-varying coefficients. \code{svc.cols} is an
  integer vector with values indicating the order of covariates specified
  in the model formula (with 1 being the intercept if specified).}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}. If \code{svc.cols} has more than one value, there should be a distinct spatial variance parameter for each spatially-varying coefficient.}

\item{phi}{a numeric value indicating the spatial decay parameter. Ignored when \code{sp = FALSE}. If \code{svc.cols} has more than one value, there should be a distinct spatial decay parameter for each spatially-varying coefficient.}

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}. If \code{svc.cols} has more than one value, there should be a distinct spatial smoothness parameter for each spatially-varying coefficient.} 

\item{x.positive}{a logical value indicating whether the simulated covariates should be simulated as random standard normal covariates (\code{x.positive = FALSE}) or restricted to positive values using a uniform distribution with lower bound 0 and upper bound 1 (\code{x.positive = TRUE}).}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times p.occ}{J x p.occ} numeric design matrix for the model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each occupancy site. Required for spatial models.}
  \item{w}{a matrix of the spatial random effect values for each site. The number of columns is determined by the \code{svc.cols} argument (the number of spatially-varying coefficients).}
  \item{psi}{a \eqn{J \times 1}{J x 1} matrix of the binomial probabilities for each site.}
  \item{y}{a length \code{J} vector of the binomial data for each site.}
  \item{X.w}{a two dimensional matrix containing the covariate effects (including an intercept) whose effects are assumed to be spatially-varying. Rows correspond to sites and columns correspond to covariate effects.}
  \item{X.re}{a numeric matrix containing the levels of any unstructured random effect included in the model. Only relevant when random effects are specified in \code{psi.RE}.}
  \item{beta.star}{a numeric vector that contains the simulated random effects for each given level of the random effects included in the model. Only relevant when random effects are included in the model.}
}

\examples{
set.seed(400)
J.x <- 10
J.y <- 10
weights <- rep(4, J.x * J.y)
beta <- c(0.5, -0.15)
svc.cols <- c(1, 2)
phi <- c(3 / .6, 3 / 0.2)
sigma.sq <- c(1.2, 0.9)
psi.RE <- list(levels = 10, 
               sigma.sq.psi = 1.2)
dat <- simBinom(J.x = J.x, J.y = J.y, weights = weights, beta = beta, 
                psi.RE = psi.RE, sp = TRUE, svc.cols = svc.cols, 
                cov.model = 'spherical', sigma.sq = sigma.sq, phi = phi)
}
