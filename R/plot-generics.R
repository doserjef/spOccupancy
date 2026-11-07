overallPlot <- function(x, param, density = TRUE, ...) {
  n.post <- x$n.post
  n.chains <- x$n.chains
  curr.samples <- vector(mode = 'list', length = n.chains) 
  indx <- 1:n.post
  # Models with a false positive detection process (true positive rate mu
  # with coefficients gamma, false positive rate omega with coefficients delta).
  fp.classes <- c('fpPGOcc', 'spFpPGOcc', 'svcFpPGOcc')
  if (param == 'beta.comm') {
    if (class(x) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc',
			'tPGOcc', 'stPGOcc', 'svcPGOcc', 'svcTPGOcc',
			'svcTPGBinom', 'svcPGBinom', fp.classes)) {
      stop("beta.comm is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$beta.comm.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'tau.sq.beta') {
    if (class(x) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc', 
			'tPGOcc', 'stPGOcc', 'svcPGOcc', 'svcTPGOcc', 
			'svcTPGBinom', 'svcPGBinom', fp.classes)) {
      stop("tau.sq.beta is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$tau.sq.beta.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'lambda') {
    if (class(x) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc', 
			'tPGOcc', 'stPGOcc', 'svcPGOcc', 'svcTPGOcc', 
			'svcTPGBinom', 'svcPGBinom', 'msPGOcc', 'spMsPGOcc', 'intMsPGOcc', 
			'tMsPGOcc', fp.classes)) {
      stop("lambda is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$lambda.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'alpha.comm') {
    if (class(x) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc', 
			'tPGOcc', 'stPGOcc', 'svcPGOcc', 'svcTPGOcc', 
			'svcTPGBinom', 'svcPGBinom', 'lfJSDM', 'sfJSDM', fp.classes)) {
      stop("alpha.comm is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$alpha.comm.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'tau.sq.alpha') {
    if (class(x) %in% c('PGOcc', 'spPGOcc', 'intPGOcc', 'spIntPGOcc', 
			'tPGOcc', 'stPGOcc', 'svcPGOcc', 'svcTPGOcc', 
			'svcTPGBinom', 'svcPGBinom', 'lfJSDM', 'sfJSDM', fp.classes)) {
      stop("tau.sq.alpha is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$tau.sq.alpha.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'theta') {
    if (!(class(x) %in% c('spPGOcc', 'spIntPGOcc', 'tPGOcc', 'stPGOcc', 
			  'svcPGOcc', 'svcPGBinom', 'svcTPGOcc', 'svcTPGBinom', 
			  'spMsPGOcc', 'sfJSDM', 'sfMsPGOcc', 'stMsPGOcc', 'svcMsPGOcc', 
			  'svcTMsPGOcc', 'tIntPGOcc', 'stIntPGOcc', 'svcTIntPGOcc',
			  'spFpPGOcc', 'svcFpPGOcc'))) {
      stop("theta is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$theta.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'beta') {
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$beta.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'alpha') {
    if (class(x) %in% c('svcPGBinom', 'svcTPGBinom', 'lfJSDM', 'sfJSDM', fp.classes)) {
      stop("alpha is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$alpha.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'sigma.sq.psi') {
    if (!x$psiRE) {
      stop("sigma.sq.psi is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$sigma.sq.psi.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'sigma.sq.p') {
    if (class(x) %in% c('svcPGBinom', 'svcTPGBinom', 'lfJSDM', 'sfJSDM', fp.classes)) {
      stop("sigma.sq.p is not a parameter in the fitted model")
    }
    if (!x$pRE) {
      stop("sigma.sq.p is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$sigma.sq.p.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'beta.star') {
    if (!x$psiRE) {
      stop("the model was not fit with any occurrence random effects (beta.star)")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$beta.star.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'alpha.star') {
    if (class(x) %in% c('svcPGBinom', 'svcTPGBinom', 'lfJSDM', 'sfJSDM', fp.classes)) {
      stop("alpha.star is not a parameter in the fitted model")
    }
    if (!x$pRE) {
      stop("the model was not fit with any detection random effects (alpha.star)")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$alpha.star.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  # False positive model parameters ---------------------------------------
  if (param == 'gamma') {
    if (!(class(x) %in% fp.classes)) {
      stop("gamma is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$gamma.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'delta') {
    if (!(class(x) %in% fp.classes)) {
      stop("delta is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$delta.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'sigma.sq.mu') {
    if (!(class(x) %in% fp.classes) || !x$muRE) {
      stop("sigma.sq.mu is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$sigma.sq.mu.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'sigma.sq.omega') {
    if (!(class(x) %in% fp.classes) || !x$omegaRE) {
      stop("sigma.sq.omega is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$sigma.sq.omega.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'gamma.star') {
    if (!(class(x) %in% fp.classes) || !x$muRE) {
      stop("the model was not fit with any true positive detection random effects (gamma.star)")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$gamma.star.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'delta.star') {
    if (!(class(x) %in% fp.classes) || !x$omegaRE) {
      stop("the model was not fit with any false positive detection random effects (delta.star)")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$delta.star.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  curr.samples <- coda::mcmc.list(curr.samples)
  plot(curr.samples, density = density)
}

plot.PGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.spPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.intPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.spIntPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.tPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.stPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcPGBinom <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcTPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcTPGBinom <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.msPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.spMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.lfJSDM <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.sfJSDM <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.lfMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.sfMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.intMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.tMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.stMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcTMsPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.tIntPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.stIntPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcTIntPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.fpPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.spFpPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcFpPGOcc <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
