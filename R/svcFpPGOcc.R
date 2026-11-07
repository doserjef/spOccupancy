svcFpPGOcc <- function(occ.formula, tp.det.formula, fp.det.formula,
                       data, inits, priors, tuning, svc.cols = 1,
                       cov.model = 'exponential', n.neighbors = 15,
                       search.type = 'cb', n.batch, batch.length,
                       accept.rate = 0.43, n.omp.threads = 1, verbose = TRUE,
                       n.report = 100, n.burn = round(.10 * n.batch * batch.length),
                       n.thin = 1, n.chains = 1, ...) {

  ptm <- proc.time()

  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing to run the model\n");
    cat("----------------------------------------\n");
  }

  # Functions ---------------------------------------------------------------
  rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
  }

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  # Call ----------------------------------------------------------------
  # Returns a call in which all of the specified arguments are
  # specified by their full names.
  cl <- match.call()
  
  # Some initial checks -------------------------------------------------
  if (missing(data)) {
    stop("data must be specified")
  }
  if (!is.list(data)) {
    stop("data must be a list")
  }
  names(data) <- tolower(names(data))
  if (missing(occ.formula)) {
    stop("occ.formula must be specified")
  }
  if (missing(tp.det.formula)) {
    stop("tp.det.formula must be specified")
  }
  if (missing(fp.det.formula)) {
    stop("fp.det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("detection data y must be specified in data")
  }
  y <- as.matrix(data$y)
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- matrix(1, dim(y)[1], 1)
    } else {
      stop("occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!is.matrix(data$occ.covs) & !is.data.frame(data$occ.covs)) {
    stop("occ.covs must be a matrix or data frame")
  }
  if (dim(y)[[1]] != dim(data$occ.covs)[[1]]){
    stop("number of sites in detection data 'y' (J = ", dim(y)[[1]],
         ") does not match number of sites in site-specific covariates 'occ.covs' (J = ", 
         dim(data$occ.covs)[[1]], ").")
  }
  
  # Check validation data
  if (!'n.val' %in% names(data)) {
    stop("human validation data n.val must be specified in data")
  }
  if (!'d' %in% names(data)) {
    stop("true positive human validation data d must be specified in data")
  }
  n.val <- as.matrix(data$n.val)
  d <- as.matrix(data$d)
  if (sum(dim(y) != dim(n.val)) > 0) {
    stop("dimension mismatch: number of sites and/or replicates does not match between y and n.val")
  }
  if (sum(dim(y) != dim(d)) > 0) {
    stop("dimension mismatch: number of sites and/or replicates does not match between y and d")
  }
  
  if (sum(is.na(data$occ.covs)) > 0) {
    stop("missing covariate values in data$occ.covs. Remove these sites from all data or impute non-missing values.")
  }
  if (!'tp.det.covs' %in% names(data)) {
    if (tp.det.formula == ~ 1) {
      if (verbose) {
        message("true positive detection covariates (tp.det.covs) not specified in data.\nAssuming intercept only true positive detection model.\n")
      }
      data$tp.det.covs <- list(int = matrix(1, dim(y)[1], dim(y)[2]))
    } else {
      stop("tp.det.covs must be specified in data for a true positive detection model with covariates")
    }
  }
  if (!is.list(data$tp.det.covs)) {
    stop("tp.det.covs must be a list of matrices, data frames, and/or vectors")
  }
  if (!'fp.det.covs' %in% names(data)) {
    if (fp.det.formula == ~ 1) {
      if (verbose) {
        message("false positive detection covariates (fp.det.covs) not specified in data.\nAssuming intercept only false positive detection model.\n")
      }
      data$fp.det.covs <- list(int = matrix(1, dim(y)[1], dim(y)[2]))
    } else {
      stop("fp.det.covs must be specified in data for a true positive detection model with covariates")
    }
  }
  if (!is.list(data$fp.det.covs)) {
    stop("fp.det.covs must be a list of matrices, data frames, and/or vectors")
  }
  if (!'coords' %in% names(data)) {
    stop("coords must be specified in data for a spatial occupancy model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("coords must be a matrix or data frame")
  }
  coords <- as.matrix(data$coords)
  # Check grid index
  if (!'grid.index' %in% names(data)) {
    if (nrow(data$coords) != nrow(data$y)) {
      stop("data$grid.index must be specified if nrow(data$coords) != nrow(data$y)")
    }
    grid.index <- 1:nrow(coords)
  } else {
    if (!is.atomic(data$grid.index) | !is.numeric(data$grid.index)) {
      stop("data$grid.index must be a numeric vector")
    }
    if (length(data$grid.index) != nrow(y)) {
      stop("length(data$grid.index) must equal the number of sites in y")
    }
    if (!all(sort(unique(data$grid.index)) == 1:nrow(coords))) {
      stop("data$grid.index must contain each integer from 1 to nrow(data$coords) at least once")
    }
    grid.index <- data$grid.index
  }
  # Check if all spatial coordinates are unique.
  unique.coords <- unique(data$coords)
  if (nrow(unique.coords) < nrow(data$coords)) {
    stop("coordinates provided in coords are not all unique. spOccupancy requires each site to have its own unique pair of spatial coordinates. This may be the result of an error in preparing the data list, or you will need to change what you consider a 'site' in order to meet this requirement. Alternatively, you can use data$grid.index to specify the SVC at a larger spatial level than the individual sites (e.g., a grid). See ?svcFpPGOcc for details.")
  }

  if (missing(n.batch)) {
    stop("must specify number of MCMC batches")
  }
  if (missing(batch.length)) {
    stop("must specify length of each MCMC batch")
  }
  n.samples <- n.batch * batch.length
  if (n.burn > n.samples) {
    stop("n.burn must be less than n.samples")
  }
  if (n.thin > n.samples) {
    stop("n.thin must be less than n.samples")
  }
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
  }
  
  # Neighbors and Ordering ----------------------------------------------
  u.search.type <- 2
  ## Order by x column. Could potentially allow this to be user defined.
  ord <- order(coords[,1])
  tmp <- lapply(ord, function (a) which(grid.index == a))
  tmp.2 <- sapply(tmp, length)
  grid.index.c <- unlist(lapply(1:length(tmp.2), function(a) rep(a, tmp.2[a]))) - 1
  long.ord <- unlist(lapply(ord, function(a) which(grid.index == a)))
  # Reorder everything to align with NN ordering
  y <- y[long.ord, , drop = FALSE]
  n.val <- n.val[long.ord, , drop = FALSE]
  d <- d[long.ord, , drop = FALSE]
  coords <- coords[ord, , drop = FALSE]
  # Occupancy covariates
  data$occ.covs <- data$occ.covs[long.ord, , drop = FALSE]
  # Detection covariates (site-level vectors or site x replicate matrices)
  for (i in 1:length(data$tp.det.covs)) {
    if (!is.null(dim(data$tp.det.covs[[i]]))) {
      data$tp.det.covs[[i]] <- data$tp.det.covs[[i]][long.ord, , drop = FALSE]
    } else {
      data$tp.det.covs[[i]] <- data$tp.det.covs[[i]][long.ord]
    }
  }
  for (i in 1:length(data$fp.det.covs)) {
    if (!is.null(dim(data$fp.det.covs[[i]]))) {
      data$fp.det.covs[[i]] <- data$fp.det.covs[[i]][long.ord, , drop = FALSE]
    } else {
      data$fp.det.covs[[i]] <- data$fp.det.covs[[i]][long.ord]
    }
  }

  # First subset all detection covariates to only use those that are included in the analysis.
  data$tp.det.covs <- data$tp.det.covs[names(data$tp.det.covs) %in% all.vars(tp.det.formula)]
  data$fp.det.covs <- data$fp.det.covs[names(data$fp.det.covs) %in% all.vars(fp.det.formula)]
  # Null model support
  if (length(data$tp.det.covs) == 0) {
    data$tp.det.covs <- list(int = matrix(1, dim(y)[1], dim(y)[2]))
  }
  if (length(data$fp.det.covs) == 0) {
    data$fp.det.covs <- list(int = matrix(1, dim(y)[1], dim(y)[2]))
  }
  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  # The likelihood is defined per survey, so expand any site-level covariates
  # (length J vectors) to J x K.max matrices before flattening like y.
  expand.site <- function(a) {
    if (is.null(dim(a)) && length(a) == nrow(y)) {
      a <- matrix(a, nrow(y), ncol(y))
    }
    unlist(c(a))
  }
  data$tp.det.covs <- data.frame(lapply(data$tp.det.covs, expand.site))
  data$fp.det.covs <- data.frame(lapply(data$fp.det.covs, expand.site))
  # Check if all detection covariates are at site level, and simplify the data
  # if necessary. This is highly unlikely in practice
  y.big <- y
  n.val.big <- n.val
  d.big <- d
  
  
  data$occ.covs <- as.data.frame(data$occ.covs)
  
  # Check whether random effects are sent in as numeric, and
  # return error if they are not. 
  # Occurrence ----------------------
  if (!is.null(findbars(occ.formula))) {
    occ.re.names <- sapply(findbars(occ.formula), all.vars)
    for (i in 1:length(occ.re.names)) {
      if (is(data$occ.covs[, occ.re.names[i]], 'factor')) {
        stop(paste("random effect variable ", occ.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$occ.covs[, occ.re.names[i]], 'character')) {
        stop(paste("random effect variable ", occ.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # TP Detection ----------------------
  if (!is.null(findbars(tp.det.formula))) {
    tp.det.re.names <- sapply(findbars(tp.det.formula), all.vars)
    for (i in 1:length(tp.det.re.names)) {
      if (is(data$tp.det.covs[, tp.det.re.names[i]], 'factor')) {
        stop(paste("random effect variable ", tp.det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$tp.det.covs[, tp.det.re.names[i]], 'character')) {
        stop(paste("random effect variable ", tp.det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # FP Detection ----------------------
  if (!is.null(findbars(fp.det.formula))) {
    fp.det.re.names <- sapply(findbars(fp.det.formula), all.vars)
    for (i in 1:length(fp.det.re.names)) {
      if (is(data$fp.det.covs[, fp.det.re.names[i]], 'factor')) {
        stop(paste("random effect variable ", fp.det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$fp.det.covs[, fp.det.re.names[i]], 'character')) {
        stop(paste("random effect variable ", fp.det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  
  # Checking missing values ---------------------------------------------
  # y -------------------------------
  y.na.test <- apply(y.big, 1, function(a) sum(!is.na(a)))
  if (sum(y.na.test == 0) > 0) {
    stop("some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
  }
  # occ.covs ------------------------
  if (sum(is.na(data$occ.covs)) != 0) {
    stop("missing values in occ.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
  }
  # tp.det.covs ------------------------
  for (i in 1:ncol(data$tp.det.covs)) {
    if (sum(is.na(data$tp.det.covs[, i])) > sum(is.na(y.big))) {
      stop("some elements in tp.det.covs have missing values where there is an observed data value in y. Please either replace the NA values in tp.det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  # Misalignment between y and det.covs
  y.missing <- which(is.na(y.big))
  tp.det.covs.missing <- lapply(data$tp.det.covs, function(a) which(is.na(a)))
  for (i in 1:length(tp.det.covs.missing)) {
    tmp.indx <- !(y.missing %in% tp.det.covs.missing[[i]])
    if (sum(tmp.indx) > 0) {
      if (i == 1 & verbose) {
        message("There are missing values in data$y with corresponding non-missing values in data$tp.det.covs.\nRemoving these site/replicate combinations for fitting the model.")
      }
      data$tp.det.covs[y.missing, i] <- NA
    }
  }
  # fp.det.covs ------------------------
  for (i in 1:ncol(data$fp.det.covs)) {
    if (sum(is.na(data$fp.det.covs[, i])) > sum(is.na(y.big))) {
      stop("some elements in fp.det.covs have missing values where there is an observed data value in y. Please either replace the NA values in fp.det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  # Misalignment between y and det.covs
  y.missing <- which(is.na(y.big))
  fp.det.covs.missing <- lapply(data$fp.det.covs, function(a) which(is.na(a)))
  for (i in 1:length(fp.det.covs.missing)) {
    tmp.indx <- !(y.missing %in% fp.det.covs.missing[[i]])
    if (sum(tmp.indx) > 0) {
      if (i == 1 & verbose) {
        message("There are missing values in data$y with corresponding non-missing values in data$fp.det.covs.\nRemoving these site/replicate combinations for fitting the model.")
      }
      data$fp.det.covs[y.missing, i] <- NA
    }
  }
  # Check missingness patterns in the validated data
  # Cells that are missing should be identical between d, n.val, and y. 
  n.val.missing <- which(is.na(n.val.big)) 
  d.missing <- which(is.na(d.big))
  if (!isTRUE(all.equal(n.val.missing, y.missing))) {
    stop("mismatch between missing values in n.val and y. These must have the same number and location of missing values")
  }
  if (!isTRUE(all.equal(n.val.missing, d.missing))) {
    stop("mismatch between missing values in n.val and d. These must have the same number and location of missing values")
  }
  # Check y, n.val, and d are non-negative integer counts with d <= n.val <= y.
  # Otherwise the Poisson/binomial likelihoods are -Inf and the sampler stalls.
  is.count <- function(a) {
    a <- a[!is.na(a)]
    all(a >= 0 & a == round(a))
  }
  if (!is.count(y.big)) {
    stop("y must contain non-negative integer counts")
  }
  if (!is.count(n.val.big)) {
    stop("n.val must contain non-negative integer counts")
  }
  if (!is.count(d.big)) {
    stop("d must contain non-negative integer counts")
  }
  if (sum(n.val.big > y.big, na.rm = TRUE) > 0) {
    stop("n.val cannot be larger than y. The number of validated signals must be less than or equal to the total number of signals")
  }
  if (sum(d.big > n.val.big, na.rm = TRUE) > 0) {
    stop("d cannot be larger than n.val. The number of validated true positives must be less than or equal to the number of validated signals")
  }

  # Formula -------------------------------------------------------------
  # Occupancy -----------------------
  if (is(occ.formula, 'formula')) {
    # NOTE: this will reorder random effect names by the number of unique 
    #       levels they have, not necessarily the order in which they are specified
    #       in the model formula. 
    tmp <- parseFormula(occ.formula, data$occ.covs)
    X <- as.matrix(tmp[[1]])
    X.re <- as.matrix(tmp[[4]])
    x.re.names <- colnames(X.re)
    x.names <- tmp[[2]]
  } else {
    stop("occ.formula is misspecified")
  }
  # Get RE level names
  re.level.names <- lapply(data$occ.covs[, x.re.names, drop = FALSE],
                           function (a) sort(unique(a)))
  
  # TP Detection ---------------------
  if (is(tp.det.formula, 'formula')) {
    tmp <- parseFormula(tp.det.formula, data$tp.det.covs)
    X.mu <- as.matrix(tmp[[1]])
    X.mu.re <- as.matrix(tmp[[4]])
    x.mu.re.names <- colnames(X.mu.re)
    x.mu.names <- tmp[[2]]
  } else {
    stop("tp.det.formula is misspecified")
  }
  mu.re.level.names <- lapply(data$tp.det.covs[, x.mu.re.names, drop = FALSE],
                             function (a) sort(unique(a)))
  # FP Detection ---------------------
  if (is(fp.det.formula, 'formula')) {
    tmp <- parseFormula(fp.det.formula, data$fp.det.covs)
    X.omega <- as.matrix(tmp[[1]])
    X.omega.re <- as.matrix(tmp[[4]])
    x.omega.re.names <- colnames(X.omega.re)
    x.omega.names <- tmp[[2]]
  } else {
    stop("fp.det.formula is misspecified")
  }
  omega.re.level.names <- lapply(data$fp.det.covs[, x.omega.re.names, drop = FALSE],
                                 function (a) sort(unique(a)))
  
  # Get basic info from inputs ------------------------------------------
  # Number of sites
  J <- nrow(y.big)
  # Number of coordinates in coordinate grid
  J.w <- nrow(coords)
  # Number of occupancy parameters
  p.occ <- ncol(X)
  # Number of occupancy random effect parameters
  p.occ.re <- ncol(X.re)
  # Number of tp detection parameters
  p.tp.det <- ncol(X.mu)
  # Number of tp detection random effect parameters
  p.tp.det.re <- ncol(X.mu.re)
  # Number of fp detection parameters
  p.fp.det <- ncol(X.omega)
  # Number of fp detection random effect parameters
  p.fp.det.re <- ncol(X.omega.re)

  # Check SVC columns -----------------------------------------------------
  if (is.character(svc.cols)) {
    # Check if all column names in svc are in occ.covs
    if (!all(svc.cols %in% x.names)) {
      missing.cols <- svc.cols[!(svc.cols %in% x.names)]
      stop(paste("variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
    }
    # Convert desired column names into the numeric column index
    svc.cols <- (1:p.occ)[x.names %in% svc.cols]
  } else if (is.numeric(svc.cols)) {
    # Check if all column indices are in 1:p.occ
    if (!all(svc.cols %in% 1:p.occ)) {
      missing.cols <- svc.cols[!(svc.cols %in% (1:p.occ))]
      stop(paste("column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
    }
  } else {
    stop("svc.cols must be a numeric or character vector")
  }
  if (any(duplicated(svc.cols))) {
    stop("svc.cols must not contain duplicate values")
  }
  p.svc <- length(svc.cols)
   
  # Number of latent occupancy random effect values
  n.occ.re <- length(unlist(apply(X.re, 2, unique)))
  n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of latent tp detection random effect values
  n.tp.det.re <- length(unlist(apply(X.mu.re, 2, unique)))
  n.tp.det.re.long <- apply(X.mu.re, 2, function(a) length(unique(a)))
  if (p.tp.det.re == 0) n.tp.det.re.long <- 0
  # Number of latent fp detection random effect values
  n.fp.det.re <- length(unlist(apply(X.omega.re, 2, unique)))
  n.fp.det.re.long <- apply(X.omega.re, 2, function(a) length(unique(a)))
  if (p.fp.det.re == 0) n.fp.det.re.long <- 0
  # Number of replicates at each site
  n.rep <- apply(y.big, 1, function(a) sum(!is.na(a)))
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- which(!is.na(y.big[j, ]))
  }
  # Max number of repeat visits
  K.max <- max(n.rep)
  # Because I like K better than n.rep
  K <- n.rep
  
  # Get indices to map z to y -------------------------------------------
  z.long.indx <- rep(1:J, dim(y.big)[2])
  z.long.indx <- z.long.indx[!is.na(c(y.big))]
  # Subtract 1 for indices in C
  z.long.indx <- z.long.indx - 1
  # y, n.val, d are stored in the following order: site, visit
  y <- c(y)
  names.long <- which(!is.na(y))
  n.val <- c(n.val)
  d <- c(d)
  # Remove missing observations when the covariate data are available but
  # there are missing detection-nondetection data. 
  if (nrow(X.mu) == length(y)) {
    X.mu <- X.mu[!is.na(y), , drop = FALSE]  
  }
  if (nrow(X.omega) == length(y)) {
    X.omega <- X.omega[!is.na(y), , drop = FALSE]
  }
  if (nrow(X.mu.re) == length(y) & p.tp.det.re > 0) {
    X.mu.re <- X.mu.re[!is.na(y), , drop = FALSE]
  }
  if (nrow(X.omega.re) == length(y) & p.fp.det.re > 0) {
    X.omega.re <- X.omega.re[!is.na(y), , drop = FALSE]
  }
  y <- y[!is.na(y)]
  n.val <- n.val[!is.na(n.val)]
  d <- d[!is.na(d)]
  # Number of data points for the y vector
  n.obs <- nrow(X.mu)
  
  # Get random effect matrices all set ----------------------------------
  if (p.occ.re > 1) {
    for (j in 2:p.occ.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }
  if (p.tp.det.re > 1) {
    for (j in 2:p.tp.det.re) {
      X.mu.re[, j] <- X.mu.re[, j] + max(X.mu.re[, j - 1]) + 1
    }
  }
  if (p.fp.det.re > 1) {
    for (j in 2:p.fp.det.re) {
      X.omega.re[, j] <- X.omega.re[, j] + max(X.omega.re[, j - 1]) + 1
    }
  }

  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))
  # beta -----------------------
  if ("beta.normal" %in% names(priors)) {
    if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
      stop("beta.normal must be a list of length 2")
    }
    mu.beta <- priors$beta.normal[[1]]
    sigma.beta <- priors$beta.normal[[2]]
    if (length(mu.beta) != p.occ & length(mu.beta) != 1) {
      if (p.occ == 1) {
        stop(paste("beta.normal[[1]] must be a vector of length ",
                   p.occ, " with elements corresponding to betas' mean", sep = ""))
      } else {
        stop(paste("beta.normal[[1]] must be a vector of length ",
                   p.occ, " or 1 with elements corresponding to betas' mean", sep = ""))
      }
    }
    if (length(sigma.beta) != p.occ & length(sigma.beta) != 1) {
      if (p.occ == 1) {
        stop(paste("beta.normal[[2]] must be a vector of length ",
                   p.occ, " with elements corresponding to betas' variance", sep = ""))
      } else {
        stop(paste("beta.normal[[2]] must be a vector of length ",
                   p.occ, " or 1 with elements corresponding to betas' variance", sep = ""))
      }
    }
    if (length(sigma.beta) != p.occ) {
      sigma.beta <- rep(sigma.beta, p.occ)
    }
    if (length(mu.beta) != p.occ) {
      mu.beta <- rep(mu.beta, p.occ)
    }
    Sigma.beta <- sigma.beta * diag(p.occ)
  } else {
    if (verbose) {
      message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.beta <- rep(0, p.occ)
    sigma.beta <- rep(2.72, p.occ)
    Sigma.beta <- diag(p.occ) * 2.72
  }
  # gamma -----------------------
  if ("gamma.normal" %in% names(priors)) {
    if (!is.list(priors$gamma.normal) | length(priors$gamma.normal) != 2) {
      stop("gamma.normal must be a list of length 2")
    }
    mu.gamma <- priors$gamma.normal[[1]]
    sigma.gamma <- priors$gamma.normal[[2]]
    if (length(mu.gamma) != p.tp.det & length(mu.gamma) != 1) {
      if (p.tp.det == 1) {
        stop(paste("gamma.normal[[1]] must be a vector of length ",
                   p.tp.det, " with elements corresponding to gammas' mean", sep = ""))
      } else {
        stop(paste("gamma.normal[[1]] must be a vector of length ",
                   p.tp.det, " or 1 with elements corresponding to gammas' mean", sep = ""))
      }
    }
    if (length(sigma.gamma) != p.tp.det & length(sigma.gamma) != 1) {
      if (p.tp.det == 1) {
        stop(paste("gamma.normal[[2]] must be a vector of length ",
                   p.tp.det, " with elements corresponding to gammas' variance", sep = ""))
      } else {
        stop(paste("gamma.normal[[2]] must be a vector of length ",
                   p.tp.det, " or 1 with elements corresponding to gammas' variance", sep = ""))
      }
    }
    if (length(sigma.gamma) != p.tp.det) {
      sigma.gamma <- rep(sigma.gamma, p.tp.det)
    }
    if (length(mu.gamma) != p.tp.det) {
      mu.gamma <- rep(mu.gamma, p.tp.det)
    }
    Sigma.gamma <- sigma.gamma * diag(p.tp.det)
  } else {
    if (verbose) {
      message("No prior specified for gamma.normal.\nSetting prior mean to 0 and prior variance to 10\n")
    }
    mu.gamma <- rep(0, p.tp.det)
    sigma.gamma <- rep(10, p.tp.det)
    Sigma.gamma <- diag(p.tp.det) * 10 
  }
  # delta -----------------------
  if ("delta.normal" %in% names(priors)) {
    if (!is.list(priors$delta.normal) | length(priors$delta.normal) != 2) {
      stop("delta.normal must be a list of length 2")
    }
    mu.delta <- priors$delta.normal[[1]]
    sigma.delta <- priors$delta.normal[[2]]
    if (length(mu.delta) != p.fp.det & length(mu.delta) != 1) {
      if (p.fp.det == 1) {
        stop(paste("delta.normal[[1]] must be a vector of length ",
                   p.fp.det, " with elements corresponding to deltas' mean", sep = ""))
      } else {
        stop(paste("delta.normal[[1]] must be a vector of length ",
                   p.fp.det, " or 1 with elements corresponding to deltas' mean", sep = ""))
      }
    }
    if (length(sigma.delta) != p.fp.det & length(sigma.delta) != 1) {
      if (p.fp.det == 1) {
        stop(paste("delta.normal[[2]] must be a vector of length ",
                   p.fp.det, " with elements corresponding to deltas' variance", sep = ""))
      } else {
        stop(paste("delta.normal[[2]] must be a vector of length ",
                   p.fp.det, " or 1 with elements corresponding to deltas' variance", sep = ""))
      }
    }
    if (length(sigma.delta) != p.fp.det) {
      sigma.delta <- rep(sigma.delta, p.fp.det)
    }
    if (length(mu.delta) != p.fp.det) {
      mu.delta <- rep(mu.delta, p.fp.det)
    }
    Sigma.delta <- sigma.delta * diag(p.fp.det)
  } else {
    if (verbose) {
      message("No prior specified for delta.normal.\nSetting prior mean to 0 and prior variance to 10\n")
    }
    mu.delta <- rep(0, p.fp.det)
    sigma.delta <- rep(10, p.fp.det)
    Sigma.delta <- diag(p.fp.det) * 10 
  }
  # sigma.sq.psi --------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.psi.ig) | length(priors$sigma.sq.psi.ig) != 2) {
        stop("sigma.sq.psi.ig must be a list of length 2")
      }
      sigma.sq.psi.a <- priors$sigma.sq.psi.ig[[1]]
      sigma.sq.psi.b <- priors$sigma.sq.psi.ig[[2]]
      if (length(sigma.sq.psi.a) != p.occ.re & length(sigma.sq.psi.a) != 1) {
        if (p.occ.re == 1) {
          stop(paste("sigma.sq.psi.ig[[1]] must be a vector of length ", 
                     p.occ.re, " with elements corresponding to sigma.sq.psis' shape", sep = ""))
        } else {
          stop(paste("sigma.sq.psi.ig[[1]] must be a vector of length ", 
                     p.occ.re, " or 1 with elements corresponding to sigma.sq.psis' shape", sep = ""))
        }
      }
      if (length(sigma.sq.psi.b) != p.occ.re & length(sigma.sq.psi.b) != 1) {
        if (p.occ.re == 1) {
          stop(paste("sigma.sq.psi.ig[[2]] must be a vector of length ", 
                     p.occ.re, " with elements corresponding to sigma.sq.psis' scale", sep = ""))
        } else {
          stop(paste("sigma.sq.psi.ig[[2]] must be a vector of length ", 
                     p.occ.re, " or 1with elements corresponding to sigma.sq.psis' scale", sep = ""))
        }
      }
      if (length(sigma.sq.psi.a) != p.occ.re) {
        sigma.sq.psi.a <- rep(sigma.sq.psi.a, p.occ.re)
      }
      if (length(sigma.sq.psi.b) != p.occ.re) {
        sigma.sq.psi.b <- rep(sigma.sq.psi.b, p.occ.re)
      }
    }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.psi.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.psi.a <- rep(0.1, p.occ.re)
      sigma.sq.psi.b <- rep(0.1, p.occ.re)
    }
  } else {
    sigma.sq.psi.a <- 0
    sigma.sq.psi.b <- 0
  }
  # sigma.sq.mu --------------------
  if (p.tp.det.re > 0) {
    if ("sigma.sq.mu.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.mu.ig) | length(priors$sigma.sq.mu.ig) != 2) {
        stop("sigma.sq.mu.ig must be a list of length 2")
      }
      sigma.sq.mu.a <- priors$sigma.sq.mu.ig[[1]]
      sigma.sq.mu.b <- priors$sigma.sq.mu.ig[[2]]
      if (length(sigma.sq.mu.a) != p.tp.det.re & length(sigma.sq.mu.a) != 1) {
        if (p.tp.det.re == 1) {
          stop(paste("sigma.sq.mu.ig[[1]] must be a vector of length ", 
                     p.tp.det.re, " with elements corresponding to sigma.sq.mu's shape", sep = ""))
        } else {
          stop(paste("sigma.sq.mu.ig[[1]] must be a vector of length ", 
                     p.tp.det.re, " or 1 with elements corresponding to sigma.sq.mu's shape", sep = ""))
        }
      }
      if (length(sigma.sq.mu.b) != p.tp.det.re & length(sigma.sq.mu.b) != 1) {
        if (p.tp.det.re == 1) {
          stop(paste("sigma.sq.mu.ig[[2]] must be a vector of length ", 
                     p.tp.det.re, " with elements corresponding to sigma.sq.mu's scale", sep = ""))
        } else {
          stop(paste("sigma.sq.mu.ig[[2]] must be a vector of length ", 
                     p.tp.det.re, " or 1 with elements corresponding to sigma.sq.mu's scale", sep = ""))
        }
      }
      if (length(sigma.sq.mu.a) != p.tp.det.re) {
        sigma.sq.mu.a <- rep(sigma.sq.mu.a, p.tp.det.re)
      }
      if (length(sigma.sq.mu.b) != p.tp.det.re) {
        sigma.sq.mu.b <- rep(sigma.sq.mu.b, p.tp.det.re)
      }
    }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.mu.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.mu.a <- rep(0.1, p.tp.det.re)
      sigma.sq.mu.b <- rep(0.1, p.tp.det.re)
    }
  } else {
    sigma.sq.mu.a <- 0
    sigma.sq.mu.b <- 0
  }
  # sigma.sq.omega --------------------
  if (p.fp.det.re > 0) {
    if ("sigma.sq.omega.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.omega.ig) | length(priors$sigma.sq.omega.ig) != 2) {
        stop("sigma.sq.omega.ig must be a list of length 2")
      }
      sigma.sq.omega.a <- priors$sigma.sq.omega.ig[[1]]
      sigma.sq.omega.b <- priors$sigma.sq.omega.ig[[2]]
      if (length(sigma.sq.omega.a) != p.fp.det.re & length(sigma.sq.omega.a) != 1) {
        if (p.fp.det.re == 1) {
          stop(paste("sigma.sq.omega.ig[[1]] must be a vector of length ", 
                     p.fp.det.re, " with elements corresponding to sigma.sq.omega's shape", sep = ""))
        } else {
          stop(paste("sigma.sq.omega.ig[[1]] must be a vector of length ", 
                     p.fp.det.re, " or 1 with elements corresponding to sigma.sq.omega's shape", sep = ""))
        }
      }
      if (length(sigma.sq.omega.b) != p.fp.det.re & length(sigma.sq.omega.b) != 1) {
        if (p.fp.det.re == 1) {
          stop(paste("sigma.sq.omega.ig[[2]] must be a vector of length ", 
                     p.fp.det.re, " with elements corresponding to sigma.sq.omega's scale", sep = ""))
        } else {
          stop(paste("sigma.sq.omega.ig[[2]] must be a vector of length ", 
                     p.fp.det.re, " or 1 with elements corresponding to sigma.sq.omega's scale", sep = ""))
        }
      }
      if (length(sigma.sq.omega.a) != p.fp.det.re) {
        sigma.sq.omega.a <- rep(sigma.sq.omega.a, p.fp.det.re)
      }
      if (length(sigma.sq.omega.b) != p.fp.det.re) {
        sigma.sq.omega.b <- rep(sigma.sq.omega.b, p.fp.det.re)
      }
    }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.omega.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.omega.a <- rep(0.1, p.fp.det.re)
      sigma.sq.omega.b <- rep(0.1, p.fp.det.re)
    }
  } else {
    sigma.sq.omega.a <- 0
    sigma.sq.omega.b <- 0
  }
  # Spatial parameters can be fixed at their initial values by setting their
  # prior to 'fixed'. phi and nu are always fixed together. Fixing applies to
  # the parameters of all SVCs.
  fixed.phi <- FALSE
  fixed.sigma.sq <- FALSE
  # Check a list prior with two elements that are each of length p.svc or 1,
  # and expand both to length p.svc.
  check.svc.prior <- function(prior, name, desc.1, desc.2) {
    if (!is.list(prior) | length(prior) != 2) {
      stop(paste(name, " must be a list of length 2", sep = ""))
    }
    a <- prior[[1]]
    b <- prior[[2]]
    if (length(a) != p.svc & length(a) != 1) {
      stop(paste(name, "[[1]] must be a vector of length ", p.svc,
                 " or 1 with elements corresponding to ", desc.1,
                 " for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(b) != p.svc & length(b) != 1) {
      stop(paste(name, "[[2]] must be a vector of length ", p.svc,
                 " or 1 with elements corresponding to ", desc.2,
                 " for each covariate with spatially-varying coefficients", sep = ""))
    }
    list(a = rep(a, length.out = p.svc), b = rep(b, length.out = p.svc))
  }
  # phi -----------------------------
  if ("phi.unif" %in% names(priors)) {
    if (priors$phi.unif[1] == 'fixed') {
      fixed.phi <- TRUE
      phi.a <- rep(1, p.svc)
      phi.b <- rep(1, p.svc)
      if ((cov.model == 'matern') & ("nu.unif" %in% names(priors))) {
        if (priors$nu.unif[1] != 'fixed') {
          message("phi is specified as fixed in priors$phi.unif but nu is not, which is not allowed. Continuing to fit the model with both phi and nu fixed.")
        }
      }
    } else {
      tmp <- check.svc.prior(priors$phi.unif, "phi.unif", "phis' lower bound",
                             "phis' upper bound")
      phi.a <- tmp$a
      phi.b <- tmp$b
    }
  } else {
    if (verbose) {
      message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
    }
    coords.D <- iDist(coords)
    phi.a <- rep(3 / max(coords.D), p.svc)
    phi.b <- rep(3 / sort(unique(c(coords.D)))[2], p.svc)
  }
  # sigma.sq -----------------------------
  # Check if both an ig and uniform prior are specified
  if (("sigma.sq.ig" %in% names(priors)) & ("sigma.sq.unif" %in% names(priors))) {
    stop("cannot specify both an IG and a uniform prior for sigma.sq")
  }
  if ("sigma.sq.ig" %in% names(priors)) { # inverse-gamma prior.
    sigma.sq.ig <- TRUE
    if (priors$sigma.sq.ig[1] == 'fixed') {
      fixed.sigma.sq <- TRUE
      sigma.sq.a <- rep(1, p.svc)
      sigma.sq.b <- rep(1, p.svc)
    } else {
      tmp <- check.svc.prior(priors$sigma.sq.ig, "sigma.sq.ig", "sigma.sqs' shape",
                             "sigma.sqs' scale")
      sigma.sq.a <- tmp$a
      sigma.sq.b <- tmp$b
    }
  } else if ("sigma.sq.unif" %in% names(priors)) { # uniform prior
    if (priors$sigma.sq.unif[1] == 'fixed') {
      sigma.sq.ig <- TRUE # This just makes the C++ side a bit easier.
      fixed.sigma.sq <- TRUE
      sigma.sq.a <- rep(1, p.svc)
      sigma.sq.b <- rep(1, p.svc)
    } else {
      sigma.sq.ig <- FALSE
      tmp <- check.svc.prior(priors$sigma.sq.unif, "sigma.sq.unif",
                             "sigma.sqs' lower bound", "sigma.sqs' upper bound")
      sigma.sq.a <- tmp$a
      sigma.sq.b <- tmp$b
    }
  } else {
    if (verbose) {
      message("No prior specified for sigma.sq.\nUsing an inverse-Gamma prior with the shape parameter set to 2 and scale parameter to 1.\n")
    }
    sigma.sq.ig <- TRUE
    sigma.sq.a <- rep(2, p.svc)
    sigma.sq.b <- rep(1, p.svc)
  }
  # nu -----------------------------
  if (cov.model == 'matern') {
    if (!"nu.unif" %in% names(priors)) {
      stop("nu.unif must be specified in priors value list")
    }
    if (priors$nu.unif[1] == 'fixed') {
      if (!fixed.phi) {
        message("nu is specified as fixed in priors$nu.unif but phi is not, which is not allowed. Continuing to fit the model with both phi and nu fixed.")
      }
      fixed.phi <- TRUE
      nu.a <- rep(1, p.svc)
      nu.b <- rep(1, p.svc)
    } else {
      tmp <- check.svc.prior(priors$nu.unif, "nu.unif", "nus' lower bound",
                             "nus' upper bound")
      nu.a <- tmp$a
      nu.b <- tmp$b
    }
  } else {
    nu.a <- rep(0, p.svc)
    nu.b <- rep(0, p.svc)
  }
  # If phi is fixed through nu.unif = 'fixed', phi's prior bounds are not used.
  if (fixed.phi) {
    phi.a <- rep(1, p.svc)
    phi.b <- rep(1, p.svc)
    if (cov.model == 'matern') {
      nu.a <- rep(1, p.svc)
      nu.b <- rep(1, p.svc)
    }
  }

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # z -------------------------------
  if ("z" %in% names(inits)) {
    z.inits <- inits$z
    if (!is.vector(z.inits)) {
      stop(paste("initial values for z must be a vector of length ",
                 J, sep = ""))
    }
    if (length(z.inits) != J) {
      stop(paste("initial values for z must be a vector of length ",
                 J, sep = ""))
    }
    if (!all(z.inits %in% c(0, 1))) {
      stop("initial values for z must be 0 or 1")
    }
    # Reorder the user supplied inits values
    z.inits <- z.inits[long.ord]
    # A validated true positive (d > 0) means the site is occupied.
    z.test <- apply(d.big, 1, function(a) sum(a, na.rm = TRUE) > 0)
    if (sum(z.inits < z.test) > 0) {
      stop("initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 at sites with at least one validated true positive (d > 0).")
    }
  } else {
    # Default is to set initial values based on wherever the species was observed, even though 
    # some of those observations may be false positives. 
    z.inits <- as.numeric(apply(y.big, 1, max, na.rm = TRUE) > 0)
    if (verbose) {
      message("z is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
  # beta -----------------------
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (length(beta.inits) != p.occ & length(beta.inits) != 1) {
      if (p.occ == 1) {
        stop(paste("initial values for beta must be of length ", p.occ,
                   sep = ""))
        
      } else {
        stop(paste("initial values for beta must be of length ", p.occ, " or 1",
                   sep = ""))
      }
    }
    if (length(beta.inits) != p.occ) {
      beta.inits <- rep(beta.inits, p.occ)
    }
  } else {
    beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
    if (verbose) {
      message('beta is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
    }
  }
  # gamma -----------------------
  if ("gamma" %in% names(inits)) {
    gamma.inits <- inits[["gamma"]]
    if (length(gamma.inits) != p.tp.det & length(gamma.inits) != 1) {
      if (p.tp.det == 1) {
        stop(paste("initial values for gamma must be of length ", p.tp.det,
                   sep = ""))
      } else {
        stop(paste("initial values for gamma must be of length ", p.tp.det, " or 1",
                   sep = ""))
      }
    }
    if (length(gamma.inits) != p.tp.det) {
      gamma.inits <- rep(gamma.inits, p.tp.det)
    }
  } else {
    gamma.inits <- rnorm(p.tp.det, mu.gamma, sqrt(sigma.gamma))
    if (verbose) {
      message("gamma is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  # delta -----------------------
  if ("delta" %in% names(inits)) {
    delta.inits <- inits[["delta"]]
    if (length(delta.inits) != p.fp.det & length(delta.inits) != 1) {
      if (p.fp.det == 1) {
        stop(paste("initial values for delta must be of length ", p.fp.det,
                   sep = ""))
      } else {
        stop(paste("initial values for delta must be of length ", p.fp.det, " or 1",
                   sep = ""))
      }
    }
    if (length(delta.inits) != p.fp.det) {
      delta.inits <- rep(delta.inits, p.fp.det)
    }
  } else {
    delta.inits <- rnorm(p.fp.det, mu.delta, sqrt(sigma.delta))
    if (verbose) {
      message("delta is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  
  # sigma.sq.psi -------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi" %in% names(inits)) {
      sigma.sq.psi.inits <- inits[["sigma.sq.psi"]]
      if (length(sigma.sq.psi.inits) != p.occ.re & length(sigma.sq.psi.inits) != 1) {
        if (p.occ.re == 1) {
          stop(paste("initial values for sigma.sq.psi must be of length ", p.occ.re, 
                     sep = ""))
        } else {
          stop(paste("initial values for sigma.sq.psi must be of length ", p.occ.re, 
                     " or 1", sep = ""))
        }
      }
      if (length(sigma.sq.psi.inits) != p.occ.re) {
        sigma.sq.psi.inits <- rep(sigma.sq.psi.inits, p.occ.re)
      }
    } else {
      sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.psi is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    beta.star.indx <- rep(0:(p.occ.re - 1), n.occ.re.long)
    beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
  } else {
    sigma.sq.psi.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
  }
  # sigma.sq.mu ------------------
  if (p.tp.det.re > 0) {
    if ("sigma.sq.mu" %in% names(inits)) {
      sigma.sq.mu.inits <- inits[["sigma.sq.mu"]]
      if (length(sigma.sq.mu.inits) != p.tp.det.re & length(sigma.sq.mu.inits) != 1) {
        if (p.tp.det.re == 1) {
          stop(paste("initial values for sigma.sq.mu must be of length ", p.tp.det.re, 
                     sep = ""))
        } else {
          stop(paste("initial values for sigma.sq.mu must be of length ", p.tp.det.re, 
                     " or 1", sep = ""))
          
        }
      }
      if (length(sigma.sq.mu.inits) != p.tp.det.re) {
        sigma.sq.mu.inits <- rep(sigma.sq.mu.inits, p.tp.det.re)
      }
    } else {
      sigma.sq.mu.inits <- runif(p.tp.det.re, 0.5, 4)
      if (verbose) {
        message("sigma.sq.mu is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    gamma.star.indx <- rep(0:(p.tp.det.re - 1), n.tp.det.re.long)
    gamma.star.inits <- rnorm(n.tp.det.re, 0, sqrt(sigma.sq.mu.inits[gamma.star.indx + 1]))
  } else {
    sigma.sq.mu.inits <- 0
    gamma.star.indx <- 0
    gamma.star.inits <- 0
  }
  # sigma.sq.omega ------------------
  if (p.fp.det.re > 0) {
    if ("sigma.sq.omega" %in% names(inits)) {
      sigma.sq.omega.inits <- inits[["sigma.sq.omega"]]
      if (length(sigma.sq.omega.inits) != p.fp.det.re & length(sigma.sq.omega.inits) != 1) {
        if (p.fp.det.re == 1) {
          stop(paste("initial values for sigma.sq.omega must be of length ", p.fp.det.re, 
                     sep = ""))
        } else {
          stop(paste("initial values for sigma.sq.omega must be of length ", p.fp.det.re, 
                     " or 1", sep = ""))
          
        }
      }
      if (length(sigma.sq.omega.inits) != p.fp.det.re) {
        sigma.sq.omega.inits <- rep(sigma.sq.omega.inits, p.fp.det.re)
      }
    } else {
      sigma.sq.omega.inits <- runif(p.fp.det.re, 0.5, 4)
      if (verbose) {
        message("sigma.sq.omega is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    delta.star.indx <- rep(0:(p.fp.det.re - 1), n.fp.det.re.long)
    delta.star.inits <- rnorm(n.fp.det.re, 0, sqrt(sigma.sq.omega.inits[delta.star.indx + 1]))
  } else {
    sigma.sq.omega.inits <- 0
    delta.star.indx <- 0
    delta.star.inits <- 0
  }
  
  # phi -----------------------------
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != p.svc & length(phi.inits) != 1) {
      stop(paste("initial values for phi must be of length ", p.svc, " or 1",
                 sep = ""))
    }
    if (length(phi.inits) != p.svc) {
      phi.inits <- rep(phi.inits, p.svc)
    }
  } else {
    if (fixed.phi) {
      stop("phi is specified as fixed in priors, so its fixed value must be supplied in inits$phi")
    }
    phi.inits <- runif(p.svc, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  # sigma.sq ------------------------
  if ("sigma.sq" %in% names(inits)) {
    sigma.sq.inits <- inits[["sigma.sq"]]
    if (length(sigma.sq.inits) != p.svc & length(sigma.sq.inits) != 1) {
      stop(paste("initial values for sigma.sq must be of length ", p.svc, " or 1",
                 sep = ""))
    }
    if (length(sigma.sq.inits) != p.svc) {
      sigma.sq.inits <- rep(sigma.sq.inits, p.svc)
    }
  } else {
    if (fixed.sigma.sq) {
      stop("sigma.sq is specified as fixed in priors, so its fixed value must be supplied in inits$sigma.sq")
    }
    if (sigma.sq.ig) {
      sigma.sq.inits <- rigamma(p.svc, sigma.sq.a, sigma.sq.b)
    } else {
      sigma.sq.inits <- runif(p.svc, sigma.sq.a, sigma.sq.b)
    }
    if (verbose) {
      message("sigma.sq is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  # w -----------------------------
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (!is.matrix(w.inits)) {
      stop(paste("initial values for w must be a matrix with dimensions ",
                 p.svc, " x ", J.w, sep = ""))
    }
    if (nrow(w.inits) != p.svc | ncol(w.inits) != J.w) {
      stop(paste("initial values for w must be a matrix with dimensions ",
                 p.svc, " x ", J.w, sep = ""))
    }
    # Reorder the user supplied inits values
    w.inits <- w.inits[, ord, drop = FALSE]
  } else {
    w.inits <- matrix(0, p.svc, J.w)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
    }
  }
  # nu ------------------------
  if (cov.model == 'matern') {
    if ("nu" %in% names(inits)) {
      nu.inits <- inits[["nu"]]
      if (length(nu.inits) != p.svc & length(nu.inits) != 1) {
        stop(paste("initial values for nu must be of length ", p.svc, " or 1",
                   sep = ""))
      }
      if (length(nu.inits) != p.svc) {
        nu.inits <- rep(nu.inits, p.svc)
      }
    } else {
      if (fixed.phi) {
        stop("nu is specified as fixed in priors, so its fixed value must be supplied in inits$nu")
      }
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      nu.inits <- runif(p.svc, nu.a, nu.b)
    }
  } else {
    nu.inits <- rep(0, p.svc)
  }

  # Should initial values be fixed --
  if ("fix" %in% names(inits)) {
    fix.inits <- inits[["fix"]]
    if ((fix.inits != TRUE) & (fix.inits != FALSE)) {
      stop(paste("inits$fix must take value TRUE or FALSE"))
    }
  } else {
    fix.inits <- FALSE
  }
  if (verbose & fix.inits & (n.chains > 1)) {
    message("Fixing initial values across all chains\n")
  }
 
  # Covariance Model ----------------------------------------------------
  # Order must match util.cpp spCor.
  cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
  if(! cov.model %in% cov.model.names){
    stop("specified cov.model '",cov.model,"' is not a valid option; choose from ",
         paste(cov.model.names, collapse=", ", sep="") ,".")}
  # Obo for cov model lookup on c side
  cov.model.indx <- which(cov.model == cov.model.names) - 1

  # Prep for SVCs ---------------------------------------------------------
  X.w <- X[, svc.cols, drop = FALSE]

  # Tuning values ---------------------------------------------------------
  # Spatial tuning values are always passed to C as blocks of length p.svc in
  # the order sigma.sq, phi, nu. Slots for parameters that are not updated
  # with Metropolis are unused.
  sigma.sq.tuning <- rep(1, p.svc)
  phi.tuning <- rep(1, p.svc)
  nu.tuning <- rep(1, p.svc)
  # Check a spatial tuning value is of length 1 or p.svc and expand it.
  check.svc.tuning <- function(val, name) {
    if (length(val) != 1 & length(val) != p.svc) {
      stop(paste(name, " tuning must be either a single value or a vector of length ",
                 p.svc, sep = ""))
    }
    rep(val, length.out = p.svc)
  }
  if (!missing(tuning)) {
    names(tuning) <- tolower(names(tuning))
    # phi ---------------------------
    if (!fixed.phi) {
      if(!"phi" %in% names(tuning)) {
        stop("phi must be specified in tuning value list")
      }
      phi.tuning <- check.svc.tuning(tuning$phi, "phi")
      if (cov.model == 'matern') {
        # nu --------------------------
        if(!"nu" %in% names(tuning)) {
          stop("nu must be specified in tuning value list")
        }
        nu.tuning <- check.svc.tuning(tuning$nu, "nu")
      }
    }
    if (!sigma.sq.ig) {
      # sigma.sq --------------------------
      if(!"sigma.sq" %in% names(tuning)) {
        stop("sigma.sq must be specified in tuning value list")
      }
      sigma.sq.tuning <- check.svc.tuning(tuning$sigma.sq, "sigma.sq")
    }
  }
  if (missing(tuning)) {
    gamma.tuning <- rep(1, p.tp.det)
    delta.tuning <- rep(1, p.fp.det)
    gamma.star.tuning <- rep(1, n.tp.det.re)
    delta.star.tuning <- rep(1, n.fp.det.re)
  } else {
    names(tuning) <- tolower(names(tuning))
    # gamma ---------------------------
    if (!"gamma" %in% names(tuning)) {
      stop("gamma must be specified in tuning value list")
    }
    gamma.tuning <- tuning$gamma
    if (length(gamma.tuning) != 1 & length(gamma.tuning) != p.tp.det) {
      stop(paste("gamma tuning must be a single value or a vector of length ",
                 p.tp.det, sep = ''))
    }
    if (length(gamma.tuning) == 1) {
      gamma.tuning <- rep(gamma.tuning, p.tp.det)
    }
    # delta ---------------------------
    if (!"delta" %in% names(tuning)) {
      stop("delta must be specified in tuning value list")
    }
    delta.tuning <- tuning$delta
    if (length(delta.tuning) != 1 & length(delta.tuning) != p.fp.det) {
      stop(paste("delta tuning must be a single value or a vector of length ",
                 p.fp.det, sep = ''))
    }
    if (length(delta.tuning) == 1) {
      delta.tuning <- rep(delta.tuning, p.fp.det)
    }
    if (p.tp.det.re > 0) {
      # gamma.star ---------------------------
      if(!"gamma.star" %in% names(tuning)) {
        stop("gamma.star must be specified in tuning value list")
      }
      gamma.star.tuning <- tuning$gamma.star
      if (length(gamma.star.tuning) != 1) {
        stop("gamma.star tuning must be a single value")
      }
      gamma.star.tuning <- rep(gamma.star.tuning, n.tp.det.re)
    } else {
      gamma.star.tuning <- NULL
    }
    if (p.fp.det.re > 0) {
      # delta.star ---------------------------
      if(!"delta.star" %in% names(tuning)) {
        stop("delta.star must be specified in tuning value list")
      }
      delta.star.tuning <- tuning$delta.star
      if (length(delta.star.tuning) != 1) {
        stop("delta.star tuning must be a single value")
      }
      delta.star.tuning <- rep(delta.star.tuning, n.fp.det.re)
    } else {
      delta.star.tuning <- NULL
    }
  }
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning, gamma.tuning,
                    delta.tuning, gamma.star.tuning, delta.star.tuning))

  # Nearest Neighbor Search ---------------------------------------------
  if(verbose){
    cat("----------------------------------------\n");
    cat("\tBuilding the neighbor list\n");
    cat("----------------------------------------\n");
  }

  search.type.names <- c("brute", "cb")

  if(!search.type %in% search.type.names){
    stop("specified search.type '",search.type,
         "' is not a valid option; choose from ",
         paste(search.type.names, collapse=", ", sep="") ,".")
  }

  storage.mode(n.neighbors) <- "integer"
  storage.mode(n.omp.threads) <- "integer"
  ## Indexes
  if(search.type == "brute"){
    indx <- mkNNIndx(coords, n.neighbors, n.omp.threads)
  } else{
    indx <- mkNNIndxCB(coords, n.neighbors, n.omp.threads)
  }

  nn.indx <- indx$nnIndx
  nn.indx.lu <- indx$nnIndxLU
  nn.indx.run.time <- indx$run.time

  storage.mode(nn.indx) <- "integer"
  storage.mode(nn.indx.lu) <- "integer"
  storage.mode(u.search.type) <- "integer"
  storage.mode(J) <- "integer"
  storage.mode(J.w) <- "integer"

  if(verbose){
    cat("----------------------------------------\n");
    cat("Building the neighbors of neighbors list\n");
    cat("----------------------------------------\n");
  }

  indx <- mkUIndx(J.w, n.neighbors, nn.indx, nn.indx.lu, u.search.type)

  u.indx <- indx$u.indx
  u.indx.lu <- indx$u.indx.lu
  ui.indx <- indx$ui.indx
  u.indx.run.time <- indx$run.time

  curr.chain <- 1
  # Set storage for all variables ---------------------------------------
  storage.mode(y) <- "double"
  storage.mode(z.inits) <- "double"
  storage.mode(n.val) <- "double"
  storage.mode(d) <- "double"
  storage.mode(X.mu) <- "double"
  storage.mode(X.omega) <- "double"
  storage.mode(X) <- "double"
  storage.mode(X.w) <- "double"
  # Whether the first column of X.mu/X.omega is an intercept. Required for the
  # joint intercept/random effect shift in the sampler.
  tp.det.int <- as.integer(x.mu.names[1] == "(Intercept)")
  fp.det.int <- as.integer(x.omega.names[1] == "(Intercept)")
  consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.tp.det, p.tp.det.re,
              n.tp.det.re, p.fp.det, p.fp.det.re, n.fp.det.re, tp.det.int,
              fp.det.int, J.w, p.svc)
  storage.mode(consts) <- "integer"
  # Each spatial prior hyperparameter is a block of length p.svc.
  sp.priors <- c(phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b)
  storage.mode(sp.priors) <- "double"
  sp.info <- c(cov.model.indx, sigma.sq.ig, fixed.phi, fixed.sigma.sq)
  storage.mode(sp.info) <- "integer"
  storage.mode(grid.index.c) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(w.inits) <- "double"
  storage.mode(u.indx) <- "integer"
  storage.mode(u.indx.lu) <- "integer"
  storage.mode(ui.indx) <- "integer"
  storage.mode(K) <- "double"
  storage.mode(beta.inits) <- "double"
  storage.mode(gamma.inits) <- "double"
  storage.mode(delta.inits) <- "double"
  storage.mode(z.long.indx) <- "integer"
  storage.mode(mu.beta) <- "double"
  storage.mode(Sigma.beta) <- "double"
  storage.mode(mu.gamma) <- "double"
  storage.mode(Sigma.gamma) <- "double"
  storage.mode(mu.delta) <- "double"
  storage.mode(Sigma.delta) <- "double"
  storage.mode(tuning.c) <- "double"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(verbose) <- "integer"
  storage.mode(n.report) <- "integer"
  chain.info <- c(curr.chain, n.chains)
  storage.mode(chain.info) <- "integer"
  n.post.samples <- length(seq(from = n.burn + 1, 
                               to = n.samples, 
                               by = as.integer(n.thin)))
  storage.mode(n.post.samples) <- "integer"
  samples.info <- c(n.burn, n.thin, n.post.samples)
  storage.mode(samples.info) <- "integer"
  # For true positive detection random effects
  storage.mode(X.mu.re) <- "integer"
  gamma.level.indx <- sort(unique(c(X.mu.re)))
  storage.mode(gamma.level.indx) <- "integer"
  storage.mode(n.tp.det.re.long) <- "integer"
  storage.mode(sigma.sq.mu.inits) <- "double"
  storage.mode(sigma.sq.mu.a) <- "double"
  storage.mode(sigma.sq.mu.b) <- "double"
  storage.mode(gamma.star.inits) <- "double"
  storage.mode(gamma.star.indx) <- "integer"
  # For false positive detection random effects
  storage.mode(X.omega.re) <- "integer"
  delta.level.indx <- sort(unique(c(X.omega.re)))
  storage.mode(delta.level.indx) <- "integer"
  storage.mode(n.fp.det.re.long) <- "integer"
  storage.mode(sigma.sq.omega.inits) <- "double"
  storage.mode(sigma.sq.omega.a) <- "double"
  storage.mode(sigma.sq.omega.b) <- "double"
  storage.mode(delta.star.inits) <- "double"
  storage.mode(delta.star.indx) <- "integer"
  # For occurrence random effects
  storage.mode(X.re) <- "integer"
  beta.level.indx <- sort(unique(c(X.re)))
  storage.mode(beta.level.indx) <- "integer"
  storage.mode(sigma.sq.psi.inits) <- "double"
  storage.mode(sigma.sq.psi.a) <- "double"
  storage.mode(sigma.sq.psi.b) <- "double"
  storage.mode(beta.star.inits) <- "double"
  storage.mode(beta.star.indx) <- "integer"
  
  # Fit the model -------------------------------------------------------
  out.tmp <- list()
  out <- list()
  for (i in 1:n.chains) {
    # Change initial values if i > 1
    if ((i > 1) & (!fix.inits)) {
      beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
      gamma.inits <- rnorm(p.tp.det, mu.gamma, sqrt(sigma.gamma))
      delta.inits <- rnorm(p.fp.det, mu.delta, sqrt(sigma.delta))
      if (p.occ.re > 0) {
        sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
        beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
      }
      if (p.tp.det.re > 0) {
        sigma.sq.mu.inits <- runif(p.tp.det.re, 0.5, 3)
        gamma.star.inits <- rnorm(n.tp.det.re, 0, sqrt(sigma.sq.mu.inits[gamma.star.indx + 1]))
      }
      if (p.fp.det.re > 0) {
        sigma.sq.omega.inits <- runif(p.fp.det.re, 0.5, 3)
        delta.star.inits <- rnorm(n.fp.det.re, 0, sqrt(sigma.sq.omega.inits[delta.star.indx + 1]))
      }
      if (!fixed.sigma.sq) {
        if (sigma.sq.ig) {
          sigma.sq.inits <- rigamma(p.svc, sigma.sq.a, sigma.sq.b)
        } else {
          sigma.sq.inits <- runif(p.svc, sigma.sq.a, sigma.sq.b)
        }
      }
      if (!fixed.phi) {
        phi.inits <- runif(p.svc, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(p.svc, nu.a, nu.b)
        }
      }
    }
    # Spatial starting values as blocks of length p.svc in the order sigma.sq, phi, nu
    theta.inits <- c(sigma.sq.inits, phi.inits, nu.inits)
    storage.mode(theta.inits) <- "double"
    storage.mode(chain.info) <- "integer"
    # Run the model in C
    out.tmp[[i]] <- .Call("svcFpPGOccNNGP", y, X, X.w, X.mu, X.omega, coords, X.re,
                          X.mu.re, X.omega.re, d, n.val, consts, n.occ.re.long,
                          n.tp.det.re.long, n.fp.det.re.long, n.neighbors, nn.indx,
                          nn.indx.lu, u.indx, u.indx.lu, ui.indx, beta.inits,
                          gamma.inits, delta.inits, sigma.sq.psi.inits,
                          sigma.sq.mu.inits, sigma.sq.omega.inits, beta.star.inits,
                          gamma.star.inits, delta.star.inits, z.inits, w.inits,
                          theta.inits, z.long.indx, beta.star.indx, beta.level.indx,
                          gamma.star.indx, gamma.level.indx, delta.star.indx,
                          delta.level.indx, mu.beta, mu.gamma, mu.delta,
                          Sigma.beta, Sigma.gamma, Sigma.delta, sp.priors,
                          sigma.sq.psi.a, sigma.sq.psi.b, sigma.sq.mu.a, sigma.sq.mu.b,
                          sigma.sq.omega.a, sigma.sq.omega.b, tuning.c, sp.info,
                          n.batch, batch.length, accept.rate, n.omp.threads, verbose,
                          n.report, samples.info, chain.info, grid.index.c)
    chain.info[1] <- chain.info[1] + 1
  } # i   
  # Calculate R-Hat ---------------
  out <- list()
  out$rhat <- list()
  if (n.chains > 1) {
    # as.vector removes the "Upper CI" when there is only 1 variable. 
    out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      mcmc(t(a$beta.samples)))), 
      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$gamma <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      mcmc(t(a$gamma.samples)))), 
      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$delta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      mcmc(t(a$delta.samples)))),
      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      mcmc(t(a$theta.samples)))),
      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
    if (p.tp.det.re > 0) {
      out$rhat$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        mcmc(t(a$sigma.sq.mu.samples)))), 
        autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
    if (p.fp.det.re > 0) {
      out$rhat$sigma.sq.omega <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        mcmc(t(a$sigma.sq.omega.samples)))), 
        autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
    if (p.occ.re > 0) {
      out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        mcmc(t(a$sigma.sq.psi.samples)))), 
        autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
  } else {
    out$rhat$beta <- rep(NA, p.occ)
    out$rhat$gamma <- rep(NA, p.tp.det)
    out$rhat$delta <- rep(NA, p.fp.det)
    out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3 * p.svc, 2 * p.svc))
    if (p.tp.det.re > 0) {
      out$rhat$sigma.sq.mu <- rep(NA, p.tp.det.re)
    }
    if (p.fp.det.re > 0) {
      out$rhat$sigma.sq.omega <- rep(NA, p.fp.det.re)
    }
    if (p.occ.re > 0) {
      out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
    }
  }
  # Put everything into MCMC objects
  out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
  colnames(out$beta.samples) <- x.names
  out$gamma.samples <- mcmc(do.call(rbind, 
                                    lapply(out.tmp, function(a) t(a$gamma.samples))))
  colnames(out$gamma.samples) <- x.mu.names
  out$delta.samples <- mcmc(do.call(rbind, 
                                    lapply(out.tmp, function(a) t(a$delta.samples))))
  colnames(out$delta.samples) <- x.omega.names
  out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
  if (cov.model != 'matern') {
    theta.names <- paste(rep(c('sigma.sq', 'phi'), each = p.svc), x.names[svc.cols], sep = '-')
  } else {
    theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = p.svc), x.names[svc.cols], sep = '-')
  }
  colnames(out$theta.samples) <- theta.names
  # Get everything back in the original order
  out$coords <- coords[order(ord), , drop = FALSE]
  out$z.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$z.samples))))
  out$z.samples <- mcmc(out$z.samples[, order(long.ord), drop = FALSE])
  # w is stored in C as (SVC, location) for each sample. Return an
  # n.post x p.svc x J.w array with locations in the original order.
  out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples,
                                         dim = c(p.svc, J.w, n.post.samples))))
  out$w.samples <- array(out$w.samples, dim = c(p.svc, J.w, n.post.samples * n.chains))
  out$w.samples <- out$w.samples[, order(ord), , drop = FALSE]
  out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
  out$psi.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$psi.samples))))
  out$psi.samples <- mcmc(out$psi.samples[, order(long.ord), drop = FALSE])
  out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
  out$like.samples <- mcmc(out$like.samples[, order(long.ord), drop = FALSE])
  if (p.occ.re > 0) {
    out$sigma.sq.psi.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.psi.samples))))
    colnames(out$sigma.sq.psi.samples) <- x.re.names
    out$beta.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
    tmp.names <- unlist(re.level.names)
    beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
    colnames(out$beta.star.samples) <- beta.star.names
    out$re.level.names <- re.level.names
  }
  if (p.tp.det.re > 0) {
    out$sigma.sq.mu.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.mu.samples))))
    colnames(out$sigma.sq.mu.samples) <- x.mu.re.names
    out$gamma.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$gamma.star.samples))))
    tmp.names <- unlist(mu.re.level.names)
    gamma.star.names <- paste(rep(x.mu.re.names, n.tp.det.re.long), tmp.names, sep = '-')
    colnames(out$gamma.star.samples) <- gamma.star.names
    out$mu.re.level.names <- mu.re.level.names
  }
  if (p.fp.det.re > 0) {
    out$sigma.sq.omega.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.omega.samples))))
    colnames(out$sigma.sq.omega.samples) <- x.omega.re.names
    out$delta.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$delta.star.samples))))
    tmp.names <- unlist(omega.re.level.names)
    delta.star.names <- paste(rep(x.omega.re.names, n.fp.det.re.long), tmp.names, sep = '-')
    colnames(out$delta.star.samples) <- delta.star.names
    out$omega.re.level.names <- omega.re.level.names
  }
  # Calculate effective sample sizes
  out$ESS <- list()
  out$ESS$beta <- effectiveSize(out$beta.samples)
  out$ESS$gamma <- effectiveSize(out$gamma.samples)
  out$ESS$delta <- effectiveSize(out$delta.samples)
  out$ESS$theta <- effectiveSize(out$theta.samples)
  if (p.tp.det.re > 0) {
    out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
  }
  if (p.fp.det.re > 0) {
    out$ESS$sigma.sq.omega <- effectiveSize(out$sigma.sq.omega.samples)
  }
  if (p.occ.re > 0) {
    out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
  }
  out$X <- X[order(long.ord), , drop = FALSE]
  out$X.re <- X.re[order(long.ord), , drop = FALSE]
  out$X.w <- X.w[order(long.ord), , drop = FALSE]
  # Observation-level design matrices hold the non-missing site/visit
  # combinations, ordered by site within visit. Expand them to J x K x p
  # arrays, put the sites back in the original order, and drop the missing
  # site/visit combinations again.
  K.full <- ncol(y.big)
  obs.keep <- !is.na(c(y.big[order(long.ord), , drop = FALSE]))
  reorder.obs <- function(a) {
    if (ncol(a) == 0) {
      return(a)
    }
    tmp <- matrix(NA, J * K.full, ncol(a))
    tmp[names.long, ] <- a
    tmp <- array(tmp, dim = c(J, K.full, ncol(a)))
    tmp <- tmp[order(long.ord), , , drop = FALSE]
    tmp <- matrix(tmp, J * K.full, ncol(a))
    tmp <- tmp[obs.keep, , drop = FALSE]
    colnames(tmp) <- colnames(a)
    tmp
  }
  out$X.mu <- reorder.obs(X.mu)
  out$X.omega <- reorder.obs(X.omega)
  out$X.mu.re <- reorder.obs(X.mu.re)
  out$X.omega.re <- reorder.obs(X.omega.re)
  out$y <- y.big[order(long.ord), , drop = FALSE]
  out$n.val <- n.val.big[order(long.ord), , drop = FALSE]
  out$d <- d.big[order(long.ord), , drop = FALSE]
  out$n.samples <- batch.length * n.batch
  out$n.neighbors <- n.neighbors
  out$cov.model.indx <- cov.model.indx
  out$svc.cols <- svc.cols
  out$type <- "NNGP"
  out$call <- cl
  out$n.post <- n.post.samples
  out$n.thin <- n.thin
  out$n.burn <- n.burn
  out$n.chains <- n.chains
  if (p.tp.det.re > 0) {
    out$muRE <- TRUE
  } else {
    out$muRE <- FALSE
  }
  if (p.fp.det.re > 0) {
    out$omegaRE <- TRUE
  } else {
    out$omegaRE <- FALSE
  }
  if (p.occ.re > 0) {
    out$psiRE <- TRUE
  } else {
    out$psiRE <- FALSE
  }

  class(out) <- "svcFpPGOcc"
  out$run.time <- proc.time() - ptm
  out
}
