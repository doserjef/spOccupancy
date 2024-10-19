svcTIntPGOcc <- function(occ.formula, det.formula, data, inits, priors, 
                         tuning, svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
                         n.neighbors = 15, search.type = 'cb', n.batch, 
                         batch.length, accept.rate = 0.43, n.omp.threads = 1, 
                         verbose = TRUE, ar1 = FALSE, n.report = 100, 
                         n.burn = round(.10 * n.batch * batch.length), 
                         n.thin = 1, n.chains = 1, ...){

  ptm <- proc.time()

  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing the data\n");
    cat("----------------------------------------\n");
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
    stop("error: data must be specified")
  }
  if (!is.list(data)) {
    stop("error: data must be a list")
  }
  names(data) <- tolower(names(data))
  if (missing(occ.formula)) {
    stop("error: occ.formula must be specified")
  }
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  if (!is.list(data$y)) {
    stop("error: y must be a list of detection-nondetection data sets")
  }
  y <- data$y
  n.data <- length(y)
  if (!is.list(det.formula)) {
    stop(paste("error: det.formula must be a list of ", n.data, " formulas", sep = ''))
  }
  for (q in 1:n.data) {
    if (length(dim(y[[q]])) != 3) {
      stop('Each individual data source in data$y must be specified as a three-dimensional array with dimensions corresponding to site, seasons, and replicate within season. Note that even if a data source is sampled only for one season or only one visit within a season, it still must be specified as a three-dimensional array')
    }
  }
  if (!'sites' %in% names(data)) {
    stop("site ids (sites) must be specified in data")
  }
  sites <- data$sites
  # Number of sites with at least one data source
  J <- length(unique(unlist(sites)))
  # Number of sites for each data set
  J.long <- sapply(y, function(a) dim(a)[[1]])
  if (!'seasons' %in% names(data)) {
    stop("seasons must be specified in data")
  }
  seasons <- data$seasons
  # Number of seasons with at least one data source
  n.years.total <- length(unique(unlist(seasons)))
  # Number of seasons for each data set
  n.years.by.data <- sapply(seasons, length)
  
  # Check occupancy covariates
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("Occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- matrix(1, J, n.years.total)
    } else {
      stop("occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!is.list(data$occ.covs)) {
    stop("error: occ.covs must be a list of matrices, data frames, and/or vectors")
  }
  # Check detection covariates
  if (!'det.covs' %in% names(data)) {
    data$det.covs <- list()
    for (i in 1:n.data) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model for each data source.\n")
      }
      det.formula.curr <- det.formula[[i]]
      if (det.formula.curr == ~ 1) {
        for (i in 1:n.data) {
          data$det.covs[[i]] <- list(int = array(1, dim = dim(y[[1]])))
        }
      } else {
        stop("error: det.covs must be specified in data for a detection model with covariates")
      }
    }
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("error: coords must be a matrix or data frame")
  }
  coords <- as.matrix(data$coords)
  # Check if all spatial coordinates are unique. 
  unique.coords <- unique(data$coords)
  if (nrow(unique.coords) < nrow(data$coords)) {
    stop("coordinates provided in coords are not all unique. spOccupancy requires each site to have its own unique pair of spatial coordinates. This may be the result of an error in preparing the data list, or you will need to change what you consider a 'site' in order to meet this requirement.") 
  }

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering. 
    coords <- coords[ord, ]
    # Occupancy covariates
    for (i in 1:length(data$occ.covs)) {
      if (!is.null(dim(data$occ.covs[[i]]))) { # Time/space varying
        data$occ.covs[[i]] <- data$occ.covs[[i]][ord, , drop = FALSE]
      } else { # Space-varying
        data$occ.covs[[i]] <- data$occ.covs[[i]][ord]
      }
    } 
    # Note that you don't need to actually reorder y or det.covs, you can just
    # reorder the site indices for integrated models.
    sites.orig <- sites
    for (i in 1:n.data) {
      for (j in 1:length(sites[[i]])) {
        sites[[i]][j] <- which(ord == sites[[i]][j])
      }
    }
  }
  # Reformat covariates ---------------------------------------------------
  # Get detection covariates in proper format for each data set
  y.big <- vector(mode = 'list', length = n.data)
  for (i in 1:n.data) {
    # First subset detection covariates to only use those that are included in the analysis.
    data$det.covs[[i]] <- data$det.covs[[i]][names(data$det.covs[[i]]) %in% all.vars(det.formula[[i]])]
    # Null model support
    if (length(data$det.covs[[i]]) == 0) {
      data$det.covs[[i]] <- list(int = array(1, dim = dim(y[[i]])))
    }
    # Turn the covariates into a data frame. Unlist is necessary for when factors
    # are supplied. 
    # Ordered by visit, year, then site. 
    data$det.covs[[i]] <- data.frame(lapply(data$det.covs[[i]], function(a) unlist(c(a))))
    # Get detection covariates in site x year x replicate format
    if (nrow(data$det.covs[[i]]) == dim(y[[i]])[1]) { # if only site-level covariates. 
      data$det.covs[[i]] <- as.data.frame(mapply(rep, data$det.covs[[i]], dim(y[[i]])[2] * dim(y[[i]])[3]))
    } else if (nrow(data$det.covs[[i]]) == dim(y[[i]])[1] * dim(y[[i]])[2]) { # if only site/year level covariates
      data$det.covs[[i]] <- as.data.frame(mapply(rep, data$det.covs[[i]], dim(y[[i]])[3]))
    }
    y.big[[i]] <- y[[i]]
  }
  # Get occurrence covariates in proper format
  # Subset covariates to only use those that are included in the analysis
  data$occ.covs <- data$occ.covs[names(data$occ.covs) %in% all.vars(occ.formula)]
  # Null model support
  if (length(data$occ.covs) == 0) {
    data$occ.covs <- list(int = matrix(1, nrow = J, ncol = n.years.total))
  }
  # Ordered by year, then site within year. 
  data$occ.covs <- data.frame(lapply(data$occ.covs, function(a) unlist(c(a))))
  # Check if only site-level covariates are included
  if (nrow(data$occ.covs) == J) {
    data$occ.covs <- as.data.frame(mapply(rep, data$occ.covs, n.years.total))
  }

  # Checking missing values ---------------------------------------------
  # y and det.covs --------------------
  for (q in 1:n.data) {
    y.na.test <- apply(y.big[[q]], 1, function(a) sum(!is.na(a)))
    if (sum(y.na.test == 0) > 0) {
      stop(paste0("some sites in data source ", q, " in y have all missing detection histories.\nRemove these sites from y and all objects in the 'data' argument if the site is not surveyed by another data source\n, then use 'predict' to obtain predictions at these locations if desired."))
    }
    for (i in 1:ncol(data$det.covs[[q]])) {
      if (sum(is.na(data$det.covs[[q]][, i])) > sum(is.na(y.big[[q]]))) {
        stop(paste0("some elements in det.covs for data set ", q, " have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.")) 
      }
    }
  }
  # occ.covs ------------------------
  if (sum(is.na(data$occ.covs)) != 0) {
    stop("error: missing values in occ.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
  }
  # Check for misalignment in y and det.covs
  for (q in 1:n.data) {
    if (det.formula[[q]] != ~ 1) {
      # Misalignment between y and det.covs
      y.missing <- which(is.na(y[[q]]))
      det.covs.missing <- lapply(data$det.covs[[q]], function(a) which(is.na(a)))
      for (i in 1:length(det.covs.missing)) {
        tmp.indx <- !(y.missing %in% det.covs.missing[[i]])
        if (sum(tmp.indx) > 0) {
          if (i == 1 & verbose) {
            message(paste0("There are missing values in data$y[[", q, "]] with corresponding non-missing values in data$det.covs[[", q, "]].\nRemoving these site/year/replicate combinations for fitting the model."))
          }
          data$det.covs[[q]][y.missing, i] <- NA
        }
      }
    }
  }
  
  # Check whether random effects are sent in as numeric, and
  # return error if they are. 
  # Occurrence ----------------------
  if (!is.null(findbars(occ.formula))) {
    occ.re.names <- sapply(findbars(occ.formula), all.vars)
    for (i in 1:length(occ.re.names)) {
      if (is(data$occ.covs[, occ.re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", occ.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$occ.covs[, occ.re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", occ.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # Detection -----------------------
  for (q in 1:n.data) {
    if (!is.null(findbars(det.formula[[q]]))) {
      det.re.names <- sapply(findbars(det.formula[[q]]), all.vars)
      for (i in 1:length(det.re.names)) {
        if (is(data$det.covs[[q]][, det.re.names[i]], 'factor')) {
          stop(paste("error: random effect variable ", det.re.names[i], " in data source ", q, " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
        } 
        if (is(data$det.covs[[q]][, det.re.names[i]], 'character')) {
          stop(paste("error: random effect variable ", det.re.names[i], "in data source ", q, " specified as character. Random effect variables must be specified as numeric.", sep = ''))
        }
      }
    }
  }

  # Check ar1 parameter ---------------------------------------------------
  if (!(ar1 %in% c(TRUE, FALSE))) {
    stop("error: ar1 must be either TRUE or FALSE")
  }
 
  # Formula -------------------------------------------------------------
  # Occupancy -----------------------
  if (is(occ.formula, 'formula')) {
    tmp <- parseFormula(occ.formula, data$occ.covs)
    X <- as.matrix(tmp[[1]])
    X.re <- as.matrix(tmp[[4]])
    x.re.names <- colnames(X.re)
    x.names <- tmp[[2]]
  } else {
    stop("error: occ.formula is misspecified")
  }
  # Get RE level names
  re.level.names <- lapply(data$occ.covs[, x.re.names, drop = FALSE],
			   function (a) sort(unique(a)))

  # Detection -----------------------
  X.p <- list()
  X.p.re <- list()
  x.p.names <- list()
  x.p.re.names <- list()
  p.re.level.names <- list()
  for (i in 1:n.data) {
    if (is(det.formula[[i]], 'formula')) {
      tmp <- parseFormula(det.formula[[i]], data$det.covs[[i]])
      X.p[[i]] <- as.matrix(tmp[[1]])
      x.p.names[[i]] <- tmp[[2]]
      if (ncol(tmp[[4]]) > 0) {
        X.p.re[[i]] <- as.matrix(tmp[[4]])
        x.p.re.names[[i]] <- colnames(X.p.re[[i]])
        p.re.level.names[[i]] <- lapply(data$det.covs[[i]][, x.p.re.names[[i]], drop = FALSE],
                                        function (a) sort(unique(a)))
      } else {
        X.p.re[[i]] <- matrix(NA, 0, 0)
        x.p.re.names[[i]] <- NULL
        p.re.level.names[[i]] <- NULL
      }
    } else {
      stop(paste("error: det.formula for data source ", i, " is misspecified", sep = ''))
    }
  }
  x.p.names <- unlist(x.p.names)
  x.p.re.names <- unlist(x.p.re.names)
  
  # Get basic info from inputs ------------------------------------------
  # J = total number of sites
  # n.years.total = total number of years
  # Number of occupancy effects
  p.occ <- ncol(X)
  # Number of occurrence random effect parameters
  p.occ.re <- ncol(X.re)
  # Number of detection random effect parameters for each data set
  p.det.re.by.data <- sapply(X.p.re, ncol)
  # Total number of detection random effects
  p.det.re <- sum(p.det.re.by.data)
  # Number of detection parameters for each data set
  p.det.long <- sapply(X.p, function(a) dim(a)[[2]])
  # Total number of detection parameters
  p.det <- sum(p.det.long)
  # Number of latent occupancy random effect values
  n.occ.re <- length(unlist(apply(X.re, 2, unique)))
  n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of levels for each detection random effect
  n.det.re.long <- unlist(sapply(X.p.re, function(a) apply(a, 2, function(b) length(unique(b)))))
  # Number of latent detection random effects for each data set
  n.det.re.by.data <- sapply(sapply(X.p.re, function(a) apply(a, 2, function(b) length(unique(b)))), sum)
  # Total number of detection random effect levels
  n.det.re <- sum(n.det.re.by.data)
  # Number of replicates at each site/year combo. This also inherently contains
  # info on which sites were sampled in which years.
  n.rep <- lapply(y, function(a1) apply(a1, c(1, 2), function(a2) sum(!is.na(a2))))
  # Max number of repeat visits for each data set
  K.long.max <- sapply(y, function(a) dim(a)[3])
  # Check SVC columns -----------------------------------------------------
  if (is.character(svc.cols)) {
    # Check if all column names in svc are in occ.covs
    if (!all(svc.cols %in% x.names)) {
      missing.cols <- svc.cols[!(svc.cols %in% x.names)]
      stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
    }
    # Convert desired column names into the numeric column index
    svc.cols <- (1:p.occ)[x.names %in% svc.cols]
    
  } else if (is.numeric(svc.cols)) {
    # Check if all column indices are in 1:p.occ
    if (!all(svc.cols %in% 1:p.occ)) {
      missing.cols <- svc.cols[!(svc.cols %in% (1:p.occ))]
      stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
    }
  }
  p.svc <- length(svc.cols)
  
  # A few more checks -----------------------------------------------------
  if (missing(n.batch)) {
    stop("error: must specify number of MCMC batches")
  }
  if (missing(batch.length)) {
    stop("error: must specify length of each MCMC batch")
  }
  n.samples <- n.batch * batch.length
  if (n.burn > n.samples) {
    stop("error: n.burn must be less than n.samples")
  }
  if (n.thin > n.samples) {
    stop("error: n.thin must be less than n.samples")
  }
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
  }
  
  # Get indices for mapping different values in Z. 
  z.site.indx <- rep(1:J, n.years.total) - 1
  z.year.indx <- rep(1:n.years.total, each = J) - 1
  # Form index that keeps track of whether or not the given site is sampled (in any data set) 
  # during a given year
  z.dat.indx <- matrix(0, J, n.years.total)
  for (q in 1:n.data) {
    for (j in 1:nrow(n.rep[[q]])) {
      for (t in 1:ncol(n.rep[[q]])) {
        if (z.dat.indx[sites[[q]][j], seasons[[q]][t]] == 0) {
          z.dat.indx[sites[[q]][j], seasons[[q]][t]] <- ifelse(n.rep[[q]][j, t] > 0, 1, 0)
        }
      }
    }
  }
  
  
  # Get indices to map z to y -------------------------------------------
  X.p.orig <- X.p
  names.long <- list()
  names.re.long <- list()
  # Remove missing observations when the covariate data are available but
  # there are missing detection-nondetection data
  for (i in 1:n.data) {
    if (nrow(X.p[[i]]) == length(y[[i]])) {
      X.p[[i]] <- X.p[[i]][!is.na(y[[i]]), , drop = FALSE]
    }
    if (nrow(X.p.re[[i]]) == length(y[[i]]) & p.det.re.by.data[i] > 0) {
      X.p.re[[i]] <- X.p.re[[i]][!is.na(y[[i]]), , drop = FALSE]
    }
    # Need these for later on
    names.long[[i]] <- which(!is.na(y[[i]]))
  }
  n.obs.long <- sapply(X.p, nrow)
  n.obs <- sum(n.obs.long)
  # z.long.indx = links each element in the psi[j, t] (J x n.year.total) matrix to 
  #               the corresponding observation in the y array. Each value is the 
  #               cell in psi[j, t]. 
  z.long.indx.r <- list()
  # Matrix indicating the index in vector format for each cell of z
  z.ind.mat <- matrix(1:(J * n.years.total), J, n.years.total)
  for (i in 1:n.data) {
    z.long.indx.r[[i]] <- rep(c(z.ind.mat[sites[[i]], seasons[[i]]]), K.long.max[i])
    z.long.indx.r[[i]] <- z.long.indx.r[[i]][!is.na(c(y[[i]]))]
  }
  z.long.indx.r <- unlist(z.long.indx.r)
  # Subtract 1 for c indices
  z.long.indx.c <- z.long.indx.r - 1
  # Get indices for WAIC calculation directly in C. 
  # n.waic is the number of site-year combinations for each data set, even if 
  # some of them aren't sampled.
  n.waic <- 0
  for (q in 1:n.data) {
    n.waic <- n.waic + J.long[q] * n.years.by.data[q]
  }
  # Links each of the total site-year combinations in a data set to the corresponding
  # cell in psi.
  waic.cell.indx <- list()
  for (q in 1:n.data) {
    waic.cell.indx[[q]] <- c(z.ind.mat[sites[[q]], seasons[[q]]])
  }
  waic.cell.indx <- unlist(waic.cell.indx) - 1
  # This should link each waic cell to a given n.obs value. 
  waic.n.obs.indx <- list()
  tmp.start <- 0
  for (i in 1:n.data) {
    tmp.vals <- rep(1:(J.long[i] * n.years.by.data[i]), K.long.max[i])
    tmp.vals <- tmp.vals[!is.na(c(y[[i]]))]
    waic.n.obs.indx[[i]] <- tmp.vals + tmp.start
    tmp.start <- tmp.start + (J.long[i] * n.years.by.data[i])
  }
  waic.n.obs.indx <- unlist(waic.n.obs.indx) - 1
  y <- unlist(y)
  y <- y[!is.na(y)]
  # Index indicating the data set associated with each data point in y
  data.indx.r <- rep(NA, n.obs)
  indx <- 1
  for (i in 1:n.data) {
    data.indx.r[indx:(indx + n.obs.long[i] - 1)] <- rep(i, n.obs.long[i])
    indx <- indx + n.obs.long[i]
  }
  data.indx.c <- data.indx.r - 1
  
  X.p.all <- matrix(NA, n.obs, max(p.det.long))
  indx <- 1
  for (i in 1:n.data) {
    X.p.all[indx:(indx + nrow(X.p[[i]]) - 1), 1:p.det.long[i]] <- X.p[[i]] 
    indx <- indx + nrow(X.p[[i]])
  }
  
  # Get random effect matrices all set ----------------------------------
  if (p.occ.re > 1) {
    for (j in 2:p.occ.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }
  # Get random effect matrices all set ----------------------------------
  if (p.occ.re > 1) {
    for (j in 2:p.occ.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }
  # Detection REs -------------------
  # Need to give a different value for each level across different random
  # effects within a given data set and across a given data set.   
  # Total number of detection random effect observations. 
  n.obs.re <- sum(sapply(X.p.re, nrow))
  curr.max <- 0
  for (i in 1:n.data) {
    if (p.det.re.by.data[i] > 0) {
      for (j in 1:p.det.re.by.data[i]) {
        X.p.re[[i]][, j] <- X.p.re[[i]][, j] + curr.max
        curr.max <- max(X.p.re[[i]]) + 1
      }
    }
  }
  # Combine all detection REs into one group. 
  X.p.re.all <- matrix(NA, n.obs, max(p.det.re.by.data))
  indx <- 1
  for (i in 1:n.data) {
    if (p.det.re.by.data[i] > 0) {
      X.p.re.all[indx:(indx + nrow(X.p.re[[i]]) - 1), 1:p.det.re.by.data[i]] <- X.p.re[[i]]
    }
    indx <- indx + nrow(X.p[[i]])
  }
  # Number of random effects for each row of X.p.re.all
  alpha.n.re.indx <- apply(X.p.re.all, 1, function(a) sum(!is.na(a)))
  

  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))
  # beta -----------------------
  if ("beta.normal" %in% names(priors)) {
    if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
      stop("error: beta.normal must be a list of length 2")
    }
    mu.beta <- priors$beta.normal[[1]]
    sigma.beta <- priors$beta.normal[[2]]
    if (length(mu.beta) != p.occ & length(mu.beta) != 1) {
      if (p.occ == 1) {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p.occ, " with elements corresponding to betas' mean", sep = ""))
      } else {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p.occ, " or 1 with elements corresponding to betas' mean", sep = ""))
      }
    }
    if (length(sigma.beta) != p.occ & length(sigma.beta) != 1) {
      if (p.occ == 1) {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
      	   p.occ, " with elements corresponding to betas' variance", sep = ""))
      } else {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
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
  # alpha -----------------------
  if ("alpha.normal" %in% names(priors)) {
    if (!is.list(priors$alpha.normal) | length(priors$alpha.normal) != 2) {
      stop("error: alpha.normal must be a list of length 2")
    }
    mu.alpha <- priors$alpha.normal[[1]]
    sigma.alpha <- priors$alpha.normal[[2]]
    if (length(mu.alpha) != n.data | !is.list(mu.alpha)) {
      stop(paste("error: alpha.normal[[1]] must be a list of length ", 
                 n.data, " with elements corresponding to alphas' mean for each data set", sep = ""))
    }
    for (q in 1:n.data) {
      if (length(mu.alpha[[q]]) != p.det.long[q] & length(mu.alpha[[q]]) != 1) {
        if (p.det.long[q] == 1) {
          stop(paste("error: prior means for alpha.normal[[1]][[", q, "]] must be a vector of length ", 
                     p.det.long[q], sep = ""))
        } else {
          stop(paste("error: prior means for alpha.normal[[1]][[", q, "]] must be a vector of length ", 
                     p.det.long[q], "or 1", sep = ""))
        }
      }
      if (length(mu.alpha[[q]]) != p.det.long[q]) {
        mu.alpha[[q]] <- rep(mu.alpha[[q]], p.det.long[q])
      }
    }
    mu.alpha <- unlist(mu.alpha)
    if (length(sigma.alpha) != n.data | !is.list(sigma.alpha)) {
      stop(paste("error: alpha.normal[[2]] must be a list of length ", 
                 n.data, " with elements corresponding to alphas' variance for each data set", sep = ""))
    }
    for (q in 1:n.data) {
      if (length(sigma.alpha[[q]]) != p.det.long[q] & length(sigma.alpha[[q]]) != 1) {
        if (p.det.long[q] == 1) {
          stop(paste("error: prior variances for alpha.normal[[2]][[", q, "]] must be a vector of length ", 
                     p.det.long[q], sep = ""))
        } else {
          stop(paste("error: prior variances for alpha.normal[[2]][[", q, "]] must be a vector of length ", 
                     p.det.long[q], " or 1", sep = ""))
          
        }
      }
      if (length(sigma.alpha[[q]]) != p.det.long[q]) {
        sigma.alpha[[q]] <- rep(sigma.alpha[[q]], p.det.long[q])
      }
    }
    sigma.alpha <- unlist(sigma.alpha)
  } else {
    if (verbose) {
      message("No prior specified for alpha.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.alpha <- rep(0, p.det)
    sigma.alpha <- rep(2.72, p.det) 
  }

  # phi -----------------------------
  if ("phi.unif" %in% names(priors)) {
    if (!is.list(priors$phi.unif) | length(priors$phi.unif) != 2) {
      stop("error: phi.unif must be a list of length 2")
    }
    phi.a <- priors$phi.unif[[1]]
    phi.b <- priors$phi.unif[[2]]
    if (length(phi.a) != p.svc & length(phi.a) != 1) {
      stop(paste("error: phi.unif[[1]] must be a vector of length ", 
                 p.svc, 
                 " or 1 with elements corresponding to phis' lower bound for each covariate with spatially-varying coefficients",
                 sep = ""))
    }
    if (length(phi.b) != p.svc & length(phi.b) != 1) {
      stop(paste("error: phi.unif[[2]] must be a vector of length ", 
                 p.svc, 
                 " or 1 with elements corresponding to phis' upper bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(phi.a) != p.svc) {
      phi.a <- rep(phi.a, p.svc)
    }
    if (length(phi.b) != p.svc) {
      phi.b <- rep(phi.b, p.svc)
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
  if (("sigma.sq.ig" %in% names(priors)) & ("sigma.sq.unif" %in% names(priors))) {
    stop("error: cannot specify both an IG and a uniform prior for sigma.sq")
  }
  if ("sigma.sq.ig" %in% names(priors)) {
    sigma.sq.ig <- TRUE
    if (!is.list(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
      stop("error: sigma.sq.ig must be a list of length 2")
    }
    sigma.sq.a <- priors$sigma.sq.ig[[1]]
    sigma.sq.b <- priors$sigma.sq.ig[[2]]
    if (length(sigma.sq.a) != p.svc & length(sigma.sq.a) != 1) {
      stop(paste("error: sigma.sq.ig[[1]] must be a vector of length ", 
                 p.svc, " or 1 with elements corresponding to sigma.sqs' shape for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(sigma.sq.b) != p.svc & length(sigma.sq.b) != 1) {
      stop(paste("error: sigma.sq.ig[[2]] must be a vector of length ", 
                 p.svc, " or 1 with elements corresponding to sigma.sqs' scale for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(sigma.sq.a) != p.svc) {
      sigma.sq.a <- rep(sigma.sq.a, p.svc)
    }
    if (length(sigma.sq.b) != p.svc) {
      sigma.sq.b <- rep(sigma.sq.b, p.svc)
    }
  } else if ("sigma.sq.unif" %in% names(priors)) { # uniform prior
    sigma.sq.ig <- FALSE
    if (!is.list(priors$sigma.sq.unif) | length(priors$sigma.sq.unif) != 2) {
      stop("error: sigma.sq.unif must be a list of length 2")
    }
    sigma.sq.a <- priors$sigma.sq.unif[[1]]
    sigma.sq.b <- priors$sigma.sq.unif[[2]]
    if (length(sigma.sq.a) != p.svc & length(sigma.sq.a) != 1) {
      stop(paste("error: sigma.sq.unif[[1]] must be a vector of length ", 
                 p.svc, " or 1 with elements corresponding to sigma.sqs' shape for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(sigma.sq.b) != p.svc & length(sigma.sq.b) != 1) {
      stop(paste("error: sigma.sq.unif[[2]] must be a vector of length ", 
                 p.svc, " or 1 with elements corresponding to sigma.sqs' scale for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(sigma.sq.a) != p.svc) {
      sigma.sq.a <- rep(sigma.sq.a, p.svc)
    }
    if (length(sigma.sq.b) != p.svc) {
      sigma.sq.b <- rep(sigma.sq.b, p.svc)
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
  if (cov.model == "matern") {
    if (!"nu.unif" %in% names(priors)) {
      stop("error: nu.unif must be specified in priors value list")
    }
    nu.a <- priors$nu.unif[[1]]
    nu.b <- priors$nu.unif[[2]]
    if (!is.list(priors$nu.unif) | length(priors$nu.unif) != 2) {
      stop("error: nu.unif must be a list of length 2")
    }
    if (length(nu.a) != p.svc & length(nu.a) != 1) {
      stop(paste("error: nu.unif[[1]] must be a vector of length ", 
                 p.svc, " or 1 with elements corresponding to nus' lower bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(nu.b) != p.svc & length(nu.b) != 1) {
      stop(paste("error: nu.unif[[2]] must be a vector of length ", 
                 p.svc, " or 1 with elements corresponding to nus' upper bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(nu.a) != p.svc) {
      nu.a <- rep(nu.a, p.svc)
    }
    if (length(nu.b) != p.svc) {
      nu.b <- rep(nu.b, p.svc)
    }
  } else {
    nu.a <- rep(0, p.svc)
    nu.b <- rep(0, p.svc)
  }
  
  # sigma.sq.psi --------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.psi.ig) | length(priors$sigma.sq.psi.ig) != 2) {
        stop("error: sigma.sq.psi.ig must be a list of length 2")
      }
      sigma.sq.psi.a <- priors$sigma.sq.psi.ig[[1]]
      sigma.sq.psi.b <- priors$sigma.sq.psi.ig[[2]]
      if (length(sigma.sq.psi.a) != p.occ.re & length(sigma.sq.psi.a) != 1) {
        if (p.occ.re == 1) {
          stop(paste("error: sigma.sq.psi.ig[[1]] must be a vector of length ", 
                     p.occ.re, " with elements corresponding to sigma.sq.psis' shape", sep = ""))
        } else {
          stop(paste("error: sigma.sq.psi.ig[[1]] must be a vector of length ", 
                     p.occ.re, " or 1 with elements corresponding to sigma.sq.psis' shape", sep = ""))
        }
      }
      if (length(sigma.sq.psi.b) != p.occ.re & length(sigma.sq.psi.b) != 1) {
        if (p.occ.re == 1) {
          stop(paste("error: sigma.sq.psi.ig[[2]] must be a vector of length ", 
                     p.occ.re, " with elements corresponding to sigma.sq.psis' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.psi.ig[[2]] must be a vector of length ", 
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
  # sigma.sq.p --------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.p.ig) | length(priors$sigma.sq.p.ig) != 2) {
        stop("error: sigma.sq.p.ig must be a list of length 2")
      }
      sigma.sq.p.a <- priors$sigma.sq.p.ig[[1]]
      sigma.sq.p.b <- priors$sigma.sq.p.ig[[2]]
      if (length(sigma.sq.p.a) != p.det.re & length(sigma.sq.p.a) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: sigma.sq.p.ig[[1]] must be a vector of length ", 
                     p.det.re, " with elements corresponding to sigma.sq.ps' shape", sep = ""))
        } else {
          stop(paste("error: sigma.sq.p.ig[[1]] must be a vector of length ", 
                     p.det.re, " or 1 with elements corresponding to sigma.sq.ps' shape", sep = ""))
        }
      }
      if (length(sigma.sq.p.b) != p.det.re & length(sigma.sq.p.b) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: sigma.sq.p.ig[[2]] must be a vector of length ", 
                     p.det.re, " with elements corresponding to sigma.sq.ps' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.p.ig[[2]] must be a vector of length ", 
                     p.det.re, " or 1 with elements corresponding to sigma.sq.ps' scale", sep = ""))
        }
      }
      if (length(sigma.sq.p.a) != p.det.re) {
        sigma.sq.p.a <- rep(sigma.sq.p.a, p.det.re)
      }
      if (length(sigma.sq.p.b) != p.det.re) {
        sigma.sq.p.b <- rep(sigma.sq.p.b, p.det.re)
      }
    }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.p.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.p.a <- rep(0.1, p.det.re)
      sigma.sq.p.b <- rep(0.1, p.det.re)
    }
  } else {
    sigma.sq.p.a <- 0
    sigma.sq.p.b <- 0
  }
  if (ar1) {
    # rho -----------------------------
    if ("rho.unif" %in% names(priors)) {
      if (!is.vector(priors$rho.unif) | !is.atomic(priors$rho.unif) | length(priors$rho.unif) != 2) {
        stop("error: rho.unif must be a vector of length 2 with elements corresponding to rho's lower and upper bounds")
      }
      rho.a <- priors$rho.unif[1]
      rho.b <- priors$rho.unif[2]
    } else {
      if (verbose) {
        message("No prior specified for rho.unif.\nSetting uniform bounds to -1 and 1.\n")
      }
      rho.a <- -1 
      rho.b <- 1 
    }
    
    # sigma.sq.t.t ----------------------
    if ("sigma.sq.t.ig" %in% names(priors)) { 
      if (!is.vector(priors$sigma.sq.t.ig) | !is.atomic(priors$sigma.sq.t.ig) | length(priors$sigma.sq.t.ig) != 2) {
        stop("error: sigma.sq.t.ig must be a vector of length 2 with elements corresponding to sigma.sq.t's shape and scale parameters")
      }
      sigma.sq.t.a <- priors$sigma.sq.t.ig[1]
      sigma.sq.t.b <- priors$sigma.sq.t.ig[2]
    } else {
      if (verbose) {
        message("No prior specified for sigma.sq.t.\nUsing an inverse-Gamma prior with the shape parameter set to 2 and scale parameter to 0.5.\n")
      }
      sigma.sq.t.a <- 2
      sigma.sq.t.b <- 0.5
    }
  } else {
    rho.a <- 0
    rho.b <- 0
    sigma.sq.t.a <- 0
    sigma.sq.t.b <- 0
  }

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # z -------------------------------
  #ORDER: stored in column-major order as a vector, with values sorted by year, then 
  #       site within year. 
  z.inits.default <- matrix(0, J, n.years.total)
  for (q in 1:n.data) {
    for (j in 1:J.long[q]) {
      for (t in 1:n.years.by.data[q]) {
        if (sum(!is.na(y.big[[q]][j, t, ])) > 0) {
          if (z.inits.default[sites[[q]][j], seasons[[q]][t]] == 0) {
            z.inits.default[sites[[q]][j], seasons[[q]][t]] <- max(y.big[[q]][j, t, ], na.rm = TRUE)
          }
        }
      }
    }
  }
  if ("z" %in% names(inits)) {
    z.inits <- inits$z
    if (!is.matrix(z.inits) & !is.data.frame(z.inits)) {
      stop(paste("error: initial values for z must be a matrix or data frame with ",
      	   J, " rows and ", n.years.total, " columns.", sep = ""))
    }
    if (nrow(z.inits) != J | ncol(z.inits) != n.years.total) {
      stop(paste("error: initial values for z must be a matrix or data frame with ",
      	   J, " rows and ", n.years.total, " columns.", sep = ""))
    }
    # Reorder the user supplied inits values for NNGP models
    if (NNGP) {
      z.inits <- z.inits[ord, ]
    }
    init.test <- sum(z.inits < z.inits.default, na.rm = TRUE)
    if (init.test > 0) {
      stop("initial values for latent occurrence (z) are invalid. Please re-specify inits$z so initial values are 1 if the species is observed at that site/year combination.")
    }
  } else {
    z.inits <- z.inits.default
    if (verbose) {
      message("z is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
  # beta -----------------------
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (length(beta.inits) != p.occ & length(beta.inits) != 1) {
      if (p.occ == 1) {
        stop(paste("error: initial values for beta must be of length ", p.occ,
      	     sep = ""))

      } else {
        stop(paste("error: initial values for beta must be of length ", p.occ, " or 1",
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
  # alpha -----------------------------
  if ("alpha" %in% names(inits)) {
    alpha.inits <- inits[["alpha"]]
    if (length(alpha.inits) != n.data | !is.list(alpha.inits)) {
      stop(paste("error: initial values for alpha must be a list of length ", n.data,
                 sep = ""))
    }
    for (q in 1:n.data) {
      if (length(alpha.inits[[q]]) != p.det.long[q] & length(alpha.inits[[q]]) != 1) {
        if (p.det.long[q] == 1) {
          stop(paste("error: initial values for alpha[[", q, "]] must be a vector of length ", 
                     p.det.long[q], sep = ""))
        } else {
          stop(paste("error: initial values for alpha[[", q, "]] must be a vector of length ", 
                     p.det.long[q], " or 1", sep = ""))
        }
      }
      if (length(alpha.inits[[q]]) != p.det.long[q]) {
        alpha.inits[[q]] <- rep(alpha.inits, p.det.long[q])
      }
    }
    alpha.inits <- unlist(alpha.inits)
  } else {
    if (verbose) {
      message("alpha is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
    alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
  }
  
  alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
  alpha.indx.c <- alpha.indx.r - 1
  # phi -------------------------------
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != p.svc & length(phi.inits) != 1) {
      stop(paste("error: initial values for phi must be of length ", p.svc, " or 1", 
                 sep = ""))
    }
    if (length(phi.inits) != p.svc) {
      phi.inits <- rep(phi.inits, p.svc)
    }
  } else {
    phi.inits <- runif(p.svc, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random values from the prior distribution\n")
    }
  }
  # sigma.sq ------------------------
  if ("sigma.sq" %in% names(inits)) {
    sigma.sq.inits <- inits[["sigma.sq"]]
    if (length(sigma.sq.inits) != p.svc & length(sigma.sq.inits) != 1) {
      stop(paste("error: initial values for sigma.sq must be of length ", p.svc,  " or 1",
                 sep = ""))
    }
    if (length(sigma.sq.inits) != p.svc) {
      sigma.sq.inits <- rep(sigma.sq.inits, p.svc)
    }
  } else {
    if (sigma.sq.ig) {
      sigma.sq.inits <- runif(p.svc, 0.1, 10)
    } else {
      sigma.sq.inits <- runif(p.svc, sigma.sq.a, sigma.sq.b)
    }
    if (verbose) {
      message("sigma.sq is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  # nu ------------------------
  if ("nu" %in% names(inits)) {
    nu.inits <- inits[["nu"]]
    if (length(nu.inits) != p.svc & length(nu.inits) != 1) {
      stop(paste("error: initial values for nu must be of length ", p.svc,  " or 1",
                 sep = ""))
    }
    if (length(nu.inits) != p.svc) {
      nu.inits <- rep(nu.inits, p.svc)
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      nu.inits <- runif(p.svc, nu.a, nu.b)
    } else {
      nu.inits <- rep(0, p.svc)
    }
  }
  # w ---------------------------------
  # Just set initial W values to 0. 
  w.inits <- rep(0, J * p.svc)
  # sigma.sq.psi -------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi" %in% names(inits)) {
      sigma.sq.psi.inits <- inits[["sigma.sq.psi"]]
      if (length(sigma.sq.psi.inits) != p.occ.re & length(sigma.sq.psi.inits) != 1) {
        if (p.occ.re == 1) {
          stop(paste("error: initial values for sigma.sq.psi must be of length ", p.occ.re, 
      	       sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.psi must be of length ", p.occ.re, 
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
  # sigma.sq.p ------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p" %in% names(inits)) {
      sigma.sq.p.inits <- inits[["sigma.sq.p"]]
      if (length(sigma.sq.p.inits) != p.det.re & length(sigma.sq.p.inits) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: initial values for sigma.sq.p must be of length ", p.det.re, 
                     sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.p must be of length ", p.det.re, 
                     " or 1", sep = ""))
          
        }
      }
      if (length(sigma.sq.p.inits) != p.det.re) {
        sigma.sq.p.inits <- rep(sigma.sq.p.inits, p.det.re)
      }
    } else {
      sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.p is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    # Keep track of which detection random effect you're on. 
    alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
    # Index that indicates the column in X.p.re.all
    alpha.col.list <- list()
    indx <- 1
    for (i in 1:n.data) {
      if (p.det.re.by.data[i] > 0) { 
        for (j in 1:p.det.re.by.data[i]) {
          if (j > 1) {
            alpha.col.list[[i]] <- c(alpha.col.list[[i]], rep(j - 1, n.det.re.long[indx]))
          } else {
            alpha.col.list[[i]] <- rep(j - 1, n.det.re.long[indx])
          }
          indx <- indx + 1
        }
      }
    }
    alpha.col.indx <- unlist(alpha.col.list)
    # Index that indicates the data source the random effect corresponds to. 
    alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
  } else {
    sigma.sq.p.inits <- 0
    alpha.star.indx <- 0
    alpha.star.inits <- 0
    alpha.col.indx <- 0
  }
  if (ar1) {
    # rho -----------------------------
    if ("rho" %in% names(inits)) {
      rho.inits <- inits[["rho"]]
      if (length(rho.inits) != 1) {
        stop("error: initial values for rho must be of length 1")
      }
    } else {
      rho.inits <- runif(1, rho.a, rho.b)
      if (verbose) {
        message("rho is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
      }
    }
    # sigma.sq.t ------------------------
    if ("sigma.sq.t" %in% names(inits)) {
      sigma.sq.t.inits <- inits[["sigma.sq.t"]]
      if (length(sigma.sq.t.inits) != 1) {
        stop("error: initial values for sigma.sq.t must be of length 1")
      }
    } else {
      sigma.sq.t.inits <- runif(1, 0.5, 10)
      if (verbose) {
        message("sigma.sq.t is not specified in initial values.\nSetting initial value to random value between 0.5 and 10\n")
      }
    }
  } else {
    rho.inits <- 0
    sigma.sq.t.inits <- 0
  }
  # Should initial values be fixed --
  if ("fix" %in% names(inits)) {
    fix.inits <- inits[["fix"]]
    if ((fix.inits != TRUE) & (fix.inits != FALSE)) {
      stop(paste("error: inits$fix must take value TRUE or FALSE"))
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
    stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ", 
         paste(cov.model.names, collapse=", ", sep="") ,".")}
  # Obo for cov model lookup on c side
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  
  # Prep for SVCs ---------------------------------------------------------
  X.w <- X[, svc.cols, drop = FALSE]

  # Get tuning values ---------------------------------------------------
  # Not accessed, but necessary to keep things in line. 
  sigma.sq.tuning <- rep(0, p.svc)
  phi.tuning <- rep(0, p.svc)
  nu.tuning <- rep(0, p.svc)
  rho.tuning <- 0
  sigma.sq.t.tuning <- 0
  if (missing(tuning)) {
    phi.tuning <- rep(1, p.svc)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, p.svc)
    }
    rho.tuning <- 1
    if (!sigma.sq.ig) {
      sigma.sq.tuning <- rep(1, p.svc)
    }
  } else {
    names(tuning) <- tolower(names(tuning))
    # phi ---------------------------
    if(!"phi" %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (length(phi.tuning) == 1) {
      phi.tuning <- rep(tuning$phi, p.svc)
    } else if (length(phi.tuning) != p.svc) {
      stop(paste("error: phi tuning must be either a single value or a vector of length ",
                 p.svc, sep = ""))
    }
    if (!sigma.sq.ig) {
      # sigma.sq --------------------------
      if(!"sigma.sq" %in% names(tuning)) {
        stop("error: sigma.sq must be specified in tuning value list")
      }
      sigma.sq.tuning <- tuning$sigma.sq
      if (length(sigma.sq.tuning) == 1) {
        sigma.sq.tuning <- rep(tuning$sigma.sq, p.svc)
      } else if (length(sigma.sq.tuning) != p.svc) {
        stop(paste("error: sigma.sq tuning must be either a single value or a vector of length ",
                   p.svc, sep = ""))
      }
    }
    if (ar1) {
      # rho ---------------------------
      if(!"rho" %in% names(tuning)) {
        stop("error: rho must be specified in tuning value list")
      }
      rho.tuning <- tuning$rho
      if (length(rho.tuning) != 1) {
        stop("error: rho tuning must be a single value")
      } 
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) == 1) {
        nu.tuning <- rep(tuning$nu, p.svc)
      } else if (length(nu.tuning) != p.svc) {
        stop(paste("error: nu tuning must be either a single value or a vector of length ",
                   p.svc, sep = ""))
      }
    }
  }
  if (ar1) {
    if (cov.model == 'matern') {
      tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning, 
                        sigma.sq.t.tuning, rho.tuning))
    } else {
      tuning.c <- log(c(sigma.sq.tuning, phi.tuning,  
                        sigma.sq.t.tuning, rho.tuning))
    }
  } else {
    tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  }
  
  curr.chain <- 1

  if (!NNGP) {
    stop("error: svcTIntPGOcc is currently only implemented for NNGP models. Please set NNGP = TRUE") 
  } else {
    # Nearest Neighbor Search ---------------------------------------------
    if(verbose){
      cat("----------------------------------------\n");
      cat("\tBuilding the neighbor list\n");
      cat("----------------------------------------\n");
    }

    search.type.names <- c("brute", "cb")
    
    if(!search.type %in% search.type.names){
      stop("error: specified search.type '",search.type,
	   "' is not a valid option; choose from ", 
	   paste(search.type.names, collapse=", ", sep="") ,".")
    }
    
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

    if(verbose){
      cat("----------------------------------------\n");
      cat("Building the neighbors of neighbors list\n");
      cat("----------------------------------------\n");
    }
    
    indx <- mkUIndx(J, n.neighbors, nn.indx, nn.indx.lu, u.search.type)
    
    u.indx <- indx$u.indx
    u.indx.lu <- indx$u.indx.lu
    ui.indx <- indx$ui.indx
    u.indx.run.time <- indx$run.time

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(z.inits) <- "double"
    storage.mode(X.p.all) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p.det.long) <- 'integer'
    storage.mode(n.obs.long) <- 'integer'
    storage.mode(J.long) <- 'integer'
    storage.mode(X.w) <- 'double'
    n.post.samples <- length(seq(from = n.burn + 1, 
                                 to = n.samples, 
                                 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"
    consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.det, p.det.re, 
                n.det.re, n.years.total, n.data, n.waic, sigma.sq.ig, ar1, 
                cov.model.indx, n.neighbors, n.batch, batch.length, 
                n.omp.threads, verbose, n.report, n.thin, n.burn, 
                n.post.samples, p.svc)
    storage.mode(consts) <- "integer"
    storage.mode(coords) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(z.long.indx.c) <- "integer"
    storage.mode(z.year.indx) <- "integer"
    storage.mode(z.dat.indx) <- "integer"
    storage.mode(data.indx.c) <- 'integer'
    storage.mode(alpha.indx.c) <- 'integer'
    storage.mode(waic.n.obs.indx) <- 'integer'
    storage.mode(waic.cell.indx) <- 'integer'
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(sigma.alpha) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(sigma.sq.a) <- "double"
    storage.mode(sigma.sq.b) <- "double"
    storage.mode(tuning.c) <- "double"
    storage.mode(accept.rate) <- "double"
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.indx) <- "integer"
    storage.mode(u.indx.lu) <- "integer"
    storage.mode(ui.indx) <- "integer"
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    # For detection random effects
    storage.mode(X.p.re.all) <- "integer"
    storage.mode(p.det.re.by.data) <- 'integer'
    alpha.level.indx <- sort(unique(c(X.p.re.all)))
    storage.mode(alpha.level.indx) <- "integer"
    storage.mode(n.det.re.long) <- "integer"
    storage.mode(sigma.sq.p.inits) <- "double"
    storage.mode(sigma.sq.p.a) <- "double"
    storage.mode(sigma.sq.p.b) <- "double"
    storage.mode(alpha.star.inits) <- "double"
    storage.mode(alpha.star.indx) <- "integer"
    storage.mode(alpha.n.re.indx) <- 'integer'
    storage.mode(alpha.col.indx) <- 'integer'
    # For occurrence random effects
    storage.mode(X.re) <- "integer"
    beta.level.indx <- sort(unique(c(X.re)))
    storage.mode(beta.level.indx) <- "integer"
    storage.mode(sigma.sq.psi.inits) <- "double"
    storage.mode(sigma.sq.psi.a) <- "double"
    storage.mode(sigma.sq.psi.b) <- "double"
    storage.mode(n.occ.re.long) <- "integer"
    storage.mode(beta.star.inits) <- "double"
    storage.mode(beta.star.indx) <- "integer"  
    # AR1 parameters
    storage.mode(ar1) <- "integer"
    ar1.vals <- c(rho.a, rho.b, sigma.sq.t.a, sigma.sq.t.b, 
                  rho.inits, sigma.sq.t.inits)
    storage.mode(ar1.vals) <- "double"
    
    out.tmp <- list()                          
    out <- list()                              
    for (i in 1:n.chains) {                  
      # Change initial values if i > 1       
      if ((i > 1)) {
        beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
        alpha.inits <- rnorm(p.det, mu.alpha, sqrt(sigma.alpha))
        if (sigma.sq.ig) {
          sigma.sq.inits <- runif(p.svc, 0.1, 10)
        } else {
            sigma.sq.inits <- runif(p.svc, sigma.sq.a, sigma.sq.b)
        }
        phi.inits <- runif(p.svc, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(p.svc, nu.a, nu.b)
        }
        if (p.det.re > 0) {
          sigma.sq.p.inits <- runif(p.det.re, 0.5, 10)
          alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
        }
        if (p.occ.re > 0) {
          sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
          beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
        }
        if (ar1) {
          ar1.vals[5] <- runif(1, rho.a, rho.b)
          ar1.vals[6] <- runif(1, 0.5, 10)	
        }
      }
      storage.mode(chain.info) <- "integer"
      out.tmp[[i]] <- .Call("svcTIntPGOccNNGP", y, X, X.w, X.p.all, coords, X.re, X.p.re.all, 
                            consts, p.det.long, J.long, n.obs.long, 
                            n.occ.re.long, n.det.re.long,
                            nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
                            beta.inits, alpha.inits, sigma.sq.psi.inits, 
                            sigma.sq.p.inits, beta.star.inits, alpha.star.inits, 
                            phi.inits, sigma.sq.inits, nu.inits,
                            w.inits, z.inits, z.long.indx.c, data.indx.c, 
                            alpha.indx.c, z.year.indx, z.dat.indx, 
                            beta.star.indx, beta.level.indx,
                            alpha.star.indx, alpha.level.indx, alpha.n.re.indx, 
                            alpha.col.indx, mu.beta, Sigma.beta, mu.alpha, sigma.alpha, 
                            phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b,
                            sigma.sq.psi.a, sigma.sq.psi.b, sigma.sq.p.a, sigma.sq.p.b,
                            ar1.vals, tuning.c, accept.rate, chain.info, 
                            waic.n.obs.indx, waic.cell.indx)
      chain.info[1] <- chain.info[1] + 1
    }

    # Calculate R-Hat ---------------
    if (cov.model != 'matern') {
      theta.names <- paste(rep(c('sigma.sq', 'phi'), each = p.svc), x.names[svc.cols], sep = '-')
    } else {
      theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = p.svc), x.names[svc.cols], sep = '-')
    } 
    if (ar1) {
      theta.names <- c(theta.names, 'sigma.sq.t', 'rho')
    }
    n.theta <- length(theta.names)
    out <- list()
    out$rhat <- list()
    if (n.chains > 1) {
      out$rhat$beta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
                                                             mcmc(t(a$beta.samples)))), 
                                   autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
      out$rhat$alpha <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
                                                              mcmc(t(a$alpha.samples)))), 
                                    autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					                                            mcmc(t(a$theta.samples)))), 
                                    autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
      if (p.det.re > 0) {
        out$rhat$sigma.sq.p <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
                                                                               mcmc(t(a$sigma.sq.p.samples)))), 
                                                     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }
      if (p.occ.re > 0) {
        out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
                                                                                 mcmc(t(a$sigma.sq.psi.samples)))), 
                                                       autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }

    } else {
      out$rhat$beta <- rep(NA, p.occ)
      out$rhat$alpha <- rep(NA, p.det)
      out$rhat$theta <- rep(NA, n.theta)
      if (p.det.re > 0) {
        out$rhat$sigma.sq.p <- rep(NA, p.det.re)
      }
      if (p.occ.re > 0) {
        out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
      }
    }

    # Put everything into an MCMC objects
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- x.names
    out$alpha.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.samples))))
    colnames(out$alpha.samples) <- x.p.names
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    colnames(out$theta.samples) <- theta.names
    if (ar1) {
      out$eta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$eta.samples))))
    } 
    # Return things back in the original order
    out$coords <- coords[order(ord), ]
    out$X <- array(X, dim = c(J, n.years.total, p.occ))
    out$X <- out$X[order(ord), , , drop = FALSE]
    dimnames(out$X)[[3]] <- x.names
    out$X.re <- array(X.re, dim = c(J, n.years.total, p.occ.re))
    out$X.re <- out$X.re[order(ord), , , drop = FALSE]
    dimnames(out$X.re)[[3]] <- x.re.names
    out$X.w <- array(X.w, dim = c(J, n.years.total, p.svc))
    out$X.w <- out$X.w[order(ord), , , drop = FALSE]
    dimnames(out$X.w)[[3]] <- x.names[svc.cols]
    # Account for case when intercept only spatial model. 
    if (p.svc == 1) {
      tmp <- do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples)))
      tmp <- tmp[, order(ord), drop = FALSE]
      out$w.samples <- array(NA, dim = c(n.post.samples * n.chains, p.svc, J))
      out$w.samples[, 1, ] <- tmp
    } else {
      out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
                                                                        dim = c(p.svc, J, n.post.samples))))
      out$w.samples <- out$w.samples[, order(ord), ]
      out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
    }
    out$z.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$z.samples, 
      								dim = c(J, n.years.total, n.post.samples))))
    out$z.samples <- out$z.samples[order(ord), , ]
    out$z.samples <- aperm(out$z.samples, c(3, 1, 2))
    out$psi.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$psi.samples, 
      								dim = c(J, n.years.total, n.post.samples))))
    out$psi.samples <- out$psi.samples[order(ord), , ]
    out$psi.samples <- aperm(out$psi.samples, c(3, 1, 2))
    out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
    # p.samples is returned as a list, where each element 
    # corresponds to a different data set. 
    out$p.samples <- do.call(rbind, lapply(out.tmp, function(a) a$p.samples))
    tmp <- list()
    indx <- 1
    for (q in 1:n.data) {
      tmp[[q]] <- array(NA, dim = c(J.long[q] * K.long.max[q] * n.years.by.data[q], 
                                    n.post.samples * n.chains))
      tmp[[q]][names.long[[q]], ] <- out$p.samples[indx:(indx + n.obs.long[q] - 1), ] 
      tmp[[q]] <- array(tmp[[q]], dim = c(J.long[q], n.years.by.data[q], 
                                          K.long.max[q], n.post.samples * n.chains))
      tmp[[q]] <- aperm(tmp[[q]], c(4, 1, 2, 3))
      indx <- indx + n.obs.long[q]
    }
    out$p.samples <- tmp
    out$y <- y.big
    out$X.p <- X.p.orig
    out$X.p.re <- X.p.re
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
    if (p.det.re > 0) {
      out$sigma.sq.p.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.p.samples))))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.star.samples))))
      tmp.names <- unlist(p.re.level.names)
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$p.re.level.names <- p.re.level.names
    }
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$alpha <- effectiveSize(out$alpha.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    if (p.det.re > 0) {
      out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
    }
    if (p.occ.re > 0) {
      out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
    }
    out$call <- cl
    out$n.samples <- batch.length * n.batch
    out$n.neighbors <- n.neighbors
    out$cov.model.indx <- cov.model.indx
    out$svc.cols <- svc.cols
    out$type <- "NNGP"
    out$sites <- sites.orig
    out$seasons <- seasons
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$ar1 <- as.logical(ar1)
    if (p.det.re > 0) {
      out$pRE <- TRUE
    } else {
      out$pRE <- FALSE
    }
    if (p.occ.re > 0) {
      out$psiRE <- TRUE
    } else {
      out$psiRE <- FALSE
    }
    out$pRELong <- ifelse(p.det.re.by.data > 0, TRUE, FALSE)
  } # NNGP
  class(out) <- "svcTIntPGOcc"
  out$run.time <- proc.time() - ptm
  out
}
