parseFormula <-  function(formula, data, intercept=TRUE, justX=FALSE){

    # Find random effect terms
    bars <- findbars(formula)
    re.terms <- NULL
    if (!is.null(bars)) {
      re.terms <- mkReTrms(bars, data)
    }

    formula <- nobars(formula)

    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data=data)
    if(missing(data)) data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    if (!intercept){
      attributes(mt)$intercept <- 0
    }


    # null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf)
    X <- as.matrix(X)         # X matrix
    xvars <- dimnames(X)[[2]] # X variable names
    xobs  <- dimnames(X)[[1]] # X observation names

    # Get random effects
    X.random <- matrix(NA, nrow(X), length(re.terms$Ztlist))
    X.re <- matrix(NA, nrow(X), length(re.terms$flist))
    # Account for RE only model when there are some missing values
    if (ncol(X.re) > 0) {
      if (nrow(X) > length(re.terms$Ztlist[[1]]@i)) {
        X.re <- matrix(NA, length(re.terms$Ztlist[[1]]@i), 
        	       length(re.terms$flist)) 
        X.random <- matrix(NA, length(re.terms$Ztlist[[1]]@i), 
        	       length(re.terms$Ztlist)) 
      }
    }
    # Get the random factors associated with each element. 
    tmp <- sub('.*\\|\\s*', "", names(re.terms$Ztlist))
    re.col.indx <- match(tmp, unique(tmp)) 
    # Get the unique instance of a random factor. 
    unique.indx <- match(unique(tmp), tmp)
    if (ncol(X.re) > 0) {
      for (j in 1:ncol(X.re)) {
        curr.indx <- unique.indx[j]
	tmp <- as.numeric(re.terms$flist[[re.col.indx[curr.indx]]])
	miss.indx <- is.na(tmp)
        X.re[, j] <- tmp[!miss.indx]
      }
      colnames(X.re) <- names(re.terms$flist)
      X.re <- X.re[, re.col.indx, drop = FALSE]
      for (j in 1:length(re.terms$Ztlist)) {
        # TODO: will need to make sure this is correct. 
        tmp <- re.terms$Ztlist[[j]]@x[re.terms$Ztlist[[j]]@p]
        X.random[, j] <- tmp[!miss.indx] 
      }
      tmp <- sapply(re.terms$cnms, function(a) a[length(a)])
      tmp.2 <- tmp
      attr(tmp.2, 'names') <- NULL
      tmp.2 <- paste(tmp.2, names(tmp), sep = '-') 
      colnames(X.random) <- tmp.2
    }
    return(list(X, xvars, xobs, X.re, X.random))
  }

