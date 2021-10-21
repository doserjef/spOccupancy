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
    X.re <- matrix(NA, nrow(X), length(re.terms$Ztlist))
    if (ncol(X.re) > 0) {
      # Support for RE only model
      if (length(re.terms$Ztlist[[1]]@i) != nrow(X)) {
        X.re <- matrix(NA, length(re.terms$Ztlist[[1]]@i), 
		       length(re.terms$Ztlist)) 
      }
      for (j in 1:ncol(X.re)) {
        X.re[, j] <- as.vector(re.terms$Ztlist[[j]]@i)
      }
      colnames(X.re) <- names(re.terms$flist)
    }
    return(list(X, xvars, xobs, X.re))
  }

