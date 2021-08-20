mkUIndx <- function(n, m, nn.indx, nn.indx.lu, search.type){
    
    n.indx <- (1+m)/2*m+(n-m-1)*m
    u.indx <- rep(0, n.indx)
    u.indx.lu <- rep(0, 2*n)
    ui.indx <- rep(0, n.indx)
    
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(nn.indx) <- "integer"
    storage.mode(u.indx) <- "integer"
    storage.mode(u.indx.lu) <- "integer"
    storage.mode(ui.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(search.type) <- "integer"
    
    ptm <- proc.time()
    
    out <- .Call("mkUIndx", n, m, nn.indx, u.indx, u.indx.lu, ui.indx, nn.indx.lu, search.type)
    
    run.time <- proc.time() - ptm
    
    list("run.time"=run.time, "u.indx"=as.integer(u.indx), "u.indx.lu"=as.integer(u.indx.lu), "ui.indx"=as.integer(ui.indx))
}


mkNNIndx <- function(coords, m, n.omp.threads=1){

    n <- nrow(coords)
    nIndx <- (1+m)/2*m+(n-m-1)*m
    nnIndx <- rep(0, nIndx)
    nnDist <- rep(0, nIndx)
    nnIndxLU <- matrix(0, n, 2)

    n <- as.integer(n)
    m <- as.integer(m)
    coords <- as.double(coords)
    nnIndx <- as.integer(nnIndx)
    nnDist <- as.double(nnDist)
    nnIndxLU <- as.integer(nnIndxLU)
    n.omp.threads <- as.integer(n.omp.threads)
    
    ptm <- proc.time()
    
    out <- .Call("mkNNIndx", n, m, coords, nnIndx, nnDist, nnIndxLU, n.omp.threads)

    run.time <- proc.time() - ptm
    
    list("run.time"=run.time, "nnIndx"=as.integer(nnIndx), "nnDist"=as.double(nnDist), "nnIndxLU"=nnIndxLU)

}

mkNNIndxCB <- function(coords, m, n.omp.threads=1){
    
    n <- nrow(coords)
    nIndx <- (1+m)/2*m+(n-m-1)*m
    nnIndx <- rep(0, nIndx)
    nnDist <- rep(0, nIndx)
    nnIndxLU <- matrix(0, n, 2)
    
    n <- as.integer(n)
    m <- as.integer(m)
    coords <- as.double(coords)
    nnIndx <- as.integer(nnIndx)
    nnDist <- as.double(nnDist)
    nnIndxLU <- as.integer(nnIndxLU)
    n.omp.threads <- as.integer(n.omp.threads)

    ptm <- proc.time()
    
    out <- .Call("mkNNIndxCB", n, m, coords, nnIndx, nnDist, nnIndxLU, n.omp.threads)

    run.time <- proc.time() - ptm
    
    list("run.time"=run.time, "nnIndx"=as.integer(nnIndx), "nnDist"=as.double(nnDist), "nnIndxLU"=nnIndxLU)
}
