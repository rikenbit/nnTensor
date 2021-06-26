GabrielNMF <- function(X, J=3, nx = 5, ny = 5, ...){
    # Argument check
    .checkGabrielNMF(X, J, nx, ny)
    # Bi-Cross-Validation
    int <- .initGabrielNMF(X, nx, ny)
    xholdouts <- int$xholdouts
    yholdouts <- int$yholdouts
    TestRecError <- int$TestRecError
    tcount <- 1
    for(x in seq(nx*ny)){
        M <- X
        M[, ] <- 1
        M[xholdouts[[x]], yholdouts[[x]]] <- 0
        out <- try(.eachGabrielNMF(X, M, J, ...))
        if(class(out) != "try-error"){
            TestRecError[tcount] <- out
        }else{
            TestRecError[tcount] <- NA
        }
        tcount <- tcount + 1
    }
    # Output
    list(TestRecError = TestRecError)
}

.eachGabrielNMF <- function(X, M, J, ...){
    # Split submatrices
    i <- which(rowSums(M) != dim(M)[2])
    j <- which(colSums(M) != dim(M)[1])
    minusi <- which(rowSums(M) == dim(M)[2])
    minusj <- which(colSums(M) == dim(M)[1])
    if(length(i) < J){
        stop("The number of rows of hold-out matrix A is too small!
            Please specify smaller nx!!")
    }
    if(length(j) < J){
        stop("The number of columns of hold-out matrix A is too small!
            Please specify smaller ny!!")
    }
    if(length(minusi) < J){
        stop("The number of rows of hold-out matrix A is too large!
            Please specify larger nx!!")
    }
    if(length(minusj) < J){
        stop("The number of columns of hold-out matrix A is too large!
            Please specify larger ny!!")
    }
    A <- as.matrix(X[i, j])
    B <- as.matrix(X[i, minusj])
    C <- as.matrix(X[minusi, j])
    D <- as.matrix(X[minusi, minusj])
    # NMF of D
    nmfD <- NMF(D, J=J, ...)
    # Prediction of U and V
    U_B <- NMF(B, initV=nmfD$V, fixV=TRUE, J=J, ...)$U
    V_C <- NMF(C, initU=nmfD$U, fixU=TRUE, J=J, ...)$V
    # Bi Cross Validation
    norm(A - U_B %*% t(V_C), "F")^2
}

.checkGabrielNMF <- function(X, J, nx, ny){
    if(!is.matrix(X)){
        stop("Please specify the X as a matrix")
    }
    if(!is.numeric(J)){
        stop("Please specify the J as a number or a numeric vector")
    }
    if(!is.numeric(nx)){
        stop("Please specify the nx as a number or a numeric vector")
    }else{
        if(!(2 <= nx) || !(nx <= nrow(X))){
            stop("Please specify the value of nx within the range 2 <= nx <= nrow(X)")
        }
    }
    if(!is.numeric(ny)){
        stop("Please specify the ny as a number or a numeric vector")
    }else{
        if(!(2 <= ny) || !(ny <= ncol(X))){
            stop("Please specify the value of ny within the range 2 <= ny <= ncol(X)")
        }
    }
}

.initGabrielNMF <- function(X, nx, ny){
    xholdouts <- list()
    yholdouts <- list()
    xdiv <- nrow(X)/nx
    ydiv <- ncol(X)/ny
    count <- 1
    for(x in seq_len(nx)){
        for(y in seq_len(ny)){
            xstart <- xdiv*(x-1) + 1
            xend <- xdiv*x
            ystart <- ydiv*(y-1) + 1
            yend <- ydiv*y
            xholdouts[[count]] <- xstart:xend
            yholdouts[[count]] <- ystart:yend
            count <- count + 1
        }
    }
    TestRecError <- rep(0, length=nx*ny)
    list(xholdouts=xholdouts, yholdouts=yholdouts, TestRecError=TestRecError)
}