GabrielNMF <- function(X, J = 3, nx = 5, ny = 5,
 algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
    "Alpha", "Beta", "PGD", "HALS", "GCD"), Alpha = 1, Beta = 2,
    eta = 1e-04, thr1 = 1e-10, thr2 = 1e-10, tol = 1e-04,
    num.iter = 100, verbose = FALSE){
    # Argument check
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
    algorithm <- match.arg(algorithm)
    if(!is.numeric(Alpha)){
        stop("Please specify the Alpha as a numeric")
    }
    if(!is.numeric(Beta)){
        stop("Please specify the Beta as a numeric")
    }
    if(!is.numeric(eta)){
        stop("Please specify the eta as a numeric")
    }
    if(!is.numeric(thr1)){
        stop("Please specify the thr1 as a numeric")
    }
    if(!is.numeric(thr2)){
        stop("Please specify the thr2 as a numeric")
    }
    if(!is.numeric(tol)){
        stop("Please specify the tol as a numeric")
    }
    if(!is.numeric(num.iter)){
        stop("Please specify the num.iter as a numeric")
    }
    if(!is.logical(verbose)){
        stop("Please specify the verbose as a logical")
    }

    # Bi-Cross-Validation
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
    tcount <- 1
    for(x in seq(nx*ny)){
        M <- X
        M[, ] <- 1
        M[xholdouts[[x]], yholdouts[[x]]] <- 0
        out <- try(.eachGabrielNMF(X, M, J, algorithm, Alpha, Beta,
            eta, thr1, thr2, tol, num.iter, verbose))
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

.eachGabrielNMF <- function(X, M, J, algorithm, Alpha, Beta, eta, thr1, thr2,
    tol, num.iter, verbose) {
    # Submatrices
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
    nmfD <- NMF(D, J=J, algorithm=algorithm, Alpha=Alpha, Beta=Beta, eta=eta,
        thr1=thr1, thr2=thr2, tol=tol, num.iter=num.iter, verbose=verbose)
    # Prediction of U and V
    U_B <- NMF(B, initV=nmfD$V, fixV=TRUE,
        J=J, algorithm=algorithm, Alpha=Alpha, Beta=Beta, eta=eta,
        thr1=thr1, thr2=thr2, tol=tol, num.iter=num.iter, verbose=verbose)$U
    V_C <- NMF(C, initU=nmfD$U, fixU=TRUE,
        J=J, algorithm=algorithm, Alpha=Alpha, Beta=Beta, eta=eta,
        thr1=thr1, thr2=thr2, tol=tol, num.iter=num.iter, verbose=verbose)$V

    # Bi Cross Validation
    norm(A - U_B %*% t(V_C), "F")^2
}
