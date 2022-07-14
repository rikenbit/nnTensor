NMTF <- function(X, pseudocount=.Machine$double.eps,
    initU=NULL, initS=NULL, initV=NULL,
    fixU=FALSE, fixS=FALSE, fixV=FALSE,
    L1_U=1e-10, L1_S=1e-10, L1_V=1e-10,
    L2_U=1e-10, L2_S=1e-10, L2_V=1e-10,
    rank = c(3, 4),
    algorithm = c("Frobenius", "KL", "IS", "ALS", "PG", "COD", "Beta"),
    Beta = 2, thr = 1e-10, num.iter = 100,
    viz = FALSE, figdir = NULL, verbose = FALSE){
    # Argument check
    algorithm <- match.arg(algorithm)
    .checkNMTF(X, pseudocount, initU, initS, initV,
    fixU, fixS, fixV, L1_U, L1_S, L1_V, L2_U, L2_S, L2_V,
    rank, Beta, thr, num.iter, viz, figdir, verbose)
    # Initizalization
    int <- .initNMTF(X, pseudocount, rank, initU, initS, initV,
    algorithm, Beta, thr, verbose)
    X <- int$X
    U <- int$U
    S <- int$S
    V <- int$V
    RecError <- int$RecError
    RelChange <- int$RelChange
    Beta <- int$Beta
    algorithm <- int$algorithm
    # Iteration
    iter <- 1
    while ((RecError[iter] > thr) && (iter <= num.iter)) {
        # Reconstruction
        X_bar <- U %*% S %*% t(V)
        pre_Error <- .recError(X, X_bar)
        # Update U
        if(!fixU){
            U <- .updateU(X, U, S, V, L1_U, L2_U, Beta, algorithm)
        }
        # Update V
        if(!fixV){
            V <- .updateV(X, U, S, V, L1_V, L2_V, Beta, algorithm)
        }
        # Update S
        if(!fixS){
            S <- .updateS(X, U, S, V, L1_S, L2_S, Beta, algorithm)
        }
        # After Update U, S, V
        iter <- iter + 1
        X_bar <- U %*% S %*% t(V)
        RecError[iter] <- .recError(X, X_bar)
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if (viz && !is.null(figdir)) {
            png(filename = paste0(figdir, "/", iter-1, ".png"))
            .multiImagePlots(list(X, X_bar, U, S, t(V)))
            dev.off()
        }
        if (viz && is.null(figdir)) {
            .multiImagePlots(list(X, X_bar, U, S, t(V)))
        }
        if (verbose) {
            cat(paste0(iter-1, " / ", num.iter, " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
        if (is.nan(RelChange[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }
    # After iteration
    if (viz && !is.null(figdir)) {
        png(filename = paste0(figdir, "/finish.png"))
        image.plot2(X_bar)
        dev.off()
        png(filename = paste0(figdir, "/original.png"))
        image.plot2(X)
        dev.off()
    }
    if (viz && is.null(figdir)) {
            .multiImagePlots(list(X, X_bar, U, S, t(V)))
    }
    names(RecError) <- c("offset", 1:(iter-1))
    names(RelChange) <- c("offset", 1:(iter-1))
    # Output
    list(U = U, S = S, V = V, rank = rank,
        RecError = RecError, RelChange = RelChange)
}

.checkNMTF <- function(X, pseudocount, initU, initS, initV,
    fixU, fixS, fixV, L1_U, L1_S, L1_V, L2_U, L2_S, L2_V,
    rank, Beta, thr, num.iter, viz, figdir, verbose){
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(pseudocount))
    if(!is.null(initU)){
        if(!identical(nrow(X), nrow(initU))){
            stop("Please specify nrow(X) and nrow(initU) are same")
        }
    }
    if(!is.null(initS)){
        if(rank[1] != nrow(initS)){
            stop("Please specify rank[1] and nrow(initS) are same")
        }
        if(rank[2] != nrow(initS)){
            stop("Please specify rank[2] and ncol(initS) are same")
        }
    }
    if(!is.null(initV)){
        if(!identical(ncol(X), nrow(initV))){
            stop("Please specify ncol(X) and nrow(initV) are same")
        }
    }
    stopifnot(is.logical(fixU))
    stopifnot(is.logical(fixS))
    stopifnot(is.logical(fixV))
    if(L1_U < 0){
        stop("Please specify the L1_U that larger than 0")
    }
    if(L1_S < 0){
        stop("Please specify the L1_S that larger than 0")
    }
    if(L1_V < 0){
        stop("Please specify the L1_V that larger than 0")
    }
    if(L2_U < 0){
        stop("Please specify the L2_U that larger than 0")
    }
    if(L2_S < 0){
        stop("Please specify the L2_S that larger than 0")
    }
    if(L2_V < 0){
        stop("Please specify the L2_V that larger than 0")
    }
    stopifnot(is.numeric(rank))
    stopifnot(is.numeric(Beta))
    stopifnot(is.numeric(thr))
    stopifnot(is.numeric(num.iter))
    stopifnot(is.logical(viz))
    if(!is.character(figdir) && !is.null(figdir)){
        stop("Please specify the figdir as a string or NULL")
    }
    stopifnot(is.logical(verbose))
}

.initNMTF <- function(X, pseudocount, rank, initU, initS, initV,
    algorithm, Beta, thr, verbose){
    X[which(X == 0)] <- pseudocount
    if(is.null(initU)){
        U <- matrix(runif(nrow(X)*rank[1]), nrow=nrow(X), ncol=rank[1])
    }else{
        U <- initU
    }
    if(is.null(initV)){
        V <- matrix(runif(ncol(X)*rank[2]), nrow=ncol(X), ncol=rank[2])
    }else{
        V <- initV
    }
    if(is.null(initS)){
        S <- matrix(runif(prod(rank)), nrow=rank[1], ncol=rank[2])
    }else{
        S <- initS
    }
    RecError = c()
    RelChange = c()
    RecError[1] <- thr * 10
    RelChange[1] <- thr * 10
    # Algorithm
    if (algorithm == "Frobenius") {
        Beta = 2
        algorithm = "Beta"
    }
    if (algorithm == "KL") {
        Beta = 1
        algorithm = "Beta"
    }
    if (algorithm == "IS") {
        Beta = 0
        algorithm = "Beta"
    }
    if (verbose) {
        cat("Iterative step is running...\n")
    }
    list(X=X, U=U, S=S, V=V,
        RecError=RecError, RelChange=RelChange,
        Beta=Beta, algorithm=algorithm)
}

.updateU <- function(X, U, S, V, L1_U, L2_U, Beta, algorithm){
    if(algorithm == "Beta"){
        VS <- V %*% t(S)
        numer <- ((U %*% S %*% t(V))^(Beta-2) * X) %*% VS
        denom <- (U %*% S %*% t(V))^(Beta-1) %*% VS
        denom <- denom + L1_U + L2_U * U
        U <- U * (numer / denom)^.rho(Beta)
    }
    if(algorithm == "ALS"){
        VS <- V %*% t(S)
        U <- .positive((X %*% V %*% t(S)) %*% ginv(crossprod(VS)))
    }
    if(algorithm == "PG"){
        Pu <- U - U / ((U %*% S %*% t(V) %*% V %*% t(S)) * (X %*% V %*% t(S)))
        numer <- sum(Pu * (U %*% S %*% t(V) %*% V %*% t(S) - X %*% V %*% t(S)))
        denom <- .trace((S %*% t(V) %*% V) %*% (t(S) %*% t(Pu) %*% Pu))
        etaU <-  numer / denom
        U <- .positive(U - etaU * Pu)
    }
    if(algorithm == "COD"){
        for(i in seq_len(ncol(U))){
            numer <- (X %*% V %*% t(S))[,i] - (U %*% S %*% t(V) %*% V %*% t(S))[,i]
            denom <- as.numeric(crossprod(V %*% S[i, ]))
            U[,i] <- U[,i] + .positive(numer / denom)
        }
    }
    U
}

.updateV <- function(X, U, S, V, L1_V, L2_V, Beta, algorithm){
    if(algorithm == "Beta"){
        SU <- t(S) %*% t(U)
        numer <- SU %*% ((U %*% S %*% t(V))^(Beta - 2) * X)
        denom <- SU %*% (U %*% S %*% t(V))^(Beta - 1)
        denom <- denom + L1_V + L2_V * t(V)
        V <- V * t((numer / denom)^.rho(Beta))
    }
    if(algorithm == "ALS"){
        US <- U %*% S
        V <- .positive((t(X) %*% U %*% S) %*% ginv(crossprod(US)))
    }
    if(algorithm == "PG"){
        Pv <- V - V / ((V %*% t(S) %*% t(U) %*% U %*% S) * (t(X) %*% U %*% S))
        numer <- sum(Pv * (V %*% t(S) %*% t(U) %*% U %*% S - t(X) %*% U %*% S))
        denom <- .trace((S %*% t(Pv) %*% Pv) %*% (t(S) %*% t(U) %*% U))
        etaV <-  numer / denom
        V <- .positive(V - etaV * Pv)
    }
    if(algorithm == "COD"){
        for(j in seq_len(ncol(V))){
            numer <- (t(X) %*% U %*% S)[,j] - (V %*% t(S) %*% t(U) %*% U %*% S)[,j]
            denom <- as.numeric(crossprod(U %*% S[,j]))
            V[,j] <- V[,j] + .positive(numer / denom)
        }
    }
    V
}

.updateS <- function(X, U, S, V, L1_S, L2_S, Beta, algorithm){
    if(algorithm == "Beta"){
        numer <- t(U) %*% ((U %*% S %*% t(V))^(Beta-2) * X) %*% V
        denom <- t(U) %*% (U %*% S %*% t(V))^(Beta-1) %*% V
        denom <- denom + L1_S + L2_S * S
        S <- S * (numer / denom)^.rho(Beta)
    }
    if(algorithm == "ALS"){
        S <- .positive(ginv(t(U) %*% U) %*% (t(U) %*% X %*% V) %*% ginv(t(V) %*% V))
    }
    if(algorithm == "PG"){
        Ps <- S - S / ((t(U) %*% U %*% S %*% t(V) %*% V) * (t(U) %*% X %*% V))
        numer <- sum(Ps * (t(U) %*% U %*% S %*% t(V) %*% V - t(U) %*% X %*% V))
        denom <- .trace((t(U) %*% U %*% Ps) %*% (t(V) %*% V %*% t(Ps)))
        etaS <-  numer / denom
        S <- .positive(S - etaS * Ps)
    }
    if(algorithm == "COD"){
        for(i in seq_len(ncol(U))){
            for(j in seq_len(ncol(V))){
                numer <- (t(U) %*% X %*% V)[i,j] - (t(U) %*% U %*% S %*% t(V) %*% V)[i,j]
                denom <- as.numeric((U[,i] %*% U[,i]) * (V[,j] %*% V[,j]))
                S[i,j] <- S[i,j] + .positive(numer / denom)
            }
        }
    }
    S
}

.trace <- function(mat){
    sum(diag(mat))
}

.multiImagePlots <- function(inputList){
    layout(rbind(1:3, 4:6))
    image.plot2(inputList[[1]], main="X")
    image.plot2(inputList[[2]], main="rec X")
    plot.new()
    image.plot2(inputList[[3]], main="U")
    image.plot2(inputList[[4]], main="S")
    image.plot2(inputList[[5]], main="t(V)")
}

image.plot2 <- function(A, ...){
    image.plot(t(A[nrow(A):1,]), ...)
}
