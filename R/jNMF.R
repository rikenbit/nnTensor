jNMF <- function(X, M=NULL, pseudocount=.Machine$double.eps,
    initW=NULL, initV=NULL, initH=NULL,
    fixW=FALSE, fixV=FALSE, fixH=FALSE,
    L1_W=1e-10, L1_V=1e-10, L1_H=1e-10,
    L2_W=1e-10, L2_V=1e-10, L2_H=1e-10,
    J = 3, w=NULL, algorithm = c("Frobenius", "KL", "IS", "PLTF"), p=1,
    thr = 1e-10, num.iter = 100,
    viz = FALSE, figdir = NULL, verbose = FALSE) {
    # Argument check
    algorithm <- match.arg(algorithm)
    chk <- .checkjNMF(X, M, pseudocount,
        initW, initV, initH, fixW, fixV, fixH, J, w,
        p, thr, num.iter, viz, figdir, verbose)
    X <- chk$X
    M <- chk$M
    pM <- chk$pM
    fixV <- chk$fixV
    fixH <- chk$fixH
    w <- chk$w
    K <- chk$K
    # Initialization of W, V, H
    int <- .initjNMF(X, initW, initV, initH, J, p, thr, algorithm, verbose)
    W <- int$W
    V <- int$V
    H <- int$H
    RecError <- int$RecError
    TrainRecError <- int$TrainRecError
    TestRecError <- int$TestRecError
    RelChange <- int$RelChange
    p <- int$p
    iter <- 1
    while ((RecError[iter] > thr) && (iter <= num.iter)) {
        # Before Update W, H_k
        X_bar <- lapply(seq_len(K), function(k){
            .recMatrix(W + V[[k]], H[[k]])
        })
        pre_Error <- sqrt(sum(unlist(lapply(seq_along(X), function(x){
            .recError(X[[x]], X_bar[[x]], notsqrt=TRUE)
        }))))
        # Update W
        if(!fixW){
            W_numer <- matrix(0, nrow=nrow(W), ncol=ncol(W))
            W_denom <- matrix(0, nrow=nrow(W), ncol=ncol(W))
            for(k in seq_len(K)){
                W_numer <- W_numer + w[k] * (pM[[k]] * X[[k]] * (pM[[k]] * X_bar[[k]])^(-p)) %*% H[[k]]
                W_denom <- W_denom + w[k] * (pM[[k]] * (W + V[[k]]) %*% t(H[[k]]))^(1-p) %*% H[[k]] + L1_W + L2_W * W
            }
            W <- .columnNorm(.positive(W * W_numer / W_denom)^.rho(p))
        }
        # Update H_k
        for(k in seq_len(K)){
            if(!fixH[k]){
                Hk_numer <- (t(pM[[k]] * X[[k]]) * t(pM[[k]] * X_bar[[k]])^(-p)) %*% (W + V[[k]])
                Hk_denom <- t(pM[[k]] * (W + V[[k]]) %*% t(H[[k]]))^(1-p) %*% W + L1_H + L2_H * H[[k]]
                H[[k]] <- H[[k]] * (Hk_numer / Hk_denom)^.rho(p)
            }
        }
        # Update V_k
        for(k in seq_len(K)){
            if(!fixV[k]){
                Vk_numer <- (pM[[k]] * X[[k]] * (pM[[k]] * X_bar[[k]])^(-p)) %*% H[[k]]
                Vk_denom <- (pM[[k]] * (W + V[[k]]) %*% t(H[[k]]))^(1-p) %*% H[[k]] + L1_V + L2_V * V[[k]]
                V[[k]] <- V[[k]] * (Vk_numer / Vk_denom)^.rho(p)
            }
        }
        # After Update W, H_k
        iter <- iter + 1
        X_bar <- lapply(seq_len(K), function(k){
            .recMatrix(W + V[[k]], H[[k]])
        })
        RecError[iter] <- sqrt(sum(unlist(lapply(seq_along(X), function(x){
            .recError(X[[x]], X_bar[[x]], notsqrt=TRUE)
        }))))
        TrainRecError[iter] <- sqrt(sum(unlist(lapply(seq_along(X), function(x){
            .recError(M[[x]]*X[[x]], M[[x]]*X_bar[[x]], notsqrt=TRUE)
        }))))
        TestRecError[iter] <- sqrt(sum(unlist(lapply(seq_along(X), function(x){
            .recError((1-M[[x]])*X[[x]], (1-M[[x]])*X_bar[[x]], notsqrt=TRUE)
        }))))
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if (viz && !is.null(figdir)) {
            png(filename = paste0(figdir, "/", iter-1, ".png"))
            .imageplot_jNMF(X, W, V, H)
            dev.off()
        }
        if (viz && is.null(figdir)) {
            .imageplot_jNMF(X, W, V, H)
        }
        if (verbose) {
            cat(paste0(iter-1, " / ", num.iter, " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
        if (is.nan(RelChange[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }
    if (viz && !is.null(figdir)) {
        png(filename = paste0(figdir, "/finish.png"))
        .imageplot_jNMF(X, W, V, H)
        dev.off()
    }
    if (viz && is.null(figdir)) {
        .imageplot_jNMF(X, W, V, H)
    }
    names(RecError) <- c("offset", seq_len(iter-1))
    names(TrainRecError) <- c("offset", seq_len(iter-1))
    names(TestRecError) <- c("offset", seq_len(iter-1))
    names(RelChange) <- c("offset", seq_len(iter-1))
    return(list(W = W, V = V, H = H,
        RecError = RecError,
        TrainRecError = TrainRecError,
        TestRecError = TestRecError,
        RelChange = RelChange))
}

.checkjNMF <- function(X, M, pseudocount, initW, initV, initH, fixW, fixV, fixH, J, w,
    p, thr, num.iter, viz, figdir, verbose){
    stopifnot(is.list(X))
    if(length(X) < 2){
        stop("input list X must have at least two datasets!")
    }
    if(!is.null(M)){
        dimX <- as.vector(unlist(lapply(X, function(x){dim(x)})))
        dimM <- as.vector(unlist(lapply(M, function(x){dim(x)})))
        if(!identical(dimX, dimM)){
            stop("Please specify the dimensions of X and M are same")
        }
    }else{
        M <- X
        for(i in seq(length(X))){
            M[[i]][] <- 1
        }
    }
    stopifnot(is.numeric(pseudocount))
    if(!is.null(initW)){
        if(!identical(nrow(X[[1]]), nrow(initW))){
            stop("Please specify nrow(X[[k]]) and nrow(W) are same")
        }
    }
    if(!is.null(initV)){
        nrowX <- as.vector(unlist(lapply(X, nrow)))
        nrowV <- as.vector(unlist(lapply(initV, nrow)))
        if(!identical(nrowX, nrowV)){
            stop("Please specify nrow(X[[k]]) and nrow(W) are same")
        }
    }
    if(!is.null(initH)){
        ncolX <- as.vector(unlist(lapply(X, ncol)))
        nrowH <- as.vector(unlist(lapply(initH, nrow)))
        if(!identical(ncolX, nrowH)){
            stop("Please specify all the ncol(initH[[k]]) are same as ncol(X[[k]]) (k=1,2,...)")
        }
    }
    stopifnot(is.logical(fixW))
    if(!is.logical(fixV)){
        if(!is.vector(fixV)){
            stop("Please specify the fixV as a logical or a logical vector such as c(TRUE, FALSE, TRUE)")
        }else{
            if(length(X) != length(fixV)){
                stop("Please specify the length of fixV same as the length of X")
            }
        }
    }else{
        fixV <- rep(fixV, length=length(X))
    }
    if(!is.logical(fixH)){
        if(!is.vector(fixH)){
            stop("Please specify the fixH as a logical or a logical vector such as c(TRUE, FALSE, TRUE)")
        }else{
            if(length(X) != length(fixH)){
                stop("Please specify the length of fixH same as the length of X")
            }
        }
    }else{
        fixH <- rep(fixH, length=length(X))
    }
    stopifnot(is.numeric(J))
    if(is.null(w)){
        w <- rep(1, length=length(X))
    }else{
        if(length(X) != length(w)){
            stop("The length of weight vector must be same as that of input list X!")
        }else{
            w <- w / sum(w)
        }
    }
    stopifnot(is.numeric(p))
    stopifnot(is.numeric(thr))
    stopifnot(is.numeric(num.iter))
    stopifnot(is.logical(viz))
    stopifnot(is.logical(verbose))
    if(!is.character(figdir) && !is.null(figdir)){
        stop("Please specify the figdir as a string or NULL")
    }
    pM <- M
    lapply(seq_along(X), function(x){
         X[[x]][which(X[[x]] == 0)] <<- pseudocount
         pM[[x]][which(pM[[x]] == 0)] <<- pseudocount
     })
    K <- length(X)
    list(X=X, M=M, pM=pM, fixV=fixV, fixH=fixH, w=w, K=K)
}

.initjNMF <- function(X, initW, initV, initH, J, p, thr, algorithm, verbose){
    if(is.null(initW)){
        W <- matrix(runif(nrow(X[[1]])*J),
            nrow=nrow(X[[1]]), ncol=J)
        W <- .columnNorm(W * W)
    }else{
        W <- initW
    }
    if(is.null(initV)){
        V <- lapply(seq_along(X), function(x){
            tmp <- matrix(runif(nrow(X[[1]])*J),
            nrow=nrow(X[[1]]), ncol=J)
            .columnNorm(tmp * tmp)
        })
    }else{
        V <- initV
    }
    if(is.null(initH)){
        H <- lapply(seq_along(X), function(x){
            tmp <- matrix(runif(ncol(X[[x]])*J),
                nrow=ncol(X[[x]]), ncol=J)
            .columnNorm(tmp * tmp)
        })
    }else{
        H <- initH
    }
    RecError = c()
    TrainRecError = c()
    TestRecError = c()
    RelChange = c()
    RecError[1] <- thr * 10
    TrainRecError[1] <- thr * 10
    TestRecError[1] <- thr * 10
    RelChange[1] <- thr * 10
    # Algorithm
    if (algorithm == "Frobenius") {
        p = 0
    }
    if (algorithm == "KL") {
        p = 1
    }
    if (algorithm == "IS") {
        p = 2
    }
    if (algorithm == "PLTF") {
        p = p
    }
    if (verbose) {
        cat("Iterative step is running...\n")
    }
    list(W=W, V=V, H=H, RecError=RecError, TrainRecError=TrainRecError, TestRecError=TestRecError, RelChange=RelChange, p=p)
}

.imageplot_jNMF <- function(X, W, V, H){
    l <- length(X)
    layout(rbind(seq_len(l), seq_len(l)+l, seq_len(l)+2*l))
    lapply(seq_along(X), function(x){
        image.plot(t(X[[x]]))
    })
    lapply(seq_along(X), function(x){
        image.plot(t(W %*% t(H[[x]])))
    })
    lapply(seq_along(X), function(x){
        image.plot(t(V[[x]] %*% t(H[[x]])))
    })
}
