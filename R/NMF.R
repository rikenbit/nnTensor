NMF <- function(X, M=NULL, pseudocount=1e-10,
    initU=NULL, initV=NULL, fixU=FALSE, fixV=FALSE,
    L1_U=1e-10, L1_V=1e-10, L2_U=1e-10, L2_V=1e-10, J = 3,
    rank.method=c("all", "ccc", "dispersion", "rss", "evar", "residuals", "sparseness.basis", "sparseness.coef", "sparseness2.basis",  "sparseness2.coef",  "norm.info.gain.basis",  "norm.info.gain.coef",  "singular",  "volume",  "condition"), runtime=30,
    algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman", "Alpha", "Beta", "PGD", "HALS", "GCD"), Alpha = 1, Beta = 2,
    eta = 1e-04, thr1 = 1e-10, thr2 = 1e-10, tol = 1e-04, num.iter = 100,
    viz = FALSE, figdir = NULL, verbose = FALSE){
    # Argument check
    rank.method <- match.arg(rank.method)
    algorithm <- match.arg(algorithm)
    chk <- .checkNMF(X, M, pseudocount, initU, initV, fixU, fixV,
        L1_U, L1_V, L2_U, L2_V, J, runtime,
        Alpha, Beta, eta, thr1, thr2, tol, num.iter, viz, figdir, verbose)
    X <- chk$X
    M <- chk$M
    pM <- chk$pM
    if(length(J) != 1){
        cat("Each rank, multiple NMF runs are performed\n")
        out1 <- .lapply_pb(J, function(j){
            # Original data
            out.original <- lapply(seq_len(runtime), function(r){
                try(.eachNMF(X, M, pM, initU, initV, fixU, fixV, L1_U, L1_V,
                    L2_U, L2_V, j, algorithm, Alpha, Beta, eta, thr1, thr2,
                    tol, num.iter, viz, figdir, verbose))
            })
            # Rand data
            out.random <- lapply(seq_len(runtime), function(r){
                try(.eachNMF(.randomize(X), M, pM, initU, initV, fixU, fixV,
                    L1_U, L1_V, L2_U, L2_V,
                    j, algorithm, Alpha, Beta, eta, thr1, thr2,
                    tol, num.iter, viz, figdir, verbose))
            })
            list(out.original=out.original, out.random=out.random)
        })
        names(out1) <- paste0("Rank", J)
        cat("Each rank estimation method\n")
        out2 <- .lapply_pb(J, function(j, other){
            out1 <- other$out1
            rank.method <- other$rank.method
            runtime <- other$runtime
            X <- other$X
            out.original <- eval(parse(text=paste0("out1$Rank", j, "$out.original")))
            out.random <- eval(parse(text=paste0("out1$Rank", j, "$out.random")))
            f <- .flist2[[rank.method]]
            f(out.original, out.random, X, runtime)
        }, other=list(out1=out1, rank.method=rank.method,
            runtime=runtime, X=X))
        names(out2) <- paste0("Rank", J)
        out <- list(U = NULL,
            V = NULL,
            J = J,
            RecError = NULL,
            TrainRecError = NULL,
            TestRecError = NULL,
            RelChange = NULL,
            Trial = out2, Runtime = runtime, RankMethod=rank.method)
        class(out) <- "NMF"
        out
    }else{
        out1 <- .eachNMF(X, M, pM, initU, initV, fixU, fixV,
            L1_U, L1_V, L2_U, L2_V,
            J, algorithm, Alpha, Beta, eta, thr1, thr2, tol, num.iter,
            viz, figdir, verbose)
        list(U = out1$U,
            V = out1$V,
            J = J,
            RecError = out1$RecError,
            TrainRecError = out1$TrainRecError,
            TestRecError = out1$TestRecError,
            RelChange = out1$RelChange,
            Trial = NULL, Runtime = NULL, RankMethod=NULL)
    }
}

.checkNMF <- function(X, M, pseudocount, initU, initV, fixU, fixV,
    L1_U, L1_V, L2_U, L2_V, J, runtime,
    Alpha, Beta, eta, thr1, thr2, tol, num.iter, viz, figdir, verbose){
    stopifnot(is.matrix(X))
    if(!is.null(M)){
        if(!identical(dim(X), dim(M))){
            stop("Please specify the dimensions of X and M are same")
        }
    }else{
        M <- X
        M[,] <- 1
    }
    stopifnot(is.numeric(pseudocount))
    if(!is.null(initU)){
        if(!identical(nrow(X), nrow(initU))){
            stop("Please specify nrow(X) and nrow(initU) are same")
        }
    }
    if(!is.null(initV)){
        if(!identical(ncol(X), nrow(initV))){
            stop("Please specify ncol(X) and nrow(initV) are same")
        }
    }
    stopifnot(is.logical(fixU))
    stopifnot(is.logical(fixV))
    if(L1_U < 0){
        stop("Please specify the L1_U that larger than 0")
    }
    if(L1_V < 0){
        stop("Please specify the L1_V that larger than 0")
    }
    if(L2_U < 0){
        stop("Please specify the L2_U that larger than 0")
    }
    if(L2_V < 0){
        stop("Please specify the L2_V that larger than 0")
    }
    stopifnot(is.numeric(J))
    stopifnot(is.numeric(runtime))
    stopifnot(is.numeric(Alpha))
    stopifnot(is.numeric(Beta))
    stopifnot(is.numeric(eta))
    stopifnot(is.numeric(thr1))
    stopifnot(is.numeric(thr2))
    stopifnot(is.numeric(tol))
    stopifnot(is.numeric(num.iter))
    stopifnot(is.logical(viz))
    if(!is.character(figdir) && !is.null(figdir)){
        stop("Please specify the figdir as a string or NULL")
    }
    stopifnot(is.logical(verbose))
    pM <- M
    X[which(X == 0)] <- pseudocount
    pM[which(pM == 0)] <- pseudocount
    list(X=X, M=M, pM=pM)
}

.eachNMF <- function(X, M, pM, initU, initV, fixU, fixV, L1_U, L1_V, L2_U, L2_V,
    J, algorithm, Alpha, Beta, eta, thr1, thr2,
    tol, num.iter, viz, figdir, verbose) {
    # Initialization of U, V
    int <- .initNMF(X, M, pM, initU, initV, J, thr1, Alpha, Beta, algorithm, verbose)
    X <- int$X
    M <- int$M
    pM <- int$pM
    U <- int$U
    V <- int$V
    RecError <- int$RecError
    TrainRecError <- int$TrainRecError
    TestRecError <- int$TestRecError
    RelChange <- int$RelChange
    Alpha <- int$Alpha
    Beta <- int$Beta
    algorithm <- int$algorithm
    iter <- 1
    while ((RecError[iter] > thr1) && (iter <= num.iter)) {
        # Update U, V
        X_bar <- .recMatrix(U, V)
        pre_Error <- .recError(X, X_bar)
        if (algorithm == "PGD") {
            if(!fixU){
                U <- .positive(U - eta * (-2 * (pM * X) %*% V + 2 * (pM * U %*% t(V)) %*% V), thr2)
            }
            if(!fixV){
                V <- .positive(V - eta * (-2 * t(pM * X) %*% U + 2 * (t(pM) * V %*% t(U)) %*% U), thr2)
            }
        }
        else if (algorithm == "Alpha") {
            if(!fixU){
                U <- U * ((((pM * X)/(pM * U %*% t(V)))^Alpha %*% V) /
                    (matrix(1, nrow = nrow(X), ncol = 1) %*%
                        colSums(V)))^(1/Alpha)
            }
            if(!fixV){
                V <- V * ((t((pM * X)/(pM * U %*% t(V)))^Alpha %*% U) /
                    (matrix(1, nrow = ncol(X), ncol = 1) %*%
                        colSums(U)))^(1/Alpha)
            }
        }
        else if (algorithm == "Beta") {
            if(!fixU){
                U <- U * (((pM * U %*% t(V))^(Beta - 2) * (pM * X)) %*% V) /
                ((pM * U %*% t(V))^(Beta - 1) %*% V + L1_U + L2_U * U)
            }
            if(!fixV){
                V <- V * (t((pM * U %*% t(V))^(Beta - 2) * (pM * X)) %*% U) /
                (t((pM * U %*% t(V))^(Beta - 1)) %*% U + L1_V + L2_V * V)
            }
        }
        else if (algorithm == "HALS") {
            A <- X %*% V
            B <- t(V) %*% V
            for (j in 1:J) {
                U[, j] <- .positive(A[, j] - U %*% B[, j] + U[,
                  j] * B[j, j], thr2)
                U[, j] <- U[, j]/norm(as.matrix(U[, j]), "F")
            }
            C <- t(X) %*% U
            D <- t(U) %*% U
            for (j in 1:J) {
                V[, j] <- .positive(C[, j] - V %*% D[, j] + V[,
                  j] * D[j, j], thr2)
            }
        }
        else if (algorithm == "GCD") {
            if(!fixU){
                Unew <- .doiter(U, V, X, tol = tol, J)
                U <- U + Unew                
            }
            if(!fixV){
                Vnew <- .doiter(V, U, t(X), tol = tol, J)
                V <- V + Vnew                
            }
        }
        else {
            stop("Please specify the appropriate algorithm\n")
        }
        # After Update U, V
        iter <- iter + 1
        X_bar <- .recMatrix(U, V)
        RecError[iter] <- .recError(X, X_bar)
        TrainRecError[iter] <- .recError(M*X, M*X_bar)
        TestRecError[iter] <- .recError((1-M)*X, (1-M)*X_bar)
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if (viz && !is.null(figdir)) {
            png(filename = paste0(figdir, "/", iter-1, ".png"))
            image.plot(X_bar)
            dev.off()
        }
        if (viz && is.null(figdir)) {
            image.plot(X_bar)
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
        image.plot(X_bar)
        dev.off()
        png(filename = paste0(figdir, "/original.png"))
        image.plot(X)
        dev.off()
    }
    if (viz && is.null(figdir)) {
        image.plot(X_bar)
    }
    names(RecError) <- c("offset", 1:(iter-1))
    names(TrainRecError) <- c("offset", 1:(iter-1))
    names(TestRecError) <- c("offset", 1:(iter-1))
    names(RelChange) <- c("offset", 1:(iter-1))

    list(U = U, V = V, RecError = RecError,
        TrainRecError = TrainRecError,
        TestRecError = TestRecError,
        RelChange = RelChange)
}

.initNMF <- function(X, M, pM, initU, initV, J, thr1, Alpha, Beta, algorithm, verbose){
    if(is.null(initU)){
        U <- matrix(runif(nrow(X) * J), nrow = nrow(X), ncol = J)
        U <- U * U
    }else{
        U <- initU
    }
    if(is.null(initV)){
        V <- matrix(runif(ncol(X) * J), nrow = ncol(X), ncol = J)
        V <- V * V
    }else{
        V <- initV
    }
    RecError = c()
    TrainRecError = c()
    TestRecError = c()
    RelChange = c()
    RecError[1] <- thr1 * 10
    TrainRecError[1] <- thr1 * 10
    TestRecError[1] <- thr1 * 10
    RelChange[1] <- thr1 * 10
    R <- NULL
    if (algorithm == "HALS") {
        R <- X - U %*% t(V)
    }
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
    if (algorithm == "Pearson") {
        Alpha = 2
        algorithm = "Alpha"
    }
    if (algorithm == "Hellinger") {
        Alpha = 0.5
        algorithm = "Alpha"
    }
    if (algorithm == "Neyman") {
        Alpha = -1
        algorithm = "Alpha"
    }
    if (verbose) {
        cat("Iterative step is running...\n")
    }
    list(X=X, M=M, pM=pM, U=U, V=V, RecError=RecError, TrainRecError=TrainRecError,
        TestRecError=TestRecError, RelChange=RelChange,
        Alpha=Alpha, Beta=Beta, algorithm=algorithm)
}


plot <- function(x, ...){
    UseMethod("plot", x, ...)
}

plot.NMF <- function(x, ...){
    rank.method <- x$RankMethod
    if(is.null(x$U)){
        if(rank.method == "all"){
            .plot.all(x)
        }else{
            .plot.each(x)
        }
    }else{
        NULL
    }
}

.out2df <- function(object, type){
    J <- object$J
    Runtime <- object$Runtime
    # for .plot.each
    if(type == "each"){
        df <- data.frame(
                Rank=.seq2(J, Runtime),
                .eachval(object, J, Runtime))
    }
    # for .plot.all
    if(type == "all"){
        df <- data.frame(
                Rank=.seq2(J, Runtime),
                .eachval2(object, .all.method, J, Runtime)
            )
    }
    df
}

.eachval <- function(object, J, Runtime){
    oval <- unlist(lapply(J, function(j){
                        eval(parse(text=paste0("object$Trial$Rank", j, "$original")))
                    }))
    rval <- unlist(lapply(J, function(j){
                        eval(parse(text=paste0("object$Trial$Rank", j, "$random")))
                    }))
    if(length(oval) == length(J)){
        oval <- .seq2(oval, Runtime)
        rval <- .seq2(rval, Runtime)
    }
    data.frame(
        Value=c(oval, rval),
        Type=.seq2(c("original", "random"),
            length(J)*Runtime))
}

.eachval2 <- function(object, Method, J, Runtime){
    oval <- unlist(lapply(Method, function(m){
        val <- unlist(lapply(J, function(j){
                            eval(parse(text=paste0("object$Trial$Rank", j, "$", m, "$original")))}))
        if(length(val) == length(J)){
            val <- .seq2(val, Runtime)
        }
        val
    }))
    rval <- unlist(lapply(Method, function(m){
        val <- unlist(lapply(J, function(j){
                            eval(parse(text=paste0("object$Trial$Rank", j, "$", m, "$random")))}))
        if(length(val) == length(J)){
            val <- .seq2(val, Runtime)
        }
        val
    }))
    data.frame(Value=c(oval, rval),
        Type=.seq2(c("original", "random"),
            length(Method)*length(J)*Runtime),
        Method=rep(.seq2(Method, length(J)*Runtime), 2))
}

.all.method <- c("ccc", "dispersion", "rss", "evar", "residuals", "sparseness.basis", "sparseness.coef", "sparseness2.basis", "sparseness2.coef", "norm.info.gain.basis", "norm.info.gain.coef", "singular", "volume", "condition")

.seq2 <- function(J, Runtime){
    unlist(lapply(J, function(j, r){
        rep(j, r)
    }, r=Runtime))
}

# 1. Line plot
.plot.each <- function(object){
    BestRank <- object$BestRank
    Rank <- NULL
    Value <- NULL
    Type <- NULL
    df <- .out2df(object, type="each")
    df <- df[which(!is.na(df$Value)), ]
    g <- ggplot(df, aes(Rank, Value, group=Type , colour=Type))
    g <- g + geom_point()
    g <- g + stat_summary(fun.y = "mean", geom="line", size = 2)  
    g <- g + stat_summary(fun.data = "mean_se", geom="errorbar")
    g <- g + geom_vline(xintercept=BestRank, linetype="dashed")
    g
}

# 2. Panel plot <2*7>
.plot.all <- function(object){
    BestRank <- object$BestRank
    Rank <- NULL
    Value <- NULL
    Type <- NULL
    df <- .out2df(object, type="all")
    df <- df[which(!is.na(df$Value)), ]
    g <- ggplot(df, aes(Rank, Value, group=Type , colour=Type))
    g <- g + geom_point()
    g <- g + stat_summary(fun.y = "mean", geom="line", size = 0.5)
    g <- g + stat_summary(fun.data = "mean_se", geom="errorbar")
    g <- g + geom_vline(xintercept=BestRank, linetype="dashed")
    g <- g + facet_wrap(~Method, scales="free", nrow=2)
    g
}
