NTF <- function(X, M=NULL, pseudocount=.Machine$double.eps,
    initA=NULL, fixA=FALSE,
    L1_A=1e-10, L2_A=1e-10, rank = 3,
    algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman", "HALS", "Alpha-HALS", "Beta-HALS", "Alpha", "Beta"),
    init = c("NMF", "ABS-SVD", "ALS", "Random"), Alpha = 1,
    Beta = 2, thr = 1e-10, num.iter = 100, viz = FALSE, figdir = NULL,
    verbose = FALSE){
    # Argument check
    algorithm <- match.arg(algorithm)
    init <- match.arg(init)
    .checkNTF(X, M, pseudocount, initA, fixA, rank, Alpha, Beta,
        thr, num.iter, viz, figdir, verbose)
    # Initialization of An
    int <- .initNTF(X, M, pseudocount, fixA, rank, init, initA, Alpha, Beta,
        algorithm, thr, verbose)
    X <- int$X
    M <- int$M
    pM <- int$pM
    M_NA <- int$M_NA
    fixA <- int$fixA
    N <- int$N
    A <- int$A
    RecError <- int$RecError
    TrainRecError <- int$TrainRecError
    TestRecError <- int$TestRecError
    RelChange <- int$RelChange
    Alpha <- int$Alpha
    Beta <- int$Beta
    algorithm <- int$algorithm
    E <- int$E
    T1 <- int$T1
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        # Before Update An
        X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
        pre_Error <- .recError(X, X_bar)
        # Fill with Machine Epsilon
        for (n in seq(N)) {
            A[[n]][which(A[[n]] < pseudocount)] <- pseudocount
            A[[n]][which(is.infinite(A[[n]]))] <- pseudocount
            A[[n]][which(is.nan(A[[n]]))] <- pseudocount
            A[[n]][which(is.nan(A[[n]]))] <- pseudocount
        }
        # Update An
        if (algorithm == "Alpha") {
            for (n in seq(N)) {
                if(!fixA[n]){
                    X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
                    A_notn <- .KhatriRao_notn(A, n)
                    A[[n]] <- A[[n]] * (t(A_notn) %*% (cs_unfold(pM*X,
                      m = n)@data/cs_unfold(pM*X_bar, m = n)@data)^Alpha)^(1/Alpha)
                    if (n != N) {
                      A[[n]] <- t(apply(A[[n]], 1, function(x) {
                        x/norm(as.matrix(x), "F")
                      }))
                    }
                }
            }
        }
        else if (algorithm == "Beta") {
            X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
            for (n in seq(N)) {
                if(!fixA[n]){
                    A_notn <- .KhatriRao_notn(A, n)
                    numer <- t(A_notn) %*% (cs_unfold(pM*X, m = n)@data/cs_unfold(pM*X_bar^(Beta - 1), m = n)@data)
                    denom <- t(A_notn) %*% cs_unfold(X_bar^(Beta),
                      m = n)@data
                    A[[n]] <- A[[n]] * (numer / (denom + L1_A + L2_A * A[[n]]))^.rho(Beta)
                    if (n != N) {
                      A[[n]] <- t(apply(A[[n]], 1, function(x) {
                        x/norm(as.matrix(x), "F")
                      }))
                    }
                }
            }
        }
        else if (algorithm == "HALS") {
            gamma <- diag((A[[N]] %*% t(A[[N]])))
            for (n in seq(N)) {
                if (n == N) {
                  gamma[] <- 1
                }
                A_notn <- .KhatriRao_notn(A, n)
                T2 = t(cs_unfold(X, m = n)@data) %*% A_notn
                T3 = T1/(A[[n]] %*% t(A[[n]]))
                for (r in seq(rank)) {
                  A[[n]][r, ] = .positive(gamma[r] * A[[n]][r,
                    ] + T2[, r] - as.vector(t(A[[n]]) %*% T3[,
                    r]))^2
                  if (n != N) {
                    A[[n]][r, ] <- A[[n]][r, ]/norm(as.matrix(A[[n]][r,
                      ]), "F")
                  }
                }
                T1 <- T3 * (A[[n]] %*% t(A[[n]]))
            }
        }
        else if (algorithm == "Alpha-HALS") {
            for (r in seq(rank)) {
                X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
                Xr <- E + X_bar
                for (n in seq(N)) {
                  if (Alpha == 0) {
                    tmp_Xr <- .positive(Xr)
                    tmp_Xr@data <- log(tmp_Xr@data)
                    numer <- tmp_Xr
                  }
                  else {
                    numer <- .positive(Xr)^Alpha
                  }
                  denom <- 0
                  not_n <- setdiff(seq(N), n)
                  for (m in not_n) {
                    numer <- ttm(numer, t(as.matrix(A[[m]][r,
                      ])), m = m)
                    if (Alpha == 0) {
                      denom <- denom + as.numeric(t(A[[m]][r,
                        ]) %*% log(A[[m]][r, ]))
                    }
                    else {
                      denom <- denom + as.numeric(t(A[[m]][r,
                        ]) %*% A[[m]][r, ]^Alpha)
                    }
                  }
                  numer <- as.vector(numer@data[])
                  if (Alpha == 0) {
                    A[[n]][r, ] <- exp(.positive(numer/denom))
                  }
                  else {
                    A[[n]][r, ] <- .positive(numer/denom)^(1/Alpha)
                  }
                  if (n != N) {
                    A[[n]][r, ] <- A[[n]][r, ]/norm(as.matrix(A[[n]][r,
                      ]), "F")
                  }
                }
                X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
                E <- Xr - X_bar
            }
        }
        else if (algorithm == "Beta-HALS") {
            for (r in seq(rank)) {
                X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
                Xr <- E + X_bar
                for (n in 1:(N - 1)) {
                  not_n <- setdiff(seq(N), n)
                  tmp_u <- Xr
                  for (m in not_n) {
                    tmp_u <- ttm(tmp_u, t(as.matrix((A[[m]][r,
                      ])))^Beta, m = m)
                  }
                  A[[n]][r, ] <- .positive(as.vector(tmp_u@data[]))
                  A[[n]][r, ] <- A[[n]][r, ]/norm(as.matrix(A[[n]][r,
                    ]), "F")
                }
                not_N <- setdiff(seq(N), N)
                numer <- Xr
                denom <- 0
                for (m in not_N) {
                  numer <- ttm(numer, t(as.matrix(A[[m]][r, ])),
                    m = m)
                  denom <- denom + as.numeric(t(A[[m]][r, ])^Beta %*%
                    A[[m]][r, ])
                }
                A[[N]][r, ] <- .positive(as.vector(numer@data[]/denom))
                X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
                E <- Xr - X_bar
            }
        }
        # Fill with Machine Epsilon
        for (n in seq(N)) {
            A[[n]][which(A[[n]] < pseudocount)] <- pseudocount
            A[[n]][which(is.infinite(A[[n]]))] <- pseudocount
            A[[n]][which(is.nan(A[[n]]))] <- pseudocount
            A[[n]][which(is.nan(A[[n]]))] <- pseudocount
        }
        # After Update U, V
        iter <- iter + 1
        X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
        RecError[iter] <- .recError(X, X_bar)
        TrainRecError[iter] <- .recError((1-M_NA+M)*X, (1-M_NA+M)*X_bar)
        TestRecError[iter] <- .recError((M_NA-M)*X, (M_NA-M)*X_bar)
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]

        if (viz && !is.null(figdir) && N == 3) {
            png(filename = paste0(figdir, "/", iter, ".png"))
            plotTensor3D(X_bar)
            dev.off()
        }
        if (viz && is.null(figdir) && N == 3) {
            plotTensor3D(X_bar)
        }
        if (verbose) {
            cat(paste0(iter-1, " / ", num.iter, " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
        if (is.nan(RelChange[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }
    if (viz && !is.null(figdir) && N == 3) {
        png(filename = paste0(figdir, "/finish.png"))
        plotTensor3D(X_bar)
        dev.off()
        png(filename = paste0(figdir, "/original.png"))
        plotTensor3D(X)
        dev.off()
    }
    if (viz && is.null(figdir) && N == 3) {
        plotTensor3D(X_bar)
    }
    names(RecError) <- c("offset", seq_len(iter-1))
    names(TrainRecError) <- c("offset", seq_len(iter-1))
    names(TestRecError) <- c("offset", seq_len(iter-1))
    names(RelChange) <- c("offset", seq_len(iter-1))

    # normalization
    S <- apply(A[[N]], 1, function(an){
        norm(as.matrix(an), "F")
    })
    A[[N]] <- A[[N]] / S

    return(list(S = S, A = A,
        RecError = RecError,
        TrainRecError = TrainRecError,
        TestRecError = TestRecError,
        RelChange = RelChange))
}

.checkNTF <- function(X, M, pseudocount, initA, fixA, rank, Alpha, Beta,
    thr, num.iter, viz, figdir, verbose){
    stopifnot(is.array(X@data))
    if(!is.null(M)){
        if(!identical(dim(X), dim(M))){
            stop("Please specify the dimensions of X and M are same")
        }
        .checkZeroNA(X, M, type="Tensor")
    }
    stopifnot(is.numeric(pseudocount))
    if(!is.null(initA)){
        dimX <- dim(X)
        ncolA <- as.vector(unlist(lapply(initA, ncol)))
        if(!identical(dimX, ncolA)){
            stop("Please specify the dimensions of X and ncol(A[[k]]) are same")
        }
    }
    if(!is.logical(fixA)){
        if(!is.vector(fixA)){
            stop("Please specify the fixA as a logical or a logical vector such as c(TRUE, FALSE, TRUE)")
        }else{
            if(length(dim(X)) != length(fixA)){
                stop("Please specify the length of fixA same as the order of X")
            }
        }
    }
    stopifnot(is.numeric(rank))
    stopifnot(is.numeric(Alpha))
    stopifnot(is.numeric(Beta))
    stopifnot(is.numeric(thr))
    stopifnot(is.numeric(num.iter))
    stopifnot(is.logical(viz))
    stopifnot(is.logical(verbose))
    if(!is.character(figdir) && !is.null(figdir)){
        stop("Please specify the figdir as a string or NULL")
    }
    if (verbose) {
        cat("Initialization step is running...\n")
    }
}

.initNTF <- function(X, M, pseudocount, fixA, rank, init, initA, Alpha, Beta,
    algorithm, thr, verbose){
    N <- length(dim(X))
    T1 <- NULL
    E <- NULL
    fixA <- rep(fixA, length=length(dim(X)))
    # NA mask
    M_NA <- X
    M_NA@data[] <- 1
    M_NA@data[which(is.na(X@data))] <- 0
    if(is.null(M)){
        M <- M_NA
    }
    pM <- M
    # Pseudo count
    X@data[which(is.na(X@data))] <- pseudocount
    X <- .pseudocount(X, pseudocount)
    pM <- .pseudocount(M, pseudocount)
    if(is.null(initA)){
        A <- list()
        length(A) <- N
        for (n in seq(N)) {
            if (init == "NMF") {
                Xn <- cs_unfold(X, m = n)@data
                res.nmf <- NMF(Xn, J = rank, algorithm = "KL")
                A[[n]] <- t(res.nmf$V)
                orderA <- order(sapply(seq(rank), function(x) {
                    norm(as.matrix(res.nmf$V[, x], "F")) * norm(as.matrix(res.nmf$U[,
                      x]), "F")
                }), decreasing = TRUE)
                if(rank != 1){
                    A[[n]] <- A[[n]][orderA, ]
                }
            } else if (init == "ABS-SVD") {
                Xn <- cs_unfold(X, m = n)@data
                res.svd <- svd(Xn)
                A[[n]] <- t(.positive(res.svd$v[, seq(rank)]))
                orderA <- order(sapply(res.svd$d[seq(rank)], function(x) {
                    norm(as.matrix(x), "F")
                }), decreasing = TRUE)
            } else if (init == "ALS") {
                Xn <- cs_unfold(X, m = n)@data
                res.svd <- svd(Xn)
                A[[n]] <- t(abs(res.svd$v[, seq(rank)]))
                orderA <- order(sapply(res.svd$d[seq(rank)], function(x) {
                    norm(as.matrix(x), "F")
                }), decreasing = TRUE)
            } else if (init == "Random") {
                A[[n]] <- matrix(runif(rank * dim(X)[n]), nrow = rank,
                    ncol = dim(X)[n])
                orderA <- order(apply(A[[n]], 1, function(x) {
                    norm(as.matrix(x), "F")
                }), decreasing = TRUE)
            }
            if (n != N) {
                A[[n]] <- t(apply(A[[n]], 1, function(x) {
                    x/norm(as.matrix(x), "F")
                }))
            }
        }
    }else{
        A <- initA
    }

    RecError = c()
    TrainRecError = c()
    TestRecError = c()
    RelChange = c()

    RecError[1] <- thr * 10
    TrainRecError[1] <- thr * 10
    TestRecError[1] <- thr * 10
    RelChange[1] <- thr * 10

    if (algorithm == "HALS") {
        T1 <- eval(parse(text=.HALSCMD8(N)))
    }
    if (algorithm == "Alpha-HALS" || algorithm == "Beta-HALS") {
        X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
        E <- X - X_bar
    }
    if (algorithm == "Frobenius") {
        Beta = 2
        algorithm = "Beta"
    }
    if (algorithm == "KL") {
        Alpha = 1
        algorithm = "Alpha"
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
    list(X=X, M=M, pM=pM, M_NA=M_NA, fixA=fixA, N=N,
        A=A, RecError=RecError, TrainRecError=TrainRecError,
        TestRecError=TestRecError, RelChange=RelChange,
        Alpha=Alpha, Beta=Beta, algorithm=algorithm, E=E, T1=T1)
}
