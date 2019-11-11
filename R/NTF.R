NTF <- function(X, rank = 3, algorithm = "KL", init = "NMF", Alpha = 1,
    Beta = 2, thr = 1e-10, num.iter = 100, viz = FALSE, figdir = NULL,
    verbose = FALSE){
    if (verbose) {
        cat("Initialization step is running...\n")
    }
    X <- .pseudocount(X)
    N <- length(dim(X))
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
        } else if (init == "ALS") {
            Xn <- cs_unfold(X, m = n)@data
            res.svd <- svd(Xn)
            A[[n]] <- .positive(res.svd$u[seq(rank), ])
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
        A[[n]] <- A[[n]][orderA, ]
        if (n != N) {
            A[[n]] <- t(apply(A[[n]], 1, function(x) {
                x/norm(as.matrix(x), "F")
            }))
        }
    }

    iter <- 1
    RecError = c()
    RelChange = c()
    X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
    RecError[1] <- .recError(X, X_bar)
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
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        # Before Update U, V
        X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
        pre_Error <- .recError(X, X_bar)
        if (algorithm == "Alpha") {
            for (n in seq(N)) {
                X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
                A_notn <- .KhatriRao_notn(A, n)
                A[[n]] <- A[[n]] * (t(A_notn) %*% (cs_unfold(X,
                  m = n)@data/cs_unfold(X_bar, m = n)@data)^Alpha)^(1/Alpha)
                if (n != N) {
                  A[[n]] <- t(apply(A[[n]], 1, function(x) {
                    x/norm(as.matrix(x), "F")
                  }))
                }
            }
        }
        else if (algorithm == "Beta") {
            X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
            for (n in seq(N)) {
                A_notn <- .KhatriRao_notn(A, n)
                numer <- t(A_notn) %*% (cs_unfold(X, m = n)@data/cs_unfold(X_bar^(Beta - 1), m = n)@data)
                denom <- t(A_notn) %*% cs_unfold(X_bar^(Beta),
                  m = n)@data
                A[[n]] <- A[[n]] * numer/denom
                if (n != N) {
                  A[[n]] <- t(apply(A[[n]], 1, function(x) {
                    x/norm(as.matrix(x), "F")
                  }))
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
        else {
            stop("Please specify the appropriate algorithm\n")
        }
        for (n in seq(N)) {
            if (any(is.infinite(A[[n]])) || any(is.nan(A[[n]]))) {
                stop("Inf or NaN is generated!\n")
            }
        }
        # After Update U, V
        iter <- iter + 1
        X_bar <- recTensor(rep(1, length = rank), A, idx=seq(N))
        RecError[iter] <- .recError(X, X_bar)
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
    names(RecError) <- c("offset", 1:(iter-1))
    names(RelChange) <- c("offset", 1:(iter-1))
    # normalization
    S <- apply(A[[N]], 1, function(an){
        norm(as.matrix(an), "F")
    })
    A[[N]] <- A[[N]] / S
    # sort
    orderS <- order(S, decreasing=TRUE)
    S <- S[orderS]
    for(n in seq(N)){
        A[[n]] <- A[[n]][orderS, ]
    }
    return(list(S = S, A = A, RecError = RecError, RelChange = RelChange))
}
