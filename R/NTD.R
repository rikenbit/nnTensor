NTD <-
function (X, rank = c(3, 3, 3), algorithm = "KL", init = "NMF",
    Alpha = 1, Beta = 2, thr = 1e-10, num.iter = 100, viz = FALSE,
    figdir = ".", verbose = FALSE)
{
    if (length(dim(X)) != length(rank)) {
        stop("Please specify the appropriate rank\n")
    }
    if (verbose) {
        cat("Initialization step is running...\n")
    }
    X <- .pseudocount(X)
    N <- length(dim(X))
    A <- list()
    length(A) <- N
    if (init == "NMF") {
        sapply(1:N, function(n) {
            Xn <- cs_unfold(X, m = n)@data
            An <- t(NMF(Xn, J = rank[n], algorithm = "KL")$V)
            A[[n]] <<- t(apply(An, 1, function(x) {
                x/norm(as.matrix(x), "F")
            }))
        })
        S <- recTensor(X, A, reverse = TRUE)
    }
    else if (init == "ALS") {
        sapply(1:N, function(n) {
            Xn <- cs_unfold(X, m = n)@data
            An <- .positive(svd(Xn)$u[1:rank[n], ])
            A[[n]] <<- t(apply(An, 1, function(x) {
                x/norm(as.matrix(x), "F")
            }))
        })
        S <- recTensor(X, A, reverse = TRUE)
    }
    else if (init == "Random") {
        A[[1]] <- matrix(runif(rank[1] * dim(X)[1]), nrow = rank[1],
            ncol = dim(X)[1])
        A[[2]] <- matrix(runif(rank[2] * dim(X)[2]), nrow = rank[2],
            ncol = dim(X)[2])
        A[[3]] <- matrix(runif(rank[3] * dim(X)[3]), nrow = rank[3],
            ncol = dim(X)[3])
        S <- recTensor(X, A, reverse = TRUE)
    }
    iter <- 1
    RecError = c() # Added in 2018.9.14
    RelChange = c() # Added in 2018.9.14
    X_bar <- recTensor(S, A)
    RecError[1] <- .recError(X, X_bar)
    RelChange[1] <- thr * 10
    if (algorithm == "HALS") {
        E <- X - recTensor(S, A)
        J_hat <- lapply(apply(expand.grid(1:rank[1], 1:rank[2],
            1:rank[3]), 1, function(x) {
            as.list(x)
        }), function(y) {
            unlist(y)
        })
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
    while ((RecError[iter] > thr) && (iter <= num.iter)) {
        # Before Update U, V
        X_bar <- recTensor(S, A)
        pre_Error <- .recError(X, X_bar)
        for (n in 1:N) {
            if (algorithm == "Alpha") {
                S_A <- t(cs_unfold(S, m = n)@data) %*% kronecker_list(sapply(rev(setdiff(1:N,
                  n)), function(x) {
                  A[[x]]
                }, simplify = FALSE))
                Xn <- cs_unfold(X, m = n)@data
                numer <- S_A %*% (Xn/t(t(A[[n]]) %*% S_A))^Alpha
                denom <- t(as.matrix(rep(1, dim(X)[n]) %*% t(rowSums(S_A))))
                A[[n]] <- A[[n]] * (numer/denom)^(1/Alpha)
            }
            else if (algorithm == "Beta") {
                S_A <- t(cs_unfold(S, m = n)@data) %*% kronecker_list(sapply(rev(setdiff(1:N,
                  n)), function(x) {
                  A[[x]]
                }, simplify = FALSE))
                Xn <- cs_unfold(X, m = n)@data
                Xn_bar <- cs_unfold(recTensor(S, A), m = n)@data
                numer <- S_A %*% (Xn * Xn_bar^(Beta - 1))
                denom <- S_A %*% (t(S_A) %*% A[[n]])^Beta
                A[[n]] <- A[[n]] * numer/denom
            }
            else if (algorithm == "HALS") {
                X_bar <- recTensor(S, A, idx = setdiff(1:N, n))
                for (jn in 1:nrow(A[[n]])) {
                  X_barkn <- .slice(X_bar, mode = n, column = jn)
                  wjn <- fnorm(X_barkn)^2
                  ajn <- .positive(A[[n]][jn, ] + .contProd(E,
                    X_barkn, mode = n)/wjn)
                  E <- E + ttm(X_barkn, as.matrix(A[[n]][jn,
                    ] - ajn), m = n)
                  A[[n]][jn, ] <- ajn/norm(as.matrix(ajn), "F")
                }
            }
            else {
                stop("Please specify the appropriate algorithm\n")
            }
        }
        if (algorithm == "Alpha") {
            S <- .positive(recTensor(X, A, reverse = TRUE))
            numer <- (X/recTensor(S, A))^Alpha
            denom <- X
            denom[, , ] <- 1
            for (n in 1:N) {
                numer <- ttm(numer, A[[n]], m = n)
                denom <- ttm(denom, A[[n]], m = n)
            }
            S <- S * (numer/denom)^(1/Alpha)
        }
        else if (algorithm == "Beta") {
            X_bar <- recTensor(S, A)
            numer <- X * X_bar^(Beta - 1)
            denom <- X_bar^Beta
            for (n in 1:N) {
                numer <- ttm(numer, A[[n]], m = n)
                denom <- ttm(denom, A[[n]], m = n)
            }
            S <- S * numer/denom
        }
        else if (algorithm == "HALS") {
            for (j_ijk in 1:length(J_hat)) {
                j1 <- J_hat[[j_ijk]][1]
                j2 <- J_hat[[j_ijk]][2]
                j3 <- J_hat[[j_ijk]][3]
                A1 <- as.matrix(A[[1]][j1, ])
                A2 <- as.matrix(A[[2]][j2, ])
                A3 <- as.matrix(A[[3]][j3, ])
                S_old <- S[j1, j2, j3]@data
                S@data[j1, j2, j3] <- .positive(as.numeric(S_old +
                  as.vector(recTensor(E, list(A1, A2, A3))@data)))
                diffS <- as.numeric(S_old - (S[j1, j2, j3])@data)
                E <- E + recTensor(diffS, list(t(A1), t(A2),
                  t(A3)))
            }
        }
        else {
            stop("Please specify the appropriate algorithm\n")
        }
        for (n in 1:N) {
            if (any(is.infinite(A[[n]])) || any(is.nan(A[[n]]))) {
                stop("Inf or NaN is generated!\n")
            }
        }
        # After Update U, V
        iter <- iter + 1
        X_bar <- recTensor(S, A)
        RecError[iter] <- .recError(X, X_bar)
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if (viz) {
            png(filename = paste0(figdir, "/", iter, ".png"))
            plotTensor3D(X_bar)
            dev.off()
        }
        if (verbose) {
            cat(paste0(iter-1, " / ", num.iter, " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
        if (is.nan(RelChange[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }
    if (viz) {
        png(filename = paste0(figdir, "/finish.png"))
        plotTensor3D(X_bar)
        dev.off()
        png(filename = paste0(figdir, "/original.png"))
        plotTensor3D(X)
        dev.off()
    }
    names(RecError) <- c("offset", 1:(iter-1))
    names(RelChange) <- c("offset", 1:(iter-1))
    return(list(S = S, A = A, RecError = RecError, RelChange = RelChange))
}
