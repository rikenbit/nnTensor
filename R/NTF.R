NTF <-
function (X, rank = 3, algorithm = "KL", init = "NMF", Alpha = 1, 
    Beta = 2, thr = 1e-10, num.iter = 100, viz = FALSE, figdir = ".", 
    verbose = FALSE) 
{
    if (verbose) {
        cat("Initialization step is running...\n")
    }
    X <- .pseudocount(X)
    N <- length(dim(X))
    A <- list()
    length(A) <- N
    for (n in 1:N) {
        if (init == "NMF") {
            Xn <- cs_unfold(X, m = n)@data
            res.nmf <- NMF(Xn, J = rank, algorithm = "KL")
            A[[n]] <- t(res.nmf$V)
            orderA <- order(sapply(1:rank, function(x) {
                norm(as.matrix(res.nmf$V[, x], "F")) * norm(as.matrix(res.nmf$U[, 
                  x]), "F")
            }), decreasing = TRUE)
        }
        else if (init == "ALS") {
            Xn <- cs_unfold(X, m = n)@data
            res.svd <- svd(Xn)
            A[[n]] <- .positive(res.svd$u[1:rank, ])
            orderA <- order(sapply(res.svd$d[1:rank], function(x) {
                norm(as.matrix(x), "F")
            }), decreasing = TRUE)
        }
        else if (init == "Random") {
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
    X_bar <- recTensor(rep(1, length = rank), A)
    Error <- .recError(X, X_bar)
    if (algorithm == "HALS") {
        T1 <- (A[[1]] %*% t(A[[1]])) * (A[[2]] %*% t(A[[2]])) * 
            (A[[3]] %*% t(A[[3]]))
    }
    if (algorithm == "Alpha-HALS" || algorithm == "Beta-HALS") {
        X_bar <- recTensor(rep(1, length = rank), A)
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
    while ((Error > thr) && (iter <= num.iter)) {
        X_bar <- recTensor(rep(1, length = rank), A)
        pre_Error <- .recError(X, X_bar)
        if (algorithm == "Alpha") {
            for (n in 1:N) {
                X_bar <- recTensor(rep(1, length = rank), A)
                tmp1 <- A[[setdiff(1:rank, n)[1]]]
                tmp2 <- A[[setdiff(1:rank, n)[2]]]
                A_nonn <- khatri_rao(t(tmp1), t(tmp2))
                A[[n]] <- A[[n]] * (t(A_nonn) %*% (cs_unfold(X, 
                  m = n)@data/cs_unfold(X_bar, m = n)@data)^Alpha)^(1/Alpha)
                if (n != N) {
                  A[[n]] <- t(apply(A[[n]], 1, function(x) {
                    x/norm(as.matrix(x), "F")
                  }))
                }
            }
        }
        else if (algorithm == "Beta") {
            X_bar <- recTensor(rep(1, length = rank), A)
            for (n in 1:N) {
                tmp1 <- A[[setdiff(1:rank, n)[1]]]
                tmp2 <- A[[setdiff(1:rank, n)[2]]]
                A_nonn <- khatri_rao(t(tmp1), t(tmp2))
                numer <- t(A_nonn) %*% (cs_unfold(X, m = n)@data/cs_unfold(X_bar^(Beta - 
                  1), m = n)@data)
                denom <- t(A_nonn) %*% cs_unfold(X_bar^(Beta), 
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
            for (n in 1:N) {
                if (n == N) {
                  gamma[] <- 1
                }
                tmp1 <- A[[setdiff(1:rank, n)[1]]]
                tmp2 <- A[[setdiff(1:rank, n)[2]]]
                A_nonn <- khatri_rao(t(tmp1), t(tmp2))
                T2 = t(cs_unfold(X, m = n)@data) %*% A_nonn
                T3 = T1/(A[[n]] %*% t(A[[n]]))
                for (r in 1:rank) {
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
            for (r in 1:rank) {
                X_bar <- recTensor(rep(1, length = rank), A)
                Xr <- E + X_bar
                for (n in 1:N) {
                  if (Alpha == 0) {
                    tmp_Xr <- .positive(Xr)
                    tmp_Xr@data[, , ] <- log(tmp_Xr@data[, , 
                      ])
                    numer <- tmp_Xr
                  }
                  else {
                    numer <- .positive(Xr)^Alpha
                  }
                  denom <- 0
                  non_n <- setdiff(1:N, n)
                  for (m in non_n) {
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
                X_bar <- recTensor(rep(1, length = rank), A)
                E <- Xr - X_bar
            }
        }
        else if (algorithm == "Beta-HALS") {
            for (r in 1:rank) {
                X_bar <- recTensor(rep(1, length = rank), A)
                Xr <- E + X_bar
                for (n in 1:(N - 1)) {
                  non_n <- setdiff(1:N, n)
                  tmp_u <- Xr
                  for (m in non_n) {
                    tmp_u <- ttm(tmp_u, t(as.matrix((A[[m]][r, 
                      ])))^Beta, m = m)
                  }
                  A[[n]][r, ] <- .positive(as.vector(tmp_u@data[]))
                  A[[n]][r, ] <- A[[n]][r, ]/norm(as.matrix(A[[n]][r, 
                    ]), "F")
                }
                non_N <- setdiff(1:N, N)
                numer <- Xr
                denom <- 0
                for (m in non_N) {
                  numer <- ttm(numer, t(as.matrix(A[[m]][r, ])), 
                    m = m)
                  denom <- denom + as.numeric(t(A[[m]][r, ])^Beta %*% 
                    A[[m]][r, ])
                }
                A[[N]][r, ] <- .positive(as.vector(numer@data[]/denom))
                X_bar <- recTensor(rep(1, length = rank), A)
                E <- Xr - X_bar
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
        X_bar <- recTensor(rep(1, length = rank), A)
        if (viz) {
            png(filename = paste0(figdir, "/", iter, ".png"))
            plotTensor3D(X_bar)
            dev.off()
        }
        Error <- .recError(X, X_bar)
        Error <- abs(pre_Error - Error)/Error
        if (is.nan(Error)) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
        if (verbose) {
            cat(paste0(iter, " / ", num.iter, " |Previous Error - Error| / Error = ", 
                Error, "\n"))
        }
        iter <- iter + 1
    }
    if (viz) {
        png(filename = paste0(figdir, "/finish.png"))
        plotTensor3D(X_bar)
        dev.off()
    }
    return(list(S = .diag(X_bar), A = A))
}
