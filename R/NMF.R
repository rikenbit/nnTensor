NMF <-
function (X, J = 3, algorithm = "Frobenius", Alpha = 1, Beta = 2,
    eta = 1e-04, thr1 = 1e-10, thr2 = 1e-10, tol = 1e-04, num.iter = 100,
    viz = FALSE, figdir = ".", verbose = FALSE)
{
    X[which(X == 0)] <- 1e-10
    U <- matrix(runif(nrow(X) * J), nrow = nrow(X), ncol = J)
    V <- matrix(runif(ncol(X) * J), nrow = ncol(X), ncol = J)
    U <- U * U
    V <- V * V
    RecError = c() # Added in 2018.9.14
    RelChange = c() # Added in 2018.9.14
    RecError[1] <- thr1 * 10
    RelChange[1] <- thr1 * 10
    iter <- 1
    if (algorithm == "HALS") {
        R <- X - U %*% t(V)
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
    while ((RecError[iter] > thr1) && (iter <= num.iter)) {
        # Before Update U, V
        X_bar <- .recMatrix(U, V)
        pre_Error <- .recError(X, X_bar)
        if (algorithm == "PGD") {
            U <- .positive(U - eta * (-2 * X %*% V + 2 * U %*%
                t(V) %*% V), thr2)
            V <- .positive(V - eta * (-2 * t(X) %*% U + 2 * V %*%
                t(U) %*% U), thr2)
        }
        else if (algorithm == "Alpha") {
            U <- U * (((X/(U %*% t(V)))^Alpha %*% V)/(matrix(1,
                nrow = nrow(X), ncol = 1) %*% colSums(V)))^(1/Alpha)
            V <- V * ((t(X/(U %*% t(V)))^Alpha %*% U)/(matrix(1,
                nrow = ncol(X), ncol = 1) %*% colSums(U)))^(1/Alpha)
        }
        else if (algorithm == "Beta") {
            U <- U * (((U %*% t(V))^(Beta - 2) * X) %*% V)/((U %*%
                t(V))^(Beta - 1) %*% V)
            V <- V * (t((U %*% t(V))^(Beta - 2) * X) %*% U)/(t((U %*%
                t(V))^(Beta - 1)) %*% U)
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
            Unew <- .doiter(U, V, X, tol = tol, J)
            U <- U + Unew
            Vnew <- .doiter(V, U, t(X), tol = tol, J)
            V <- V + Vnew
        }
        else {
            stop("Please specify the appropriate algorithm\n")
        }
        # After Update U, V
        iter <- iter + 1
        X_bar <- .recMatrix(U, V)
        RecError[iter] <- .recError(X, X_bar)
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if (viz) {
            png(filename = paste0(figdir, "/", iter-1, ".png"))
            image.plot(X_bar)
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
        image.plot(X_bar)
        dev.off()
        png(filename = paste0(figdir, "/original.png"))
        image.plot(X)
        dev.off()
    }
    names(RecError) <- c("offset", 1:(iter-1))
    names(RelChange) <- c("offset", 1:(iter-1))
    return(list(U = U, V = V, RecError = RecError,
        RelChange = RelChange))
}
