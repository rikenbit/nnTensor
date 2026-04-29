PSDTF <- function(X, rank=3, initH=NULL, initV=NULL,
    thr=1e-10, num.iter=100, verbose=FALSE) {

    ######################################
    # Argument Check
    ######################################
    .checkPSDTF(X, rank, initH, initV, thr, num.iter, verbose)

    ######################################
    # Initialization
    ######################################
    int <- .initPSDTF(X, rank, initH, initV)
    M <- int$M
    N <- int$N
    K <- rank
    H <- int$H
    V <- int$V
    X_slices <- int$X_slices

    RecError <- c()
    RelChange <- c()
    RecError[1] <- thr * 10
    RelChange[1] <- thr * 10

    ######################################
    # Iteration
    ######################################
    for (iter in seq_len(num.iter)) {
        pre_Error <- RecError[length(RecError)]

        # Compute Y_n = sum_k H[k,n] * V[[k]] and their inverses
        Y_list <- vector("list", N)
        Y_inv <- vector("list", N)
        for (n in seq_len(N)) {
            Y_n <- matrix(0, nrow = M, ncol = M)
            for (k in seq_len(K)) {
                Y_n <- Y_n + H[k, n] * V[[k]]
            }
            Y_list[[n]] <- Y_n
            Y_inv[[n]] <- solve(Y_n)
        }

        # Update H (Eq. 31)
        for (k in seq_len(K)) {
            for (n in seq_len(N)) {
                YiV <- Y_inv[[n]] %*% V[[k]]
                numer_h <- sum(diag(YiV %*% Y_inv[[n]] %*% X_slices[[n]]))
                denom_h <- sum(diag(YiV))
                H[k, n] <- H[k, n] * sqrt(numer_h / denom_h)
            }
        }

        # Recompute Y_inv after H update
        for (n in seq_len(N)) {
            Y_n <- matrix(0, nrow = M, ncol = M)
            for (k in seq_len(K)) {
                Y_n <- Y_n + H[k, n] * V[[k]]
            }
            Y_list[[n]] <- Y_n
            Y_inv[[n]] <- solve(Y_n)
        }

        # Update V (Eq. 34)
        for (k in seq_len(K)) {
            P_k <- matrix(0, nrow = M, ncol = M)
            Q_k <- matrix(0, nrow = M, ncol = M)
            for (n in seq_len(N)) {
                P_k <- P_k + H[k, n] * Y_inv[[n]]
                Q_k <- Q_k + H[k, n] * Y_inv[[n]] %*% X_slices[[n]] %*% Y_inv[[n]]
            }
            L_k <- t(chol(Q_k))
            inner <- t(L_k) %*% V[[k]] %*% P_k %*% V[[k]] %*% L_k
            inner_inv_sqrt <- .matrix_inv_sqrt(inner)
            V[[k]] <- V[[k]] %*% L_k %*% inner_inv_sqrt %*% t(L_k) %*% V[[k]]
            # Scale normalization: tr(V_k) = s, H_kn <- s * H_kn
            s <- sum(diag(V[[k]]))
            if (s > 0) {
                V[[k]] <- V[[k]] / s
                H[k, ] <- H[k, ] * s
            }
        }

        # Compute LogDet divergence
        total_err <- 0
        for (n in seq_len(N)) {
            Y_n <- matrix(0, nrow = M, ncol = M)
            for (k in seq_len(K)) {
                Y_n <- Y_n + H[k, n] * V[[k]]
            }
            total_err <- total_err + .logDetDivergence(X_slices[[n]], Y_n)
        }
        RecError <- c(RecError, total_err)
        rel <- abs(pre_Error - total_err) / max(abs(total_err), .Machine$double.eps)
        RelChange <- c(RelChange, rel)

        if (verbose) {
            cat(paste0(iter, " / ", num.iter,
                " |Previous Error - Error| / Error = ", rel, "\n"))
        }

        if (is.nan(total_err)) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
        if (iter > 1 && rel < thr) {
            break
        }
    }

    names(RecError) <- c("offset", seq_len(length(RecError) - 1))
    names(RelChange) <- c("offset", seq_len(length(RelChange) - 1))

    list(H = H, V = V, RecError = RecError, RelChange = RelChange)
}
