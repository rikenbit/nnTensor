ZITNMF <- function(X, Z=NULL, pseudocount=1e-10,
    initF=NULL, initA=NULL, initW=0.5, fixF=FALSE, fixA=FALSE,
    initializer = c("NMF", "ALS", "Random"),
    J=3, Beta=2, phi=1, thr=1e-10, num.iter=100, verbose=FALSE) {

    ######################################
    # Argument Check
    ######################################
    .checkZITNMF(X, Z, pseudocount, initF, initA, fixF, fixA,
        J, Beta, phi, thr, num.iter, verbose)
    initializer <- match.arg(initializer)

    ######################################
    # Initialization
    ######################################
    int <- .initZITNMF(X, Z, pseudocount, initF, initA, initW,
        J, Beta, initializer)
    X <- int$X
    E <- int$E
    Fmat <- int$Fmat
    A <- int$A
    w <- int$w
    BetaDivergenceHistory <- int$BetaDivergenceHistory
    BetaDivergenceRelativeChanges <- int$BetaDivergenceChanges

    ######################################
    # Iteration
    ######################################
    X_estimate <- tcrossprod(Fmat, A)

    for (iter in seq_len(num.iter)) {
        if (!is.null(thr) && iter > 2) {
            if (BetaDivergenceRelativeChanges[iter - 1] > thr) {
                break
            }
        }

        # Update Z
        Z <- w / (w + (1 - w) * .dtweedie_beta_representation(
            y = as.numeric(X) * 0, mu = as.numeric(X_estimate),
            phi = as.numeric(phi), beta = as.numeric(Beta)))
        Z <- Z * (X == 0)

        # Update w
        epsilon <- .Machine$double.eps
        w <- mean(Z)
        if (w == 0) w <- epsilon

        # Update F
        if (!fixF) {
            numerF <- ((E - Z) * X * tcrossprod(Fmat, A)^(Beta - 2)) %*% A
            denomF <- ((E - Z) * tcrossprod(Fmat, A)^(Beta - 1)) %*% A
            Fmat <- Fmat * (numerF / denomF)^.rho(Beta)
        }

        # Update A
        if (!fixA) {
            numerA <- crossprod((E - Z) * X * tcrossprod(Fmat, A)^(Beta - 2), Fmat)
            denomA <- crossprod((E - Z) * tcrossprod(Fmat, A)^(Beta - 1), Fmat)
            A <- A * (numerA / denomA)^.rho(Beta)
        }

        # After Update
        X_estimate <- tcrossprod(Fmat, A)

        # Calculate Beta-divergence
        BetaDivergenceHistory[iter] <- .betaDivergence(X,
            (E - Z) * X_estimate, Beta)
        if (iter > 1) {
            BetaDivergenceRelativeChanges[iter] <-
                (BetaDivergenceHistory[iter] - BetaDivergenceHistory[iter - 1]) /
                BetaDivergenceHistory[iter - 1]
        } else {
            BetaDivergenceRelativeChanges[iter] <- 0
        }

        if (verbose) {
            cat(paste0(iter - 1, " / ", num.iter, " BetaDivergence = ",
                BetaDivergenceHistory[iter], "\n"))
        }

        if (is.infinite(BetaDivergenceHistory[iter])) {
            stop("Inf is generated. Please run again or change the parameters.\n")
        }
        if (is.nan(BetaDivergenceHistory[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }

    list(F = Fmat, A = A, Z = Z, w = w,
        RecError = BetaDivergenceHistory,
        RelChange = BetaDivergenceRelativeChanges)
}
