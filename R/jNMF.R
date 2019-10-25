jNMF <- function(X, J = 3, w=NULL, algorithm = "KL", p=1,
    thr = 1e-10, num.iter = 100,
    viz = FALSE, figdir = NULL, verbose = FALSE) {
    if(!is.list(X)){
        stop("input X must be specified as a list!")
    }
    if(length(X) < 2){
        stop("input list X must have at least two datasets!")
    }
    if(is.null(w)){
        w <- rep(1, length=length(X))
    }else{
        if(length(X) != length(w)){
            stop("The length of weight vector must be same as that of input list X!")
        }else{
            w <- w / sum(w)
        }
    }
    K <- length(X)
    lapply(seq_along(X), function(x){
        X[[x]][which(X[[x]] == 0)] <<- 1e-10
    })

    # Initialization
    W <- matrix(runif(nrow(X[[1]])*J),
        nrow=nrow(X[[1]]), ncol=J)
    W <- .columnNorm(W * W)

    V <- lapply(seq_along(X), function(x){
        tmp <- matrix(runif(nrow(X[[1]])*J),
        nrow=nrow(X[[1]]), ncol=J)
        .columnNorm(tmp * tmp)
    })

    H <- lapply(seq_along(X), function(x){
        tmp <- matrix(runif(ncol(X[[x]])*J),
            nrow=ncol(X[[x]]), ncol=J)
        .columnNorm(tmp * tmp)
    })
    RecError = c()
    RelChange = c()
    RecError[1] <- thr * 10
    RelChange[1] <- thr * 10
    iter <- 1
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
    while ((RecError[iter] > thr) && (iter <= num.iter)) {
        # Before Update W, H_k
        X_bar <- lapply(seq_len(K), function(k){
            .recMatrix(W + V[[k]], H[[k]])
        })
        pre_Error <- sqrt(sum(unlist(lapply(seq_along(X), function(x){
            .recError(X[[x]], X_bar[[x]], notsqrt=TRUE)
        }))))

        # Update W
        W_numer <- matrix(0, nrow=nrow(W), ncol=ncol(W))
        W_denom <- matrix(0, nrow=nrow(W), ncol=ncol(W))
        for(k in seq_len(K)){
            W_numer <- W_numer + w[k] * (X[[k]] * (X_bar[[k]])^(-p)) %*% H[[k]]
            W_denom <- W_denom + w[k] * ((W + V[[k]]) %*% t(H[[k]]))^(1-p) %*% H[[k]]
        }
        W <- .columnNorm(.positive(W * W_numer / W_denom))

        # Update H_k
        for(k in seq_len(K)){
            Hk_numer <- (t(X[[k]]) * t(X_bar[[k]])^(-p)) %*% (W + V[[k]])
            Hk_denom <- (H[[k]] %*% t(W + V[[k]]))^(1-p) %*% W
            H[[k]] <- H[[k]] * Hk_numer / Hk_denom
        }

        # Update V_k
        for(k in seq_len(K)){
            Vk_numer <- (X[[k]] * (X_bar[[k]])^(-p)) %*% H[[k]]
            Vk_denom <- ((W + V[[k]]) %*% t(H[[k]]))^(1-p) %*% H[[k]]
            V[[k]] <- V[[k]] * Vk_numer / Vk_denom
        }

        # After Update W, H_k
        iter <- iter + 1
        X_bar <- lapply(seq_len(K), function(k){
            .recMatrix(W + V[[k]], H[[k]])
        })
        RecError[iter] <- sqrt(sum(unlist(lapply(seq_along(X), function(x){
            .recError(X[[x]], X_bar[[x]], notsqrt=TRUE)
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
    names(RecError) <- c("offset", 1:(iter-1))
    names(RelChange) <- c("offset", 1:(iter-1))
    return(list(W = W, V = V, H = H, RecError = RecError,
        RelChange = RelChange))
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
