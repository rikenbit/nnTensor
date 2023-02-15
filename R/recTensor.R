recTensor <- function(S = NULL, A = NULL, idx = seq_along(dim(S)), reverse = FALSE){
    if (is(S)[1] == "Tensor") {
        X_bar <- S
        if(length(dim(S)) > 3){
            idx <- seq_along(A)
        }
    }
    else if ("numeric" %in% is(S) && length(S) == 1) {
        X_bar <- rand_tensor(rep(1, length=length(idx)))
        X_bar@data[] <- S
    }
    else if ("numeric" %in% is(S) && length(S) != 1) {
        X_bar <- rand_tensor(sapply(A, nrow))
        X_bar@data[] <- 0
        for (n in seq_along(S)) {
            cmd = paste0("X_bar@data[",
                paste(rep("n", length=length(A)), collapse=","),
                "] <- S[n]")
            eval(parse(text=cmd))
        }
    }
    if (!reverse) {
        for (n in idx) {
            X_bar <- ttm(X_bar, t(A[[n]]), m = n)
        }
    }
    else {
        for (n in idx) {
            X_bar <- ttm(X_bar, A[[n]], m = n)
        }
    }
    return(X_bar)
}
