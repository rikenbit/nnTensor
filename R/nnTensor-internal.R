.columnNorm <- function(X){
    X_norm <- apply(X, 2, function(x){
        norm(as.matrix(x), "F")
    })
    t(t(X) / X_norm)
}

.recError <- function (X = NULL, X_bar = NULL, notsqrt = FALSE){
    if (is(X)[1] == "matrix" && is(X_bar)[1] == "matrix") {
        v <- as.vector(X_bar - X)
    }
    else if (is(X)[1] == "Tensor" && is(X_bar)[1] == "Tensor") {
        v <- vec(X_bar - X)
    }
    if(notsqrt){
        sum(v * v)
    }else{
        sqrt(sum(v * v))
    }
}

.positive <- function(X, thr = 1e-10){
    if (is(X)[1] == "matrix") {
        X[which(X < thr)] <- thr
    }
    else if (is(X)[1] == "Tensor") {
        X@data[which(X@data < thr)] <- thr
    }
    else if ("numeric" %in% is(X) && length(X) != 1) {
        X[which(X < thr)] <- thr
    }
    else if ("numeric" %in% is(X) && length(X) == 1) {
        X <- max(X, thr)
    }
    X
}

.recMatrix <- function(U = NULL, V = NULL){
    if (is(U)[1] != "matrix" || is(V)[1] != "matrix") {
        stop("Please specify the appropriate U and V\n")
    }
    return(U %*% t(V))
}

.argmaxj <- function(D){
    colmax <- apply(D, 2, max)
    which(colmax == max(colmax))
}

.doiter <- function(U, V, X, tol = 1e-04, J){
    Unew <- matrix(0, nrow = nrow(X), ncol = J)
    G <- U %*% t(V) %*% V - X %*% V
    VV <- matrix(1, nrow = nrow(X), ncol = 1) %*% diag(t(V) %*%
        V)
    S <- .positive(U - G/VV, 0) - U
    D <- -G * S - (VV * S * S)/2
    qi <- .argmaxj(D)[1]
    pinit <- max(D[, qi])
    for (i in 1:nrow(X)) {
        iter2 <- 1
        while ((D[i, qi] < tol * pinit) && (iter2 < J^2)) {
            s <- S[i, qi]
            Unew[i, qi] <- Unew[i, qi] + s
            G[i, ] <- G[i, ] + s * VV[qi, ]
            S[i, ] <- .positive(U[i, ] - (G[i, ]/VV[1, ]), 0) -
                U[i, ]
            D[i, ] <- -G[i, ] * S[i, ] - (VV[1, ] * S[i, ] *
                S[i, ])/2
            qi <- .argmaxj(D)[1]
            iter2 <- iter2 + 1
        }
    }
    Unew
}

.insertNULL <- function(rank, Iposition, N){
    out <- rep(0, length=N)
    out[setdiff(1:N, Iposition)] <- rank
    out
}

.pseudocount <- function(X, pseudocount = 1e-10){
    X@data[which(X@data == 0)] <- pseudocount
    X
}

.KhatriRao_notn <- function(A, n){
    idx <- setdiff(seq_len(length(A)), n)
    out <- t(A[[idx[1]]])
    for(notn in setdiff(idx, idx[1])){
        out <- khatri_rao(out, t(A[[notn]]))
    }
    out
}

.CPCore <- function(A){
    ranks <- unique(unlist(lapply(A, nrow)))
    if(length(ranks) >= 2){
        stop("The number of rows of all factor matrices must be same")
    }
    norms <- lapply(A, function(a){

    })
}

.slice <- function(X, mode = 1, column = 1){
    N <- length(dim(X))
    notmode <- setdiff(seq(N), mode)
    d <- dim(X)[notmode]
    modes <- rep(1, length=N)
    modes[notmode] <- d
    out <- rand_tensor(modes)
    tmp1 <- rep("", length=N)
    tmp2 <- rep("", length=N)
    tmp1[mode] <- 1
    tmp2[mode] <- "column"
    cmd <- paste0(
        "out[",
        paste0(tmp1, collapse=","),
        "] <- X[",
        paste0(tmp2, collapse=","),
        "]"
    )
    out
}

.contProd <- function(A, B, mode = 1){
    l <- dim(A)
    N <- length(l)
    out <- rep(0, length=l[mode])
    tmp1 <- rep("", length=N)
    tmp2 <- rep("", length=N)
    tmp1[mode] <- "i"
    tmp2[mode] <- 1
    cmd <- paste0(
        "out[i] <- sum(A[",
        paste(tmp1, collapse=","),
        "]@data * B[",
        paste(tmp2, collapse=","),
        "]@data)"
    )
    for(i in seq(l[mode])){
        eval(parse(text=cmd))
    }
    out
}

.HALSCMD1 <- function(N){
    paste0(
        "J_hat <- lapply(apply(expand.grid(",
        paste0("1:rank[", seq(N), "]", collapse=","),
        "), 1, function(x){",
        "as.list(x)",
        "}), function(y){",
        "unlist(y)",
        "})")
}

.HALSCMD2 <- function(N){
    paste0("j", seq(N),
        " <- J_hat[[j_ijk]][", seq(N), "]", collapse=";")
}

.HALSCMD3 <- function(N){
    paste0("A", seq(N), " <- as.matrix(A[[", seq(N),
        "]][j", seq(N), ", ])", collapse=";")
}

.HALSCMD4 <- function(N){
    paste0("S_old <- S[",
        paste0("j", seq(N), collapse=","), "]@data")
}

.HALSCMD5 <- function(N){
    paste0("S@data[",
        paste0("j", seq(N), collapse=","),
        "] <- .positive(as.numeric(S_old + as.vector(",
        "recTensor(S=E, A=list(",
        paste0("A", seq(N), collapse=","),
        "))@data)))"
        )
}

.HALSCMD6 <- function(N){
    paste0("diffS <- as.numeric(S_old - S[",
        paste0("j", seq(N), collapse=","),
        "]@data)"
        )
}

.HALSCMD7 <- function(N){
    paste0("E <- E + recTensor(S=diffS, A=list(",
        paste0("t(A", seq(N), ")", collapse=","),
        "), idx=modes)")
}

.HALSCMD8 <- function(N){
    tmp <- paste0("(A[[", seq(N), "]] %*% t(A[[", seq(N), "]]))")
    paste(tmp, collapse=" * ")
}

