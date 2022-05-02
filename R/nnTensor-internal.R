.lapply_pb <- function(X, FUN, ...)
{
 env <- environment()
 pb_Total <- length(X)
 counter <- 0
 pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

 # wrapper around FUN
 wrapper <- function(...){
   curVal <- get("counter", envir = env)
   assign("counter", curVal +1 ,envir=env)
   setTxtProgressBar(get("pb", envir=env), curVal +1)
   FUN(...)
 }
 res <- lapply(X, wrapper, ...)
 close(pb)
 res
}

.binalize <- function(V){
    t(apply(V, 1, function(x){
        max.pos <- which(max(x)[1] == x)
        x[max.pos] <- 1L
        x[setdiff(seq_along(x), max.pos)] <- 0L
        x
    }))
}

.consensus <- function(V){
    if(ncol(V) == 1){
        matrix(1, nrow=nrow(V), ncol=nrow(V))
    }else{
        B <- .binalize(V)
        B %*% t(B)
    }
}

.averageC <- function(out, runtime){
    C <- matrix(0, nrow=nrow(out[[1]]$V), ncol=nrow(out[[1]]$V))
    for(i in seq_len(runtime)){
        C <- C + .consensus(out[[i]]$V)
    }
    C / runtime
}

.cophcor <- function(C){
    if(all(C == 1)){
        1
    }else{
        d.consensus <- as.dist(1 - C)
        hc <- hclust(d.consensus, method="average")
        d.coph <- cophenetic(hc)
        cor(d.consensus, d.coph, method="pearson")
    }
}

.ccc <- function(out.original, out.random, X, runtime){
    tmp1 <- .cophcor(.averageC(out.original, runtime))
    tmp2 <- .cophcor(.averageC(out.random, runtime))
    list(original=tmp1, random=tmp2)
}

.dispval <- function(C){
    sum(4 * (C - 0.5)^2) / nrow(C)^2
}

.dispersion <- function(out.original, out.random, X, runtime){
    tmp1 <- .dispval(.averageC(out.original, runtime))
    tmp2 <- .dispval(.averageC(out.random, runtime))
    list(original=tmp1, random=tmp2)
}

# Residual Sum of Squares: RecError^2
.rss <- function(out.original, out.random, X, runtime){
    tmp1 <- unlist(lapply(out.original, function(x){
        rev(x$RecError)[1]^2
    }))
    tmp2 <- unlist(lapply(out.random, function(x){
        rev(x$RecError)[1]^2
    }))
    list(original=tmp1, random=tmp2)
}
# Explained Variance
.evar <- function(out.original, out.random, X, runtime){
    tmp1 <- unlist(lapply(out.original, function(x){
        1 - rev(x$RecError)[1]^2 / sum(X)^2
    }))
    tmp2 <- unlist(lapply(out.random, function(x){
        1 - rev(x$RecError)[1]^2 / sum(X)^2
    }))
    list(original=tmp1, random=tmp2)
}
# sqrt(RSS) = RecError
.residuals <- function(out.original, out.random, X, runtime){
    tmp1 <- unlist(lapply(out.original, function(x){
        rev(x$RecError)[1]
    }))
    tmp2 <- unlist(lapply(out.random, function(x){
        rev(x$RecError)[1]
    }))
    list(original=tmp1, random=tmp2)
}

.L1 <- function(x){
    norm(as.matrix(x), "1")
}
.L2 <- function(x){
    norm(as.matrix(x), "2")
}

.sparseness <- function(x){
    nc <- ncol(x)
    nr <- nrow(x)
    mean(
        unlist(lapply(nc, function(n){
        (sqrt(nr) - .L1(x[,n]) / .L2(x[,n])) / (sqrt(nr) - 1)}))
    )
}

.spval.basis <- function(out){
    unlist(lapply(out, function(x){
        .sparseness(x$U)
    }))
}

.sparseness.basis <- function(out.original, out.random, X, runtime){
    tmp1 <- .spval.basis(out.original)
    tmp2 <- .spval.basis(out.random)
    list(original=tmp1, random=tmp2)
}

.spval.coef <- function(out){
    unlist(lapply(out, function(x){
        .sparseness(x$V)
    }))
}

.sparseness.coef <- function(out.original, out.random, X, runtime){
    tmp1 <- .spval.coef(out.original)
    tmp2 <- .spval.coef(out.random)
    list(original=tmp1, random=tmp2)
}

.F <- function(x){
    norm(as.matrix(x), "2")^2
}

.sparseness2 <- function(x){
    nc <- ncol(x)
    nr <- nrow(x)
    if(nc == 1){
        NA
    }else{
        a <- 1 / sqrt(nc)
        sp <- mean(unlist(lapply(nc, function(n){.F(x[,n])})))
        (sp - a) / (1 - a)        
    }
}

.spval2.basis <- function(out){
    unlist(lapply(out, function(x){
        .sparseness2(x$U)
    }))
}
.sparseness2.basis <- function(out.original, out.random, X, runtime){
    tmp1 <- .spval2.basis(out.original)
    tmp2 <- .spval2.basis(out.random)
    list(original=tmp1, random=tmp2)
}

.spval2.coef <- function(out){
    unlist(lapply(out, function(x){
        .sparseness2(x$V)
    }))
}

.sparseness2.coef <- function(out.original, out.random, X, runtime){
    tmp1 <- .spval2.coef(out.original)
    tmp2 <- .spval2.coef(out.random)
    list(original=tmp1, random=tmp2)
}

.norm.info.gain <- function(x){
    nc <- ncol(x)
    if(nc == 1){
        NA
    }else{
        1 - 1 / (nc*log2(nc)) * (- sum(x * log2(x)))
    }
}

.nig.basis <- function(out){
    unlist(lapply(out, function(x){
        .norm.info.gain(x$U)
    }))
}

.norm.info.gain.basis <- function(out.original, out.random, X, runtime){
    tmp1 <- .nig.basis(out.original)
    tmp2 <- .nig.basis(out.random)
    list(original=tmp1, random=tmp2)
}

.nig.coef <- function(out){
    unlist(lapply(out, function(x){
        .norm.info.gain(x$V)
    }))
}

.norm.info.gain.coef <- function(out.original, out.random, X, runtime){
    tmp1 <- .nig.coef(out.original)
    tmp2 <- .nig.coef(out.random)
    list(original=tmp1, random=tmp2)
}

.singular <- function(out.original, out.random, X, runtime){
    nc <- ncol((out.original[[1]]$V))
    d1 <- svd(X)$d
    tmp1 <- rev(cumsum(d1[seq(nc)]))[1] / sum(d1)
    tmp1 <- rep(tmp1, length(out.original))

    tmp2 <- unlist(lapply(seq_len(runtime), function(x, misc){
        X <- misc$X
        nc <- misc$nc
        d1 <- svd(.randomize(X))$d
        rev(cumsum(d1[seq(nc)]))[1] / sum(d1)
    }, misc=list(X=X, nc=nc)))
    list(original=tmp1, random=tmp2)
}

.volval <- function(out){
    unlist(lapply(out, function(x){
        det(t(x$V) %*% x$V)
    }))
}

.volume <- function(out.original, out.random, X, runtime){
    tmp1 <- .volval(out.original)
    tmp2 <- .volval(out.random)
    list(original=tmp1, random=tmp2)
}

.conval <- function(out, runtime){
    Vs <- matrix(0, nrow(out[[1]]$V), ncol(out[[1]]$V))
    for(i in seq_len(runtime)){
        vnorm <- apply(out[[i]]$V, 2,
            function(v){norm(as.matrix(v), "F")})
        unorm <- apply(out[[i]]$U, 2,
            function(u){norm(as.matrix(u), "F")})
        Vs <- Vs + out[[i]]$V[,
            order(vnorm*unorm, decreasing=TRUE)]
    }
    Vs <- Vs / runtime
    res.svd <- svd(Vs)
    lambda <- res.svd$d^2
    lambda[1]/rev(lambda)[1]
}

.condition <- function(out.original, out.random, X, runtime){
    tmp1 <- .conval(out.original, runtime)
    tmp2 <- .conval(out.random, runtime)
    list(original=tmp1, random=tmp2)
}

.all <- function(out.original, out.random, X, runtime){
    cat("ccc\n")
    ccc = .ccc(out.original, out.random, X, runtime)
    cat("dispersion\n")
    dispersion = .dispersion(out.original, out.random, X, runtime)
    cat("rss\n")
    rss = .rss(out.original, out.random, X, runtime)
    cat("evar\n")
    evar = .evar(out.original, out.random, X, runtime)
    cat("residuals\n")
    residuals = .residuals(out.original, out.random, X, runtime)
    cat("sparseness.basis\n")
    sparseness.basis = .sparseness.basis(out.original, out.random, X, runtime)
    cat("sparseness.coef\n")
    sparseness.coef = .sparseness.coef(out.original, out.random, X, runtime)
    cat("sparseness2.basis\n")
    sparseness2.basis = .sparseness2.basis(out.original, out.random, X, runtime)
    cat("sparseness2.coef\n")
    sparseness2.coef = .sparseness2.coef(out.original, out.random, X, runtime)
    cat("norm.info.gain.basis\n")
    norm.info.gain.basis = .norm.info.gain.basis(out.original, out.random, X, runtime)
    cat("norm.info.gain.coef\n")
    norm.info.gain.coef = .norm.info.gain.coef(out.original, out.random, X, runtime)
    cat("singular\n")
    singular = .singular(out.original, out.random, X, runtime)
    cat("volume\n")
    volume = .volume(out.original, out.random, X, runtime)
    cat("condition\n")
    condition = .condition(out.original, out.random, X, runtime)
    list(ccc=ccc, dispersion=dispersion,
        rss=rss, evar=evar, residuals=residuals,
        sparseness.basis=sparseness.basis, sparseness.coef=sparseness.coef,
        sparseness2.basis=sparseness2.basis, sparseness2.coef=sparseness2.coef,
        norm.info.gain.basis=norm.info.gain.basis,
        norm.info.gain.coef=norm.info.gain.coef,
        singular=singular, volume=volume, condition=condition
    )
}

.flist2 <- list(
    "all" = .all,
    "ccc" = .ccc,
    "dispersion" = .dispersion,
    "rss" = .rss,
    "evar" = .evar,
    "residuals" = .residuals,
    "sparseness.basis" = .sparseness.basis,
    "sparseness.coef" = .sparseness.coef,
    "sparseness2.basis" = .sparseness2.basis,
    "sparseness2.coef" = .sparseness2.coef,
    "norm.info.gain.basis" = .norm.info.gain.basis,
    "norm.info.gain.coef" = .norm.info.gain.coef,
    "singular" = .singular,
    "volume" = .volume,
    "condition" = .condition
)

.randomize <- function(x){
    nr <- nrow(x)
    nc <- ncol(x)
    matrix(sample(as.vector(x)), nrow=nr, ncol=nc)
}

.errorPosition <- function(recerror){
    runtime <- length(recerror)
    unlist(lapply(seq_len(runtime), function(x, recerror){
            if(inherits(recerror[x], "try-error")){
                x
            }
        }, recerror=recerror))
}

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