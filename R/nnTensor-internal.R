.recError <-
function (X = NULL, X_bar = NULL) 
{
    if (is(X)[1] == "matrix" && is(X_bar)[1] == "matrix") {
        v <- as.vector(X_bar - X)
    }
    else if (is(X)[1] == "Tensor" && is(X_bar)[1] == "Tensor") {
        v <- vec(X_bar - X)
    }
    sqrt(sum(v * v))
}
.positive <-
function (X, thr = 1e-10) 
{
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
.recMatrix <-
function (U = NULL, V = NULL) 
{
    if (is(U)[1] != "matrix" || is(V)[1] != "matrix") {
        stop("Please specify the appropriate U and V\n")
    }
    return(U %*% t(V))
}
.argmaxj <-
function (D) 
{
    colmax <- apply(D, 2, max)
    which(colmax == max(colmax))
}
.doiter <-
function (U, V, X, tol = 1e-04, J) 
{
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
.pseudocount <-
function (X, pseudocount = 1e-10) 
{
    d <- dim(X)
    for (i in 1:d[1]) {
        Xi <- X@data[i, , ]
        Xi[which(Xi == 0)] <- pseudocount
        X@data[i, , ] <- Xi
    }
    X
}
.diag <-
function (S) 
{
    if (dim(S)[1] != dim(S)[2] || dim(S)[2] != dim(S)[3]) {
        stop("Symmetric Tensor is required!")
    }
    out <- rep(0, length = dim(S)[1])
    for (x in 1:length(out)) {
        out[x] <- S@data[x, x, x]
    }
    out
}
.slice <-
function (X, mode = 1, column = 1) 
{
    if (mode == 1) {
        d <- dim(X[column, , ])
        out <- rand_tensor(modes = c(1, d[1:2]))
        out[1, , ] <- X[column, , ]
    }
    else if (mode == 2) {
        d <- dim(X[, column, ])
        out <- rand_tensor(modes = c(d[1], 1, d[2]))
        out[, 1, ] <- X[, column, ]
    }
    else if (mode == 3) {
        d <- dim(X[, , column])
        out <- rand_tensor(modes = c(d[1:2], 1))
        out[, , 1] <- X[, , column]
    }
    else {
        stop("Wrong mode!\n")
    }
    out
}
.contProd <-
function (A, B, mode = 1) 
{
    l <- dim(A)
    if (mode == 1) {
        out <- rep(0, length = l[1])
        for (i in 1:l[1]) {
            out[i] <- sum(A[i, , ]@data * B[1, , ]@data)
        }
    }
    else if (mode == 2) {
        out <- rep(0, length = l[2])
        for (i in 1:l[2]) {
            out[i] <- sum(A[, i, ]@data * B[, 1, ]@data)
        }
    }
    else if (mode == 3) {
        out <- rep(0, length = l[3])
        for (i in 1:l[3]) {
            out[i] <- sum(A[, , i]@data * B[, , 1]@data)
        }
    }
    else {
        stop("Wrong mode!\n")
    }
    out
}
