# Generate toy PSD tensor
set.seed(123)
M <- 4
N <- 10
K <- 2
V_true <- list()
for (k in 1:K) {
    A <- matrix(runif(M * M), M, M)
    V_true[[k]] <- A %*% t(A)
}
H_true <- matrix(runif(K * N), K, N)
X <- array(0, dim = c(M, M, N))
for (n in 1:N) {
    for (k in 1:K) {
        X[, , n] <- X[, , n] + H_true[k, n] * V_true[[k]]
    }
}

out <- PSDTF(X, rank = K, num.iter = 5)

expect_equivalent(length(out), 4)
expect_identical(names(out), c("H", "V", "RecError", "RelChange"))

expect_identical(is.matrix(out$H), TRUE)
expect_equal(dim(out$H), c(K, N))

expect_identical(is.list(out$V), TRUE)
expect_equal(length(out$V), K)
expect_equal(dim(out$V[[1]]), c(M, M))

expect_identical(is.numeric(out$RecError), TRUE)
expect_identical(is.numeric(out$RelChange), TRUE)

# Check V_k are PSD (eigenvalues >= 0)
for (k in 1:K) {
    eigvals <- eigen(out$V[[k]], symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigvals >= -1e-10))
}
