X <- toyModel("NMF")
M <- X
M[,] <- 0
M[,] <- rbinom(length(M), 1, 0.5)

out1 <- NMF(X, M=M, J=3, algorithm="Frobenius", num.iter=2)
out2 <- NMF(X, M=M, J=3, algorithm="KL", num.iter=2)
out3 <- NMF(X, M=M, J=3, algorithm="IS", num.iter=2)
out4 <- NMF(X, M=M, J=3, algorithm="Pearson", num.iter=2)
out5 <- NMF(X, M=M, J=3, algorithm="Hellinger", num.iter=2)
out6 <- NMF(X, M=M, J=3, algorithm="Neyman", num.iter=2)
out7 <- NMF(X, M=M, J=3, algorithm="Alpha", num.iter=2)
out8 <- NMF(X, M=M, J=3, algorithm="Beta", num.iter=2)
out9 <- NMF(X, M=M, J=3, algorithm="PGD", num.iter=2)

expect_equivalent(length(out1), 10)
expect_equivalent(length(out2), 10)
expect_equivalent(length(out3), 10)
expect_equivalent(length(out4), 10)
expect_equivalent(length(out5), 10)
expect_equivalent(length(out6), 10)
expect_equivalent(length(out7), 10)
expect_equivalent(length(out8), 10)
expect_equivalent(length(out9), 10)
