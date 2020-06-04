#
# 3-order tensor
#
X <- toyModel("Tucker")
M <- X
M@data[] <- rbinom(length(X@data), 1, 0.5)

out1_1 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out1_2 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTD(X, M=M, rank=c(1,2,3), algorithm="KL", num.iter=2)
out3 <- NTD(X, M=M, rank=c(1,2,3), algorithm="IS", num.iter=2)
out4 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Pearson", num.iter=2)
out5 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Hellinger", num.iter=2)
out6 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Neyman", num.iter=2)
out7 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Alpha", num.iter=2)
out8 <- NTD(X, M=M, rank=c(1,2,3), algorithm="Beta", num.iter=2)

expect_equivalent(length(out1_1), 6)
expect_equivalent(length(out1_2), 6)
expect_equivalent(length(out1_3), 6)
expect_equivalent(length(out2), 6)
expect_equivalent(length(out3), 6)
expect_equivalent(length(out4), 6)
expect_equivalent(length(out5), 6)
expect_equivalent(length(out6), 6)
expect_equivalent(length(out7), 6)
expect_equivalent(length(out8), 6)
