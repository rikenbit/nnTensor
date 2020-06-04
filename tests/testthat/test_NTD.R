#
# 3-order tensor
#
X <- toyModel("Tucker")

out1_1 <- NTD(X, rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out1_2 <- NTD(X, rank=c(1,2,3), algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTD(X, rank=c(1,2,3), algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTD(X, rank=c(1,2,3), algorithm="KL", num.iter=2)
out3 <- NTD(X, rank=c(1,2,3), algorithm="IS", num.iter=2)
out4 <- NTD(X, rank=c(1,2,3), algorithm="Pearson", num.iter=2)
out5 <- NTD(X, rank=c(1,2,3), algorithm="Hellinger", num.iter=2)
out6 <- NTD(X, rank=c(1,2,3), algorithm="Neyman", num.iter=2)
out7 <- NTD(X, rank=c(1,2,3), algorithm="Alpha", num.iter=2)
out8 <- NTD(X, rank=c(1,2,3), algorithm="Beta", num.iter=2)
out9 <- NTD(X, rank=c(1,2,3), algorithm="HALS", num.iter=2)

out_NTD2_1 <- NTD(X, rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)
out_NTD2_2 <- NTD(X, rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)
out_NTD2_3 <- NTD(X, rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)

out_NTD1_1 <- NTD(X, rank=3, modes=1, algorithm="Frobenius", num.iter=2)
out_NTD1_2 <- NTD(X, rank=4, modes=2, algorithm="Frobenius", num.iter=2)
out_NTD1_3 <- NTD(X, rank=5, modes=3, algorithm="Frobenius", num.iter=2)

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
expect_equivalent(length(out9), 6)

expect_equivalent(length(out_NTD2_1), 6)
expect_equivalent(length(out_NTD2_2), 6)
expect_equivalent(length(out_NTD2_3), 6)

expect_equivalent(length(out_NTD1_1), 6)
expect_equivalent(length(out_NTD1_2), 6)
expect_equivalent(length(out_NTD1_3), 6)

#
# 4-order tensor
#
# library("rTensor")
# XX <- array(0, dim=c(50,50,50,4))
# XX[,,,1] <- X@data
# XX[,,,2] <- X@data * runif(length(X@data))
# XX[,,,3] <- X@data * 1.2 * runif(length(X@data))
# XX[,,,4] <- X@data * 3 * runif(length(X@data))
# XX <- as.tensor(XX)

# out_NTD4_1_1 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Frobenius", num.iter=2)
# out_NTD4_1_2 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Frobenius", init="ALS", num.iter=2)
# out_NTD4_1_3 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Frobenius", init="Random", num.iter=2)
# out_NTD4_2 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="KL", num.iter=2)
# out_NTD4_3 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="IS", num.iter=2)
# out_NTD4_4 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Pearson", num.iter=2)
# out_NTD4_5 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Hellinger", num.iter=2)
# out_NTD4_6 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Neyman", num.iter=2)
# out_NTD4_7 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Alpha", num.iter=2)
# out_NTD4_8 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="Beta", num.iter=2)
# out_NTD4_9 <- NTD(XX, rank=c(6,5,4,3), modes=1:4, algorithm="HALS", num.iter=2)

# out_NTD4_3_1 <- NTD(XX, rank=c(1,2,3), modes=2:4, algorithm="Frobenius", num.iter=2)
# out_NTD4_3_2 <- NTD(XX, rank=c(1,2,3), modes=c(1,3:4), algorithm="Frobenius", num.iter=2)
# out_NTD4_3_3 <- NTD(XX, rank=c(1,2,3), modes=c(1:2,4), algorithm="Frobenius", num.iter=2)
# out_NTD4_3_4 <- NTD(XX, rank=c(1,2,3), modes=1:3, algorithm="Frobenius", num.iter=2)

# out_NTD4_2_1 <- NTD(XX, rank=c(3,3), modes=3:4, algorithm="Frobenius", num.iter=2)
# out_NTD4_2_2 <- NTD(XX, rank=c(3,3), modes=c(2,4), algorithm="Frobenius", num.iter=2)
# out_NTD4_2_3 <- NTD(XX, rank=c(3,3), modes=2:3, algorithm="Frobenius", num.iter=2)
# out_NTD4_2_4 <- NTD(XX, rank=c(3,3), modes=c(1,4), algorithm="Frobenius", num.iter=2)
# out_NTD4_2_5 <- NTD(XX, rank=c(3,3), modes=c(1,3), algorithm="Frobenius", num.iter=2)
# out_NTD4_2_6 <- NTD(XX, rank=c(3,3), modes=c(1,2), algorithm="Frobenius", num.iter=2)

# out_NTD4_1_1 <- NTD(XX, rank=3, modes=1, algorithm="Frobenius", num.iter=2)
# out_NTD4_1_2 <- NTD(XX, rank=3, modes=2, algorithm="Frobenius", num.iter=2)
# out_NTD4_1_3 <- NTD(XX, rank=3, modes=3, algorithm="Frobenius", num.iter=2)
# # Too heavy
# # out_NTD4_1_4 <- NTD(XX, rank=3, modes=4, algorithm="Frobenius", num.iter=2)

# expect_equivalent(length(out_NTD4_1_1), 6)
# expect_equivalent(length(out_NTD4_1_2), 6)
# expect_equivalent(length(out_NTD4_1_3), 6)
# expect_equivalent(length(out_NTD4_2), 6)
# expect_equivalent(length(out_NTD4_3), 6)
# expect_equivalent(length(out_NTD4_4), 6)
# expect_equivalent(length(out_NTD4_5), 6)
# expect_equivalent(length(out_NTD4_6), 6)
# expect_equivalent(length(out_NTD4_7), 6)
# expect_equivalent(length(out_NTD4_8), 6)
# expect_equivalent(length(out_NTD4_9), 6)

# expect_equivalent(length(out_NTD4_3_1), 6)
# expect_equivalent(length(out_NTD4_3_2), 6)
# expect_equivalent(length(out_NTD4_3_3), 6)
# expect_equivalent(length(out_NTD4_3_4), 6)

# expect_equivalent(length(out_NTD4_2_1), 6)
# expect_equivalent(length(out_NTD4_2_2), 6)
# expect_equivalent(length(out_NTD4_2_3), 6)
# expect_equivalent(length(out_NTD4_2_4), 6)
# expect_equivalent(length(out_NTD4_2_5), 6)
# expect_equivalent(length(out_NTD4_2_6), 6)

# expect_equivalent(length(out_NTD4_1_1), 6)
# expect_equivalent(length(out_NTD4_1_2), 6)
# expect_equivalent(length(out_NTD4_1_3), 6)
# # Too heavy
# # expect_equivalent(length(out_NTD4_1_4), 6)
