X <- toyModel("NMF")
X[sample(seq(length(X)), 0.1*length(X))] <- NA

out1 <- NMTF(X, rank=c(3,4), algorithm="Frobenius", num.iter=2)
out2 <- NMTF(X, rank=c(3,4), algorithm="KL", num.iter=2)
out3 <- NMTF(X, rank=c(3,4), algorithm="IS", num.iter=2)
out4 <- NMTF(X, rank=c(3,4), algorithm="ALS", num.iter=2)
out5 <- NMTF(X, rank=c(3,4), algorithm="PG", num.iter=2)
out6 <- NMTF(X, rank=c(3,4), algorithm="COD", num.iter=2)
out7 <- NMTF(X, rank=c(3,4), algorithm="Beta", Beta=-1, num.iter=2)

expect_equivalent(length(out1), 8)
expect_equivalent(length(out2), 8)
expect_equivalent(length(out3), 8)
expect_equivalent(length(out4), 8)
expect_equivalent(length(out5), 8)
expect_equivalent(length(out6), 8)
expect_equivalent(length(out7), 8)
