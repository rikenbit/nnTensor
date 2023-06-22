X <- toyModel("NMF")

#
# initU
#
initU <- matrix(runif(nrow(X)*3),
	nrow=nrow(X), ncol=3)

out1 <- NMTF(X, initU=initU, rank=c(3,4), algorithm="Frobenius", num.iter=2)
out2 <- NMTF(X, initU=initU, rank=c(3,4), algorithm="KL", num.iter=2)
out3 <- NMTF(X, initU=initU, rank=c(3,4), algorithm="IS", num.iter=2)
out4 <- NMTF(X, initU=initU, rank=c(3,4), algorithm="ALS", num.iter=2)
out5 <- NMTF(X, initU=initU, rank=c(3,4), algorithm="PG", num.iter=2)
out6 <- NMTF(X, initU=initU, rank=c(3,4), algorithm="COD", num.iter=2)
out7 <- NMTF(X, initU=initU, rank=c(3,4), algorithm="Beta", Beta=-1, num.iter=2)

expect_equivalent(length(out1), 6)
expect_equivalent(length(out2), 6)
expect_equivalent(length(out3), 6)
expect_equivalent(length(out4), 6)
expect_equivalent(length(out5), 6)
expect_equivalent(length(out6), 6)
expect_equivalent(length(out7), 6)

#
# initS
#
initS <- matrix(runif(3*4), nrow=3, ncol=4)

out8 <- NMTF(X, initS=initS, rank=c(3,4), algorithm="Frobenius", num.iter=2)
out9 <- NMTF(X, initS=initS, rank=c(3,4), algorithm="KL", num.iter=2)
out10 <- NMTF(X, initS=initS, rank=c(3,4), algorithm="IS", num.iter=2)
out11 <- NMTF(X, initS=initS, rank=c(3,4), algorithm="ALS", num.iter=2)
out12 <- NMTF(X, initS=initS, rank=c(3,4), algorithm="PG", num.iter=2)
out13 <- NMTF(X, initS=initS, rank=c(3,4), algorithm="COD", num.iter=2)
out14 <- NMTF(X, initS=initS, rank=c(3,4), algorithm="Beta", Beta=-1, num.iter=2)

expect_equivalent(length(out8), 6)
expect_equivalent(length(out9), 6)
expect_equivalent(length(out10), 6)
expect_equivalent(length(out11), 6)
expect_equivalent(length(out12), 6)
expect_equivalent(length(out13), 6)
expect_equivalent(length(out14), 6)

#
# initV
#
initV <- matrix(runif(ncol(X)*4), nrow=ncol(X), ncol=4)

out15 <- NMTF(X, initV=initV, rank=c(3,4), algorithm="Frobenius", num.iter=2)
out16 <- NMTF(X, initV=initV, rank=c(3,4), algorithm="KL", num.iter=2)
out17 <- NMTF(X, initV=initV, rank=c(3,4), algorithm="IS", num.iter=2)
out18 <- NMTF(X, initV=initV, rank=c(3,4), algorithm="ALS", num.iter=2)
out19 <- NMTF(X, initV=initV, rank=c(3,4), algorithm="PG", num.iter=2)
out20 <- NMTF(X, initV=initV, rank=c(3,4), algorithm="COD", num.iter=2)
out21 <- NMTF(X, initV=initV, rank=c(3,4), algorithm="Beta", Beta=-1, num.iter=2)

expect_equivalent(length(out15), 6)
expect_equivalent(length(out16), 6)
expect_equivalent(length(out17), 6)
expect_equivalent(length(out18), 6)
expect_equivalent(length(out19), 6)
expect_equivalent(length(out20), 6)
expect_equivalent(length(out21), 6)
