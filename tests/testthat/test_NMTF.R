X <- toyModel("NMF")

out1 <- NMTF(X, rank=c(3,4), algorithm="Frobenius", num.iter=2)
out2 <- NMTF(X, rank=c(3,4), algorithm="KL", num.iter=2)
out3 <- NMTF(X, rank=c(3,4), algorithm="IS", num.iter=2)
out4 <- NMTF(X, rank=c(3,4), algorithm="ALS", num.iter=2)
out5 <- NMTF(X, rank=c(3,4), algorithm="PG", num.iter=2)
out6 <- NMTF(X, rank=c(3,4), algorithm="COD", num.iter=2)
out7 <- NMTF(X, rank=c(3,4), algorithm="Beta", Beta=-1, num.iter=2)

expect_equivalent(length(out1), 6)
expect_equivalent(length(out2), 6)
expect_equivalent(length(out3), 6)
expect_equivalent(length(out4), 6)
expect_equivalent(length(out5), 6)
expect_equivalent(length(out6), 6)
expect_equivalent(length(out7), 6)

# out4 <- NMTF(X, rank=c(4,4), algorithm="ALS", num.iter=100, viz=TRUE)

# initU <- out4$U
# initS <- out4$S
# initV <- out4$V

# outU <- NMTF(X, rank=c(4,4), algorithm="ALS", initU=initU, fixU=TRUE, num.iter=30, viz=TRUE)
# outS <- NMTF(X, rank=c(4,4), algorithm="ALS", initS=initS, fixS=TRUE, num.iter=30, viz=TRUE)
# outV <- NMTF(X, rank=c(4,4), algorithm="ALS", initV=initV, fixV=TRUE, num.iter=30, viz=TRUE)

# outUS <- NMTF(X, rank=c(4,4), algorithm="ALS", initU=initU, initS=initS, fixU=TRUE, fixS=TRUE, num.iter=30, viz=TRUE)
# outUV <- NMTF(X, rank=c(4,4), algorithm="ALS", initU=initU, initV=initV, fixU=TRUE, fixV=TRUE, num.iter=30, viz=TRUE)
# outSV <- NMTF(X, rank=c(4,4), algorithm="ALS", initS=initS, initV=initV, fixS=TRUE, fixV=TRUE, num.iter=30, viz=TRUE)
