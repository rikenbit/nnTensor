X <- toyModel("NMF")

#
# fixU
#
initU <- matrix(runif(nrow(X)*3),
	nrow=nrow(X), ncol=3)

out1 <- NMTF(X, initU=initU, fixU=TRUE,
	rank=c(3,4), algorithm="Frobenius", num.iter=2)

expect_equivalent(out1$U, initU)

#
# fixV
#
initV <- matrix(runif(ncol(X)*4),
	nrow=ncol(X), ncol=4)

out2 <- NMTF(X, initV=initV, fixV=TRUE,
	rank=c(3,4), algorithm="Frobenius", num.iter=2)

expect_equivalent(out2$V, initV)
