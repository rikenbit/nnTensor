X <- toyModel("NMF")

#
# fixU
#
initU <- matrix(runif(nrow(X)*3),
	nrow=nrow(X), ncol=3)

out1 <- NMF(X, initU=initU, fixU=TRUE,
	J=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out1$U, initU)

#
# fixV
#
initV <- matrix(runif(ncol(X)*3),
	nrow=ncol(X), ncol=3)

out2 <- NMF(X, initV=initV, fixV=TRUE,
	J=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out2$V, initV)
