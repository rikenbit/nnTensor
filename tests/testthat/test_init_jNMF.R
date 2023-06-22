X <- toyModel("siNMF_Hard")

#
# initW
#
initW <- matrix(runif(nrow(X[[1]])*3),
	nrow=nrow(X[[1]]), ncol=3)
out1 <- jNMF(X, initW=initW,
	J=3, algorithm="Frobenius", num.iter=2)
out2 <- jNMF(X, initW=initW,
	J=3, algorithm="KL", num.iter=2)
out3 <- jNMF(X, initW=initW,
	J=3, algorithm="IS", num.iter=2)
out4 <- jNMF(X, initW=initW,
	J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out1), 7)
expect_equivalent(length(out2), 7)
expect_equivalent(length(out3), 7)
expect_equivalent(length(out4), 7)

#
# initH
#
initH <- list(
	H1=matrix(runif(ncol(X[[1]])*3),
		nrow=ncol(X[[1]]), ncol=3),
	H2=matrix(runif(ncol(X[[2]])*3),
		nrow=ncol(X[[2]]), ncol=3),
	H3=matrix(runif(ncol(X[[3]])*3),
		nrow=ncol(X[[3]]), ncol=3)
)
out5 <- jNMF(X, initH=initH,
	J=3, algorithm="Frobenius", num.iter=2)
out6 <- jNMF(X, initH=initH,
	J=3, algorithm="KL", num.iter=2)
out7 <- jNMF(X, initH=initH,
	J=3, algorithm="IS", num.iter=2)
out8 <- jNMF(X, initH=initH,
	J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out5), 7)
expect_equivalent(length(out6), 7)
expect_equivalent(length(out7), 7)
expect_equivalent(length(out8), 7)

#
# initV
#
initV <- list(
	V1=matrix(runif(nrow(X[[1]])*3),
		nrow=nrow(X[[1]]), ncol=3),
	V2=matrix(runif(nrow(X[[2]])*3),
	nrow=nrow(X[[2]]), ncol=3),
	V3=matrix(runif(nrow(X[[3]])*3),
	nrow=nrow(X[[3]]), ncol=3)
	)
out9 <- jNMF(X, initV=initV,
	J=3, algorithm="Frobenius", num.iter=2)
out10 <- jNMF(X, initV=initV,
	J=3, algorithm="KL", num.iter=2)
out11 <- jNMF(X, initV=initV,
	J=3, algorithm="IS", num.iter=2)
out12 <- jNMF(X, initV=initV,
	J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out9), 7)
expect_equivalent(length(out10), 7)
expect_equivalent(length(out11), 7)
expect_equivalent(length(out12), 7)
