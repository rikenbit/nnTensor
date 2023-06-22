X <- toyModel("siNMF_Easy")

#
# initW
#
initW <- matrix(runif(nrow(X[[1]])*3),
	nrow=nrow(X[[1]]), ncol=3)

out1 <- siNMF(X, initW=initW,
	J=3, algorithm="Frobenius", num.iter=2)
out2 <- siNMF(X, initW=initW,
	J=3, algorithm="KL", num.iter=2)
out3 <- siNMF(X, initW=initW,
	J=3, algorithm="IS", num.iter=2)
out4 <- siNMF(X, initW=initW,
	J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out1), 6)
expect_equivalent(length(out2), 6)
expect_equivalent(length(out3), 6)
expect_equivalent(length(out4), 6)


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

out5 <- siNMF(X, initH=initH,
	J=3, algorithm="Frobenius", num.iter=2)
out6 <- siNMF(X, initH=initH,
	J=3, algorithm="KL", num.iter=2)
out7 <- siNMF(X, initH=initH,
	J=3, algorithm="IS", num.iter=2)
out8 <- siNMF(X, initH=initH,
	J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out5), 6)
expect_equivalent(length(out6), 6)
expect_equivalent(length(out7), 6)
expect_equivalent(length(out8), 6)
