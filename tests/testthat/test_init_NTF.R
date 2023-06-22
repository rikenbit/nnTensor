#
# 3-order tensor
#
X <- toyModel("CP")

initA <- list(
	A1=matrix(runif(3*dim(X)[1]),
		nrow=3, ncol=dim(X)[1]),
	A2=matrix(runif(3*dim(X)[2]),
		nrow=3, ncol=dim(X)[2]),
	A3=matrix(runif(3*dim(X)[3]),
		nrow=3, ncol=dim(X)[3])
	)

out1_1 <- NTF(X, initA=initA,
	rank=3, algorithm="Frobenius", num.iter=2)
out1_2 <- NTF(X, initA=initA,
	rank=3, algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTF(X, initA=initA,
	rank=3, algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTF(X, initA=initA,
	rank=3, algorithm="KL", num.iter=2)
out3 <- NTF(X, initA=initA,
	rank=3, algorithm="IS", num.iter=2)
out4 <- NTF(X, initA=initA,
	rank=3, algorithm="Pearson", num.iter=2)
out5 <- NTF(X, initA=initA,
	rank=3, algorithm="Hellinger", num.iter=2)
out6 <- NTF(X, initA=initA,
	rank=3, algorithm="Neyman", num.iter=2)
out7 <- NTF(X, initA=initA,
	rank=3, algorithm="Alpha", num.iter=2)
out8 <- NTF(X, initA=initA,
	rank=3, algorithm="Beta", num.iter=2)
out9 <- NTF(X, initA=initA,
	rank=3, algorithm="HALS", num.iter=2)
out10 <- NTF(X, initA=initA,
	rank=3, algorithm="Alpha-HALS", num.iter=2)
out11 <- NTF(X, initA=initA,
	rank=3, algorithm="Beta-HALS", num.iter=2)

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
expect_equivalent(length(out10), 6)
expect_equivalent(length(out11), 6)
