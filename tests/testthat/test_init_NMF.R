X <- toyModel("NMF")

#
# initU
#
initU <- matrix(runif(nrow(X)*3),
	nrow=nrow(X), ncol=3)

out1 <- NMF(X, initU=initU,
	J=3, algorithm="Frobenius", num.iter=2)
out2 <- NMF(X, initU=initU,
	J=3, algorithm="KL", num.iter=2)
out3 <- NMF(X, initU=initU,
	J=3, algorithm="IS", num.iter=2)
out4 <- NMF(X, initU=initU,
	J=3, algorithm="Pearson", num.iter=2)
out5 <- NMF(X, initU=initU,
	J=3, algorithm="Hellinger", num.iter=2)
out6 <- NMF(X, initU=initU,
	J=3, algorithm="Neyman", num.iter=2)
out7 <- NMF(X, initU=initU,
	J=3, algorithm="Alpha", num.iter=2)
out8 <- NMF(X, initU=initU,
	J=3, algorithm="Beta", num.iter=2)
out9 <- NMF(X, initU=initU,
	J=3, algorithm="HALS", num.iter=2)
out10 <- NMF(X, initU=initU,
	J=3, algorithm="PGD", num.iter=2)
out11 <- NMF(X, initU=initU,
	J=3, algorithm="GCD", num.iter=2)
out12 <- NMF(X, initU=initU,
	J=3, algorithm="Projected", num.iter=2)
out13 <- NMF(X, initU=initU,
	J=3, algorithm="NHR", num.iter=2)
out14 <- NMF(X, initU=initU,
	J=3, algorithm="DTPP", num.iter=2)
out15 <- NMF(X, initU=initU,
	J=3, algorithm="Orthogonal", num.iter=2)
out16 <- NMF(X, initU=initU,
	J=3, algorithm="OrthReg", num.iter=2)

expect_equivalent(length(out1), 10)
expect_equivalent(length(out2), 10)
expect_equivalent(length(out3), 10)
expect_equivalent(length(out4), 10)
expect_equivalent(length(out5), 10)
expect_equivalent(length(out6), 10)
expect_equivalent(length(out7), 10)
expect_equivalent(length(out8), 10)
expect_equivalent(length(out9), 10)
expect_equivalent(length(out10), 10)
expect_equivalent(length(out11), 10)
expect_equivalent(length(out12), 10)
expect_equivalent(length(out13), 10)
expect_equivalent(length(out14), 10)
expect_equivalent(length(out15), 10)
expect_equivalent(length(out16), 10)

#
# initV
#
initV <- matrix(runif(ncol(X)*3),
	nrow=ncol(X), ncol=3)

out12 <- NMF(X, initV=initV,
	J=3, algorithm="Frobenius", num.iter=2)
out13 <- NMF(X, initV=initV,
	J=3, algorithm="KL", num.iter=2)
out14 <- NMF(X, initV=initV,
	J=3, algorithm="IS", num.iter=2)
out15 <- NMF(X, initV=initV,
	J=3, algorithm="Pearson", num.iter=2)
out16 <- NMF(X, initV=initV,
	J=3, algorithm="Hellinger", num.iter=2)
out17 <- NMF(X, initV=initV,
	J=3, algorithm="Neyman", num.iter=2)
out18 <- NMF(X, initV=initV,
	J=3, algorithm="Alpha", num.iter=2)
out19 <- NMF(X, initV=initV,
	J=3, algorithm="Beta", num.iter=2)
out20 <- NMF(X, initV=initV,
	J=3, algorithm="HALS", num.iter=2)
out21 <- NMF(X, initV=initV,
	J=3, algorithm="PGD", num.iter=2)
out22 <- NMF(X, initV=initV,
	J=3, algorithm="GCD", num.iter=2)
out23 <- NMF(X, initV=initV,
	J=3, algorithm="Projected", num.iter=2)
out24 <- NMF(X, initV=initV,
	J=3, algorithm="NHR", num.iter=2)
out25 <- NMF(X, initV=initV,
	J=3, algorithm="DTPP", num.iter=2)
out26 <- NMF(X, initV=initV,
	J=3, algorithm="Orthogonal", num.iter=2)
out27 <- NMF(X, initV=initV,
	J=3, algorithm="OrthReg", num.iter=2)

expect_equivalent(length(out12), 10)
expect_equivalent(length(out13), 10)
expect_equivalent(length(out14), 10)
expect_equivalent(length(out15), 10)
expect_equivalent(length(out16), 10)
expect_equivalent(length(out17), 10)
expect_equivalent(length(out18), 10)
expect_equivalent(length(out19), 10)
expect_equivalent(length(out20), 10)
expect_equivalent(length(out21), 10)
expect_equivalent(length(out22), 10)
expect_equivalent(length(out23), 10)
expect_equivalent(length(out24), 10)
expect_equivalent(length(out25), 10)
expect_equivalent(length(out26), 10)
expect_equivalent(length(out27), 10)
