#
# 3-order tensor
#
X <- toyModel("Tucker")

#
# initS (NTD3)
#
initS <- as.tensor(array(runif(1*2*3), dim=c(1,2,3)))

out1_1 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out1_2 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="KL", num.iter=2)
out3 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="IS", num.iter=2)
out4 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Pearson", num.iter=2)
out5 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Hellinger", num.iter=2)
out6 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Neyman", num.iter=2)
out7 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Alpha", num.iter=2)
out8 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="Beta", num.iter=2)
out9 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="HALS", num.iter=2)
out10 <- NTD(X, initS=initS,
	rank=c(1,2,3), algorithm="NMF",
	nmf.algorithm="Frobenius", num.iter=3, num.iter2=2)

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

#
# initS (NTD2)
#
initS2_1 <- as.tensor(
	array(runif(2*3*dim(X)[3]), dim=c(2,3,dim(X)[3])))
initS2_2 <- as.tensor(
	array(runif(dim(X)[1]*3*4), dim=c(dim(X)[1],3,4)))
initS2_3 <- as.tensor(
	array(runif(4*dim(X)[2]*6), dim=c(4,dim(X)[2],6)))

out_NTD2_1 <- NTD(X, initS=initS2_1,
	rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)
out_NTD2_2 <- NTD(X, initS=initS2_2,
	rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)
out_NTD2_3 <- NTD(X, initS=initS2_3,
	rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)

expect_equivalent(length(out_NTD2_1), 6)
expect_equivalent(length(out_NTD2_2), 6)
expect_equivalent(length(out_NTD2_3), 6)


#
# initS (NTD1)
#
initS1_1 <- as.tensor(
	array(runif(3*dim(X)[2]*dim(X)[3]), dim=c(3,dim(X)[2],dim(X)[3])))
initS1_2 <- as.tensor(
	array(runif(dim(X)[1]*4*dim(X)[3]), dim=c(dim(X)[1],4,dim(X)[3])))
initS1_3 <- as.tensor(
	array(runif(dim(X)[1]*dim(X)[2]*5), dim=c(dim(X)[1],dim(X)[1],5)))

out_NTD1_1 <- NTD(X, initS=initS1_1,
	rank=3, modes=1, algorithm="Frobenius", num.iter=2)
out_NTD1_2 <- NTD(X, initS=initS1_2,
	rank=4, modes=2, algorithm="Frobenius", num.iter=2)
out_NTD1_3 <- NTD(X, initS=initS1_3,
	rank=5, modes=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(length(out_NTD1_1), 6)
expect_equivalent(length(out_NTD1_2), 6)
expect_equivalent(length(out_NTD1_3), 6)


#
# initA (NTD3)
#
initA <- list(
	A1=matrix(runif(1*dim(X)[1]), nrow=1, ncol=dim(X)[1]),
	A2=matrix(runif(2*dim(X)[1]), nrow=2, ncol=dim(X)[1]),
	A3=matrix(runif(3*dim(X)[1]), nrow=3, ncol=dim(X)[1])
	)

out1_1 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out1_2 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="KL", num.iter=2)
out3 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="IS", num.iter=2)
out4 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Pearson", num.iter=2)
out5 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Hellinger", num.iter=2)
out6 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Neyman", num.iter=2)
out7 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Alpha", num.iter=2)
out8 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="Beta", num.iter=2)
out9 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="HALS", num.iter=2)
out10 <- NTD(X, initA=initA,
	rank=c(1,2,3), algorithm="NMF",
	nmf.algorithm="Frobenius", num.iter=3, num.iter2=2)

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

#
# initA (NTD2)
#
initA2_1 <- list(
	A1=matrix(runif(2*dim(X)[1]), nrow=2, ncol=dim(X)[1]),
	A2=matrix(runif(3*dim(X)[2]), nrow=3, ncol=dim(X)[2]),
	I1=diag(dim(X)[3])
	)
initA2_2 <- list(
	I1=diag(dim(X)[1]),
	A2=matrix(runif(3*dim(X)[2]), nrow=3, ncol=dim(X)[2]),
	A3=matrix(runif(4*dim(X)[3]), nrow=4, ncol=dim(X)[3])
	)
initA2_3 <- list(
	A1=matrix(runif(4*dim(X)[1]), nrow=4, ncol=dim(X)[1]),
	I2=diag(dim(X)[2]),
	A3=matrix(runif(6*dim(X)[3]), nrow=6, ncol=dim(X)[3])
	)

out_NTD2_1 <- NTD(X, initA=initA2_1,
	rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)
out_NTD2_2 <- NTD(X, initA=initA2_2,
	rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)
out_NTD2_3 <- NTD(X, initA=initA2_3,
	rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)

expect_equivalent(length(out_NTD2_1), 6)
expect_equivalent(length(out_NTD2_2), 6)
expect_equivalent(length(out_NTD2_3), 6)


#
# initA (NTD1)
#
initA1_1 <- list(
	A1=matrix(runif(3*dim(X)[1]), nrow=3, ncol=dim(X)[1]),
	I1=diag(dim(X)[2]),
	I2=diag(dim(X)[3])
	)
initA1_2 <- list(
	I1=diag(dim(X)[1]),
	A2=matrix(runif(4*dim(X)[2]), nrow=4, ncol=dim(X)[2]),
	I2=diag(dim(X)[3])
	)
initA1_3 <- list(
	I1=diag(dim(X)[1]),
	I2=diag(dim(X)[2]),
	A3=matrix(runif(5*dim(X)[3]), nrow=5, ncol=dim(X)[3])
	)

out_NTD1_1 <- NTD(X, initA=initA1_1,
	rank=3, modes=1, algorithm="Frobenius", num.iter=2)
out_NTD1_2 <- NTD(X, initA=initA1_2,
	rank=4, modes=2, algorithm="Frobenius", num.iter=2)
out_NTD1_3 <- NTD(X, initA=initA1_3,
	rank=5, modes=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(length(out_NTD1_1), 6)
expect_equivalent(length(out_NTD1_2), 6)
expect_equivalent(length(out_NTD1_3), 6)
