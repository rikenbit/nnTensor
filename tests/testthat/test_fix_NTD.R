#
# 3-order tensor
#
X <- toyModel("Tucker")

#
# fixS (NTD3)
#
initS <- as.tensor(array(runif(1*2*3), dim=c(1,2,3)))

out1_1 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out1_2 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="KL", num.iter=2)
out3 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="IS", num.iter=2)
out4 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Pearson", num.iter=2)
out5 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Hellinger", num.iter=2)
out6 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Neyman", num.iter=2)
out7 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Alpha", num.iter=2)
out8 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="Beta", num.iter=2)
out9 <- NTD(X, initS=initS, fixS=TRUE,
	rank=c(1,2,3), algorithm="HALS", num.iter=2)

expect_equivalent(out1_1$S, initS)
expect_equivalent(out1_2$S, initS)
expect_equivalent(out1_3$S, initS)
expect_equivalent(out2$S, initS)
expect_equivalent(out3$S, initS)
expect_equivalent(out4$S, initS)
expect_equivalent(out5$S, initS)
expect_equivalent(out6$S, initS)
expect_equivalent(out7$S, initS)
expect_equivalent(out8$S, initS)
expect_equivalent(out9$S, initS)

#
# initS (NTD2)
#
initS2_1 <- as.tensor(
	array(runif(2*3*dim(X)[3]), dim=c(2,3,dim(X)[3])))
initS2_2 <- as.tensor(
	array(runif(dim(X)[1]*3*4), dim=c(dim(X)[1],3,4)))
initS2_3 <- as.tensor(
	array(runif(4*dim(X)[2]*6), dim=c(4,dim(X)[2],6)))

out_NTD2_1 <- NTD(X, initS=initS2_1, fixS=TRUE,
	rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)
out_NTD2_2 <- NTD(X, initS=initS2_2, fixS=TRUE,
	rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)
out_NTD2_3 <- NTD(X, initS=initS2_3, fixS=TRUE,
	rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)

expect_equivalent(out_NTD2_1$S, initS2_1)
expect_equivalent(out_NTD2_2$S, initS2_2)
expect_equivalent(out_NTD2_3$S, initS2_3)

#
# initS (NTD1)
#
initS1_1 <- as.tensor(
	array(runif(3*dim(X)[2]*dim(X)[3]), dim=c(3,dim(X)[2],dim(X)[3])))
initS1_2 <- as.tensor(
	array(runif(dim(X)[1]*4*dim(X)[3]), dim=c(dim(X)[1],4,dim(X)[3])))
initS1_3 <- as.tensor(
	array(runif(dim(X)[1]*dim(X)[2]*5), dim=c(dim(X)[1],dim(X)[1],5)))

out_NTD1_1 <- NTD(X, initS=initS1_1, fixS=TRUE,
	rank=3, modes=1, algorithm="Frobenius", num.iter=2)
out_NTD1_2 <- NTD(X, initS=initS1_2, fixS=TRUE,
	rank=4, modes=2, algorithm="Frobenius", num.iter=2)
out_NTD1_3 <- NTD(X, initS=initS1_3, fixS=TRUE,
	rank=5, modes=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out_NTD1_1$S, initS1_1)
expect_equivalent(out_NTD1_2$S, initS1_2)
expect_equivalent(out_NTD1_3$S, initS1_3)


#
# initA (NTD3)
#
initA <- list(
	A1=matrix(runif(1*dim(X)[1]), nrow=1, ncol=dim(X)[1]),
	A2=matrix(runif(2*dim(X)[2]), nrow=2, ncol=dim(X)[2]),
	A3=matrix(runif(3*dim(X)[3]), nrow=3, ncol=dim(X)[3])
	)
# Normalization
initA$A1 <- t(apply(initA$A1, 1, function(x){
	x / norm(as.matrix(x), "F")
}))
initA$A2 <- t(apply(initA$A2, 1, function(x){
	x / norm(as.matrix(x), "F")
}))
initA$A3 <- t(apply(initA$A3, 1, function(x){
	x / norm(as.matrix(x), "F")
}))

out_NTD3_1 <- NTD(X, initA=initA, fixA=c(TRUE, TRUE, TRUE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out_NTD3_2 <- NTD(X, initA=initA, fixA=c(FALSE, FALSE, FALSE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out_NTD3_3 <- NTD(X, initA=initA, fixA=c(TRUE, TRUE, FALSE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out_NTD3_4 <- NTD(X, initA=initA, fixA=c(TRUE, FALSE, TRUE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out_NTD3_5 <- NTD(X, initA=initA, fixA=c(FALSE, TRUE, TRUE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out_NTD3_6 <- NTD(X, initA=initA, fixA=c(FALSE, FALSE, TRUE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out_NTD3_7 <- NTD(X, initA=initA, fixA=c(FALSE, TRUE, FALSE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out_NTD3_8 <- NTD(X, initA=initA, fixA=c(TRUE, FALSE, FALSE),
	rank=c(1,2,3), algorithm="Frobenius", num.iter=2)

expect_equivalent(out_NTD3_1$A$A1, initA$A1)
expect_equivalent(out_NTD3_1$A$A2, initA$A2)
expect_equivalent(out_NTD3_1$A$A3, initA$A3)

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
	I1=diag(dim(X)[2]),
	A3=matrix(runif(6*dim(X)[3]), nrow=6, ncol=dim(X)[3])
	)

# Normalization
initA2_1$A1 <- t(apply(initA2_1$A1, 1, function(x){
	x / norm(as.matrix(x), "F")
}))
initA2_1$A2 <- t(apply(initA2_1$A2, 1, function(x){
	x / norm(as.matrix(x), "F")
}))
initA2_2$A2 <- t(apply(initA2_2$A2, 1, function(x){
	x / norm(as.matrix(x), "F")
}))
initA2_2$A3 <- t(apply(initA2_2$A3, 1, function(x){
	x / norm(as.matrix(x), "F")
}))
initA2_3$A1 <- t(apply(initA2_3$A1, 1, function(x){
	x / norm(as.matrix(x), "F")
}))
initA2_3$A3 <- t(apply(initA2_3$A3, 1, function(x){
	x / norm(as.matrix(x), "F")
}))

out_NTD2_1_1 <- NTD(X, initA=initA2_1, fixA=c(TRUE, TRUE),
	rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)
out_NTD2_1_2 <- NTD(X, initA=initA2_1, fixA=c(TRUE, FALSE),
	rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)
out_NTD2_1_3 <- NTD(X, initA=initA2_1, fixA=c(FALSE, TRUE),
	rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)
out_NTD2_1_4 <- NTD(X, initA=initA2_1, fixA=c(FALSE, FALSE),
	rank=c(2,3), modes=1:2, algorithm="Frobenius", num.iter=2)

out_NTD2_2_1 <- NTD(X, initA=initA2_2, fixA=c(TRUE, TRUE),
	rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)
out_NTD2_2_2 <- NTD(X, initA=initA2_2, fixA=c(TRUE, FALSE),
	rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)
out_NTD2_2_3 <- NTD(X, initA=initA2_2, fixA=c(FALSE, TRUE),
	rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)
out_NTD2_2_4 <- NTD(X, initA=initA2_2, fixA=c(FALSE, FALSE),
	rank=c(3,4), modes=2:3, algorithm="Frobenius", num.iter=2)

out_NTD2_3_1 <- NTD(X, initA=initA2_3, fixA=c(TRUE, TRUE),
	rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)
out_NTD2_3_2 <- NTD(X, initA=initA2_3, fixA=c(TRUE, FALSE),
	rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)
out_NTD2_3_3 <- NTD(X, initA=initA2_3, fixA=c(FALSE, TRUE),
	rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)
out_NTD2_3_4 <- NTD(X, initA=initA2_3, fixA=c(FALSE, FALSE),
	rank=c(4,6), modes=c(1,3), algorithm="Frobenius", num.iter=2)

expect_equivalent(out_NTD2_1_1$A$A1, initA2_1$A1)
expect_equivalent(out_NTD2_1_1$A$A2, initA2_1$A2)
expect_equivalent(out_NTD2_1_2$A$A1, initA2_1$A1)
expect_equivalent(out_NTD2_1_3$A$A2, initA2_1$A2)

expect_equivalent(out_NTD2_2_1$A$A2, initA2_2$A2)
expect_equivalent(out_NTD2_2_1$A$A3, initA2_2$A3)
expect_equivalent(out_NTD2_2_2$A$A2, initA2_2$A2)
expect_equivalent(out_NTD2_2_3$A$A3, initA2_2$A3)

expect_equivalent(out_NTD2_3_1$A$A1, initA2_3$A1)
expect_equivalent(out_NTD2_3_1$A$A3, initA2_3$A3)
expect_equivalent(out_NTD2_3_2$A$A1, initA2_3$A1)
expect_equivalent(out_NTD2_3_3$A$A3, initA2_3$A3)



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
# Normalizatison
initA1_1$A1 <- t(apply(initA1_1$A1, 1, function(x){
	x/norm(as.matrix(x), "F")
}))
initA1_2$A2 <- t(apply(initA1_2$A2, 1, function(x){
	x/norm(as.matrix(x), "F")
}))
initA1_3$A3 <- t(apply(initA1_3$A3, 1, function(x){
	x/norm(as.matrix(x), "F")
}))

out_NTD1_1_1 <- NTD(X, initA=initA1_1, fixA=TRUE,
	rank=3, modes=1, algorithm="Frobenius", num.iter=2)
out_NTD1_1_2 <- NTD(X, initA=initA1_1, fixA=FALSE,
	rank=3, modes=1, algorithm="Frobenius", num.iter=2)

out_NTD1_2_1 <- NTD(X, initA=initA1_2, fixA=TRUE,
	rank=4, modes=2, algorithm="Frobenius", num.iter=2)
out_NTD1_2_2 <- NTD(X, initA=initA1_2, fixA=FALSE,
	rank=4, modes=2, algorithm="Frobenius", num.iter=2)

out_NTD1_3_1 <- NTD(X, initA=initA1_3, fixA=TRUE,
	rank=5, modes=3, algorithm="Frobenius", num.iter=2)
out_NTD1_3_2 <- NTD(X, initA=initA1_3, fixA=FALSE,
	rank=5, modes=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out_NTD1_1_1$A$A1, initA1_1$A1)
expect_equivalent(out_NTD1_2_1$A$A2, initA1_2$A2)
expect_equivalent(out_NTD1_3_1$A$A3, initA1_3$A3)
