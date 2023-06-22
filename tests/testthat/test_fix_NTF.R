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

# Normalization
initA$A1 <- t(apply(initA$A1, 1, function(a){
	a / norm(as.matrix(a), "F")
}))
initA$A2 <- t(apply(initA$A2, 1, function(a){
	a / norm(as.matrix(a), "F")
}))
initA$A3 <- t(apply(initA$A3, 1, function(a){
	a / norm(as.matrix(a), "F")
}))

out1 <- NTF(X, initA=initA, fixA=c(TRUE,TRUE,TRUE),
	rank=3, algorithm="Frobenius", num.iter=2)
out2 <- NTF(X, initA=initA, fixA=c(FALSE,FALSE,FALSE),
	rank=3, algorithm="Frobenius", num.iter=2)
out3 <- NTF(X, initA=initA, fixA=c(TRUE,TRUE,FALSE),
	rank=3, algorithm="Frobenius", num.iter=2)
out4 <- NTF(X, initA=initA, fixA=c(TRUE,FALSE,TRUE),
	rank=3, algorithm="Frobenius", num.iter=2)
out5 <- NTF(X, initA=initA, fixA=c(FALSE,TRUE,TRUE),
	rank=3, algorithm="Frobenius", num.iter=2)
out6 <- NTF(X, initA=initA, fixA=c(FALSE,FALSE,TRUE),
	rank=3, algorithm="Frobenius", num.iter=2)
out7 <- NTF(X, initA=initA, fixA=c(FALSE,TRUE,FALSE),
	rank=3, algorithm="Frobenius", num.iter=2)
out8 <- NTF(X, initA=initA, fixA=c(TRUE,FALSE,FALSE),
	rank=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out1$A$A1, initA$A1)
expect_equivalent(out1$A$A2, initA$A2)
expect_equivalent(out1$A$A3, initA$A3)

expect_equivalent(out3$A$A1, initA$A1)
expect_equivalent(out3$A$A2, initA$A2)

expect_equivalent(out4$A$A1, initA$A1)
expect_equivalent(out4$A$A3, initA$A3)

expect_equivalent(out5$A$A2, initA$A2)
expect_equivalent(out5$A$A3, initA$A3)

expect_equivalent(out6$A$A3, initA$A3)
expect_equivalent(out7$A$A2, initA$A2)
expect_equivalent(out8$A$A1, initA$A1)
