X <- toyModel("siNMF_Hard")

#
# fixW
#
initW <- matrix(runif(nrow(X[[1]])*3),
	nrow=nrow(X[[1]]), ncol=3)

out1 <- jNMF(X, initW=initW, fixW=TRUE,
	J=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out1$W, initW)

#
# fixH
#
initH <- list(
	H1=matrix(runif(ncol(X[[1]])*3),
		nrow=ncol(X[[1]]), ncol=3),
	H2=matrix(runif(ncol(X[[2]])*3),
		nrow=ncol(X[[2]]), ncol=3),
	H3=matrix(runif(ncol(X[[3]])*3),
		nrow=ncol(X[[3]]), ncol=3)
)

out2 <- jNMF(X, initH=initH, fixH=c(TRUE, TRUE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out3 <- jNMF(X, initH=initH, fixH=c(FALSE, FALSE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out4 <- jNMF(X, initH=initH, fixH=c(TRUE, TRUE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out5 <- jNMF(X, initH=initH, fixH=c(TRUE, FALSE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out6 <- jNMF(X, initH=initH, fixH=c(FALSE, TRUE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out7 <- jNMF(X, initH=initH, fixH=c(FALSE, FALSE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out8 <- jNMF(X, initH=initH, fixH=c(FALSE, TRUE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out9 <- jNMF(X, initH=initH, fixH=c(TRUE, FALSE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out2$H$H1, initH$H1)
expect_equivalent(out2$H$H2, initH$H2)
expect_equivalent(out2$H$H3, initH$H3)

expect_equivalent(out4$H$H1, initH$H1)
expect_equivalent(out4$H$H2, initH$H2)

expect_equivalent(out5$H$H1, initH$H1)
expect_equivalent(out5$H$H3, initH$H3)

expect_equivalent(out6$H$H2, initH$H2)
expect_equivalent(out6$H$H3, initH$H3)

expect_equivalent(out7$H$H3, initH$H3)
expect_equivalent(out8$H$H2, initH$H2)
expect_equivalent(out9$H$H1, initH$H1)

#
# fixV
#
initV <- list(
	V1=matrix(runif(nrow(X[[1]])*3),
		nrow=nrow(X[[1]]), ncol=3),
	V2=matrix(runif(nrow(X[[2]])*3),
	nrow=nrow(X[[2]]), ncol=3),
	V3=matrix(runif(nrow(X[[3]])*3),
	nrow=nrow(X[[3]]), ncol=3)
	)

out10 <- jNMF(X, initV=initV, fixV=c(TRUE, TRUE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out11 <- jNMF(X, initV=initV, fixV=c(FALSE, FALSE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out12 <- jNMF(X, initV=initV, fixV=c(TRUE, TRUE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out13 <- jNMF(X, initV=initV, fixV=c(TRUE, FALSE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out14 <- jNMF(X, initV=initV, fixV=c(FALSE, TRUE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out15 <- jNMF(X, initV=initV, fixV=c(FALSE, FALSE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out16 <- jNMF(X, initV=initV, fixV=c(FALSE, TRUE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out17 <- jNMF(X, initV=initV, fixV=c(TRUE, FALSE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)

expect_equivalent(out10$V$V1, initV$V1)
expect_equivalent(out10$V$V2, initV$V2)
expect_equivalent(out10$V$V3, initV$V3)

expect_equivalent(out12$V$V1, initV$V1)
expect_equivalent(out12$V$V2, initV$V2)

expect_equivalent(out13$V$V1, initV$V1)
expect_equivalent(out13$V$V3, initV$V3)

expect_equivalent(out14$V$V2, initV$V2)
expect_equivalent(out14$V$V3, initV$V3)

expect_equivalent(out15$V$V3, initV$V3)
expect_equivalent(out16$V$V2, initV$V2)
expect_equivalent(out17$V$V1, initV$V1)
