X <- toyModel("siNMF_Easy")

#
# fixW
#
initW <- matrix(runif(nrow(X[[1]])*3),
	nrow=nrow(X[[1]]), ncol=3)

out1 <- siNMF(X, initW=initW, fixW=TRUE,
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

out2 <- siNMF(X, initH=initH, fixH=c(TRUE, TRUE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out3 <- siNMF(X, initH=initH, fixH=c(FALSE, FALSE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out4 <- siNMF(X, initH=initH, fixH=c(TRUE, TRUE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out5 <- siNMF(X, initH=initH, fixH=c(TRUE, FALSE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out6 <- siNMF(X, initH=initH, fixH=c(FALSE, TRUE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out7 <- siNMF(X, initH=initH, fixH=c(FALSE, FALSE, TRUE),
	J=3, algorithm="Frobenius", num.iter=2)
out8 <- siNMF(X, initH=initH, fixH=c(FALSE, TRUE, FALSE),
	J=3, algorithm="Frobenius", num.iter=2)
out9 <- siNMF(X, initH=initH, fixH=c(TRUE, FALSE, FALSE),
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
