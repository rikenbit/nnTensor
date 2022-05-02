#
# 3-order tensor
#
X <- toyModel("Tucker")
M1 <- X
M2 <- X
M1@data[] <- rbinom(length(X@data), 1, 0.5)
M2@data[] <- rbinom(length(X@data), 1, 0.5)



out1_1 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out1_2 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="KL", num.iter=2)
out3 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="IS", num.iter=2)
out4 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Pearson", num.iter=2)
out5 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Hellinger", num.iter=2)
out6 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Neyman", num.iter=2)
out7 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Alpha", num.iter=2)
out8 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="Beta", num.iter=2)
out9 <- NTD(X, M=M1, rank=c(1,2,3), algorithm="NMF",
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



out10_1 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Frobenius", num.iter=2)
out10_2 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Frobenius", init="ALS", num.iter=2)
out10_3 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Frobenius", init="Random", num.iter=2)
out11 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="KL", num.iter=2)
out12 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="IS", num.iter=2)
out13 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Pearson", num.iter=2)
out14 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Hellinger", num.iter=2)
out15 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Neyman", num.iter=2)
out16 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Alpha", num.iter=2)
out17 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="Beta", num.iter=2)
out18 <- NTD(X, M=M2, rank=c(1,2,3), algorithm="NMF",
	nmf.algorithm="Frobenius", num.iter=3, num.iter2=2)

expect_true(rev(out1_1$TestRecError)[1] != rev(out10_1$TestRecError)[1])
expect_true(rev(out1_2$TestRecError)[1] != rev(out10_2$TestRecError)[1])
expect_true(rev(out1_3$TestRecError)[1] != rev(out10_3$TestRecError)[1])
expect_true(rev(out2$TestRecError)[1] != rev(out11$TestRecError)[1])
expect_true(rev(out3$TestRecError)[1] != rev(out12$TestRecError)[1])
expect_true(rev(out4$TestRecError)[1] != rev(out13$TestRecError)[1])
expect_true(rev(out5$TestRecError)[1] != rev(out14$TestRecError)[1])
expect_true(rev(out6$TestRecError)[1] != rev(out15$TestRecError)[1])
expect_true(rev(out7$TestRecError)[1] != rev(out16$TestRecError)[1])
expect_true(rev(out8$TestRecError)[1] != rev(out17$TestRecError)[1])
expect_true(rev(out9$TestRecError)[1] != rev(out18$TestRecError)[1])



out19_1 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Frobenius", num.iter=2)
out19_2 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Frobenius", init="ALS", num.iter=2)
out19_3 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Frobenius", init="Random", num.iter=2)
out20 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="KL", num.iter=2)
out21 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="IS", num.iter=2)
out22 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Pearson", num.iter=2)
out23 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Hellinger", num.iter=2)
out24 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Neyman", num.iter=2)
out25 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Alpha", num.iter=2)
out26 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="Beta", num.iter=2)
out27 <- NTD(X, M=M1, rank=c(2,3,1), algorithm="NMF",
	nmf.algorithm="Frobenius", num.iter=3, num.iter2=2)

expect_true(rev(out1_1$TestRecError)[1] != rev(out19_1$TestRecError)[1])
expect_true(rev(out1_2$TestRecError)[1] != rev(out19_2$TestRecError)[1])
expect_true(rev(out1_3$TestRecError)[1] != rev(out19_3$TestRecError)[1])
expect_true(rev(out2$TestRecError)[1] != rev(out20$TestRecError)[1])
expect_true(rev(out3$TestRecError)[1] != rev(out21$TestRecError)[1])
expect_true(rev(out4$TestRecError)[1] != rev(out22$TestRecError)[1])
expect_true(rev(out5$TestRecError)[1] != rev(out23$TestRecError)[1])
expect_true(rev(out6$TestRecError)[1] != rev(out24$TestRecError)[1])
expect_true(rev(out7$TestRecError)[1] != rev(out25$TestRecError)[1])
expect_true(rev(out8$TestRecError)[1] != rev(out26$TestRecError)[1])
expect_true(rev(out9$TestRecError)[1] != rev(out27$TestRecError)[1])
