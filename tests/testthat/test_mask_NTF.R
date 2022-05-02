#
# 3-order tensor
#
X <- toyModel("CP")
M1 <- X
M2 <- X
M1@data[] <- rbinom(length(X@data), 1, 0.5)
M2@data[] <- rbinom(length(X@data), 1, 0.5)

out1_1 <- NTF(X, M=M1, rank=3, algorithm="Frobenius", num.iter=2)
out1_2 <- NTF(X, M=M1, rank=3, algorithm="Frobenius", init="ALS", num.iter=2)
out1_3 <- NTF(X, M=M1, rank=3, algorithm="Frobenius", init="Random", num.iter=2)
out2 <- NTF(X, M=M1, rank=3, algorithm="KL", num.iter=2)
out3 <- NTF(X, M=M1, rank=3, algorithm="IS", num.iter=2)
out4 <- NTF(X, M=M1, rank=3, algorithm="Pearson", num.iter=2)
out5 <- NTF(X, M=M1, rank=3, algorithm="Hellinger", num.iter=2)
out6 <- NTF(X, M=M1, rank=3, algorithm="Neyman", num.iter=2)
out7 <- NTF(X, M=M1, rank=3, algorithm="Alpha", num.iter=2)
out8 <- NTF(X, M=M1, rank=3, algorithm="Beta", num.iter=2)

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



out9_1 <- NTF(X, M=M2, rank=3, algorithm="Frobenius", num.iter=2)
out9_2 <- NTF(X, M=M2, rank=3, algorithm="Frobenius", init="ALS", num.iter=2)
out9_3 <- NTF(X, M=M2, rank=3, algorithm="Frobenius", init="Random", num.iter=2)
out10 <- NTF(X, M=M2, rank=3, algorithm="KL", num.iter=2)
out11 <- NTF(X, M=M2, rank=3, algorithm="IS", num.iter=2)
out12 <- NTF(X, M=M2, rank=3, algorithm="Pearson", num.iter=2)
out13 <- NTF(X, M=M2, rank=3, algorithm="Hellinger", num.iter=2)
out14 <- NTF(X, M=M2, rank=3, algorithm="Neyman", num.iter=2)
out15 <- NTF(X, M=M2, rank=3, algorithm="Alpha", num.iter=2)
out16 <- NTF(X, M=M2, rank=3, algorithm="Beta", num.iter=2)

expect_true(rev(out1_1$TestRecError)[1] != rev(out9_1$TestRecError)[1])
expect_true(rev(out1_2$TestRecError)[1] != rev(out9_2$TestRecError)[1])
expect_true(rev(out1_3$TestRecError)[1] != rev(out9_3$TestRecError)[1])
expect_true(rev(out2$TestRecError)[1] != rev(out10$TestRecError)[1])
expect_true(rev(out3$TestRecError)[1] != rev(out11$TestRecError)[1])
expect_true(rev(out4$TestRecError)[1] != rev(out12$TestRecError)[1])
expect_true(rev(out5$TestRecError)[1] != rev(out13$TestRecError)[1])
expect_true(rev(out6$TestRecError)[1] != rev(out14$TestRecError)[1])
expect_true(rev(out7$TestRecError)[1] != rev(out15$TestRecError)[1])
expect_true(rev(out8$TestRecError)[1] != rev(out16$TestRecError)[1])



out17_1 <- NTF(X, M=M1, rank=4, algorithm="Frobenius", num.iter=2)
out17_2 <- NTF(X, M=M1, rank=4, algorithm="Frobenius", init="ALS", num.iter=2)
out17_3 <- NTF(X, M=M1, rank=4, algorithm="Frobenius", init="Random", num.iter=2)
out18 <- NTF(X, M=M1, rank=4, algorithm="KL", num.iter=2)
out19 <- NTF(X, M=M1, rank=4, algorithm="IS", num.iter=2)
out20 <- NTF(X, M=M1, rank=4, algorithm="Pearson", num.iter=2)
out21 <- NTF(X, M=M1, rank=4, algorithm="Hellinger", num.iter=2)
out22 <- NTF(X, M=M1, rank=4, algorithm="Neyman", num.iter=2)
out23 <- NTF(X, M=M1, rank=4, algorithm="Alpha", num.iter=2)
out24 <- NTF(X, M=M1, rank=4, algorithm="Beta", num.iter=2)

expect_true(rev(out1_1$TestRecError)[1] != rev(out17_1$TestRecError)[1])
# expect_true(rev(out1_2$TestRecError)[1] != rev(out17_2$TestRecError)[1])
expect_true(rev(out1_3$TestRecError)[1] != rev(out17_3$TestRecError)[1])
expect_true(rev(out2$TestRecError)[1] != rev(out18$TestRecError)[1])
expect_true(rev(out3$TestRecError)[1] != rev(out19$TestRecError)[1])
expect_true(rev(out4$TestRecError)[1] != rev(out20$TestRecError)[1])
expect_true(rev(out5$TestRecError)[1] != rev(out21$TestRecError)[1])
expect_true(rev(out6$TestRecError)[1] != rev(out22$TestRecError)[1])
expect_true(rev(out7$TestRecError)[1] != rev(out23$TestRecError)[1])
expect_true(rev(out8$TestRecError)[1] != rev(out24$TestRecError)[1])
