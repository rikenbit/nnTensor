X <- toyModel("siNMF_Easy")

out1 <- siNMF(X, J=3, algorithm="Frobenius", num.iter=2)
out2 <- siNMF(X, J=3, algorithm="KL", num.iter=2)
out3 <- siNMF(X, J=3, algorithm="IS", num.iter=2)
out4 <- siNMF(X, J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out1), 4)
expect_equivalent(length(out2), 4)
expect_equivalent(length(out3), 4)
expect_equivalent(length(out4), 4)
