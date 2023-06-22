X <- toyModel("siNMF_Easy")
X[[1]][sample(seq(length(X[[1]])), 0.1*length(X[[1]]))] <- NA
X[[2]][sample(seq(length(X[[2]])), 0.1*length(X[[2]]))] <- NA
X[[3]][sample(seq(length(X[[3]])), 0.1*length(X[[3]]))] <- NA

out1 <- siNMF(X, J=3, algorithm="Frobenius", num.iter=2)
out2 <- siNMF(X, J=3, algorithm="KL", num.iter=2)
out3 <- siNMF(X, J=3, algorithm="IS", num.iter=2)
out4 <- siNMF(X, J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out1), 6)
expect_equivalent(length(out2), 6)
expect_equivalent(length(out3), 6)
expect_equivalent(length(out4), 6)
