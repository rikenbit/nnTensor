X <- toyModel("siNMF_Easy")
M <- X
for(i in seq(length(X))){
	M[[i]][] <- 0
	M[[i]][] <- rbinom(length(M[[i]]), 1, 0.5)
}

out1 <- siNMF(X, M=M, J=3, algorithm="Frobenius", num.iter=2)
out2 <- siNMF(X, M=M, J=3, algorithm="KL", num.iter=2)
out3 <- siNMF(X, M=M, J=3, algorithm="IS", num.iter=2)
out4 <- siNMF(X, M=M, J=3, algorithm="PLTF", p=1, num.iter=2)

expect_equivalent(length(out1), 6)
expect_equivalent(length(out2), 6)
expect_equivalent(length(out3), 6)
expect_equivalent(length(out4), 6)
