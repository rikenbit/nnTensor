X <- toyModel("NMF")
J <- 3

out1 <- ZITNMF(X, J=J, Beta=2, num.iter=2, initializer="Random")
out2 <- ZITNMF(X, J=J, Beta=1, num.iter=2, initializer="Random")

expect_equivalent(length(out1), 6)
expect_equivalent(length(out2), 6)

expect_identical(names(out1), c("F", "A", "Z", "w", "RecError", "RelChange"))
expect_identical(names(out2), c("F", "A", "Z", "w", "RecError", "RelChange"))

expect_identical(is.matrix(out1$F), TRUE)
expect_identical(is.matrix(out1$A), TRUE)
expect_identical(is.matrix(out1$Z), TRUE)

expect_equal(dim(out1$F), c(nrow(X), J))
expect_equal(dim(out1$A), c(ncol(X), J))
expect_equal(dim(out1$Z), dim(X))

expect_equal(length(out1$w), 1)

expect_identical(is.vector(out1$RecError), TRUE)
expect_identical(is.vector(out1$RelChange), TRUE)
