out1 <- toyModel(model="NMF")
out2 <- toyModel(model="CP")
out3 <- toyModel(model="Tucker")

expect_equivalent(length(dim(out1)), 2)
expect_equivalent(length(dim(out2)), 3)
expect_equivalent(length(dim(out3)), 3)
