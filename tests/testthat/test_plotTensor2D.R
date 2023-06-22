tmp <- tempdir()
out <- toyModel(model="NMF")
filename <- paste0(tmp, "/tmp.png")

png(filename=filename)
plotTensor2D(out)
dev.off()

expect_true(file.exists(filename))
