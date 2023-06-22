tmp <- tempdir()
out <- toyModel(model="CP")
filename <- paste0(tmp, "/tmp.png")

png(filename=filename)
plotTensor3D(out)
dev.off()

expect_true(file.exists(filename))
