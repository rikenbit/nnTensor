library("testthat")
library("nnTensor")

options(testthat.use_colours = FALSE)

test_file("testthat/test_toyModel.R")
test_file("testthat/test_NMF.R")
test_file("testthat/test_siNMF.R")
test_file("testthat/test_jNMF.R")
test_file("testthat/test_NTF.R")
test_file("testthat/test_NTD.R")
test_file("testthat/test_plotTensor3D.R")
test_file("testthat/test_recTensor.R")
