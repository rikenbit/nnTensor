library("testthat")
library("nnTensor")
library("rTensor")

options(testthat.use_colours = FALSE)

# Basic usage
test_file("testthat/test_toyModel.R")
test_file("testthat/test_NMF.R")
test_file("testthat/test_NMTF.R")
test_file("testthat/test_siNMF.R")
test_file("testthat/test_jNMF.R")
test_file("testthat/test_NTF.R")
test_file("testthat/test_NTD.R")
test_file("testthat/test_plotTensor2D.R")
test_file("testthat/test_plotTensor3D.R")
test_file("testthat/test_recTensor.R")
test_file("testthat/test_kFoldMaskTensor.R")

# # Gabriel-type rank estimation for NMF
# test_file("testthat/test_GabrielNMF.R")

# # Rank estimation for NMF
# test_file("testthat/test_NMF_rankestimation.R")

# # Mask matrix/tensor for imputation
# test_file("testthat/test_mask_NMF.R")
# test_file("testthat/test_mask_siNMF.R")
# test_file("testthat/test_mask_jNMF.R")
# test_file("testthat/test_mask_NTF.R")
# test_file("testthat/test_mask_NTD.R")

# # Initialization
# test_file("testthat/test_init_NMF.R")
# test_file("testthat/test_init_siNMF.R")
# test_file("testthat/test_init_jNMF.R")
# test_file("testthat/test_init_NTF.R")
# test_file("testthat/test_init_NTD.R")

# # Fix options for transfer learning
# test_file("testthat/test_fix_NMF.R")
# test_file("testthat/test_fix_siNMF.R")
# test_file("testthat/test_fix_jNMF.R")
# test_file("testthat/test_fix_NTF.R")
# test_file("testthat/test_fix_NTD.R")
