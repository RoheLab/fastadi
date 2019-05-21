context("test-adaptive-impute-citation")

# check that results agree with dense computation
#
# test_that("sparse citation agrees with dense computation", {
#
#   set.seed(27)
#
#   M <- rsparsematrix(8, 12, nnz = 30)
#
#   M * upper.tri(M)
#
#   sparse <- sparse_adaptive_initialize(M, 5)
#   dense <- dense_adaptive_initialize(M, 5)
#
#   expect_true(equal_svds(sparse, dense))
#
#   expect_equal(2 * 2, 4)
# })
