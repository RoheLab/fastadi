context("masked approximation")

test_that("masked_approximation", {

  library(Matrix)
  library(RSpectra)

  set.seed(27)

  # matrix to build an SVD from
  A <- rsparsematrix(30, 40, nnz = 50)

  s <- svds(A, 5)

  # expand the SVD only at observed values by hand

  Z <- s$u %*% diag(s$d) %*% t(s$v)  # full SVD
  Y <- A != 0                        # observed indicator

  expected <- Z * Y                  # project SVD onto observed matrix

  A_triplet <- as(A, "dgTMatrix")

  impl_result <- masked_approximation_impl(s$u, s$d, s$v,
                                           A_triplet@i, A_triplet@j)

  wrapped_result <- masked_approximation(s, A)

  expect_equal(impl_result, expected)
  expect_equal(wrapped_result, expected)
})
