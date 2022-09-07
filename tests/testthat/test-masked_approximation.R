library(testthat)

skip_if_not_installed("Matrix")
skip_if_not_installed("RSpectra")

library(Matrix)
library(RSpectra)

test_that("masked_approximation_impl (small sparse matrix)", {

  set.seed(27)

  # matrix to build an SVD from
  A <- rsparsematrix(30, 40, nnz = 50)

  s <- svds(A, 5)

  # expand the SVD only at observed values by hand

  Z <- s$u %*% diag(s$d) %*% t(s$v)  # full SVD
  Y <- A != 0                        # observed indicator

  expected <- Z * Y                  # project SVD onto observed matrix

  A_triplet <- as(A, "TsparseMatrix")

  impl_result <- masked_approximation_impl(s$u, tcrossprod(s$v, diag(s$d)),
                                           A_triplet@i, A_triplet@j)

  expect_equal(impl_result, expected)
})

test_that("masked_approximation_impl (4 billion+ element sparse matrix)", {

  # need to check that this will work for sparse matrices with more
  # elements than 32-bit integers, which will happen often for
  # large sparse social networks

  # this is related to https://github.com/RcppCore/RcppArmadillo/pull/90
  # and https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define and came up during the
  # web of science applied work

  set.seed(27)

  n <- sqrt(.Machine$integer.max) + 1000

  # matrix to build an SVD from
  A <- rsparsematrix(n, n, nnz = 2000)

  # make sure the
  stopifnot(prod(dim(A)) > .Machine$integer.max)

  s <- svds(A, 5)

  A_triplet <- as(A, "TsparseMatrix")

  expect_silent(
    masked_approximation_impl(
      s$u, tcrossprod(s$v, diag(s$d)),
      A_triplet@i, A_triplet@j
    )
  )
})
