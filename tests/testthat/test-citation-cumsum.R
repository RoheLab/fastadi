library(Matrix)
library(RSpectra)
library(testthat)

# M must be square of you'll get dimension errors
M <- rsparsematrix(11, 11, nnz = 20)
x <- rnorm(11)
s <- svds(M, 3)

test_that("p_u_zx_impl", {

  set.seed(27)

  # NOTE: this currently works with square and wide matrices
  # it could probably be extended to long matrices but that would
  # take more work

  Y <- as(upper.tri(M), "CsparseMatrix")

  expected <- drop((s$u %*% diag(s$d) %*% t(s$v) * Y) %*% x)
  impl_result <- p_u_zx_impl(s$u, s$d, s$v, x, 1L)

  expect_equal(
    drop(impl_result),
    drop(expected)
  )
})

test_that("p_u_tilde_zx_impl", {
  # lower triangular non-zero mask
  L <- M & lower.tri(M)

  # project
  lower_expected <- (s$u %*% diag(s$d) %*% t(s$v) * L) %*% x

  mask <- as(L, "TsparseMatrix")
  lower_impl <- p_u_tilde_zx_impl(s$u, s$d, s$v, mask@i, mask@j, x, 1L)

  expect_equal(
    drop(lower_impl),
    drop(lower_expected)
  )
})

test_that("p_u_ztx_impl", {

  Y <- as(upper.tri(M), "CsparseMatrix")

  expected <- drop(t(s$u %*% diag(s$d) %*% t(s$v) * Y) %*% x)
  impl_result <- p_u_ztx_impl(s$u, s$d, s$v, x, 1L)

  expect_equal(
    drop(impl_result),
    drop(expected)
  )
})

test_that("p_u_tilde_ztx_impl", {
  # lower triangular non-zero mask
  L <- M & lower.tri(M)

  # project
  lower_expected <- t(s$u %*% diag(s$d) %*% t(s$v) * L) %*% x

  mask <- as(L, "TsparseMatrix")
  lower_impl <- p_u_tilde_ztx_impl(s$u, s$d, s$v, mask@i, mask@j, x, 1L)

  expect_equal(
    drop(lower_impl),
    drop(lower_expected)
  )
})

rm(M, x, s)
