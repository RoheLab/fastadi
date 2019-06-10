context("test-svd_frob_inner_prod")

test_that("svd_frob_inner_prod", {

  library(Matrix)
  library(RSpectra)
  library(testthat)

  set.seed(27)

  n <- 11
  r <- 3

  A <- rsparsematrix(n, n, nnz = 20)
  s_A <- svds(A, r)

  B <- rsparsematrix(n, n, nnz = 20)
  s_B <- svds(B, r)

  U_A <- s_A$u
  DVt_A <- diag(s_A$d) %*% t(s_A$v)

  U_B <- s_B$u
  DVt_B <- diag(s_B$d) %*% t(s_B$v)

  Z_A <- U_A %*% DVt_A
  Z_B <- U_B %*% DVt_B

  # inner product check

  frob_expected <- sum(Z_A * Z_B)
  frob_impl <- svd_frob_inner_prod_impl(s_A$u, s_A$d, s_A$v, s_B$u, s_B$d, s_B$v)

  expect_equal(frob_impl, frob_expected)

  # frobenius norm check

  norm_expected <- sum(Z_A^2)
  norm_impl <- svd_frob_inner_prod_impl(s_A$u, s_A$d, s_A$v, s_A$u, s_A$d, s_A$v)

  expect_equal(norm_impl, norm_expected)

  # absolute diff in norms check

  expected_diff <- sum((Z_A - Z_B)^2)
  impl_diff <- svd_frob_inner_prod_impl(s_A$u, s_A$d, s_A$v, s_A$u, s_A$d, s_A$v) +
    svd_frob_inner_prod_impl(s_B$u, s_B$d, s_B$v, s_B$u, s_B$d, s_B$v) - 2 *
    svd_frob_inner_prod_impl(s_A$u, s_A$d, s_A$v, s_B$u, s_B$d, s_B$v)

  expect_equal(impl_diff, expected_diff)

  # relative diff in norms check

  expected_rel_diff <- sum((Z_A - Z_B)^2) / sum(Z_B^2)

  rel_diff <- relative_f_norm_change(s_A, s_B)
  expect_equal(rel_diff, expected_rel_diff)
})
