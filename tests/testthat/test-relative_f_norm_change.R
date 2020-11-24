library(Matrix)
library(RSpectra)

test_that("relative_f_norm_change", {

  set.seed(27)

  A <- rsparsematrix(30, 40, nnz = 50)
  B <- rsparsematrix(30, 40, nnz = 50)

  s_new <- svds(A, 5)
  s <- svds(B, 5)

  # dense calculation
  Z_A <- s_new$u %*% diag(s_new$d) %*% t(s_new$v)
  Z_B <- s$u %*% diag(s$d) %*% t(s$v)

  expected_result <- norm(Z_A - Z_B, type = "F")^2 / norm(Z_B, type = "F")^2

  # Armadillo calculation
  impl_result <- relative_f_norm_change_impl(s_new$u, s_new$d, s_new$v,
                                             s$u, s$d, s$v, 1L)

  # wrapped Armadillo calculation
  wrapped_result <- relative_f_norm_change(s_new, s)

  expect_equal(impl_result, expected_result)
  expect_equal(wrapped_result, expected_result)
})
