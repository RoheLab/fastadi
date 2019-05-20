context("citation helpers")

test_that("p_omega_f_norm_ut", {

  M <- Matrix(
    rbind(
      c(0, 0, 3, 1, 0),
      c(3, 0, 0, 8, 0),
      c(0, -1, 0, 0, 0),
      c(0, 0, 0, 0, 0),
      c(0, 2, 0, 0, 0),
      c(5, 0, 7, 0, 4)
    )
  )

  # NOTE: not treating nodes as citing themselves!
  Y <- rbind(
    c(0, 1, 1, 1, 1),
    c(1, 0, 1, 1, 1),
    c(0, 1, 0, 1, 1),
    c(0, 0, 0, 0, 1),
    c(0, 1, 0, 0, 0),
    c(1, 0, 1, 0, 1)
  )

  s <- svds(M, 2)

  # triplet form
  # compressed column matrix form even better but don't
  # understand the format
  Y <- as(Y, "TsparseMatrix")

  Z <- s$u %*% diag(s$d) %*% t(s$v)
  approximation <- Z * Y

  expected <- sum(approximation@x^2)

  impl_result <- p_omega_f_norm_ut_impl(s$u, s$d, s$v, Y@i, Y@j)
  wrapped_result <- p_omega_f_norm_ut(s, Y)

  expect_equal(impl_result, expected)
  expect_equal(wrapped_result, expected)
})


test_that("Ax_citation", {

  library(Matrix)
  library(RSpectra)

  set.seed(27)

  # matrix to build an SVD from
  A <- rsparsematrix(30, 40, nnz = 50)

  s <- svds(A, 5)

  # vector to multiply by
  x <- rnorm(40)

  # build a residual matrix with matching missingness
  R <- A
  R@x <- rnorm(50)

  expect_false(isTRUE(all.equal(A, R)))
  expect_equal(A != 0, R != 0)

  # expand the SVD only at observed values by hand

  Z <- s$u %*% diag(s$d) %*% t(s$v)  # full SVD
  Y <- A != 0                        # observed indicator

  # add the upper triangle back to Y since we're in the citation
  # setting
  Y <- Y | upper.tri(Y)

  obs_Z <- Z * Y                     # project SVD onto observed matrix

  A0 <- drop0(A)
  expected <- A0 %*% x - obs_Z %*% x + Z %*% x

  A <- as(A, "TsparseMatrix")
  impl_result <- p_omega_zx_impl(s$u, s$d, s$v, A@i, A@j, x)

  expect_equal(
    drop(impl_result),
    drop(obs_Z %*% x)
  )

  args <- list(u = s$u, d = s$d, v = s$v, M = A)
  wrapped_result <- Ax_citation(x, args)

  expect_equal(wrapped_result, expected)
})

test_that("Atx_citation", {

  library(Matrix)
  library(RSpectra)

  set.seed(27)

  # matrix to build an SVD from
  A <- rsparsematrix(30, 40, nnz = 50)

  s <- svds(A, 5)

  # vector to multiply by
  x <- rnorm(30)

  # build a residual matrix with matching missingness
  R <- A
  R@x <- rnorm(50)

  expect_false(isTRUE(all.equal(A, R)))
  expect_equal(A != 0, R != 0)

  # expand the SVD only at observed values by hand

  Z <- s$u %*% diag(s$d) %*% t(s$v)  # full SVD
  Y <- A != 0                        # observed indicator

  # add the upper triangle back to Y since we're in the citation
  # setting
  Y <- Y | upper.tri(Y)

  obs_Zt <- t(Z * Y)                 # project SVD onto observed matrix

  A0 <- drop0(A)
  expected <- t(A0) %*% x - obs_Zt %*% x + t(Z) %*% x

  A <- as(A, "TsparseMatrix")
  impl_result <- p_omega_ztx_impl(s$u, s$d, s$v, A@i, A@j, x)

  expect_equal(
    drop(impl_result),
    drop(obs_Zt %*% x)
  )

  args <- list(u = s$u, d = s$d, v = s$v, M = A)
  wrapped_result <- Atx_citation(x, args)

  expect_equal(wrapped_result, expected)
})
