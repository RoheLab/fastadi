library(Matrix)
library(RSpectra)

test_that("p_omega_f_norm_ut", {

  # important: make sure to test when there are elements on the diagonal
  # in practice this shouldn't happen, but we still want correct
  # calculations when this is the case. must be diagonal!
  M <- Matrix(
    rbind(
      c(4, 0, 3, 1, 0),
      c(3, 0, 0, 8, 0),
      c(0, -1, 0, 0, 0),
      c(0, 2, 0, 0, 0),
      c(5, 0, 7, 0, 4)
    )
  )

  # NOTE: not treating nodes as citing themselves!
  Y <- rbind(
    c(1, 1, 1, 1, 1),
    c(1, 0, 1, 1, 1),
    c(0, 1, 0, 1, 1),
    c(0, 1, 0, 0, 1),
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

  impl_result <- p_omega_f_norm_ut_impl(s$u, s$d, s$v, Y@i, Y@j, 1L)
  wrapped_result <- p_omega_f_norm_ut(s, Y)

  expect_equal(impl_result, expected)
  expect_equal(wrapped_result, expected)
})

test_that("citation svd", {

  set.seed(27)

  # matrix to build an SVD from
  M <- rsparsematrix(40, 40, nnz = 50)
  x <- rnorm(40)

  r <- 5
  s <- svds(M, r)

  # expand the SVD only at observed values by hand

  Z <- s$u %*% diag(s$d) %*% t(s$v)  # full SVD
  Y <- M != 0                        # observed indicator

  # add the upper triangle back to Y since we're in the citation
  # setting
  Y <- Y | upper.tri(Y)

  # KEY! M_tilde = P_\Omega(M) + P_\Omega^\perp(Z_t) x
  M_tilde <- M + as(!Y, "CsparseMatrix") * Z

  ### SVD: reference calculation

  expected_svd <- svd(M_tilde, r, r)

  ### SVD: fancy calculation

  args <- list(u = s$u, d = s$d, v = s$v, M = M)

  s_new <- svds(
    Ax_citation,
    k = r,
    Atrans = Atx_citation,
    dim = dim(M),
    args = args
  )

  expect_true(equal_svds(s_new, expected_svd))

  ### alpha: reference calculation

  expected_alpha <- mean(expected_svd$d[(r+1):length(expected_svd$d)]^2)

  ### alpha: fancy calculation

  f_norm_M <- sum(M@x^2)

  # p_omega_f_norm_ut() expects r singular values where r is the rank
  # of the decomposition, but svd() always returns *all* singular values
  # be sure to use s, corresponding to Z_(t-1) here
  clean_svd <- s
  clean_svd$d <- clean_svd$d[1:r]

  # s is the SVD of M
  M_tilde_f_norm <- f_norm_M + sum(s$d^2) - p_omega_f_norm_ut(clean_svd, M)

  # make sure we got the norm of M tilde right
  expect_equal(M_tilde_f_norm, sum(M_tilde^2))

  # as a general rule we want to find the average of the remaining
  # singular values. but in this case, we are working with a wide
  # matrix instead of a tall matrix, so it is impossible to calculate
  # d = ncol(M) singular values. this serves as a brief reminder of
  # that:

  num_s_values <- min(dim(M))

  fancy_alpha <- (M_tilde_f_norm - sum(s_new$d^2)) / (num_s_values - r)

  expect_equal(
    fancy_alpha,
    expected_alpha
  )
})
