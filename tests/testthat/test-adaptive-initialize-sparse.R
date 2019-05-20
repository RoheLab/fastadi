context("sparse adaptive initialize")

test_that("Mx()", {
  M <- rsparsematrix(8, 12, nnz = 30)
  x <- rnorm(12)
  p <- 0.3

  expected <- (t(M) %*% M / p^2 - (1 - p) * diag(diag(t(M) %*% M))) %*% x
  Mx_result <- Mx(x, args = list(M = M, p = p))

  expect_equal(
    drop(expected),
    drop(Mx_result)
  )
})

test_that("agrees with dense implementation", {

  set.seed(27)

  M <- rsparsematrix(8, 12, nnz = 30)

  sparse <- sparse_adaptive_initialize(M, 5)
  dense <- dense_adaptive_initialize(M, 5)

  expect_true(equal_svds(sparse, dense))
})

