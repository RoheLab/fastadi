library(RSpectra)

equal_svds <- function(s, s2) {

  # svd() always gives you all the singular values, but we only
  # want to compare the first k
  k <- ncol(s$u)

  # the term sign(s$u) * sign(s2$u) performs a sign correction
  # see https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
  # for details

  # isTRUE because output of all.equal is not a boolean, it's something
  # weird when the inputs aren't equal

  # check.attributes = FALSE because different matrix classes doesn't matter
  # and vector versus single column matrix doesn't matter

  u_ok <- isTRUE(
    all.equal(s$u, s2$u * sign(s$u) * sign(s2$u),
              check.attributes = FALSE,
              tolerance = 1e-6)
  )

  v_ok <- isTRUE(
    all.equal(s$v, s2$v * sign(s$v) * sign(s2$v),
              check.attributes = FALSE,
              tolerance = 1e-6)
  )

  d_ok <- isTRUE(
    all.equal(s$d[1:k], s2$d[1:k],
              check.attributes = FALSE,
              tolerance = 1e-6)
  )

  u_ok && v_ok && d_ok
}

check_svd_like <- function(s, check_orthonormal = FALSE) {

  # must have u, v, d
  expect_true(all(c("u", "d", "v") %in% names(s)))

  # original data is n x d

  expect_equal(ncol(s$u), length(s$d))
  expect_equal(length(s$d), ncol(s$v))

  # U is n x r
  # d is length(r)
  # V is r by d

  # if check orthonormal, make sure that inverse and transpose
  # are the same
}

check_svd_like_matches_matrix <- function(s, A) {
  expect_equal(nrow(s$u), nrow(A))
  expect_equal(nrow(s$v), ncol(A))
}

test_that("helpers work", {
  set.seed(27)

  M <- rsparsematrix(8, 12, nnz = 30)
  k <- 4

  s <- svd(M, k, k)
  s2 <- svds(M, k, k)

  expect_silent(equal_svds(s, s2))

  # errors since svd() returns all the singular values no matter what
  # expect_silent(check_svd_like(s))
  expect_silent(check_svd_like(s2))

  expect_silent(check_svd_like_matches_matrix(s, M))
  expect_silent(check_svd_like_matches_matrix(s2, M))
})

