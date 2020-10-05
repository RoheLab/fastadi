

test_that("citation impute svd initialization example", {

  set.seed(284)
  n <- 100
  A <- rsparsematrix(n, n, 0.1, rand.x = NULL)

  expect_warning(
    mf <- citation_impute(A, rank = 3L, max_iter = 3L)
  )

})

test_that("citation impute adaptive initialization example", {
  set.seed(284)
  n <- 100
  A <- rsparsematrix(n, n, 0.1, rand.x = NULL)

  expect_warning(
    mf2 <- citation_impute(
      A,
      rank = 3L,
      max_iter = 3L,
      initialization = "adaptive-initialize"
    )
  )
})


