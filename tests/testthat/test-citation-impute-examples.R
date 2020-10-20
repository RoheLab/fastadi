

test_that("citation impute svd initialization example", {

  set.seed(284)
  n <- 100
  A <- rsparsematrix(n, n, 0.1, rand.x = NULL)

  expect_error(
    citation_impute(A, rank = 2L, max_iter = 1L),
    regexp = NA
  )

  expect_error(
    citation_impute(A, rank = 3L, max_iter = 3L),
    regexp = NA
  )

})

test_that("citation impute adaptive initialization example", {
  set.seed(284)
  n <- 100
  A <- rsparsematrix(n, n, 0.1, rand.x = NULL)

  expect_warning(
    citation_impute(
      A,
      rank = 3L,
      max_iter = 3L,
      initialization = "adaptive-initialize"
    )
  )
})


test_that("citation impute approximate initialization example", {

  set.seed(284)
  n <- 100
  A <- rsparsematrix(n, n, 0.1, rand.x = NULL)

  expect_error(
    citation_impute(
      A,
      rank = 3L,
      max_iter = 3L,
      initialization = "approximate"
    ),
    regexp = "Must specify `additional` when using approximate initialization."
  )

  expect_error(
    citation_impute(
      A,
      rank = 3L,
      max_iter = 3L,
      initialization = "approximate",
      additional = n + 10
    )
  )

  expect_warning(
    citation_impute(
      A,
      rank = 3L,
      max_iter = 3L,
      initialization = "approximate",
      additional = 5
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})

