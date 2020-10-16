
set.seed(284)
n <- 100
A <- rsparsematrix(n, n, 0.1, rand.x = NULL)

test_that("citation impute check_interval = NULL", {
  expect_warning(
    citation_impute(
      A,
      rank = 2L,
      max_iter = 5L,
      check_interval = NULL
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})

test_that("citation impute check_interval = 1L", {
  expect_warning(
    citation_impute(
      A,
      rank = 2L,
      max_iter = 5L,
      check_interval = 1L
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )

})

test_that("citation impute check_interval > max_iter", {
  expect_warning(
    citation_impute(
      A,
      rank = 2L,
      max_iter = 5L,
      check_interval = 6L
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})

rm(n, A)