test_that("adaptive impute check_interval = NULL", {

  expect_warning(
    mf <- adaptive_impute(
      ml100k,
      rank = 2L,
      max_iter = 5L,
      check_interval = NULL
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})

test_that("adaptive impute check_interval = 1L", {

  expect_warning(
    mf <- adaptive_impute(
      ml100k,
      rank = 2L,
      max_iter = 5L,
      check_interval = 1L
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )

})

test_that("adaptive impute check_interval > max_iter", {

  expect_warning(
    mf <- adaptive_impute(
      ml100k,
      rank = 2L,
      max_iter = 5L,
      check_interval = 6L
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )

})