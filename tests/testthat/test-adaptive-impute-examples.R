test_that("adaptive impute svd initialization", {

  expect_warning(
    mf <- adaptive_impute(ml100k, rank = 2L, max_iter = 1L),
    regexp = "Reached maximum allowed iterations. Returning early."
  )

  expect_warning(
    mf <- adaptive_impute(ml100k, rank = 3L, max_iter = 3L),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})

test_that("adaptive impute adaptive initialization example", {

  expect_warning(
    mf2 <- adaptive_impute(
      ml100k,
      rank = 3L,
      max_iter = 3L,
      initialization = "adaptive-initialize"
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})

test_that("adaptive impute approximate initialization example", {

  expect_error(
    adaptive_impute(
      ml100k,
      rank = 3L,
      max_iter = 3L,
      initialization = "approximate"
    ),
    regexp = "Must specify `additional` when using approximate initialization."
  )

  expect_warning(
    adaptive_impute(
      ml100k,
      rank = 3L,
      max_iter = 3L,
      initialization = "approximate",
      additional = 10
    ),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})
