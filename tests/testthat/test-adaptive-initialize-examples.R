
test_that("adaptive initialization example", {

  expect_error(
    adaptive_initialize(ml100k, rank = 2L),
    regexp = NA
  )

  expect_error(
    adaptive_initialize(ml100k, rank = 2L, alpha_method = "approximate"),
    regexp = "Must specify an integer value >= 1 for `additional`."
  )

  expect_error(
    adaptive_initialize(
      ml100k,
      rank = 2L,
      alpha_method = "approximate",
      additional = 25
    ),
    regexp = NA
  )
})
