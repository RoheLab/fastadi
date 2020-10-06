
test_that("adaptive initialization example", {

  expect_error(
    adaptive_initialize(ml100k, rank = 2L),
    regexp = NA
  )

  expect_error(
    adaptive_initialize(ml100k, rank = 5L),
    regexp = NA
  )
})
