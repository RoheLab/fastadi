library(LRMF3)

skip_if_not_installed("vsp")

library(vsp)

test_that("fastadi, no transformations", {

  expect_warning(
    mf <- adaptive_impute(ml100k, rank = 3, max_iter = 5),
    regexp = "Reached maximum allowed iterations. Returning early."
  )

  skip_on_cran()

  expect_silent(
    fa <- vsp(mf)
  )
})

test_that("fastadi, with scaling", {

  scaler <- RegularizedLaplacian(ml100k)
  M <- invertiforms::transform(scaler, ml100k)

  expect_warning(
    mf <- adaptive_impute(M, rank = 3, max_iter = 5),
    regexp = "Reached maximum allowed iterations. Returning early."
  )

  skip_on_cran()

  expect_silent(
    vsp(mf, scaler = scaler, rescale = TRUE)
  )
})
