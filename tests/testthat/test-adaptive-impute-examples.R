test_that("adaptive impute svd initialization", {

  expect_warning(
    mf <- adaptive_impute(ml100k, rank = 3L, max_iter = 3L)
  )

  preds <- predict(mf, ml100k)

  R <- resid(mf, ml100k)
  norm(R, type = "F") / nnzero(ml100k)
})

test_that("adaptive impute adaptive initialization example", {

  expect_warning(
    mf2 <- adaptive_impute(
      ml100k,
      rank = 3L,
      max_iter = 3L,
      initialization = "adaptive-initialize"
    )
  )


  R2 <- resid(mf2, ml100k)
  norm(R2, type = "F") / nnzero(ml100k)
})