library(Matrix)

test_that("sijia's negative alpha example", {

  D <- read.csv("genetreePairwiseDistance_primary.csv")
  C <- exp(-D)
  Csub <- C[-(1:3),-(1:3)]
  Cimcomp <- as.matrix(Csub)
  diag(Cimcomp) <- 0
  matSparse <- as(Cimcomp, "sparseMatrix")

  expect_warning(
    adaptive_impute(matSparse, rank = 3L, max_iter = 5L),
    regexp = "Reached maximum allowed iterations. Returning early."
  )
})
