X <- ml100k
p_hat <- nnzero(ml100k) / prod(dim(ml100k))
k <- 10

test_that("svds.SigmaP", {

  sigma_p <- SigmaP(X, p_hat)
  s <- svds(sigma_p, k = k)

  XtX <- crossprod(X)
  sigma_p2 <- XtX - (1 - p_hat) * Diagonal(n = ncol(X), x = diag(XtX))

  s_expected <- svds(sigma_p2, k = k)

  # checking other attributes is a PITA, looks fine for now
  expect_equal(s$d, s_expected$d)
})


test_that("svds.SigmaT", {

  sigma_p <- SigmaT(X, p_hat)
  s <- svds(sigma_p, k = k)

  XXt <- tcrossprod(X)
  sigma_t2 <- XXt - (1 - p_hat) * Diagonal(n = nrow(X), x = diag(XXt))

  s_expected <- svds(sigma_t2, k = k)

  # checking other attributes is a PITA, looks fine for now
  expect_equal(s$d, s_expected$d)
})

rm(X, p_hat, k)