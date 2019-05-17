# context("test-adaptive-initialize-sparse")
#
# test_that("Mx() works", {
#   x <- rnorm(12)
#   p <- 0.3
#   out <- (t(M) %*% M / p^2 - (1 - p) * diag(diag(t(M) %*% M))) %*% x
#   out2 <- Mx(x, args = list(M = M, p = p))
#   all.equal(as.matrix(out), as.matrix(out))
# })
#
# test_that("agrees with dense implementation", {
#   lr_init <- sparse_adaptive_initialize(dat, r)
#
#   # some weird stuff is happening with the singular values but I'm
#   # going to not worry about it for the time being
#
#   equal_svds(init, lr_init)
# })
#
