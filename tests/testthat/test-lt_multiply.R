# context("test-lt_multiply")
#
# n <- 7
# d <- 5
# k <- 3
#
# M <- matrix(rep(0, n * d), n, d)
# L <- M
# L[lower.tri(L)] <- 1
# L
#
# U <- matrix(1:(n * k), n, k)
# V <- matrix(rev(1:(d * k)), d, k)
#
# x <- 1:d
#
# # R solution
# (L * U %*% t(V)) %*% x
# lt_multiply(U, V, x)
#
#
