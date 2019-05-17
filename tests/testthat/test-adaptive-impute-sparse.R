# context("test-adaptive-impute-sparse")
#
# set.seed(17)
#
# M <- rsparsematrix(8, 12, nnz = 30)
# s <- svds(M, 5)
#
# y <- as(M, "lgCMatrix")
#
# Z <- s$u %*% diag(s$d) %*% t(s$v)
#
# all.equal(
#   svd_perp(s, M),
#   Z * y
# )
#
#
# set.seed(17)
# r <- 5
#
# M <- rsparsematrix(8, 12, nnz = 30)
# y <- as(M, "lgCMatrix")
#
# s <- svds(M, r)
# Z <- s$u %*% diag(s$d) %*% t(s$v)
#
# M_tilde <- M + Z * (1 - y)  # dense!
#
# Z_perp <- svd_perp(s, M)
# sum_singular_squared <- sum(M@x^2) + sum(s$d^2) - sum(Z_perp@x^2)
#
# all.equal(
#   sum(svd(M_tilde)$d^2),
#   sum_singular_squared
# )
#
# set.seed(17)
# r <- 5
#
# M <- rsparsematrix(8, 12, nnz = 30)
# y <- as(M, "lgCMatrix")
#
# s <- svds(M, r)
# Z <- s$u %*% diag(s$d) %*% t(s$v)
#
# M_tilde <- M + Z * (1 - y)  # dense!
#
# svd_M_tilde <- svds(M_tilde, r)
# svd_M_tilde
#
# Ax <- function(x, args) {
#   drop(M_tilde %*% x)
# }
#
# Atx <- function(x, args) {
#   drop(t(M_tilde) %*% x)
# }
#
# # is eigs_sym() with a two-sided multiply faster?
# args <- list(u = s$u, d = s$d, v = s$v, m = M)
# test1 <- svds(Ax, k = r, Atrans = Atx, dim = dim(M), args = args)
#
# test1
# svd_M_tilde
#
# all.equal(
#   svd_M_tilde,
#   test1
# )
#
#
#
# R <- M - svd_perp(s, M)  # residual matrix
# args <- list(u = s$u, d = s$d, v = s$v, R = R)
#
# # is eigs_sym() with a two-sided multiply faster?
# test2 <- svds(Ax, k = r, Atrans = Atx, dim = dim(M), args = args)
#
# all.equal(
#   svd_M_tilde,
#   test2
# )
#
# # mask as a pair list
# # L and Z / svd are both n x d matrices
# # x is a d x 1 matrix / vector
# masked_svd_times_x <- function(s, mask, x) {
#
#   stopifnot(inherits(mask, "lgTMatrix"))
#
#   u <- s$u
#   d <- s$d
#   v <- s$v
#
#   zx <- numeric(nrow(u))
#
#   # lgTMatrix uses zero based indexing, add one
#   row <- mask@i + 1
#   col <- mask@j + 1
#
#   # need to loop over index of indexes
#   # double looping over i and j here feels intuitive
#   # but is incorrect
#   for (idx in seq_along(row)) {
#     i <- row[idx]
#     j <- col[idx]
#
#     z_ij <- sum(u[i, ] * d * v[j, ])
#     zx[i] <- zx[i] + x[j] * z_ij
#   }
#
#   zx
# }
#
# # how to calculate just one element of the reconstructed
# # data using the SVD
#
# i <- 6
# j <- 4
#
# sum(s$u[i, ] * s$d * s$v[j, ])
# Z[i, j]
#
# # the whole masked matrix multiply
#
# Z <- s$u %*% diag(s$d) %*% t(s$v)
# out <- drop((Z * Y) %*% x)
#
# # check that we did this right
# all.equal(
#   masked_svd_times_x(s, Y, x),
#   out
# )
