low_rank_adaptive_impute <- function(M, r, epsilon = 1e-03) {

  # coerce M to sparse matrix such that we use sparse operations
  M <- as(M, "dgCMatrix")

  # low rank svd-like object, s ~ Z_1
  s <- low_rank_adaptive_initialization(M, r)  # line 1
  delta <- Inf

  norm_M <- sum(M^2)

  while (delta > epsilon) {
    s_new <- low_rank_threshold(s, M, norm_M, r)
    delta <- low_rank_relative_change(s, s_new, M)
    s <- s_new
  }

  s
}

low_rank_threshold <- function(s, M, norm_M, r) {

  d <- ncol(M)

  # update s: lines 4 and 5
  # take the SVD of M-tilde
  s <- update_svd(s, M)
  MtM <- norm_M + sum(s$d^2) - sum(expand_nnzero_only(s, M)^2)
  alpha <- (sum(MtM) - sum(s$d^2)) / (d - r)  # line 6
  s$d <- sqrt(s$d^2 - alpha)  # line 7. will explode if negative

  # NOTE: skip explicit computation of line 8
  s
}

low_rank_relative_change <- function(s, s_new, M) {
  # pass
}

# s is a matrix defined in terms of it's svd
# G is a sparse matrix
# compute only elements of U %*% diag(d) %*% t(V) only on non-zero elements of G
# G and U %*% t(V) must have same dimensions
expand_nnzero_only <- function(s, G) {

  # note: must be dgTMatrix to get column indexes j larger
  mt <- as(G, "dgTMatrix")

  # the indices for which we want to compute the matrix multiplication
  # turn zero based indices into one based indices
  i <- mt@i + 1
  j <- mt@j + 1

  # gets rows and columns of U and V to multiply, then multiply
  ud <- s$u %*% diag(s$d)
  left <- ud[i, ]
  right <- s$v[j, ]

  # compute inner products to get elements of U %*% t(V)
  uv <- rowSums(left * right)

  # NOTE: specify dimensions just in case
  sparseMatrix(i = i, j = j, x = uv, dims = dim(G))
}
