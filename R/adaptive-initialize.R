#' AdaptiveInitialize (sparse)
#'
#' A reference implemention of the `AdaptiveInitialize` using sparse
#' matrix computations.
#'
#' @param M A sparse Matrix created by one of the [Matrix] pkg
#'   constructors.
#' @param r Desired rank to use in the low rank approximation.
#'
#' @return A list with elements:
#'
#'   - `u`: Left-singular-ish vectors
#'   - `d`: Singular-ish values
#'   - `v`: Right-singular-ish vectors
#'
#' @export
#'
#' @examples
#'
#' library(Matrix)
#'
#' set.seed(27)
#'
#' # create a random 8 x 12 sparse matrix with 30 nonzero entries
#' M <- rsparsematrix(8, 12, nnz = 30)
#'
#' s <- sparse_adaptive_initialize(M, 5)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
sparse_adaptive_initialize <- function(M, r,
                                       additional = min(100, ncol(M) / 4)) {

  M <- as(M, "CsparseMatrix")

  n <- nrow(M)
  d <- ncol(M)

  # as.numeric() promotes us from 32 bit integers to 53 bit doubles
  # for some reason this happens by default in the console but not
  # in compiled code
  p_hat <- nnzero(M) / (as.numeric(n) * d)  # line 1

  # now that we have counted explicitly observed zeros
  # we can make the explicit zeros implicit for more
  # efficiency in the initializer
  M <- drop0(M)

  # NOTE: skip explicit computation of line 2
  # NOTE: skip explicit computation of line 3

  # since sigma_p and sigma_t are symmetric the eigenvectors
  # and left and right singular vectors are all the same
  # and it is slightly nicer to user the eigs_sym() interface
  # here than svds(). however, we need to be careful about
  # the singular values. the singular values will be the
  # absolute values of the eigenvalues

  args <- list(M = M, p = p_hat)

  additional <- min(d - r - 1, additional)

  # next we run into the issue that computing the entire
  # SVD of sigma_p will almost always be computational
  # infeasible. instead of computing all the remaining
  # singular values after the first r singular values,
  # we only compute `additional` many more, and assume
  # the uncalculate singular values are zero. this leads
  # to slightly lower value of alpha and less truncation
  # than in the AdapativeInitialize algorithm itself.
  # users can specify the `additional` argument if they
  # better information about the rank of their matrix

  # need the minimum to keep from asking for more singular values
  # of sigma_p than exist. note that the eigenvalues are
  svd_p <- svds(Mx, r + additional, dim = c(d, d), Atrans = Mx,
                args = args)

  svd_t <- svds(Mtx, r, dim = c(n, n), Atrans = Mtx, args = args)

  v_hat <- svd_p$v[, 1:r]  # line 4
  u_hat <- svd_t$u         # line 5

  # approximate calculation of line 6
  alpha <- sum(svd_p$d[r + 1:additional]) / (d - r)

  lambda_hat <- sqrt(svd_p$d[1:r] - alpha) / p_hat  # line 7

  svd_M <- svds(M, r)

  # diag(crossprod()) patterns
  v_sign <- crossprod(rep(1, d), svd_M$v * v_hat)
  u_sign <- crossprod(rep(1, n), svd_M$u * u_hat)
  s_hat <- drop(sign(v_sign * u_sign))  # line 8

  # make the sign adjustment to v_hat so we don't have
  # to carry s_hat around with use. multiplies each
  # *row* of v_hat (i.e. column v_hat^T) by the corresponding
  # element of s_hat
  v_hat <- sweep(v_hat, 2, s_hat, "*")

  list(u = u_hat, d = lambda_hat, v = v_hat)
}

# args are `M` and `p`
# M is n x d, x is a vector in R^d. M^T M is then d x d.
# computes M^T M x - (1 - p) diag(M^T M) x
# using tricks from before to calculate diagonals of cross-products
Mx <- function(x, args) {
  drop(
    crossprod(args$M, args$M %*% x) -
      (1 - args$p) * Diagonal(ncol(args$M), colSums(args$M^2)) %*% x
  )
}

# args are `M` and `p`
# M is n x d, x is a vector in R^d. M M^T is then n x n
# computes M^T M x - (1 - p) diag(M^T M) x
# using tricks from before to calculate diagonals of cross-products
Mtx <- function(x, args) {
  drop(
    args$M %*% crossprod(args$M, x) -
      (1 - args$p) * Diagonal(nrow(args$M), rowSums(args$M^2)) %*% x
  )
}

