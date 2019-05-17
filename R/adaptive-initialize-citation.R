#' AdaptiveInitialize (citations)
#'
#' A reference implemention of the `AdaptiveInitialize` using sparse
#' matrix computations. For citation matrices where missing values in
#' the upper triangle are taken to be *explicitly observed* observed
#' zeros, as opposed to missing values.
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
#' s <- citation_adaptive_initialize(M, 5)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
citation_adaptive_initialize <- function(M, r) {

  M <- as(M, "dgCMatrix")

  p_hat <- nnzero(M) / prod(dim(M))  # line 1

  # NOTE: skip explicit computation of line 2
  # NOTE: skip explicit computation of line 3

  eig_p <- eigen_helper(M, r)
  eig_t <- eigen_helper(t(M), r)

  v_hat <- eig_p$vectors  # line 4
  u_hat <- eig_t$vectors  # line 5

  d <- ncol(M)
  n <- nrow(M)

  # NOTE: alpha is again incorrect since we work with eigenvalues
  # rather than singular values here
  sum_eigen_values <- sum(M@x^2) / p^2 - (1 - p) * sum(colSums(M^2))
  alpha <- (sum_eigen_values - sum(eig_p$values)) / (d - r)  # line 6

  lambda_hat <- sqrt(eig_p$values - alpha) / p_hat  # line 7

  # TODO: Karl had another sign computation here that he said was faster
  # but it wasn't documented anywhere, so I'm going with what was in the
  # paper

  svd_M <- svds(M, r)

  # v_hat is d by r
  v_sign <- crossprod(rep(1, d), svd_M$v * v_hat)
  u_sign <- crossprod(rep(1, n), svd_M$u * u_hat)
  s_hat <- c(sign(v_sign * u_sign))  # line 8

  lambda_hat <- lambda_hat * s_hat

  list(u = u_hat, d = lambda_hat, v = v_hat)
}

# Take the eigendecomposition of t(M) %*% M - (1 - p) * diag(t(M) %*% M)
# using sparse computations only
eigen_helper <- function(M, r) {
  eigs_sym(
    Mx, r,
    n = ncol(M),
    args = list(
      M = M,
      p = nnzero(M) / prod(dim(M))
    )
  )
}

# compute (t(M) %*% M / p^2 - (1 - p) * diag(diag(t(M) %*% M))) %*% x
# using sparse operations
Mx <- function(x, args) {
  drop(
    crossprod(args$M, args$M %*% x) / args$p^2 -
      (1 - args$p) * Diagonal(ncol(args$M), colSums(args$M^2)) %*% x
  )
}
