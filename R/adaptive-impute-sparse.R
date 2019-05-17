#' AdaptiveImpute (sparse)
#'
#' A reference implemention of the `AdaptiveImpute` algorithm using sparse
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
#' s <- sparse_adaptive_impute(M, 5)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
sparse_adaptive_impute <- function(M, r, epsilon = 1e-7) {

  # coerce M to sparse matrix such that we use sparse operations
  M <- as(M, "dgCMatrix")

  # low rank svd-like object, s ~ Z_1
  s <- low_rank_adaptive_initialize(M, r)  # line 1
  delta <- Inf
  d <- ncol(M)
  norm_M <- sum(M@x^2)

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    R <- M - masked_approximation(s, M)  # residual matrix
    args <- list(u = s$u, d = s$d, v = s$v, R = R)

    s_new <- svds(Ax, k = r, Atrans = Atx, dim = dim(M), args = args)

    MtM <- norm_M + sum(s_new$d^2) - sum(masked_approximation(s_new, M)^2)

    alpha <- (sum(MtM) - sum(s_new$d^2)) / (d - r)  # line 6

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    # NOTE: skip explicit computation of line 8
    delta <- relative_f_norm_change(s_new, s)

    s <- s_new

    print(glue::glue("delta: {round(delta, 8)}, alpha: {round(alpha, 3)}"))
  }

  s
}

relative_f_norm_change <- function(s_new, s) {
  # TODO: don't do the dense calculation here

  Z_new <- s_new$u %*% diag(s_new$d) %*% t(s_new$v)
  Z <- s$u %*% diag(s$d) %*% t(s$v)
  sum((Z_new - Z)^2) / sum(Z^2)
}

Ax <- function(x, args) {
  drop(args$R %*% x + args$u %*% diag(args$d) %*% crossprod(args$v, x))
}

Atx <- function(x, args) {
  drop(t(args$R) %*% x + args$v %*% diag(args$d) %*% crossprod(args$u, x))
}

