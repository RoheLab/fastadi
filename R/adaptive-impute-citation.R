#' AdaptiveImpute (sparse)
#'
#' A reference implemention of the `AdaptiveImpute` algorithm using sparse
#' matrix computations. For citation matrices where missing values in
#' the upper triangle are taken to be *explicitly observed* observed
#' zeros, as opposed to missing values.
#'
#' @param M A sparse Matrix created by one of the [Matrix] pkg
#'   constructors.
#' @param r Desired rank to use in the low rank approximation.
#' @param epsilon Tolerance, measured in terms of relative change
#'   in Frobenius norm of the full imputed matrix.
#'
#' @return A list with elements:
#'
#'   - `u`: Left singular-ish vectors
#'   - `d`: Singular-ish values
#'   - `v`: Right singular-ish vectors
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
#' s <- citation_adaptive_impute(M, 5)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
citation_adaptive_impute <- function(M, r, epsilon = 1e-7) {

  # NOTE: no differences are necessary from the sparse
  # adaptive_initialize since M just has more zeros

  s <- sparse_adaptive_initialize(M, r)  # line 1
  delta <- Inf
  d <- ncol(M)
  f_norm_M <- sum(M@x^2)

  # coerce M to sparse matrix such that we use sparse operations
  M <- as(M, "CsparseMatrix")

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    args <- list(u = s$u, d = s$d, v = s$v, M = M)

    s_new <- svds(
      Ax_citation,
      k = r,
      Atrans = Atx_citation,
      dim = dim(M),
      args = args
    )

    M_tilde_f_norm <- f_norm_M + sum(s_new$d^2) -
      p_omega_f_norm_ut(s_new, M)

    alpha <- (M_tilde_f_norm - sum(s_new$d^2)) / (d - r)  # line 6

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    # NOTE: skip explicit computation of line 8

    delta <- relative_f_norm_change(s_new, s)

    s <- s_new

    print(glue::glue("delta: {round(delta, 8)}, alpha: {round(alpha, 3)}"))
  }

  s
}

# mask needs to be a sparse matrix stored as triplets
p_omega_f_norm_ut <- function(s, mask) {
  mask <- as(mask, "TsparseMatrix")
  p_omega_f_norm_ut_impl(s$u, s$d, s$v, mask@i, mask@j)
}

# KEY: must return a vector! test this otherwise you will be sad!
Ax_citation <- function(x, args) {
  mask <- as(args$M, "TsparseMatrix")

  out <- drop0(args$M) %*% x -
    p_omega_zx_impl(args$u, args$d, args$v, mask@i, mask@j, x) +
    args$u %*% diag(args$d) %*% crossprod(args$v, x)

  drop(out)
}

Atx_citation <- function(x, args) {

  mask <- as(args$M, "TsparseMatrix")

  out <- t(drop0(args$M)) %*% x -
    p_omega_ztx_impl(args$u, args$d, args$v, mask@i, mask@j, x) +
    args$v %*% diag(args$d) %*% crossprod(args$u, x)

  drop(out)
}

