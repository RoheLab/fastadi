#' Upper triangular SVD (sparse)
#'
#' A reference implemention of the `AdaptiveImpute` algorithm using sparse
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
#' M <- rsparsematrix(12, 12, nnz = 30)
#'
#' s <- citation_svds(M, 5, max_iter = 20)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
citation_svds <- function(M, r) {

  # only works for square matrices
  stopifnot(ncol(M) == nrow(M))

  svds(
    Ax_citation_svds,
    k = r,
    Atrans = Atx_citation_svds,
    dim = dim(M),
    args = args
  )
}

# mask needs to be a sparse matrix stored as triplets
p_omega_f_norm_ut <- function(s, mask) {
  mask <- as(mask, "TsparseMatrix")
  p_omega_f_norm_ut_impl(s$u, s$d, s$v, mask@i, mask@j)
}

# KEY: must return a vector! test this otherwise you will be sad!
Ax_citation_svds <- function(x, args) {
  mask <- as(args$M, "TsparseMatrix")

  # out <- drop0(args$M) %*% x
  # out <- out - p_omega_zx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  # out <- out + args$u %*% diag(args$d) %*% crossprod(args$v, x)

  out <- drop0(args$M) %*% x
  out <- out - p_u_zx_impl(args$u, args$d, args$v, x)
  out <- out - p_u_tilde_zx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  out <- out + args$u %*% diag(args$d) %*% crossprod(args$v, x)

  drop(out)
}

Atx_citation_svds <- function(x, args) {

  mask <- as(args$M, "TsparseMatrix")

  # out <- t(drop0(args$M)) %*% x
  # out <- out - p_omega_ztx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  # out <- out + args$v %*% diag(args$d) %*% crossprod(args$u, x)

  out <- t(drop0(args$M)) %*% x
  out <- out - p_u_ztx_impl(args$u, args$d, args$v, x)
  out <- out - p_u_tilde_ztx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  out <- out + args$v %*% diag(args$d) %*% crossprod(args$u, x)

  drop(out)
}

