#' AdaptiveInitialize (dense)
#'
#' A reference implemention of the `AdaptiveInitialize` using dense
#' matrix computations. Designed to be readable rather than performant.
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
#' s <- dense_adaptive_initialize(M, 5)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
dense_adaptive_initialize <- function(M, r) {

  # TODO: ignores observed zeros!
  p_hat <- nnzero(M) / prod(dim(M))  # line 1

  MtM <- crossprod(M)
  MMt <- tcrossprod(M)

  # need to divide by p^2 from Cho et al 2016 to get the "right"
  # singular values / singular values on a comparable scale

  # both of these matrices are symmetric, but not necessarily positive
  # this has important implications for the SVD / eigendecomp relationship

  sigma_p <- MtM / p_hat^2 - (1 - p_hat) * diag(diag(MtM))  # line 2
  sigma_t <- MMt / p_hat^2 - (1 - p_hat) * diag(diag(MMt))  # line 3

  # crossprod() and tcrossprod() return dsCMatrix objects,
  # sparse matrix objects that know they are symmetric

  # unfortunately, RSpectra doesn't support dsCMatrix objects,
  # but does support dgCMatrix objects, a class representing sparse
  # but not symmetric matrices

  # support for dsCMatrix objects in RSpectra is on the way,
  # which will eliminate the need for the following coercions.
  # see: https://github.com/yixuan/RSpectra/issues/15

  sigma_p <- as(sigma_p, "dgCMatrix")
  sigma_t <- as(sigma_t, "dgCMatrix")

  # TODO: is eigs_sym() faster?
  svd_p <- svds(sigma_p, r)
  svd_t <- svds(sigma_t, r)

  v_hat <- svd_p$v  # line 4
  u_hat <- svd_t$u  # line 5

  n <- nrow(M)
  d <- ncol(M)

  # NOTE: alpha is incorrect due to singular values and eigenvalues
  # being different when sigma_p is not positive

  alpha <- (sum(diag(sigma_p)) - sum(svd_p$d)) / (d - r)  # line 6
  lambda_hat <- sqrt(svd_p$d - alpha) / p_hat             # line 7

  svd_M <- svds(M, r)

  v_sign <- crossprod(rep(1, d), svd_M$v * v_hat)
  u_sign <- crossprod(rep(1, n), svd_M$u * u_hat)
  s_hat <- drop(sign(v_sign * u_sign))

  lambda_hat <- lambda_hat * s_hat  # line 8

  list(u = u_hat, d = lambda_hat, v = v_hat)
}
