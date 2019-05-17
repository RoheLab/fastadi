#' AdaptiveImpute (dense)
#'
#' A reference implemention of the `AdaptiveImpute` algorithm using dense
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
#' s <- dense_adaptive_impute(M, 5)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
dense_adaptive_impute <- function(M, r, epsilon = 1e-7) {

  s <- dense_adaptive_initialize(M, r)
  Z <- s$u %*% diag(s$d) %*% t(s$v)  # line 1
  delta <- Inf

  while (delta > epsilon) {

    y <- as(M, "lgCMatrix")  # indicator if entry of M observed
    M_tilde <- M + Z * (1 - y)  # line 3

    svd_M <- svds(M_tilde, r)

    u_hat <- svd_M$u  # line 4
    v_hat <- svd_M$v  # line 5

    d <- ncol(M)

    alpha <- (sum(M_tilde^2) - sum(svd_M$d^2)) / (d - r)  # line 6

    lambda_hat <- sqrt(svd_M$d^2 - alpha)  # line 7

    Z_new <- u_hat %*% diag(lambda_hat) %*% t(v_hat)

    delta <- sum((Z_new - Z)^2) / sum(Z^2)
    Z <- Z_new

    print(glue::glue("delta: {round(delta, 8)}, alpha: {round(alpha, 3)}"))
  }

  list(u = u_hat, d = lambda_hat, v = v_hat)
}
