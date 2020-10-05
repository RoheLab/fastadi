#' AdaptiveInitialize
#'
#' An implementation of the `AdaptiveInitialize` algorithm for
#' matrix imputation for sparse matrices. At the moment the implementation
#' is only suitable for small matrices with on the order of thousands
#' of rows and columns at most.
#'
#' @param X A sparse matrix of `sparseMatrix` class. Explicit (observed)
#'   zeroes in `X` can be dropped for
#'
#' @param rank Desired rank (integer) to use in the low rank approximation.
#'
#' @param p_hat The portion of `X` that is observed. Defaults to `NULL`,
#'   in which case `p_hat` is set to the number of observed elements of
#'   `X`. Useful when you have many explicit zeros that w
#'
#' @param ... Ignored.
#'
#' @return A low rank matrix factorization represented by an
#'   `LRMF` object. See [LRMF3::svd_like()] for details.
#'
#' @export
#'
#' @examples
#'
#' mf <- adaptive_initialize(ml100k, rank = 3L)
#' mf
#'
#' # build a rank-5 approximation only for
#' # observed elements of ml100k
#' preds <- predict(mf, ml100k)
#'
adaptive_initialize <- function(
  X,
  rank,
  p_hat = NULL,
  ...
) {

  ellipsis::check_dots_used()

  # avoid issues with svds() not supporting strictly binary Matrix classes
  X <- X * 1

  rank <- as.integer(rank)

  if (rank <= 2)
    stop("`rank` must be an integer >= 2L.", call. = FALSE)

  UseMethod("adaptive_initialize")
}

#' @export
adaptive_initialize.default <- function(
  X,
  rank,
  p_hat = NULL,
  ...) {

  stop(
    glue("No `adaptive_initialize` method for objects of class {class(X)}."),
    call. = FALSE
  )
}

#' @export
#' @rdname adaptive_impute
adaptive_initialize.sparseMatrix <- function(
  X,
  rank,
  p_hat = NULL,
  ...
) {

  log_info("Beginning AdaptiveInitialize.")

  if (is.null(p_hat)) {
    p_hat <- nnzero(X) / prod(dim(X))  # line 1
  }

  log_info(
    glue(
      "p_hat = {p_hat}, non-zero entries = {nnzero(X)}, ",
      "total entries = {prod(as.numeric(dim(X)))}"
    )
  )

  XtX <- crossprod(X)
  XXt <- tcrossprod(X)

  # need to divide by p^2 from Cho et al 2016 to get the "right"
  # singular values / singular values on a comparable scale

  # both of these matrices are symmetric, but not necessarily positive
  # this has important implications for the SVD / eigendecomp relationship

  sigma_p <- XtX - (1 - p_hat) * diag(diag(XtX))  # line 2
  sigma_t <- XXt - (1 - p_hat) * diag(diag(XXt))  # line 3

  # these computations are dense
  svd_p <- svds(sigma_p, rank)
  svd_t <- svds(sigma_t, rank)

  v_hat <- svd_p$v  # line 4
  u_hat <- svd_t$u  # line 5

  n <- nrow(X)
  d <- ncol(X)

  log_info("Starting nuclear norm computation (slow step).")

  # this next step is very slow and should be sped up when we have
  svd_p_full <- svd(sigma_p)

  alpha <- (sum(svd_p_full$d) - sum(svd_p$d)) / (d - rank)  # line 6
  lambda_hat <- sqrt(svd_p$d - alpha) / p_hat               # line 7

  assert_alpha_positive(alpha)

  log_info(glue("Computation complete, alpha = {alpha}"))
  log_info(glue("lambda_hat = ", paste(lambda_hat, collapse = ", ")))

  svd_X <- svds(X, rank)

  v_sign <- crossprod(rep(1, d), svd_X$v * v_hat)
  u_sign <- crossprod(rep(1, n), svd_X$u * u_hat)
  s_hat <- drop(sign(v_sign * u_sign))  # line 8

  # make the sign adjustment to v_hat so we don't have
  # to carry s_hat around with use. multiplies each
  # *row* of v_hat (i.e. column v_hat^T) by the corresponding
  # element of s_hat
  v_hat <- sweep(v_hat, 2, s_hat, "*")

  adaptive_imputation(
    u = u_hat,
    d = lambda_hat,
    v = v_hat,
    rank = rank,
    alpha = alpha,
    ...
  )
}

