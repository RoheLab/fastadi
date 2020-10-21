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
#'   Must be at least `2L` and at most the rank of `X`.
#'
#' @param ... Ignored.
#'
#' @param p_hat The portion of `X` that is observed. Defaults to `NULL`,
#'   in which case `p_hat` is set to the number of observed elements of
#'   `X`. Primarily for internal use in [citation_impute()] or
#'   advanced users.
#'
#' @param alpha_method Either `"exact"` or `"approximate"`, defaulting to
#'   `"exact"`. `"exact"` is computationally expensive and requires taking
#'   a complete SVD of matrix of size `nrow(X)` x `nrow(X)`, and matches
#'   the `AdaptiveInitialize` algorithm exactly. `"approximate"`
#'   departs from the `AdaptiveInitialization` algorithm to compute
#'   a truncated SVD of rank `rank` + `additional` instead of a complete
#'   SVD. This reduces computational burden, but the resulting estimates
#'   of singular-ish values will not be penalized as much as in the
#'   `AdaptiveInitialize` algorithm.
#'
#' @param additional Ignored except when `alpha_method = "approximate"`
#'   in which case it controls the precise of the approximation to `alpha`.
#'   The approximate computation of `alpha` will always understand `alpha`,
#'   but the approximation will be better for larger values of `additional`.
#'   We recommend making `additional` as large as computationally tolerable.
#'
#' @return A low rank matrix factorization represented by an
#'   [adaptive_imputation()] object.
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
  ...,
  p_hat = NULL,
  alpha_method = c("exact", "approximate"),
  additional = NULL
) {

  ellipsis::check_dots_used()

  # avoid issues with svds() not supporting strictly binary Matrix classes
  X <- X * 1

  rank <- as.integer(rank)

  if (rank <= 1)
    stop("`rank` must be an integer >= 2L.", call. = FALSE)

  UseMethod("adaptive_initialize")
}

#' @export
adaptive_initialize.default <- function(
  X,
  rank,
  ...,
  p_hat = NULL,
  alpha_method = c("exact", "approximate"),
  additional = NULL) {

  stop(
    glue("No `adaptive_initialize` method for objects of class {class(X)}."),
    call. = FALSE
  )
}

#' @export
#' @rdname adaptive_initialize
adaptive_initialize.sparseMatrix <- function(
  X,
  rank,
  ...,
  p_hat = NULL,
  alpha_method = c("exact", "approximate"),
  additional = NULL) {

  alpha_method <- match.arg(alpha_method)

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

  # need to divide by p^2 from Cho et al 2016 to get the "right"
  # singular values / singular values on a comparable scale

  n <- nrow(X)
  d <- ncol(X)

  sigma_t <- SigmaT(X, p_hat)

  # (near) dense computation, but truncated
  svd_t <- svds(sigma_t, rank)

  log_info(
    glue(
      "Calculating nuclear norm (slow step). Using {alpha_method} method."
    )
  )

  if (alpha_method == "exact") {

    XtX <- crossprod(X)
    sigma_p <- XtX - (1 - p_hat) * Diagonal(n = d, x = diag(XtX))  # line 2

    # full (near) dense svd -- slow
    svd_p <- svd(sigma_p)
    alpha <- (sum(svd_p$d) - sum(svd_p$d[1:rank])) / (d - rank)  # line 6

  } else if (alpha_method == "approximate") {

    if (is.null(additional))
      stop(
        "Must specify an integer value >= 1 for `additional`.",
        call. = FALSE
      )

    if (additional < 1)
      stop("`additional` must be >= 1L.", call = FALSE)

    if (additional + rank >= d)
      stop("`additional` + `rank` must be less than `ncol(X)`.", call. = FALSE)

    log_info(
      glue(
        "Approximating alpha using rank {rank + additional} approximation."
      )
    )

    sigma_p <- SigmaP(X, p_hat) # line 2, implicitly
    svd_p <- svds(sigma_p, k = rank + additional, nu = rank, nv = rank)

  } else {
    stop("This should not happen.", call. = FALSE)
  }

  v_hat <- svd_p$v[, 1:rank, drop = FALSE]  # line 4
  u_hat <- svd_t$u  # line 5

  alpha <- (sum(svd_p$d) - sum(svd_p$d[1:rank])) / (d - rank)  # line 6

  assert_alpha_positive(alpha)

  lambda_hat <- sqrt(svd_p$d[1:rank] - alpha) / p_hat               # line 7

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
    alpha = alpha,
    ...
  )
}

