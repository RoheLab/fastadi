#' AdaptiveImpute
#'
#' An implementation of the `AdaptiveImpute` algorithm for matrix completion
#' for sparse matrices.
#'
#' @param X A sparse matrix of [Matrix::sparseMatrix()] class.
#'
#' @param rank Desired rank (integer) to use in the low rank approximation.
#'   Must be at least `2L` and at most the rank of `X`. Note that the rank
#'   of `X` is typically unobserved and computations may be unstable or
#'   even fail when `rank` is near or exceeds this threshold.
#'
#' @param ... Unused additional arguments.
#'
#' @param initialization How to initialize the low rank approximation.
#'   Options are:
#'
#'   - `"svd"` (default). In the initialization step, this treats
#'     unobserved values as zeroes.
#'
#'   - `"adaptive-initialize"`. In the initialization step, this treats
#'     unobserved values as actually unobserved. However, the current
#'     `AdaptiveInitialize` implementation relies on dense matrix
#'     computations that are only suitable for relatively small matrices.
#'
#'   Note that initialization matters as `AdaptiveImpute` optimizes
#'   a non-convex objective. The current theory shows that initializing
#'   with `AdaptiveInitialize` leads to a consistent estimator, but it
#'   isn't know if this is the case for SVD initialization. Empirically
#'   we have found that SVD initialization works well nonetheless.
#'
#' @param epsilon Convergence criteria, measured in terms of relative change
#'   in Frobenius norm of the full imputed matrix. Defaults to `1e-7`.
#'
#' @param max_iter Maximum number of iterations to perform (integer). Defaults
#'   to `200L`. In practice 10 or so iterations will get you a decent
#'   approximation to use in exploratory analysis, and and 50-100 will get
#'   you most of the way to convergence. Must be at least `1L`.
#'
#' @param check_interval Integer specifying how often to perform convergence
#'   checks. Defaults to `1L`. In practice, check for convergence requires
#'   a norm calculation that is expensive for large matrices and decreasing
#'   the frequency of convergence checks will reduce computation time. Can
#'   also be set to `NULL`, which case `max_iter` iterations of the algorithm
#'   will occur with no possibility of stopping due to small relative change
#'   in the imputed matrix. In this case `delta` will be reported as `Inf`.
#'
#' @return A low rank matrix factorization represented by an
#'   [adaptive_imputation()] object.
#'
#' @export
#'
#' @examples
#'
#' mf <- adaptive_impute(ml100k, rank = 3L, max_iter = 20L)
#' # mf
#'
#' # build a rank-5 approximation only for
#' # observed elements of ml100k
#'
#' preds <- predict(mf, ml100k)
#'
#' # estimate the in-sample reconstruction mse
#'
#' R <- resid(mf, ml100k)
#' norm(R, type = "F") / nnzero(ml100k)
#'
#'
#' mf2 <- adaptive_impute(
#'   ml100k,
#'   rank = 3L,
#'   max_iter = 20L,
#'   initialization = "adaptive-initialize"
#' )
#'
#' # mf2
#'
#' R2 <- resid(mf2, ml100k)
#' norm(R2, type = "F") / nnzero(ml100k)
#'
adaptive_impute <- function(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize"),
  max_iter = 200L,
  check_interval = 1L,
  epsilon = 1e-7
) {

  ellipsis::check_dots_used()

  rank <- as.integer(rank)

  if (length(rank) > 1)
    stop(
      "`rank` must be an integer vector with a single element.",
      call. = FALSE
    )

  if (length(max_iter) > 1)
    stop(
      "`max_iter` must be an integer vector with a single element.",
      call. = FALSE
    )

  if (!is.null(check_interval) && length(check_interval) > 1)
    stop(
      "`check_interval` must be a single integer, or NULL.",
      call. = FALSE
    )

  if (rank <= 1 || rank >= min(nrow(X), ncol(X)))
    stop(
      "rank must satisfy 1 < rank < min(nrow(X), ncol(X)).",
      call. = FALSE
    )

  if (max_iter < 1)
    stop("`max_iter` must be an integer >= 1L.", call. = FALSE)

  if (!is.null(check_interval) && check_interval < 1)
    stop("`check_interval` must be an integer >= 1L, or NULL.", call. = FALSE)

  UseMethod("adaptive_impute")
}

#' @export
adaptive_impute.default <- function(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize"),
  max_iter = 200L,
  check_interval = 1L,
  epsilon = 1e-7) {

  stop(
    glue("No `adaptive_impute` method for objects of class {class(X)}."),
    call. = FALSE
  )
}

#' @export
#' @rdname adaptive_impute
adaptive_impute.sparseMatrix <- function(
  X,
  rank,
  initialization = c("svd", "adaptive-initialize"),
  ...
) {

  initialization <- match.arg(initialization)

  log_info(glue("Use {initialization} initialization."))

  if (initialization == "svd") {
    s <- svds(X, rank)
    mf <- as_svd_like(s)
  } else if (initialization == "adaptive-initialize") {
    mf <- adaptive_initialize(X, rank)
  } else {
    stop("This should not happen.", call. = FALSE)
  }

  log_info("Done initializing.")

  adaptive_impute.LRMF(mf, X, ...)
}

#' @export
#' @rdname adaptive_impute
adaptive_impute.LRMF <- function(
  X,
  rank,
  ...,
  epsilon = 1e-7,
  max_iter = 200L,
  check_interval = 1L
) {

  log_info(glue("Beginning AdaptiveImpute (max {max_iter} iterations)."))
  log_info(glue("Checking convergence every {check_interval} iteration(s)."))

  # first argument is the svd_like object, second is the data
  # do some renaming here

  s <- X
  X <- rank
  rank <- s$rank

  ### ITERATION STAGE

  delta <- Inf
  d <- ncol(X)
  norm_M <- norm(X, type = "F")^2
  iter <- 1L

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    R <- X - masked_approximation(s, X)  # residual matrix
    args <- list(u = s$u, d = s$d, v = s$v, R = R)

    s_new <- svds(Ax, k = rank, Atrans = Atx, dim = dim(X), args = args)

    # if there's an alpha < 0 bug it almost certainly comes from these
    # next two lines of code
    MtM <- norm_M + sum(s$d^2) - sum(masked_approximation(s, X)^2)
    alpha <- (MtM - sum(s_new$d^2)) / (d - rank)  # line 6

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    log_debug(glue("lambda_hat = ", paste(s_new$d, collapse = ", ")))

    # save a little bit on computation and only check for
    # # convergence intermittently
    if (!is.null(check_interval) && iter %% check_interval == 0) {
      log_debug("Computing relative change in Frobenius norm.")
      delta <- relative_f_norm_change(s_new, s)
    }

    s <- s_new

    log_info(
      glue(
        "Iter {iter} complete. ",
        "delta = {if (!is.null(check_interval)) round(delta, 8) else Inf}, ",
        "alpha = {round(alpha, 3)}"
      )
    )

    assert_alpha_positive(alpha)

    iter <- iter + 1

    if (iter > max_iter) {
      warning(
        "\nReached maximum allowed iterations. Returning early.\n",
        call. = FALSE
      )
      break
    }

  }

  adaptive_imputation(
    u = s$u,
    d = s$d,
    v = s$v,
    alpha = alpha
  )
}

