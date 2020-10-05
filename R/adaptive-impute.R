#' AdaptiveImpute
#'
#' An implementation of the `AdaptiveImpute` algorithm for matrix completion
#' for sparse matrices.
#'
#' @param X A sparse matrix of [Matrix::sparseMatrix()] class.
#'
#' @param rank Desired rank (integer) to use in the low rank approximation.
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
#'   you most of the way to convergence.
#'
#' @param check_interval Integer specifying how often to perform convergence
#'   checks. Defaults to `1L`. In practice, check for convergence requires
#'   a norm calculation that is expensive for large matrices and decreasing
#'   the frequency of convergence checks will reduce computation time.
#'
#' @return A low rank matrix factorization represented by an
#'   `LRMF` object. See [LRMF3::svd_like()] for details.
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

  if (rank <= 2)
    stop("`rank` must be an integer >= 2L.", call. = FALSE)

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
  check_interval = 1L,
  verbose = FALSE,
  p_hat = NULL
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

    MtM <- norm_M + sum(s$d^2) - sum(masked_approximation(s_new, X)^2)
    alpha <- (MtM - sum(s_new$d^2)) / (d - rank)  # line 6

    assert_alpha_positive(alpha)

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    log_debug(glue("lambda_hat = ", paste(s_new$d, collapse = ", ")))

    # NOTE: skip explicit computation of line 8
    delta <- relative_f_norm_change(s_new, s)

    # save a little bit on computation
    if (iter %% check_interval == 0) {
      # log_info("Finding relative change in Frobenius norm.")
      delta <- relative_f_norm_change(s_new, s)
    }

    s <- s_new

    if (iter %% check_interval == 0)
      log_info(
        glue(
          "Iter {iter} complete. ",
          "delta = {round(delta, 8)}, ",
          "alpha = {round(alpha, 3)}"
        )
      )

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

