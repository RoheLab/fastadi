#' CitationImpute
#'
#' An implementation of the `AdaptiveImpute` algorithm using efficient
#' sparse matrix computations, specialized for the case when missing
#' values in the upper triangle are taken to be *explicitly observed*
#' zeros, as opposed to missing values. This is primarily
#' useful for spectral decompositions of adjacency matrices of graphs
#' with (near) tree structure, such as citation networks.
#'
#' @details If OpenMP is available, `citation_impute` will automatically
#'   use `getOption("Ncpus", 1L)` OpenMP threads to parallelize some
#'   key computations. Note that some computations are performed with
#'   the Armadillo C++ linear algebra library and may also be parallelized
#'   dependent on your BLAS and LAPACK installations and configurations.
#'
#' @param X A *square* sparse matrix of [Matrix::sparseMatrix()] class.
#'   Implicit zeros in the upper triangle of this matrix are considered
#'   observed and predictions on these elements contribute to the
#'   objective function minimized by `AdaptiveImpute`.
#'
#' @inherit adaptive_impute params return
#' @export
#' @include masked_approximation
#'
#' @examples
#'
#' # create a (binary) square sparse matrix to demonstrate on
#'
#' library(logger)
#'
#' log_threshold(INFO)
#'
#' set.seed(887)
#'
#' n <- 10000
#' A <- rsparsematrix(n, n, 0.1, rand.x = NULL) * 1
#' A <- as(triu(A), "generalMatrix")
#'
#' mf <- citation_impute(A, rank = 50, max_iter = 10L, check_interval = NULL)
#' mf
#'
#'
#' mf2 <- citation_impute2(A, rank = 50L, max_iter = 10L, check_interval = NULL)
#' mf2
#'
#'
citation_impute2 <- function(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize", "approximate"),
  max_iter = 200L,
  check_interval = 1L,
  epsilon = 1e-7,
  additional = NULL
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

  UseMethod("citation_impute2")
}

#' @export
citation_impute2.default <- function(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize", "approximate"),
  max_iter = 200L,
  check_interval = 1L,
  epsilon = 1e-7,
  additional = NULL) {

  stop(
    glue("No `citation_impute` method for objects of class {class(X)}."),
    call. = FALSE
  )
}

#' @export
#' @rdname citation_impute
citation_impute2.sparseMatrix <- function(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize", "approximate"),
  additional = NULL
) {

  initialization <- match.arg(initialization)

  # *explicitly* observed elements of X

  obs_upper <- nnzero(triu(X))
  obs_lower <- nnzero(tril(X, k = -1))
  obs_total <- nnzero(X)

  log_info(
    glue(
      "Matrix has {obs_total} non-zero elements, {obs_upper} in the upper ",
      "triangle (including the diagonal), and {obs_lower} in the strict ",
      "lower triangle."
    )
  )

  log_info(glue("Using {initialization} initialization."))

  if (initialization == "svd") {
    # multiply by one to coerce to type that svds can handle,
    # svds doesn't like binary matrices
    s <- svds(X * 1, rank)
    mf <- as_svd_like(s)
  } else if (initialization == "adaptive-initialize") {

    # *total* observed elements of X, including entries in the
    # upper triangle that are implicitly observed

    n <- ncol(X)  # recall that X is square
    implicit_total <- (n - 1) * (n - 2) / 2 + obs_lower
    p_hat <- implicit_total / prod(dim(X))

    mf <- adaptive_initialize(X * 1, rank, p_hat = p_hat)
  } else if (initialization == "approximate") {

    if (is.null(additional))
      stop(
        "Must specify `additional` when using approximate initialization.",
        call. = FALSE
      )

    # *total* observed elements of X, including entries in the
    # upper triangle that are implicitly observed

    n <- ncol(X)  # recall that X is square
    implicit_total <- (n - 1) * (n - 2) / 2 + obs_lower
    p_hat <- implicit_total / prod(dim(X))

    mf <- adaptive_initialize(
      X * 1, rank = rank,
      p_hat = p_hat,
      alpha_method = "approximate",
      additional = additional
    )

  } else {
    stop("This should not happen.", call. = FALSE)
  }

  log_info("Done initializing.")

  citation_impute2.LRMF(mf, X, ...)
}

#' @export
#' @rdname citation_impute
citation_impute2.LRMF <- function(
  X,
  rank,
  ...,
  epsilon = 1e-7,
  max_iter = 200L,
  check_interval = 1L
) {

  log_info(glue("Beginning AdaptiveImpute (max {max_iter} iterations)."))

  if (!is.null(check_interval))
    log_info(glue("Checking convergence every {check_interval} iteration(s)."))

  # first argument is the svd_like object, second is the data
  # do some renaming here

  s <- X
  X <- rank
  rank <- s$rank

  ### ITERATION STAGE

  delta <- Inf
  d <- ncol(X)
  norm_X <- norm(X, type = "F")^2
  iter <- 1L

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    log_debug(glue("Constructing Atilde({iter})..."))
    m <- makeCitationEstimate(X, s$u, s$d, s$v)
    log_debug(glue("Constructing Atilde({iter})... done."))

    log_debug(glue("Computing singular value decomposition of Atilde({iter})..."))
    s_new <- svds(
      left,
      k = rank,
      Atrans = right,
      dim = dim(A),
      args = m
    )
    log_debug(glue("Computing singular value decomposition of Atilde({iter})... done."))

    log_debug(glue("Computing alpha({iter})..."))
    X_tilde_f_norm <- norm_X + sum(s$d^2) -
      p_omega_z(m)
    alpha <- (X_tilde_f_norm - sum(s_new$d^2)) / (d - rank)  # line 6
    log_debug(glue("Computing alpha({iter})... done."))

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    # save a little bit on computation and only check for
    # convergence intermittently
    if (!is.null(check_interval) && iter %% check_interval == 0) {
      log_debug("Computing relative change in Frobenius norm...")
      delta <- relative_f_norm_change(s_new, s)
      log_debug("Computing relative change in Frobenius norm... done.")
    }

    s <- s_new

    log_info(
      glue(
        "Iter {iter} complete. ",
        "delta = {if (!is.null(check_interval)) delta else Inf}, ",
        "alpha = {alpha}"
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
    alpha = alpha,
    ...
  )

}
