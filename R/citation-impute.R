#' CitationImpute
#'
#' An implementation of the `AdaptiveImpute` algorithm using efficient
#' sparse matrix computations, specialized for the case when missing
#' values in the upper triangle are taken to be *explicitly observed*
#' zeros, as opposed to missing values. This is primarily
#' useful for spectral decompositions of adjacency matrices of graphs
#' with (near) tree structure, such as citation networks.
#'
#' @param X A *square* sparse matrix of [Matrix::sparseMatrix()] class.
#'   Implicit zeros in the upper triangle of this matrix are considered
#'   observed and predictions on these elements contribute to the
#'   objective function minimized by `AdaptiveImpute`.
#'
#' @inherit adaptive_impute params return
#' @export
#'
#' @examples
#'
#' # create a (binary) square sparse matrix to demonstrate on
#'
#' n <- 100
#' A <- rsparsematrix(n, n, 0.1, rand.x = NULL)
#'
#' mf <- citation_impute(A, rank = 3L, max_iter = 10L)
#'
#' mf2 <- citation_impute(
#'   A,
#'   rank = 3L,
#'   max_iter = 10L,
#'   initialization = "adaptive-initialize"
#' )
#'
citation_impute <- function(
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

  if (length(check_interval) > 1)
    stop(
      "`check_interval` must be an integer vector with a single element.",
      call. = FALSE
    )

  if (rank <= 2)
    stop("`rank` must be an integer >= 2L.", call. = FALSE)

  if (max_iter <= 2)
    stop("`max_iter` must be an integer >= 2L.", call. = FALSE)

  if (check_interval < 1)
    stop("`check_interval` must be an integer >= 1L.", call. = FALSE)

  UseMethod("citation_impute")
}

#' @export
citation_impute.default <- function(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize"),
  max_iter = 200L,
  check_interval = 1L,
  epsilon = 1e-7) {

  stop(
    glue("No `citation_impute` method for objects of class {class(X)}."),
    call. = FALSE
  )
}

#' @export
#' @rdname citation_impute
citation_impute.sparseMatrix <- function(
  X,
  rank,
  initialization = c("svd", "adaptive-initialize"),
  ...
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

  log_info(glue("Use {initialization} initialization."))

  if (initialization == "svd") {
    # multiply by one to coerce to type that svds can handle,
    # svds doesn't like binary matrices
    s <- svds(X * 1, rank)
    mf <- as_svd_like(s)
  } else if (initialization == "adaptive-initialize") {

    # *total* observed elements of X, including entries in the
    # upper triangle that are implicitly observed

    n <- ncol(X)  # recall that X is square
    implicit_total <- n * (n - 1) / 2 + obs_lower
    p_hat <- implicit_total / prod(dim(X))

    mf <- adaptive_initialize(X * 1, rank, p_hat = p_hat)
  } else {
    stop("This should not happen.", call. = FALSE)
  }

  log_info("Done initializing.")

  citation_impute.LRMF(mf, X, ...)
}

#' @export
#' @rdname citation_impute
citation_impute.LRMF <- function(
  X,
  rank,
  ...,
  epsilon = 1e-7,
  max_iter = 200L,
  check_interval = 1L,
  verbose = FALSE
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
  norm_X <- norm(X, type = "F")^2
  iter <- 1L

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    args <- list(u = s$u, d = s$d, v = s$v, M = X)

    s_new <- svds(
      Ax_citation,
      k = rank,
      Atrans = Atx_citation,
      dim = dim(X),
      args = args
    )

    X_tilde_f_norm <- norm_X + sum(s$d^2) -
      p_omega_f_norm_ut(s_new, X)

    alpha <- (X_tilde_f_norm - sum(s_new$d^2)) / (d - rank)  # line 6

    assert_alpha_positive(alpha)

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    # NOTE: skip explicit computation of line 8
    delta <- relative_f_norm_change(s_new, s)

    # save a little bit on computation
    if (iter %% check_interval == 0) {
      log_info("Finding relative change in Frobenius norm.")
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
    alpha = alpha,
    ...
  )

}
