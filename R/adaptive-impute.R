#' AdaptiveImpute
#'
#' An implementation of the `AdaptiveImpute` algorithm for matrix completion
#' for sparse matrices.
#'
#' @param X A sparse matrix of `sparseMatrix` class.
#'
#' @param rank Desired rank (integer) to use in the low rank approximation.
#'
#' @param initialization TODO.
#'
#' @param epsilon Convergence criteria, measured in terms of relative change
#'   in Frobenius norm of the full imputed matrix. Defaults to `1e-7`.
#'
#' @param max_iter Maximum number of iterations to perform (integer). Defaults
#'   to `200L`. In practice 10 or so iterations will get you a decent
#'   approximation to use in exploratory analysis, and and 50-100 will get
#'   you most of the way to convergence.
#'
#' @param verbose TODO.
#'
#' @return A `low_rank_matrix_factorization` object.
#'
#' @export
#'
#' @examples
#'
#' mf <- adaptive_impute(ml100k, rank = 5L, max_iter = 10L, verbose = TRUE)
#' mf
#'
#' # build a rank-5 approximation only for
#' # observed elements of ml100k
#' masked_approximation(mf, ml100k)
#'
adaptive_impute <- function(
  X,
  rank,
  initialization = c("svd", "adaptive-initialize"),
  epsilon = 1e-7,
  max_iter = 200L,
  check_in = 1L,
  verbose = FALSE
) {

  # TODO: use parallel matrix multiplications from rsparse
  # where possible

  ### INPUT VALIDATION

  if (!inherits(X, "sparseMatrix")) {
    stop(
      glue("`X` must be a `sparseMatrix, not a {class(X)} object."),
      call. = FALSE
    )
  }

  rank <- as.integer(rank)

  if (rank <= 2)
    stop("`rank` must be an integer >= 2L.", call. = FALSE)

  initialization <- match.arg(initialization)

  ### INITIALIZATION

  if (verbose) {
    message(glue("Initialization with {initialization} strategy."))
  }

  if (initialization == "svd") {
    s <- svds(X, rank)
  } else if (initialization == "adaptive-initialize") {
    s <- adaptive_initialize(X, rank)  # line 1
  } else {
    stop("This should not happen.", call. = FALSE)
  }

  if (verbose) {
    message("Done initializing, beginning iterations.")
  }

  ### ITERATION STAGE

  delta <- Inf
  d <- ncol(X)
  norm_M <- norm(X, type = "F")^2
  iter <- 0L

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    R <- X - masked_approximation(s, X)  # residual matrix
    args <- list(u = s$u, d = s$d, v = s$v, R = R)

    s_new <- svds(Ax, k = rank, Atrans = Atx, dim = dim(X), args = args)

    MtM <- norm_M + sum(s_new$d^2) - sum(masked_approximation(s_new, X)^2)

    alpha <- (sum(MtM) - sum(s_new$d^2)) / (d - rank)  # line 6

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    # NOTE: skip explicit computation of line 8
    delta <- relative_f_norm_change(s_new, s)

    if (check_in == 1)
      message(Sys.time(), " Finding relative change in Frobenius norm.")

    # save a little bit on computation
    if (iter %% check_in == 0)
      delta <- relative_f_norm_change(s_new, s)

    s <- s_new

    if (verbose) {
      glue("delta: {round(delta, 8)}, alpha: {round(alpha, 3)}")
    }

    if (iter %% check_in == 0)
      message(
        Sys.time(),
        glue::glue(
          " Iter {iter} complete. ",
          "delta = {round(delta, 8)}, ",
          "alpha = {round(alpha, 3)}"
        )
      )

    iter <- iter + 1

    if (iter > max_iter) {
      warning(
        "\nReached maximum allowed iterations. Returning early.",
        call. = FALSE
      )
      break
    }

  }

  new_low_rank_matrix_factorization(
    U = s$u,
    d = s$d,
    V = s$v,
    rank = rank,
    alpha = alpha
  )
}

Ax <- function(x, args) {
  drop(args$R %*% x + args$u %*% diag(args$d) %*% crossprod(args$v, x))
}

Atx <- function(x, args) {
  drop(t(args$R) %*% x + args$v %*% diag(args$d) %*% crossprod(args$u, x))
}





