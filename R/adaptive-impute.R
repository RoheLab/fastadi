#' AdaptiveImpute
#'
#' An implementation of the `AdaptiveImpute` algorithm for matrix completion
#' for sparse matrices.
#'
#' @param X A sparse matrix of `sparseMatrix` class.
#'
#' @param rank Desired rank (integer) to use in the low rank approximation.
#'
#' @param epsilon Convergence criteria, measured in terms of relative change
#'   in Frobenius norm of the full imputed matrix. Defaults to `1e-7`.
#'
#' @param max_iter Maximum number of iterations to perform (integer). Defaults
#'   to `200L`. In practice 10 or so iterations will get you a decent
#'   approximation to use in exploratory analysis, and and 50-100 will get
#'   you most of the way to convergence.
#'
#' @return A `low_rank_matrix_factorization` object.
#'
#' @export
#'
#' @examples
#'
#' library(Matrix)
#'
#' mf <- adaptive_impute(ml100k, rank = 5L)
#' mf
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
adaptive_impute <- function(
  X,
  rank,
  initialization = c("svd", "adaptive-initialize"),
  epsilon = 1e-7,
  max_iter = 200L,
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

  if (initialization == "svd") {
    s <- svds(X, rank)
  } else if (initialization == "adaptive-initialize") {
    s <- adaptive_initialize(X, rank)  # line 1
  } else {
    stop("This should not happen.", call. = FALSE)
  }

  ### ITERATION STAGE

  # coerce M to sparse matrix such that we use sparse operations
  M <- as(X, "dgCMatrix")

  delta <- Inf
  d <- ncol(X)
  norm_M <- sum(X@x^2)

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    R <- X - masked_approximation(s, X)  # residual matrix
    args <- list(u = s$u, d = s$d, v = s$v, R = R)

    s_new <- svds(Ax, k = , Atrans = Atx, dim = dim(X), args = args)

    MtM <- norm_M + sum(s_new$d^2) - sum(masked_approximation(s_new, X)^2)

    alpha <- (sum(MtM) - sum(s_new$d^2)) / (d - rank)  # line 6

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    # NOTE: skip explicit computation of line 8
    delta <- relative_f_norm_change(s_new, s)

    s <- s_new

    if (verbose) {
      glue("delta: {round(delta, 8)}, alpha: {round(alpha, 3)}")
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

### standard variant

Ax <- function(x, args) {
  drop(args$R %*% x + args$u %*% diag(args$d) %*% crossprod(args$v, x))
}

Atx <- function(x, args) {
  drop(t(args$R) %*% x + args$v %*% diag(args$d) %*% crossprod(args$u, x))
}

### entirely observed upper triangle variant

#' AdaptiveImpute (sparse)
#'
#' A reference implemention of the `AdaptiveImpute` algorithm using sparse
#' matrix computations. For citation matrices where missing values in
#' the upper triangle are taken to be *explicitly observed* observed
#' zeros, as opposed to missing values.
#'
#' @param M A sparse Matrix created by one of the [Matrix] pkg
#'   constructors.
#' @param r Desired rank to use in the low rank approximation.
#' @param epsilon Tolerance, measured in terms of relative change
#'   in Frobenius norm of the full imputed matrix.
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
#' s <- citation_adaptive_impute(M, 5, max_iter = 20)
#' s
#'
#' # reconstruct a rank-5 approximation of M
#' s$u %*% diag(s$d) %*% t(s$v)
#'
#' # build a rank-5 approximation to M only for
#' # observed elements of M
#' masked_approximation(s, M)
#'
citation_adaptive_impute <- function(M, r, epsilon = 1e-7, max_iter = 200,
                                     additional = 100, check_in = 50) {

  # only works for square matrices
  stopifnot(ncol(M) == nrow(M))

  # NOTE: no differences are necessary from the sparse
  # adaptive_initialize since M just has more zeros

  message(Sys.time(), " Beginning AdaptiveInitialize step.")

  # this line can run into class issues
  # s <- svds(M, r)

  s <- sparse_adaptive_initialize(M, r, additional)  # line 1

  message(Sys.time(), " AdaptiveInitialize step complete.")

  delta <- Inf
  d <- ncol(M)
  f_norm_M <- sum(M@x^2)

  # coerce M to sparse matrix such that we use sparse operations
  M <- as(M * 1, "CsparseMatrix")

  iter <- 0

  while (delta > epsilon) {

    # update s: lines 4 and 5
    # take the SVD of M-tilde

    if (check_in == 1)
      message(Sys.time(), " Taking SVD.")

    args <- list(u = s$u, d = s$d, v = s$v, M = M)

    s_new <- svds(
      Ax_citation,
      k = r,
      Atrans = Atx_citation,
      dim = dim(M),
      args = args
    )

    if (check_in == 1)
      message(Sys.time(), " Finding average of remaining singular values.")

    M_tilde_f_norm <- f_norm_M + sum(s$d^2) -
      p_omega_f_norm_ut(s_new, M)

    alpha <- (M_tilde_f_norm - sum(s_new$d^2)) / (d - r)  # line 6

    s_new$d <- sqrt(s_new$d^2 - alpha)  # line 7

    # NOTE: skip explicit computation of line 8

    if (check_in == 1)
      message(Sys.time(), " Finding relative change in Frobenius norm.")

    # save a little bit on computation
    if (iter %% check_in == 0)
      delta <- relative_f_norm_change(s_new, s)

    s <- s_new

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

  s
}

# mask needs to be a sparse matrix stored as triplets
p_omega_f_norm_ut <- function(s, mask) {
  mask <- as(mask, "TsparseMatrix")
  p_omega_f_norm_ut_impl(s$u, s$d, s$v, mask@i, mask@j)
}

# KEY: must return a vector! test this otherwise you will be sad!
Ax_citation <- function(x, args) {
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

Atx_citation <- function(x, args) {

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





