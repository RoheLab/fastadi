#' Create an Adaptive Imputation object
#'
#' `adaptive_imputation` objects are a subclass of
#' [LRMF3::svd_like()], with an additional field `alpha`.
#'
#' @param u A *matrix* "left singular-ish" vectors.
#'
#' @param d A *vector* of "singular-ish" values.
#'
#' @param v A *matrix* of "right singular-ish" vectors.
#'
#' @param alpha Value of `alpha` after final iteration.
#'
#' @param ... Optional additional items to pass to the constructor.
#'
#' @return An `adaptive_imputation` object.
#'
adaptive_imputation <- function(u, d, v, alpha, ...) {

  ai <- svd_like(
    u = u,
    d = d,
    v = v,
    subclasses = "adaptive_imputation",
    alpha = alpha,
    ...
  )

  validate_adaptive_imputation(ai)
}

new_adaptive_imputation <- function(u, d, v, rank, alpha, ...) {
  svd_like(
    u = u,
    d = d,
    v = v,
    rank = rank,
    subclasses = "adaptive_imputation",
    alpha = alpha,
    ...
  )
}

validate_adaptive_imputation <- function(ai) {

  if (is.null(ai$alpha)) {
    stop(
      "Must have `alpha` field in adaptive imputation object.",
      call. = FALSE
      )
  }

  if (!is.numeric(ai$alpha) || length(ai$alpha) != 1) {
    stop(
      "`alpha` must be a numeric vector of length 1.",
      call. = FALSE
    )
  }

  ai
}

#' @importFrom LRMF3 dim_and_class
#' @method print adaptive_imputation
#' @export
print.adaptive_imputation <- function(x, ...) {
  cat("\nAdaptively Imputed Low Rank Matrix Factorization\n")
  cat("------------------------------------------------\n\n")

  cat(glue("Rank: {x$rank}\n\n", .trim = FALSE))

  cat(glue("Rows: {nrow(x$u)}\n", .trim = FALSE))
  cat(glue("Cols: {nrow(x$v)}\n\n", .trim = FALSE))

  cat(glue("d[rank]: {round(x$d[x$rank], 3)}\n", .trim = FALSE))
  cat(glue("alpha:   {round(x$alpha, 3)}\n\n", .trim = FALSE))

  cat("Components\n\n")

  cat("u:", dim_and_class(x$u), "\n")
  cat("d:", dim_and_class(x$d), "\n")
  cat("v:", dim_and_class(x$v), "\n")
}
