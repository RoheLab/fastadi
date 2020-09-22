#' @export
new_adaptive_imputation <- function(
  u, d, v, rank, alpha, subclasses = NULL, ...) {

  object <- list(
    u = u,
    d = d,
    v = v,
    rank = rank,
    alpha = alpha,
    ...
  )

  class(object) <- c(subclasses, "adaptive_imputation", "LRMF")
  object
}

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

  dim_and_class <- function(x) {
    if (is.vector(x))
      paste0(length(x), "      [", class(x)[1], "]")
    else
      # is a matrix
      paste0(nrow(x), " x ", ncol(x), " [", class(x)[1], "]")
  }

  cat("u:", dim_and_class(x$u), "\n")
  cat("d:", dim_and_class(x$d), "\n")
  cat("v:", dim_and_class(x$v), "\n")
}
