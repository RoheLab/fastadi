new_low_rank_matrix_factorization <- function(
  U, d, V, rank, alpha) {

  object <- list(
    U = U,
    d = d,
    V = V,
    rank = rank,
    alpha =  alpha
  )

  class(object) <- "low_rank_matrix_factorization"
  object
}

#' @export
print.low_rank_matrix_factorization <- function(x, ...) {
  cat("Low Rank Matrix Factorization\n\n")


  cat(glue("Rank: {x$rank}"), sep = "\n\n")

  cat(glue("Rows: {nrow(x$U)}"), sep = "\n")
  cat(glue("Cols: {nrow(x$V)}"), sep = "\n\n")

  cat(glue("d[rank]: {x$d[x$rank]}"), sep = "\n")
  cat(glue("alpha:   {x$alpha}"), sep = "\n")

  cat("\nPre-Processing\n\n")

  cat(" - Centered: FALSE \n")
  cat(" - Scaled:   FALSE \n\n")

  cat("Components\n\n")

  dim_and_class <- function(x) {
    if (is.vector(x))
      paste0(length(x), "      [", class(x)[1], "]")
    else
      # is a matrix
      paste0(nrow(x), " x ", ncol(x), " [", class(x)[1], "]")
  }

  cat("U:", dim_and_class(x$U), "\n")
  cat("d:", dim_and_class(x$d), "\n")
  cat("V:", dim_and_class(x$V), "\n")
}
