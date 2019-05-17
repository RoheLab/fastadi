#' Expand an SVD only at observed values of a sparse matrix
#'
#' TODO: describe what it looks like for dimensions to match up between
#' `s` and `mask`. See `vignette("sparse-computations")` for mathematical
#' details.
#'
#' @param s An svd-like list with elements `u`, `d`, `v`, where `u` and `v`
#'   are low-rank matrices and `d` is a vector.
#' @param mask A sparse matrix that can be coerced into triplet form.
#'
#'
#' @return A reconstruction of `mask`A [Matrix::sparseMatrix()] of the same
#' @export
#'
#' @examples
#'
#' # U %*% diag(d) %*% t(V)
#'
masked_approximation <- function(s, mask) {

  # note: must be dgTMatrix to get column indexes j larger
  # what if we used dlTMatrix here?
  m <- as(mask, "dgTMatrix")

  # the indices for which we want to compute the matrix multiplication
  # turn zero based indices into one based indices
  i <- m@i + 1
  j <- m@j + 1

  # gets rows and columns of U and V to multiply, then multiply
  ud <- s$u %*% diag(s$d)
  left <- ud[i, ]
  right <- s$v[j, ]

  # compute inner products to get elements of U %*% t(V)
  uv <- rowSums(left * right)

  # NOTE: specify dimensions just in case
  sparseMatrix(i = i, j = j, x = uv, dims = dim(mask))
}


