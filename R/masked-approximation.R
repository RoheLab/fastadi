masked_approximation <- function(s, mask) {
  mask <- methods::as(mask, "TsparseMatrix")
  masked_approximation_impl(s$u, s$v %*% diag(s$d), mask@i, mask@j)
}