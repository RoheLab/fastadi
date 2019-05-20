
# slightly nicer wrappers for the C++ implementations

masked_approximation <- function(s, mask) {
  mask <- as(mask, "dgTMatrix")
  masked_approximation_impl(s$u, s$d, s$v, mask@i, mask@j)
}

masked_svd_times_x <- function(s, mask, x) {
  drop(masked_svd_times_x_impl(s$u, s$d, s$v, mask@i, mask@j, x))
}

relative_f_norm_change <- function(s_new, s) {
  relative_f_norm_change_impl(s_new$u, s_new$d, s_new$v, s$u, s$d, s$v)
}
