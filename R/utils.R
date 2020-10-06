masked_svd_times_x <- function(s, mask, x) {
  drop(masked_svd_times_x_impl(s$u, s$d, s$v, mask@i, mask@j, x))
}

relative_f_norm_change <- function(s_new, s) {

  # expressed into terms of frobenius inner products of low rank
  # SVDs

  num <- sum(s_new$d^2) + sum(s$d^2) -
    2 * svd_frob_inner_prod_impl(s_new$u, s_new$d, s_new$v, s$u, s$d, s$v)

  denom <- sum(s$d^2)

  num / denom
}

assert_alpha_positive <- function(alpha) {
  if (alpha < 0)
    stop(
      "Negative alpha value. This should never happen. Please report ",
      "this issue at https://github.com/RoheLab/fastadi/issues/, ideally ",
      "as a reproducible example (see https://reprex.tidyverse.org/).",
      call. = FALSE
    )
}
