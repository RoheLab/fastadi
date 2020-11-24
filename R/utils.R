relative_f_norm_change <- function(s_new, s) {
  relative_f_norm_change_impl(
    s_new$u, s_new$d, s_new$v,
    s$u, s$d, s$v,
    getOption("Ncpus", 1L)
  )
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


### adaptive impute svd variant

## TODO: how to parallelize?

Ax <- function(x, args) {
  drop(args$R %*% x + args$u %*% diag(args$d) %*% crossprod(args$v, x))
}

Atx <- function(x, args) {
  drop(t(args$R) %*% x + args$v %*% diag(args$d) %*% crossprod(args$u, x))
}

### citation svd variant

### citation impute svd variant

# mask needs to be a sparse matrix stored as triplets
p_omega_f_norm_ut <- function(s, mask) {
  mask <- as(mask, "TsparseMatrix")
  p_omega_f_norm_ut_impl(s$u, s$d, s$v, mask@i, mask@j, getOption("Ncpus", 1L))
}

# KEY: must return a vector! test this otherwise you will be sad!
Ax_citation <- function(x, args) {
  mask <- as(args$M, "TsparseMatrix")

  out <- args$M %*% x
  out <- out - p_u_zx_impl(args$u, args$d, args$v, x, getOption("Ncpus", 1L))
  out <- out - p_u_tilde_zx_impl(args$u, args$d, args$v, mask@i, mask@j, x, getOption("Ncpus", 1L))
  out <- out + args$u %*% diag(args$d) %*% crossprod(args$v, x)

  drop(out)
}

Atx_citation <- function(x, args) {

  mask <- as(args$M, "TsparseMatrix")

  out <- t(args$M) %*% x
  out <- out - p_u_ztx_impl(args$u, args$d, args$v, x, getOption("Ncpus", 1L))
  out <- out - p_u_tilde_ztx_impl(args$u, args$d, args$v, mask@i, mask@j, x, getOption("Ncpus", 1L))
  out <- out + args$v %*% diag(args$d) %*% crossprod(args$u, x)

  drop(out)
}
