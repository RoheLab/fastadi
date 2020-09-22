
### adaptive impute svd variant

Ax <- function(x, args) {
  drop(args$R %*% x + args$u %*% diag(args$d) %*% crossprod(args$v, x))
}

Atx <- function(x, args) {
  drop(t(args$R) %*% x + args$v %*% diag(args$d) %*% crossprod(args$u, x))
}

### citation svd variant

### citation impute svd variant

citation_svds <- function(X, rank) {
  svds(
    Ax_citation_svds,
    k = rank,
    Atrans = Atx_citation_svds,
    dim = dim(X),
    args = args
  )
}

# mask needs to be a sparse matrix stored as triplets
p_omega_f_norm_ut <- function(s, mask) {
  mask <- as(mask, "TsparseMatrix")
  p_omega_f_norm_ut_impl(s$u, s$d, s$v, mask@i, mask@j)
}

# KEY: must return a vector! test this otherwise you will be sad!
Ax_citation <- function(x, args) {
  mask <- as(args$M, "TsparseMatrix")

  out <- drop0(args$M) %*% x
  out <- out - p_u_zx_impl(args$u, args$d, args$v, x)
  out <- out - p_u_tilde_zx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  out <- out + args$u %*% diag(args$d) %*% crossprod(args$v, x)

  drop(out)
}

Atx_citation <- function(x, args) {

  mask <- as(args$M, "TsparseMatrix")

  out <- t(drop0(args$M)) %*% x
  out <- out - p_u_ztx_impl(args$u, args$d, args$v, x)
  out <- out - p_u_tilde_ztx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  out <- out + args$v %*% diag(args$d) %*% crossprod(args$u, x)

  drop(out)
}
