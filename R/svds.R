#
# this file contains some helpers to compute the SVD of Sigma_p
# and Sigma_t as defined in the `AdaptiveInitialize` algorithm
#
# the tl;dr here is that we truncated svds of two matrices:
#
# XtX <- crossprod(X)
# XXt <- tcrossprod(X)
#
# sigma_p <- XtX - (1 - p_hat) * Diagonal(n = ncol(X), x = diag(XtX))  # line 2
# sigma_t <- XXt - (1 - p_hat) * Diagonal(n = nrow(X), x = diag(XXt))  # line 3
#
# and at first it seems like we should be able to use RSpectra::svds() since
# X is sparse, such that sigma_p and sigma_t will be sparse. in practice,
# sigma_p and sigma_t become sufficiently dense during the
#
# XtX <- crossprod(X)
# XXt <- tcrossprod(X)
#
# operations that blindly applying RSpectra::svds() is no longer plausible.
# in one example, where X had 1.4 M non-zero elements, I found XXt had
# 100+ M non-zero elements, which takes the problem from something that you
# can do on a laptop to a pretty hard task
#
# anyway, this is a tiny bit of OOP to make these SVD computations easy
# and readable, mostly copy-pasting and slightly modifying code from
# RoheLab/sparseLRMatrix and using the standard "find an efficient matrix
# multiplication" strategy together with RSpectra
#

setClass(
  Class = "SigmaP",
  slots = c(
    X = "sparseMatrix",
    p_hat = "numeric"
  )
)

setClass(
  Class = "SigmaT",
  slots = c(
    X = "sparseMatrix",
    p_hat = "numeric"
  )
)

SigmaP <- function(X, p_hat) {
  methods::new(Class = "SigmaP", X = X, p_hat = p_hat)
}

SigmaT <- function(X, p_hat) {
  methods::new(Class = "SigmaT", X = X, p_hat = p_hat)
}

# should do a setValidity() method down the line if these get
# a lot of use

SigmaPx <- function(x, args) {
  out <- crossprod(args$X, args$X %*% x) - args$D %*% x
  drop(out)
}

SigmaTx <- function(x, args) {
  out <- args$X %*% crossprod(args$X, x) - args$D %*% x
  drop(out)
}

#' @export
svds.SigmaP <- function(A, k, nu = k, nv = k, opts = list(), ...) {

  # maaaaaybe worth doing something more efficient here, i don't know
  XtX <- crossprod(A@X)

  args <- list(
    X = A@X,
    D = (1 - A@p_hat) * Diagonal(n = ncol(A@X), x = diag(XtX))
  )

  RSpectra::svds(
    SigmaPx,
    k,
    nu = nu,
    nv = nv,
    opts = opts,
    Atrans = SigmaPx,  # symmetric!
    dim = dim(XtX),
    args = args,
    ...
  )
}

#' @export
svds.SigmaT <- function(A, k, nu = k, nv = k, opts = list(), ...) {

  # maaaaaybe worth doing something more efficient here, i don't know
  XXt <- tcrossprod(A@X)

  args <- list(
    X = A@X,
    D = (1 - A@p_hat) * Diagonal(n = nrow(A@X), x = diag(XXt))
  )

  RSpectra::svds(
    SigmaTx,
    k,
    nu = nu,
    nv = nv,
    opts = opts,
    Atrans = SigmaTx, # also symmetric!
    dim = dim(XXt),
    args = args,
    ...
  )
}
