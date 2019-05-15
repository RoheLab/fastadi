library(Matrix)
library(RSpectra)

fast_adi <- function(M, r, tol = 1e-03) {

  n <- nrow(M)
  d <- ncol(M)
  M <- as(M, "dgCMatrix")
  p.hat <- nnzero(M) / prod(dim(M))

  s <- fastInitial(M, r)

  frobMp <- sum(M@x^2)


  resids <- graphMult(s$u, s$v %*% diag(s$d), M, resids = TRUE)@x

  error <- Inf

  while (error > tol) {
    oldResids <- resids

    s <- threshSVD(s, M, frobMp, r)

    resids <- graphMult(s$u, s$v %*% diag(s$d), M, resids = TRUE)@x

    error <- sum((oldResids - resids)^2) / sum(resids^2)

  }
  s
}


#######################################################################################
# svdSUV and its subfunctions Ax and Atx
# are source files for taking the svd of
# M = S + U%*%V^T
# where S is sparse with nonzeros that match the Data, it is the residuals:
# S = Data - U%*%V^T  (where Data is sparse and this operation is only on the nonzeros)
# Thank you to Yixuan Qiu for help.

svdSUV <- function(U, V, Data) {
  # U and V are skinny matrices, probably the parameters.
  # Data is a sparse data matrix.
  # returns svd of
  # M = S + U%*%V^T
  # where S is sparse with nonzeros that match the Data
  # S = Data - U%*%V^T  (where Data is sparse and this operation is only on the nonzeros)

  S <- graphMult(U, V, Data, resids = TRUE)

  args <- list(S = S, U = U, V = V)

  ei <- eigs_sym(f, ncol(U), n = ncol(S), args = args)

  u <- as.matrix(Ax(ei$vectors, args))

  u <- apply(u, 2, function(x) x / sqrt(sum(x^2)))

  list(u = u, v = ei$vectors, d = sqrt(ei$values))
}

# Has the effect of A * x
Ax <- function(x, args) {
  drop(args$S %*% x + args$U %*% crossprod(args$V, x))
}

# Has the effect of A' * x
Atx <- function(x, args) {
  drop(t(args$S) %*% x + args$V %*% crossprod(args$U, x))
}

f <- function(x, args) {
  Atx(Ax(x, args), args)
}

####################################################
# threshSVD is the key internal function for fastAdi

threshSVD <- function(so, M, frobMp, r) {

  # uses sparse+lowrank svd
  # computes thresholding parameters
  # returns the thresholded SVD
  # M.p is the data
  # frobMp only needs to be computed once.
  # It is the squared frobeinus norm of M.p (frobMp = sum(M.p@x^2))
  d <- ncol(M)

  s <- svdSUV(so$u, so$v %*% diag(so$d), M)

  MTM <- frobMp + sum(so$d^2) - sum(graphMult(so$u, so$v %*% diag(so$d), M, resids = FALSE)@x^2)

  tau <- (sum(MTM) - sum(s$d^2)) / (d - r)

  s$d <- sqrt(positive(s$d^2 - tau))

  s
}

##########################################################################
# graphMult is a subfunction in several pieces for computing the residuals
# S = G - U%*%t(V)
# but only computing it on the nonzero elements of G.

# elementwise multiplication but only on the observed entries
graphMult <- function(U, V, G, resids) {
  # U and V are matrices with U %*% V^T defined.
  # G is a sparse matrix
  # if resids = F, this function computes the elements of U%*%V^T on *only* non-zero elements of G.
  # if resids = T, this function computes the elements of G - U%*%V^T on *only* non-zero elements of G.

  mt <- as(G, "dgTMatrix")
  # the indices for which we want to compute the matrix multiplication:
  i <- mt@i + 1
  j <- mt@j + 1

  # the rows of U and the columns of V for which we want to compute the multiplication:
  left <- U[i, ]
  right <- V[j, ] #  Note:  there is most certainly a more memory efficient way to do this.

  # the inner products to compute the elements of the U%*%V:
  uv <- rowSums(left * right)

  if (resids)
    sparseMatrix(i = i, j = j, x = mt@x - uv)
  else
    sparseMatrix(i = i, j = j, x = uv)
}

######################################################################

# initialization


fastInitial <- function(M.p, r) {
  ei <- eigInit(M.p, r)

  p <- nnzero(M.p) / prod(dim(M.p))
  tau <- (sum(M.p@x^2) / p - sum(ei$values)) / (ncol(M.p) - r)
  s <- list(u = c(), d = c(), v = c())

  # TODO: why positive?
  s$d <- sqrt(positive(ei$values - tau))
  s$v <- ei$vectors
  s$u <- eigInit(t(M.p), r)$vectors

  B <- svd(t(s$u) %*% (M.p %*% s$v))
  s$u <- s$u %*% B$u
  s$v <- s$v %*% B$v
  s
}

# initialization helpers

positive <- function(x) {
  x[x < 0] <- .001
  x
}

# NOTE: divide by p^2 some numeric / aesthetic trick

# Has the effect of  [ t(M)%*%M - (1-p)*diag(t(M)%*%M) ] *x
Mx <- function(x, args) {
  drop(
    crossprod(args$M, args$M %*% x) / args$p^2 -
      (1 - args$p) * Diagonal(ncol(args$M), colSums(args$M^2)) %*% x / args$p^2)
}


#' given sparse M compute eigen of
#' t(M)%*%M/p^2 - (1-p)*diag(t(M)%*%M)
#' where p is the mean nnz of M.
eigInit <- function(M, r) {
  eigs_sym(
    Mx, r,
    n = ncol(M),
    args = list(
      M = M,
      p = nnzero(M) / prod(dim(M))
    )
  )
}






set.seed(27)

n <- 500
d <- 100
r <- 5

A <- matrix(runif(n * r, -5, 5), n, r)
B <- matrix(runif(d * r, -5, 5), d, r)
M0 <- A %*% t(B)

err <- matrix(rnorm(n * d), n, d)
Mf <- M0 + err

p <- 0.3
y <- matrix(rbinom(n * d, 1, p), n, d)
dat <- Mf * y

filled <- fast_adi(dat, r)
filled


