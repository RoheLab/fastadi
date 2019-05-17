set.seed(17)
M <- rsparsematrix(8, 12, nnz = 30) # small example, not very sparse

# number of singular vectors to compute
k <- 4

s <- svd(M, k, k)
s2 <- svds(M, k, k)

# irritating: svd() always gives you all the singular values even if you
# only request the first K singular vectors
s$u %*% diag(s$d[1:k]) %*% t(s$v)

# based on the flip_signs function of
# https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
equal_svds <- function(s, s2) {

  # svd() always gives you all the singular values, but we only
  # want to compare the first k
  k <- ncol(s$u)

  # the term sign(s$u) * sign(s2$u) performs a sign correction

  # isTRUE because output of all.equal is not a boolean, it's something
  # weird when the inputs aren't equal. lol why

  u_ok <- isTRUE(
    all.equal(s$u, s2$u * sign(s$u) * sign(s2$u), check.attributes = FALSE)
  )

  v_ok <- isTRUE(
    all.equal(s$v, s2$v * sign(s$v) * sign(s2$v), check.attributes = FALSE)
  )

  d_ok <- isTRUE(all.equal(s$d[1:k], s2$d[1:k], check.attributes = FALSE))

  u_ok && v_ok && d_ok
}

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

init <- adaptive_initialize(dat, r)
filled <- adaptive_impute(dat, r)



is_svd_like <- function(s, check_orthonormal = FALSE) {

  # must have u, v, d

  # original data is n x d

  # U is n x r
  # d is length(r)
  # V is r by d

  # if check orthonormal, make sure that inverse and transpose
  # are the same
}

#' Check that dimensions of a svd like object agree with matrix dimensions
#'
#' @param s An SVD
#' @param A
#'
#' @return
#' @export
#'
#' @examples
svd_like_matches_matrix_dims <- function(s, A) {

}
