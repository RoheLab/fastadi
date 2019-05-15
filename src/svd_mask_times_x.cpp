#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec masked_svd_times_x_cpp(
    const arma::mat& U,
    const arma::vec& d,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col,
    const arma::vec& x) {

  int i, j;
  double z_ij;

  arma::vec zx = arma::zeros<arma::vec>(U.n_rows);

  for (int idx; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // % does elementwise multiplication in Armadillo
    // accu() gives the sum of elements of resulting vector
    z_ij = accu(U.row(i) % d % V.row(j));
    zx(i) += x(j) * z_ij;
  }

  return zx;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
library(Matrix)
library(RSpectra)

M <- Matrix(
  rbind(
    c(0, 0, 3, 1, 0),
    c(3, 0, 0, 8, 0),
    c(0, -1, 0, 0, 0),
    c(0, 0, 0, 0, 0),
    c(0, 2, 0, 0, 0),
    c(5, 0, 7, 0, 4)
  )
)

Y <- rbind(
  c(1, 1, 1, 1, 1),
  c(1, 1, 1, 1, 1),
  c(0, 1, 1, 1, 1),
  c(0, 0, 1, 1, 1),
  c(0, 1, 0, 1, 1),
  c(1, 0, 1, 0, 1)
)

s <- svds(M, 2)

Y <- as(Y, "CsparseMatrix")

# triplet form
# compressed column matrix form even better but don't
# understand the format
Y <- as(Y, "lgCMatrix")
Y <- as(Y, "lgTMatrix")

Y

x <- rnorm(5)

# want to calculate
Z <- s$u %*% diag(s$d) %*% t(s$v)
out <- drop((Z * Y) %*% x)
out

# mask as a pair list
# L and Z / svd are both n x d matrices
# x is a d x 1 matrix / vector
masked_svd_times_x <- function(s, mask, x) {

  stopifnot(inherits(mask, "lgTMatrix"))

  u <- s$u
  d <- s$d
  v <- s$v

  zx <- numeric(nrow(u))

  # lgTMatrix uses zero based indexing, add one
  row <- mask@i + 1
  col <- mask@j + 1

  # need to loop over index of indexes
  # double looping over i and j here feels intuitive
  # but is incorrect
  for (idx in seq_along(row)) {
    i <- row[idx]
    j <- col[idx]

    z_ij <- sum(u[i, ] * d * v[j, ])
    zx[i] <- zx[i] + x[j] * z_ij
  }

  zx
}

# how to calculate just one element of the reconstructed
# data using the SVD

i <- 6
j <- 4

sum(s$u[i, ] * s$d * s$v[j, ])
Z[i, j]

# the whole masked matrix multiply

Z <- s$u %*% diag(s$d) %*% t(s$v)
out <- drop((Z * Y) %*% x)

# check that we did this right
all.equal(
  masked_svd_times_x(s, Y, x),
  out
)
*/
