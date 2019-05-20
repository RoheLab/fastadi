#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec masked_svd_times_x_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col,
    const arma::vec& x) {

  int i, j;
  double z_ij;

  arma::vec zx = zeros<vec>(U.n_rows);

  for (int idx = 0; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // % does elementwise multiplication in Armadillo
    // accu() gives the sum of elements of resulting vector
    z_ij = arma::accu(U.row(i) % d % V.row(j));

    zx(i) += x(j) * z_ij;
  }

  return zx;
}
