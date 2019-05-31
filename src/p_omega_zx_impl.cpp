#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec p_omega_zx_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col,
    const arma::vec& x) {

  // first add the observed elements on the lower triangle

  int i, j;
  double z_ij;

  arma::vec zx = zeros<vec>(U.n_rows);

  for (int idx = 0; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // only elements of the lower triangle
    if (i >= j) {
      z_ij = arma::accu(U.row(i) % d % V.row(j));
      zx(i) += x(j) * z_ij;
    }
  }

  // second add all the elements of the upper triangle

  // i indexes row, j indexes column
  for (int i = 0; i < U.n_rows; i++) {
    for (int j = i + 1; j < V.n_rows; j++) {
      z_ij = arma::accu(U.row(i) % d % V.row(j));
      zx(i) += x(j) * z_ij;
    }
  }

  return zx;
}
