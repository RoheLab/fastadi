#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double p_omega_f_norm_ut_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col) {

  // first add the observed elements on the lower triangle and diagonal

  int i, j;
  double z_ij, total = 0;
  arma::mat DVt = diagmat(d) * V.t();

  for (int idx = 0; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // only elements of the lower triangle
    // KEY: include the diagonal, but only *explicitly* observed elements
    if (i >= j) {
      z_ij = arma::dot(U.row(i), DVt.col(j));
      total += pow(z_ij, 2);
    }
  }

  // second add all the elements of the upper triangle
  int n = U.n_rows;

  arma::rowvec Z_i_trunc;

  for (int i = 0; i < n - 1; i++) {
    Z_i_trunc = U.row(i) * DVt.cols(i + 1, n - 1);
    total += dot(Z_i_trunc, Z_i_trunc);
  }

  return total;
}
