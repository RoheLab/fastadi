#include <RcppArmadillo.h>

using namespace arma;

//' @param row Row indices of original mask before transposition
//' @param col Col indices of original mask before transposition
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec p_omega_ztx_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col,
    const arma::vec& x) {

  // first add the observed elements on the upper triangle
  // since we're working with Z_t^T now!

  int i, j;
  double z_ij;

  arma::vec ztx = zeros<vec>(V.n_rows);

  for (int idx = 0; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // only elements of the (strict) upper triangle
    if (i > j) {
      z_ij = arma::accu(U.row(i) % d % V.row(j));
      ztx(j) += x(i) * z_ij;
    }
  }

  // second add all the elements of the lower triangle

  // i indexes row, j indexes column
  for (int i = 0; i < U.n_rows; i++) {

    // j = i + 1 ensures we don't add elements from the diagonal
    for (int j = i + 1; j < V.n_rows; j++) {
      z_ij = arma::accu(U.row(i) % d % V.row(j));

      // switch i, j from at_citation setup
      ztx(j) += x(i) * z_ij;
    }
  }

  return ztx;
}
