#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec p_u_tilde_ztx_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col,
    const arma::vec& x,
    const int num_threads) {

  omp_set_num_threads(num_threads);

  // first add the observed elements on the lower triangle

  int i, j;
  double z_ij;

  arma::vec ztx = zeros<vec>(V.n_rows);

  for (int idx = 0; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // only elements of the lower triangle + diagonal!
    if (i >= j) {
      z_ij = accu(U.row(i) % d % V.row(j));
      ztx(j) += x(i) * z_ij;
    }
  }

  return ztx;
}
