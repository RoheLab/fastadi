#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double relative_f_norm_change_impl(
    const arma::mat& new_U,
    const arma::rowvec& new_d,
    const arma::mat& new_V,
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V) {

  double new_z_ij, z_ij;

  double z_norm_sq = 0;
  double diff_norm_sq = 0;

  // i indexes row, j indexes column
  for (int i = 0; i < U.n_rows; i++) {
    for (int j = 0; j < V.n_rows; j++) {

      new_z_ij = arma::accu(new_U.row(i) % new_d % new_V.row(j));
      z_ij = arma::accu(U.row(i) % d % V.row(j));

      diff_norm_sq += pow(new_z_ij - z_ij, 2);
    }
  }

  // this is unhappy
  z_norm_sq = sum(pow(d, 2));

  return diff_norm_sq / z_norm_sq;
}
