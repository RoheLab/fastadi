#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double relative_f_norm_change_impl(
    const arma::mat& new_U,
    const arma::rowvec& new_d,
    const arma::mat& new_V,
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const int num_threads) {

  omp_set_num_threads(num_threads);

  arma::mat new_DVt = diagmat(new_d) * new_V.t();
  arma::mat DVt = diagmat(d) * V.t();
  
  double diff_norm_sq = 0;

  // expand a single row of new_Z and Z at a time
  // to avoid hitting memory limits
  
  // http://jakascorner.com/blog/2016/06/omp-data-sharing-attributes.html
  #pragma omp parallel for reduction(+: diff_norm_sq)
  for (int i = 0; i < U.n_rows; i++) {

    arma::rowvec new_Z_i, Z_i, diff;

    new_Z_i = new_U.row(i) * new_DVt;
    Z_i = U.row(i) * DVt;

    diff = new_Z_i - Z_i;
    diff_norm_sq += dot(diff, diff);
  }

  // exploit connection between frobenius norm
  // and singular values
  double z_norm_sq = accu(pow(d, 2));

  return diff_norm_sq / z_norm_sq;
}
