#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec p_u_ztx_impl(
    const arma::mat& U,
    const arma::vec& d,
    const arma::mat& V,
    const arma::vec& x,
    const int num_threads) {

  omp_set_num_threads(num_threads);

  // just UD at this point
  arma::mat S = U * diagmat(d);

  // multiply rows by x to obtain W
  #pragma omp parallel for
  for (int i = 0; i < S.n_rows; i++) {
    S.row(i) *= x(i);
  }

  // cumulative summation step

  S = cumsum(S); // column-wise by default
  S.insert_rows(0, 1); // add 1 row of zeros at zeroth row index
  S.shed_row(S.n_rows - 1);

  // do the matrix-vector multiplication
  arma::vec ztx = zeros<vec>(V.n_rows);

  #pragma omp parallel for
  for (int i = 0; i < V.n_rows; i++) {
    ztx(i) = dot(S.row(i), V.row(i));
  }

  return ztx;
}
