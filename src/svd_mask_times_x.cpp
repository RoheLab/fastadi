#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec masked_svd_times_x_impl(
    const mat& U,
    const rowvec& d,
    const mat& V,
    const vec& row,
    const vec& col,
    const vec& x) {

  int i, j;
  double z_ij;

  vec zx = zeros<vec>(U.n_rows);

  for (int idx = 0; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // % does elementwise multiplication in Armadillo
    // accu() gives the sum of elements of resulting vector
    z_ij = accu(U.row(i) % d % V.row(j));

    zx(i) += x(j) * z_ij;
  }

  return zx;
}
