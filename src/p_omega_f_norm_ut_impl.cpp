#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double p_omega_f_norm_ut_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col,
    const int num_threads) {

  omp_set_num_threads(num_threads);

  // first add the observed elements on the lower triangle and diagonal

  arma::mat DVt = diagmat(d) * V.t();
  double f_norm_sq = 0;

  #pragma omp parallel for reduction(+: f_norm_sq)
  for (int idx = 0; idx < row.n_elem; idx++) {

    int i, j;
    double z_ij;

    i = row(idx);
    j = col(idx);

    // only elements of the lower triangle
    // KEY: include the diagonal, but only *explicitly* observed elements
    if (i >= j) {
      z_ij = dot(U.row(i), DVt.col(j));
      f_norm_sq += pow(z_ij, 2);
    }
  }

  // second add all the elements of the upper triangle

  int r = U.n_cols;

  #pragma omp parallel for reduction(+: f_norm_sq)
  for (int l = 0; l < r; l++) {
    for (int q = 0; q < r; q++) {

      arma::vec U_lq;
      arma::rowvec V_lq, V_lq_tri;

      U_lq = U.col(l) % U.col(q);
      V_lq = DVt.row(l) % DVt.row(q);
      V_lq_tri = sum(V_lq) - cumsum(V_lq);

      f_norm_sq += dot(U_lq, V_lq_tri);
    }
  }

  return f_norm_sq;
}
