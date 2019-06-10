#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double svd_frob_inner_prod_impl(
    const arma::mat& new_U,
    const arma::rowvec& new_d,
    const arma::mat& new_V,
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V) {

  // assumes dimensions of new and old SVD match up, both internally
  // and between the SVDs
  int r = new_U.n_cols;

  arma::mat new_DVt = diagmat(new_d) * new_V.t();
  arma::mat DVt = diagmat(d) * V.t();

  arma::vec U_lq;
  arma::rowvec V_lq;

  double frob_inner_prod;

  for (int l = 0; l < r; l++) {
    for (int q = 0; q < r; q++) {
      U_lq = new_U.col(l) % U.col(q);
      V_lq = new_DVt.row(l) % DVt.row(q);

      frob_inner_prod += sum(U_lq) * sum(V_lq);
    }
  }

  return frob_inner_prod;
}
