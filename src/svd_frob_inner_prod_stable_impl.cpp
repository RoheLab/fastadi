#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double svd_frob_inner_prod_stable_impl(
    arma::mat& new_U,
    arma::rowvec& new_d,
    arma::mat& new_V,
    arma::mat& U,
    arma::rowvec& d,
    arma::mat& V) {

  // issue: want to avoid overflow and underflow.
  // solution: scale by maximum
  // https://timvieira.github.io/blog/post/2014/11/10/numerically-stable-p-norms/

  double new_alpha = 2; // new_U.max() * new_d.max() * new_V.max();
  double alpha = 2; // U.max() * d.max() * V.max();

  // TODO: does this mutate the original R matrix object??
  new_U = new_U / new_alpha;
  new_d = new_d; // / new_alpha;
  new_V = new_V; // / new_alpha;

  U = U / alpha;
  d = d; /// alpha;
  V = V; // / alpha;

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

  // don't forget to rescale at the end
  return alpha * new_alpha * frob_inner_prod;
}
