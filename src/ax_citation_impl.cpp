#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec ax_citation_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& x) {

  double z_ij;
  arma::vec zx = zeros<vec>(U.n_rows);

  // i indexes row, j indexes column
  for (int i = 0; i < U.n_rows; i++) {
    for (int j = i + 1; j < V.n_rows; j++) {
      z_ij = arma::accu(U.row(i) % d % V.row(j));
      zx(i) += x(j) * z_ij;
    }
  }

  return zx;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec atx_citation_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const arma::vec& x) {

  double zt_ij;
  arma::vec ztx = zeros<vec>(V.n_rows);

  // i indexes row, j indexes column
  for (int i = 0; i < U.n_rows; i++) {
    for (int j = i + 1; j < V.n_rows; j++) {
      zt_ij = arma::accu(U.row(j) % d % V.row(i));
      ztx(i) += x(j) * zt_ij;
    }
  }

  return ztx;
}
