#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec p_u_zx_impl(
    const arma::mat& U,
    const arma::vec& d,
    const arma::mat& V,
    const arma::vec& x,
    const int num_threads) {

  omp_set_num_threads(num_threads);

  // just DVt at this point
  arma::mat W = diagmat(d) * V.t();

  
  // multiply columns by x to obtain W
  #pragma omp parallel for
  for (int j = 0; j < W.n_cols; j++) {
    W.col(j) *= x(j);
  }

  // cumulative summation step

  // the farthest right column should should always be all zero
  // add a column of zeros to the right side of matrix
  W.insert_cols(W.n_cols, 1);

  // perform the cumulative summation from right to left
  // skip the rightmost two columns, which are already fine
  // and the leftmost column, which we will drop
  for (int j = W.n_cols - 3; j > 0; j--) {
    W.col(j) += W.col(j + 1);
  }

  // drop the leftmost column. now we have W_tilde in full
  W.shed_col(0);

  // do the matrix-vector multiplication
  arma::vec zx = zeros<vec>(U.n_rows);

  #pragma omp parallel for
  for (int j = 0; j < U.n_rows; j++) {
    zx(j) = dot(U.row(j), W.col(j));
  }

  return zx;
}
