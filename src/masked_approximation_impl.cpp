#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' Expand an SVD only at observed values of a sparse matrix
//'
//' TODO: describe what it looks like for dimensions to match up between
//' `s` and `mask`. See `vignette("sparse-computations")` for mathematical
//' details.
//'
//' @param U Low-rank matrix of left singular-ish vectors.
//' @param V Low-rank matrix of right singular-ish vectors.
//' @param row Zero-based row indices of observed elements.
//' @param col Zero-based col indices of observed elements.
//'
//' @details The idea is to populate `U`, `d` and `V` with using the
//'   elements of an SVD-like list. You can generate `row` and `col`
//'   most easily from a sparse masking Matrix (Matrix package),
//'   coercing to triplet format, and extracting `mask@i` for `row`
//'   and `mask@j` for column.
//'
//' @return A sparse matrix representing the low-rank reconstruction
//'   from `U`, `d` and `V`, only at the index pairs indicated by
//'   `row` and `col`.
//'
// [[Rcpp::export]]
arma::sp_mat masked_approximation_impl(
    const arma::mat& U,
    const arma::mat& V,
    const arma::vec& row,
    const arma::vec& col) {

  int i, j;
  arma::sp_mat reconstruction = sp_mat(U.n_rows, V.n_rows);

  for (int idx = 0; idx < row.n_elem; idx++) {

    i = row(idx);
    j = col(idx);

    // % does elementwise multiplication in Armadillo
    // accu() gives the sum of elements of resulting vector
    reconstruction(i, j) = arma::accu(U.row(i) % V.row(j));
  }

  return reconstruction;
}
