#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List AdaptiveInitialize(const arma::sp_mat& M, const int r) {

  // coerce to double to avoid integer division
  double p_hat = static_cast<double>(M.n_nonzero) / (M.n_cols * M.n_rows);

  arma::sp_mat MtM = M.t() * M;
  arma::sp_mat MMt = M * M.t();

  arma::sp_mat sigma_p = MtM / pow(p_hat, 2) - (1 - p_hat) * diagmat(MtM);
  arma::sp_mat sigma_t = MMt / pow(p_hat, 2) - (1 - p_hat) * diagmat(MMt);

  arma::mat U_p, V_p, U_t, V_t;
  arma::vec s_p, s_t;

  arma::svds(U_p, s_p, V_p, sigma_p, r);
  arma::svds(U_t, s_t, V_t, sigma_t, r);

  // TODO: this is still the eigenvalue calculation
  double alpha = (sum(sigma_p.diag()) - sum(s_p)) / (M.n_cols - r);

  arma::vec lambda_hat = sqrt(s_p - alpha) / p_hat;

  arma::mat U_m, V_m;
  arma::vec s_m;

  arma::svds(U_m, s_m, V_m, M, r);

  // sum(A % B) finds the diag(A^T B) / the diagonal of the cross product
  // sum() does a *column-wise sum*, % an elementwise multiplication
  arma::rowvec u_sign = arma::sign(arma::sum(U_m % U_t));
  arma::rowvec v_sign = arma::sign(arma::sum(V_m % V_p));

  arma::rowvec s_hat = u_sign % v_sign;
  lambda_hat = lambda_hat % arma::conv_to< arma::vec >::from(s_hat);

  return Rcpp::List::create(Rcpp::Named("u") = U_t,
                            Rcpp::Named("d") = lambda_hat,
                            Rcpp::Named("v") = V_p);
}
