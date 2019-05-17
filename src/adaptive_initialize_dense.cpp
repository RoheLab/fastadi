#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List AdaptiveInitialize(const sp_mat& M, const int r) {

  // coerce to double to avoid integer division
  double p_hat = static_cast<double>(M.n_nonzero) / (M.n_cols * M.n_rows);

  sp_mat MtM = M.t() * M;
  sp_mat MMt = M * M.t();

  sp_mat sigma_p = MtM / pow(p_hat, 2) - (1 - p_hat) * diagmat(MtM);
  sp_mat sigma_t = MMt / pow(p_hat, 2) - (1 - p_hat) * diagmat(MMt);

  mat U_p, V_p, U_t, V_t;
  vec s_p, s_t;

  svds(U_p, s_p, V_p, sigma_p, r);
  svds(U_t, s_t, V_t, sigma_t, r);

  // TODO: this is still the eigenvalue calculation
  double alpha = (sum(sigma_p.diag()) - sum(s_p)) / (M.n_cols - r);

  vec lambda_hat = sqrt(s_p - alpha) / p_hat;

  mat U_m, V_m;
  vec s_m;

  svds(U_m, s_m, V_m, M, r);

  // sum(A % B) finds the diag(A^T B) / the diagonal of the cross product
  // sum() does a *column-wise sum*, % an elementwise multiplication
  rowvec u_sign = sign(sum(U_m % U_t));
  rowvec v_sign = sign(sum(V_m % V_p));

  rowvec s_hat = u_sign % v_sign;
  lambda_hat = lambda_hat % conv_to< vec >::from(s_hat);

  return Rcpp::List::create(Rcpp::Named("u") = U_t,
                            Rcpp::Named("d") = lambda_hat,
                            Rcpp::Named("v") = V_p);
}
