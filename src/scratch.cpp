#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

class CitationEstimate {
  public:
    arma::mat U, V;
    arma::rowvec d;
    arma::sp_mat A;
    CitationEstimate(arma::sp_mat, arma::mat, arma::rowvec, arma::mat);
    arma::vec left_multiply(arma::vec);
    arma::vec right_multiply(arma::vec);
    double p_omega_z();
};

CitationEstimate::CitationEstimate(arma::sp_mat citations, arma::mat U1, arma::rowvec d1, arma::mat V1) {
  U = U1;
  d = d1;
  V = V1;
  A = citations;
}

arma::vec CitationEstimate::left_multiply(arma::vec x) {

  //
  // want to compute A_tilde x where
  //
  // A_tilde = P_{omega tilde} (A) - P_L (Z) - P_U (Z) + Z, where Z = U d V^t
  //
  // see Section 3.3.1 Feasible Implementation of Hayes and Rohe (2024)
  //

  arma::vec p1 = A * x; //  P_{omega tilde} (A) x

  // see p_u_zx_impl(args$u, args$d, args$v, x, getOption("Ncpus", 1L))
  // for initial implementation
  arma::vec p2 =  arma::zeros<arma::vec>(U.n_rows); //  P_L (Z)

  // just DVt at this point
  arma::mat W = diagmat(d) * V.t();

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

  for (int j = 0; j < U.n_rows; j++) {
    p2(j) = dot(U.row(j), W.col(j));
  }

  // see p_u_tilde_zx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  // for initial implementation
  arma::vec p3 = arma::zeros<arma::vec>(U.n_rows); // P_U (Z)

  arma::sp_mat::const_iterator it     = A.begin();
  arma::sp_mat::const_iterator it_end = A.end();

  int i, j;
  double z_ij;

  for(; it != it_end; ++it) {

    i = it.row(); // row index of current item
    j = it.col(); // col index of current item

    if (i >= j) {
      z_ij = accu(U.row(i) % d % V.row(j));
      p3(i) += x(j) * z_ij;
    }
  }

  arma::vec p4 =  U * diagmat(d) * V.t() * x; //  Z x
  return p1 - p2 - p3 + p4;
}

arma::vec CitationEstimate::right_multiply(arma::vec x) {

  //
  // want to compute x A_tilde where
  //
  // A_tilde = P_{omega tilde} (A) - P_L (Z) - P_U (Z) + Z, where Z = U d V^t
  //
  // we compute Atilde^T x instead
  //
  // see Section 3.3.1 Feasible Implementation of Hayes and Rohe (2024)
  //

  arma::vec p1 = A.t() * x; //  P_{omega tilde} (A)^T x

  // see p_u_ztx_impl(args$u, args$d, args$v, x, getOption("Ncpus", 1L))
  // for initial implementation
  arma::vec p2 =  arma::zeros<arma::vec>(V.n_rows); //  P_L (Z)^T

  // just UD at this point
  arma::mat S = U * diagmat(d);

  // multiply rows by x to obtain W
  for (int i = 0; i < S.n_rows; i++) {
    S.row(i) *= x(i);
  }

  // cumulative summation step

  S = cumsum(S); // column-wise by default
  S.insert_rows(0, 1); // add 1 row of zeros at zeroth row index
  S.shed_row(S.n_rows - 1);

  for (int i = 0; i < V.n_rows; i++) {
    p2(i) = dot(S.row(i), V.row(i));
  }

  // see p_u_tilde_ztx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
  // for initial implementation
  arma::vec p3 = arma::zeros<arma::vec>(V.n_rows); // P_U (Z)

  arma::sp_mat::const_iterator it     = A.begin();
  arma::sp_mat::const_iterator it_end = A.end();

  int i, j;
  double z_ij;

  for(; it != it_end; ++it) {

    i = it.row(); // row index of current item
    j = it.col(); // col index of current item

    // only elements of lower triangle and the diagonal
    if (i >= j) {
      z_ij = accu(U.row(i) % d % V.row(j));
      p3(j) += x(i) * z_ij;
    }
  }

  arma::vec p4 =  V * diagmat(d) * U.t() * x; //  Z^T x
  return p1 - p2 - p3 + p4;
}


double CitationEstimate::p_omega_z() {

  // computes two of the terms in proposition 3.2
  //
  // p_L(Z) and P_U (Z) terms

  arma::mat DVt = diagmat(d) * V.t();
  double f_norm_sq = 0;

  arma::sp_mat::const_iterator it     = A.begin();
  arma::sp_mat::const_iterator it_end = A.end();

  int i, j;
  double z_ij;

  for(; it != it_end; ++it) {

    i = it.row(); // row index of current item
    j = it.col(); // col index of current item

    // only elements of lower triangle and the diagonal
    if (i >= j) {
      z_ij = dot(U.row(i), DVt.col(j));
      f_norm_sq += pow(z_ij, 2);
    }
  }

  // second add all the elements of the upper triangle

  int r = U.n_cols;

  // #pragma omp parallel for reduction(+: f_norm_sq)
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

// [[Rcpp::export]]
Rcpp::XPtr<CitationEstimate> makeCitationEstimate(arma::sp_mat citations, arma::mat U, arma::rowvec d, arma::mat V) {
  CitationEstimate* p = new CitationEstimate(citations, U, d, V);
  Rcpp::XPtr<CitationEstimate> ptr(p);
  return ptr;
}

// [[Rcpp::export]]
arma::vec left(arma::vec x, Rcpp::XPtr<CitationEstimate> A) {
  return A.get() -> left_multiply(x);
}

// [[Rcpp::export]]
arma::vec right(arma::vec x, Rcpp::XPtr<CitationEstimate> A) {
  return A.get() -> right_multiply(x);
}

// [[Rcpp::export]]
double p_omega_z(Rcpp::XPtr<CitationEstimate> A) {
  return A.get() -> p_omega_z();
}