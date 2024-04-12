// usage example: https://github.com/zdebruine/RcppML/blob/5449a5b479908f40f56cf911f11e0a7e156d207f/src/RcppFunctions.cpp#L9

#include "SparseMatrix.h"

using namespace Rcpp;

class CitationEstimate {
  public:
    NumericMatrix U, D, V;
    SparseMatrix A;
    CitationEstimate(SparseMatrix, NumericMatrix, NumericVector, NumericMatrix);
    NumericVector left_multiply(NumericVector);
};

CitationEstimate::CitationEstimate(SparseMatrix A1, NumericMatrix U1, NumericVector d1, NumericMatrix V1) {
  U = U1;

  NumericMatrix D1(d1.size(), d1.size());

  for (int i = 0; i < d1.size(); i++) {
    D1(i, i) = d1(i);
  }

  D = D1;
  V = V1;
  A = A1;
}

NumericVector CitationEstimate::left_multiply(NumericVector x) {

  //
  // want to compute A_tilde x where
  //
  // A_tilde = P_{omega tilde} (A) - P_L (Z) - P_U (Z) + Z, where Z = U d V^t
  //
  // see Section 3.3.1 Feasible Implementation of Hayes and Rohe (2024)
  //

  NumericVector p1(x.size());  // A * x; //  P_{omega tilde} (A) x
  NumericVector p2(x.size());
  NumericVector p3(x.size());
  NumericVector p4(x.size());

  int i;

  // A * x is not implemented for SparseMatrix classes, so we implement
  // it here ourselves

  for (size_t col = 0; col < A.cols(); ++col) {
    for (SparseMatrix::InnerIterator it(A, col); it; ++it) {
      i = it.row();
      p1(i) += it.value() * x(col);
    }
  }

//   // see p_u_zx_impl(args$u, args$d, args$v, x, getOption("Ncpus", 1L))
//   // for initial implementation
//   arma::vec p2 =  arma::zeros<arma::vec>(U.n_rows); //  P_L (Z)
//
//   // just DVt at this point
//   arma::mat W = diagmat(d) * V.t();
//
//   for (int j = 0; j < W.n_cols; j++) {
//     W.col(j) *= x(j);
//   }
//
//   // cumulative summation step
//
//   // the farthest right column should should always be all zero
//   // add a column of zeros to the right side of matrix
//   W.insert_cols(W.n_cols, 1);
//
//   // perform the cumulative summation from right to left
//   // skip the rightmost two columns, which are already fine
//   // and the leftmost column, which we will drop
//   for (int j = W.n_cols - 3; j > 0; j--) {
//     W.col(j) += W.col(j + 1);
//   }
//
//   // drop the leftmost column. now we have W_tilde in full
//   W.shed_col(0);
//
//   for (int j = 0; j < U.n_rows; j++) {
//     p2(j) = dot(U.row(j), W.col(j));
//   }
//
//   // see p_u_tilde_zx_impl(args$u, args$d, args$v, mask@i, mask@j, x)
//   // for initial implementation
//   arma::vec p3 = arma::zeros<arma::vec>(U.n_rows); // P_U (Z)
//
//   arma::sp_mat::const_iterator it     = A.begin();
//   arma::sp_mat::const_iterator it_end = A.end();
//
//   int i, j;
//   double z_ij;
//
//   for(; it != it_end; ++it) {
//
//     i = it.row(); // row index of current item
//     j = it.col(); // col index of current item
//
//     if (i >= j) {
//       z_ij = accu(U.row(i) % d % V.row(j));
//       p3(i) += x(j) * z_ij;
//     }
//   }
// //
//

  Rcout << "D is" << std::endl << D << std::endl;
  Rcout << "p4.size() is" <<  p4.size() << std::endl;
  Rcout << "p4 is" << std::endl << p4 << std::endl;
  Rcout << "x.size() is" <<  x.size() << std::endl;
  Rcout << "x is" << std::endl << x << std::endl;
//
  p4 = x;

  Rcout << "p4 is now" << std::endl << p4 << std::endl;
  Rcout << "p4.size() is now" << p4.size() << std::endl;

  // p4 = U * D * Rcpp::transpose(V) * x;
//   return p1 - p2 - p3 + p4;
  return p1;
}

// [[Rcpp::export]]
Rcpp::XPtr<CitationEstimate> makeCitationEstimate(SparseMatrix A, NumericMatrix U, NumericVector d, NumericMatrix V) {
  CitationEstimate* p = new CitationEstimate(A, U, d, V);
  Rcpp::XPtr<CitationEstimate> ptr(p);
  return ptr;
}

// [[Rcpp::export]]
NumericVector left(NumericVector x, Rcpp::XPtr<CitationEstimate> A) {
  return A.get() -> left_multiply(x);
}

/***R
library(RSpectra)
library(Matrix)

set.seed(111)

n <- 100

A <- rsparsematrix(n, n, 0.1, rand.x = NULL) * 1
k <- 5

s <- svds(A, k)

m <- makeCitationEstimate(A, s$u, s$d, s$v)
m

x <- rnorm(n)

# s$u %*% diag(s$d) %*% t(s$v) %*% x

l1 <- left(x, m)
l1

*/
