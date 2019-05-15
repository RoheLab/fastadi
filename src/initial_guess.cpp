// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
using namespace Eigen;

// implements Algorithm 1 from Cho et al 2018
// [[Rcpp::export]]
VectorXd initial_guess(const MatrixXd & M, const int r)
{

  // where a particular observation was observed
  Eigen::MatrixXd y = M == 0;
  int d = M.cols();

  // estimated sparsity in M
  double p = y.mean();

  MatrixXd MtM = M.transpose() * M;
  MatrixXd MMt = M * M.tranpose();

  MatrixXd Sigma_p = MtM - (1 - p) * MtM.diag();
  MatrixXd Sigma_tp = MMt - (1 - p) * MMt.diag();

  // what is the best svd algorithm to use here?
  JacobiSVD<MatrixXd> svd_p(Sigma_p, ComputeThinU | ComputeThinV);
  JacobiSVD<MatrixXd> svd_pt(Sigma_p, ComputeThinU | ComputeThinV);

  VectorXd lambda = svd_p.singularValues();

  double alpha = lambda.segment(r + 1, d).sum() / (d - r);

  VectorXd tau = VectorXd::Zero(r);
  VectorXd lambda_hat = VectorXd::Zero(r);

  tau = lambda.segment(1, r) - sqrt(lambda.segment(1, r) - alpha) / p;
  tau = lambda.segment(1, r) - sqrt(lambda.segment(1, r) - alpha) / p;

  MatrixXd Cp = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
  MatrixXd diff = Cp - C;
  cout << "diff:\n" << diff.array().abs().sum() << "\n";
  return 0;


  - (1 - p) * )

  // U is n x k
  // V is d by k
  // x is d by 1 as a vector

  // must have:
  // U.cols() == V.cols();
  // x.size() == V.rows();

  int n = U.rows();
  int d = V.rows();

  VectorXd y = Eigen::VectorXd::Zero(n);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < std::min(d, i); j++)
    {
      y(i) += U.row(i).dot(V.row(j)) * x(j);
    }
  }

  return y;
}

/*** R
n <- 7
d <- 5
k <- 3

M <- matrix(rep(0, n * d), n, d)
L <- M
L[lower.tri(L)] <- 1
L

U <- matrix(1:(n * k), n, k)
V <- matrix(rev(1:(d * k)), d, k)

x <- 1:d

# R solution
(L * U %*% t(V)) %*% x
lt_multiply(U, V, x)
*/

