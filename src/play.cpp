// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
using namespace Eigen;

// [[Rcpp::export]]
VectorXd lt_multiply(const MatrixXd & U, const MatrixXd & V,const VectorXd & x)
{

    // U is n x k
    // V is d by k
    // x is d by 1 as a vector

    // must have:
    // U.cols() == V.cols();
    // x.size() == V.rows();

    int n = U.rows();
    int d = V.rows();

    VectorXd y = VectorXd::Zero(n);

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

