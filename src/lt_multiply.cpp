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

