#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

struct Wworker : public Worker
{
  const arma::vec &x;
  arma::mat &W;

  Wworker(const arma::vec &x, arma::mat &W) : x(x), W(W) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    for (int j = begin; j < end; j++)
    {
      W.col(j) *= x(j);
    }
  }
};


struct zxworker : public Worker
{
  const arma::mat & U;
  const arma::mat & W;
  arma::vec zx; 


  zxworker(const arma::mat &U, const arma::mat &W) : U(U), W(W), zx(zeros<vec>(U.n_rows)) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    for (int j = begin; j < end; j++)
    {
      zx(j) = dot(U.row(j), W.col(j));
    }
  }
};




// [[Rcpp::export]]
arma::vec p_u_zx_impl_cpp(
    const arma::mat &U,
    const arma::vec &d,
    const arma::mat &V,
    const arma::vec &x)
{


  // just DVt at this point
  arma::mat W = diagmat(d) * V.t();

// multiply columns by x to obtain W
  Wworker Ww(x,W);
  parallelFor(0,W.n_cols, Ww);

  // cumulative summation step

  // the farthest right column should should always be all zero
  // add a column of zeros to the right side of matrix
  W.insert_cols(W.n_cols, 1);

  // perform the cumulative summation from right to left
  // skip the rightmost two columns, which are already fine
  // and the leftmost column, which we will drop
  for (int j = W.n_cols - 3; j > 0; j--)
  {
    W.col(j) += W.col(j + 1);
  }

  // drop the leftmost column. now we have W_tilde in full
  W.shed_col(0);

  // do the matrix-vector multiplication
  

  zxworker zxw(U,W);
  parallelFor(0,U.n_rows,zxw);

  return(zxw.zx);
}

