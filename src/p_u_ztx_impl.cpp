#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

struct getSworker : public Worker
{
  const arma::vec & x;
  arma::mat & S;

  getSworker(const arma::vec & x,
             arma::mat & S) : x(x), S(S) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (int i = begin; i < end; i++)
    {
      S.row(i) *= x(i);
    }
  }
};

struct getztxworker: public Worker{
  const arma::mat & S;
  const arma::mat & V;

  arma::vec ztx;
  getztxworker(const arma::mat & S,
  const arma::mat & V) : S(S), V(V), ztx(zeros<vec>(V.n_rows)) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    for (int i = begin; i < end; i++)
    {
      ztx(i) = dot(S.row(i), V.row(i));
    }
  }
};



// [[Rcpp::export]]
arma::vec p_u_ztx_impl_cpp(
    const arma::mat &U,
    const arma::vec &d,
    const arma::mat &V,
    const arma::vec &x)
{


  // just UD at this point
  arma::mat S = U * diagmat(d);

// multiply rows by x to obtain W
  getSworker Swork(x,S);
  parallelFor(0,S.n_rows,Swork);

  // cumulative summation step

  S = cumsum(S);       // column-wise by default
  S.insert_rows(0, 1); // add 1 row of zeros at zeroth row index
  S.shed_row(S.n_rows - 1);

  // do the matrix-vector multiplication
  getztxworker ztxwork(S,V);
  parallelFor(0, V.n_rows, ztxwork);
  return(ztxwork.ztx);
}