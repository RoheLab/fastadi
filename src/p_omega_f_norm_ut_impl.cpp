#include <RcppArmadillo.h>

#include <RcppParallel.h>
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

// first part of fnorm_sq
struct f_norm_sq_Worker1 : public Worker
{
  // input
  const arma::mat &U;
  const arma::mat &DVt;
  const arma::vec &row;
  const arma::vec &col;
  // output
  double f_norm_sq;
  // constructors
  f_norm_sq_Worker1(const arma::mat &U,
                    const arma::mat &DVt,
                    const arma::vec &row,
                    const arma::vec &col) : U(U), DVt(DVt), row(row), col(col), f_norm_sq(0) {}

  f_norm_sq_Worker1(const f_norm_sq_Worker1 &theworker, Split) : U(theworker.U), DVt(theworker.DVt), row(theworker.row), col(theworker.col), f_norm_sq(0) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (int idx = begin; idx < end; idx++)
    {
      int i, j;
      double z_ij;

      i = row(idx);
      j = col(idx);
      if (i >= j)
      {
        z_ij = dot(U.row(i), DVt.col(j));
        f_norm_sq += pow(z_ij, 2);
      }
    }
  }

  void join(const f_norm_sq_Worker1 &theworker)
  {
    f_norm_sq += theworker.f_norm_sq;
  }
};

// this is for the second part of fnorm_sq
struct f_norm_sq_Worker2 : public Worker
{
  // input
  const arma::mat &U;
  const arma::mat &DVt;
  const int &r;
  // output
  double f_norm_sq;
  // constructors
  f_norm_sq_Worker2(const arma::mat &U,
                    const arma::mat &DVt,
                    const int &r) : U(U), DVt(DVt), r(r), f_norm_sq(0) {}

  f_norm_sq_Worker2(const f_norm_sq_Worker2 &theworker, Split) : U(theworker.U), DVt(theworker.DVt), r(theworker.r), f_norm_sq(0) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (int l = begin; l < end; l++)
    {
      for (int q = 0; q < r; q++)
      {

        arma::vec U_lq;
        arma::rowvec V_lq, V_lq_tri;

        U_lq = U.col(l) % U.col(q);
        V_lq = DVt.row(l) % DVt.row(q);
        V_lq_tri = sum(V_lq) - cumsum(V_lq);

        f_norm_sq += dot(U_lq, V_lq_tri);
      }
    }
  }

  void join(const f_norm_sq_Worker2 &theworker)
  {
    f_norm_sq += theworker.f_norm_sq;
  }
};

// [[Rcpp::export]]
double p_omega_f_norm_ut_impl_cpp(
    const arma::mat &U,
    const arma::rowvec &d,
    const arma::mat &V,
    const arma::vec &row,
    const arma::vec &col)
{


  // first add the observed elements on the lower triangle and diagonal

  arma::mat DVt = diagmat(d) * V.t();
  double f_norm_sq = 0;

  f_norm_sq_Worker1 worker1(U, DVt, row, col);
  parallelReduce(0, row.n_elem, worker1);
  f_norm_sq += worker1.f_norm_sq;

  // second add all the elements of the upper triangle

  int r = U.n_cols;

  f_norm_sq_Worker2 worker2(U, DVt, r);
  parallelReduce(0, r, worker2);
  f_norm_sq += worker2.f_norm_sq;

  return f_norm_sq;
}