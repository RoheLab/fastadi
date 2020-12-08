#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

struct diff_norm_sq_worker : public Worker
{
  const arma::mat &U;
  const arma::mat &DVt;
  const arma::mat &new_U;
  const arma::mat &new_DVt;
  double diff_norm_sq;

  diff_norm_sq_worker(const arma::mat &U,
                      const arma::mat &DVt,
                      const arma::mat &new_U,
                      const arma::mat &new_DVt) : U(U), DVt(DVt), new_U(new_U), new_DVt(new_DVt), diff_norm_sq(0) {}

  diff_norm_sq_worker(const diff_norm_sq_worker &theworker,
                      Split) : U(theworker.U), DVt(theworker.DVt), new_U(theworker.new_U), new_DVt(theworker.new_DVt), diff_norm_sq(0) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (int i = begin; i < end; i++)
    {

      arma::rowvec new_Z_i, Z_i, diff;

      new_Z_i = new_U.row(i) * new_DVt;
      Z_i = U.row(i) * DVt;

      diff = new_Z_i - Z_i;
      diff_norm_sq += dot(diff, diff);
    }
  }

  void join(const diff_norm_sq_worker &theworker)
  {
    diff_norm_sq += theworker.diff_norm_sq;
  }

};



// [[Rcpp::export]]
double relative_f_norm_change_impl_cpp(
    const arma::mat &new_U,
    const arma::rowvec &new_d,
    const arma::mat &new_V,
    const arma::mat &U,
    const arma::rowvec &d,
    const arma::mat &V)
{


  arma::mat new_DVt = diagmat(new_d) * new_V.t();
  arma::mat DVt = diagmat(d) * V.t();


// expand a single row of new_Z and Z at a time
// to avoid hitting memory limits

  diff_norm_sq_worker diff_norm(U,DVt,new_U,new_DVt);
  parallelReduce(0,U.n_rows, diff_norm);

  // exploit connection between frobenius norm
  // and singular values
  double z_norm_sq = accu(pow(d, 2));

  return (diff_norm.diff_norm_sq / z_norm_sq);
}


