#include <RcppArmadillo.h>


class LRMatrix {
  public:
    arma::mat U, V;
    arma::rowvec d;
    LRMatrix(arma::mat, arma::rowvec, arma::mat);
    arma::vec left_multiply(arma::vec);
};

// constructor
LRMatrix::LRMatrix(arma::mat U, arma::rowvec d, arma::mat V) {
  U = U;
  d = d;
  V = V;
}

// method
arma::vec LRMatrix::left_multiply(arma::vec x) {
  return x;
}

class SparseLRMatrix {
  public:
    SparseLRMatrix( arma::sp_mat sparse,  ) : value(v) {}
    double square() {return value*value;}
  private:
    double value;
};


class Rectangle {
  int width, height;
public:
  void set_values (int,int);
  int area() {return width*height;}
};

void Rectangle::set_values (int x, int y) {
  width = x;
  height = y;
}

int main () {
  Rectangle rect;
  rect.set_values (3,4);
  cout << "area: " << rect.area();
  return 0;
}



area: 12

// [[Rcpp::export]]
Rcpp::XPtr<Double> makeDouble(double x) {
  Double* pd = new Double(x);
  Rcpp::XPtr<Double> ptr(pd);
  return ptr;
}

// [[Rcpp::export]]
double squareDouble(Rcpp::XPtr<Double> x) {
  return x.get()->square();
}

/***R
(d2 <- makeDouble(5.4))
squareDouble(d2)

m = 100
n = 20
k = 5
set.seed(111)
A = matrix(rnorm(m * n), m)

svds(A, k)
svds(t(A), k, nu = 0, nv = 3)

## Sparse matrices
library(Matrix)
A[sample(m * n, m * n / 2)] = 0
Asp1 = as(A, "dgCMatrix")
Asp2 = as(A, "dgRMatrix")

svds(Asp1, k)
svds(Asp2, k, nu = 0, nv = 0)

## Function interface
Af = function(x, args)
{
  as.numeric(args %*% x)
}

Atf = function(x, args)
{
  as.numeric(crossprod(args, x))
}

svds(Af, k, Atrans = Atf, dim = c(m, n), args = Asp1)

# s_new <- svds(
#   Ax_citation,
#   k = rank,
#   Atrans = Atx_citation,
#   dim = dim(X),
#   args = args
# )

*/
