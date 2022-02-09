
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastadi

<!-- badges: start -->

[![R-CMD-check](https://github.com/RoheLab/fastadi/workflows/R-CMD-check/badge.svg)](https://github.com/RoheLab/fastadi/actions)
[![Codecov test
coverage](https://codecov.io/gh/RoheLab/fastadi/branch/main/graph/badge.svg)](https://app.codecov.io/gh/RoheLab/fastadi?branch=main)
<!-- badges: end -->

`fastadi` implements the `AdaptiveImpute` matrix completion algorithm.
`fastadi` is a self-tuning alternative to algorithms such as
`SoftImpute` (implemented in the
[`softImpute`](https://cran.r-project.org/package=softImpute) package),
truncated SVD, maximum margin matrix factorization, and weighted
regularized matrix factorization (implemented in the
[`rsparse`](https://github.com/rexyai/rsparse) package). In simulations
`fastadi` often outperforms `softImpute` by a small margin.

You may find `fastadi` useful if you are developing embeddings for
sparsely observed data, if you are working in natural language
processing, or building a recommendation system.

## Installation

You can install the released version from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("fastadi")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("RoheLab/fastadi")
```

## Example usage

Here we embed users and items in the MovieLens 100K dataset.

``` r
library(fastadi)
#> Loading required package: LRMF3
#> Loading required package: Matrix
mf <- adaptive_impute(ml100k, rank = 3L, max_iter = 5L)
#> Warning: 
#> Reached maximum allowed iterations. Returning early.
```

``` r
mf
#> 
#> Adaptively Imputed Low Rank Matrix Factorization
#> ------------------------------------------------
#> 
#> Rank: 3
#> 
#> Rows: 943
#> Cols: 1682
#> 
#> d[rank]: 467.486
#> alpha:   144.663
#> 
#> Components
#> 
#> u: 943 x 3 [matrix] 
#> d: 3      [numeric] 
#> v: 1682 x 3 [matrix]
```

Note that the vignettes are currently scratch work for reference by the
developers and are not yet ready for general consumption.

## References

1.  Cho, Juhee, Donggyu Kim, and Karl Rohe. “Asymptotic Theory for
    Estimating the Singular Vectors and Values of a Partially-Observed
    Low Rank Matrix with Noise.” Statistica Sinica, 2018.
    <https://doi.org/10.5705/ss.202016.0205>.

2.  ———. “Intelligent Initialization and Adaptive Thresholding for
    Iterative Matrix Completion: Some Statistical and Algorithmic Theory
    for Adaptive-Impute.” Journal of Computational and Graphical
    Statistics 28, no. 2 (April 3, 2019): 323–33.
    <https://doi.org/10.1080/10618600.2018.1518238>.

You can find the original implementation accompanying these papers
[here](https://github.com/chojuhee/hello-world).
