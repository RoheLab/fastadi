
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastadi

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
status](https://www.r-pkg.org/badges/version/fastadi)](https://CRAN.R-project.org/package=fastadi)
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

At the moment `fastadi` is a relatively low-level package, performing
matrix factorization and returning a low rank approximation `U D V'` of
an original data matrix. It is up to you to manipulate the matrices `U`,
`D` and `V` appropriately for your application; if you are interested in
a higher level interface for a task you are working on, please file an
issue\!

## Installation

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
#> Loading required package: Matrix
#> Loading required package: LRMF3

mf <- adaptive_impute(ml100k, rank = 3L, max_iter = 5L)
#> INFO [2020-09-22 17:07:15] Use svd initialization.
#> INFO [2020-09-22 17:07:15] Done initializing.
#> INFO [2020-09-22 17:07:15] Beginning AdaptiveImpute (max 5 iterations).
#> INFO [2020-09-22 17:07:15] Checking convergence every 1 iteration(s).
#> INFO [2020-09-22 17:07:15] Iter 1 complete. delta = 0.22312337, alpha = 184.71
#> INFO [2020-09-22 17:07:16] Iter 2 complete. delta = 0.05159124, alpha = 154.411
#> INFO [2020-09-22 17:07:16] Iter 3 complete. delta = 0.02125919, alpha = 135.61
#> INFO [2020-09-22 17:07:17] Iter 4 complete. delta = 0.01114668, alpha = 122.308
#> INFO [2020-09-22 17:07:17] Iter 5 complete. delta = 0.00669206, alpha = 112.354
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
#> d[rank]: 467.762
#> alpha:   112.354
#> 
#> Components
#> 
#> u: 943 x 3 [matrix] 
#> d: 3      [numeric] 
#> v: 1682 x 3 [matrix]
```

## References

1.  Cho, J., Kim, D. & Rohe, K. Asymptotic Theory for Estimating the
    Singular Vectors and Values of a Partially-observed Low Rank Matrix
    with Noise. arXiv:1508.05431 \[stat\] (2015).
    <http://arxiv.org/abs/1508.05431>

2.  Cho, J., Kim, D. & Rohe, K. Intelligent Initialization and Adaptive
    Thresholding for Iterative Matrix Completion; Some Statistical and
    Algorithmic Theory for Adaptive-Impute. Journal of Computational and
    Graphical Statistics 1–26 (2018)
    <doi:10.1080/10618600.2018.1518238>.
    <https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2018.1518238>

You can find the original implementation accompanying these papers
[here](https://github.com/chojuhee/hello-world).
