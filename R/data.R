#' MovieLens 100K dataset
#'
#' Standard benchmarking dataset for recommendation systems.
#' 100k movie ratings on 1682 movies by 943 users. Each user
#' has rated at least 20 movies.
#'
#' Stored as a `Matrix::dgCMatrix` object, which is a sparse
#' matrix. Each row corresponds to a user and each column to
#' a movie.
#'
#' @references
#'
#'   F. Maxwell Harper and Joseph A. Konstan. 2015. The MovieLens Datasets:
#'   History and Context. ACM Transactions on Interactive Intelligent
#'   Systems (TiiS) 5, 4, Article 19 (December 2015), 19 pages.
#'   DOI=http://dx.doi.org/10.1145/2827872
#'
#'   \url{https://grouplens.org/datasets/movielens/100k/}
#'
"ml100k"
