#' @keywords internal
#'
#' @importFrom glue glue
#' @import LRMF3
#' @import Matrix
#' @importFrom logger log_info log_debug
#' @importFrom methods as
#' @importFrom RSpectra svds
#' @importFrom stats sd
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp evalCpp
#' @useDynLib fastadi, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("fastadi", libpath)
}
