#ifndef _FASTADI_TYPES
#define _FASTADI_TYPES

#include "SparseMatrix.h"

// see section 2.5 of https://cran.r-project.org/web//packages/Rcpp/vignettes/Rcpp-attributes.pdf
// and also https://stackoverflow.com/questions/69831212/non-intrusive-extension-of-rcpp-with-rcpptraitsexporter

// we have to forward declare any classes wrapped in Rcpp::XPtr that should
// up as parts of type signatures. without this, sourceCpp() will work
// but the package won't compile correctly
class CitationEstimate;

#endif
