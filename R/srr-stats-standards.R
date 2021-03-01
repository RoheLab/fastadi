#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.3} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#'
#' @srrstats {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
#'
#' @srrstats {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#'
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#'
#' @srrstatsTODO {G1.1} *Statistical Software should document whether the algorithm(s) it implements are:* - *The first implementation of a novel algorithm*; or - *The first implementation within **R** of an algorithm which has previously been implemented in other languages or contexts*; or - *An improvement on other implementations of similar algorithms in **R***.
#'
#' @srrstatsTODO {G1.2} *All statistical terminology should be clarified and unambiguously defined.*
#'
#' @srrstatsTODO {G1.2a} *All internal (non-exported) functions should also be documented in standard [`roxygen2`](https://roxygen2.r-lib.org/) format, along with a final `@noRd` tag to suppress automatic generation of `.Rd` files.*
#'
#' @srrstatsTODO {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#'
#' @srrstatsTODO {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
#'
#' @srrstatsTODO {G2.3} *For univariate character input:*
#'
#' @srrstatsTODO {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
#'
#' @srrstatsTODO {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.*
#'
#' @srrstatsTODO {G2.13} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
#'
#' @srrstatsTODO {G3.0} *Statistical software should never compare floating point numbers for equality. All numeric equality comparisons should either ensure that they are made between integers, or use appropriate tolerances for approximate equality.*
#'
#' @srrstatsTODO {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#'
#' @srrstatsTODO {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#'
#' @srrstatsTODO {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#'
#' @srrstatsTODO {G5.3} *For functions which are expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values, the absence of any such values in return objects should be explicitly tested.*
#'
#' @srrstatsTODO {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as [RStata](https://github.com/lbraglia/RStata)).*
#'
#' @srrstatsTODO {G5.4a} *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#'
#' @srrstatsTODO {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
#'
#' @srrstatsTODO {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
#'
#' @srrstatsTODO {G5.5} *Correctness tests should be run with a fixed random seed*
#'
#' @srrstatsTODO {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#'
#' @srrstatsTODO {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#'
#' @srrstatsTODO {G5.6b} *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G4.10-4.12, below).*
#'
#' @srrstatsTODO {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
#'
#' @srrstatsTODO {G5.8} **Edge condition tests** *to test that these conditions produce expected behaviour such as clear warnings or errors when confronted with data with extreme properties including but not limited to:*
#'
#' @srrstatsTODO {G5.8a} *Zero-length data*
#'
#' @srrstatsTODO {G5.8b} *Data of unsupported types (e.g., character or complex numbers in for functions designed only for numeric data)*
#'
#' @srrstatsTODO {G5.8c} *Data with all-`NA` fields or columns or all identical fields or columns*
#'
#' @srrstatsTODO {G5.8d} *Data outside the scope of the algorithm (for example, data with more fields (columns) than observations (rows) for some regression algorithms)*
#'
#' @srrstatsTODO {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#'
#' @srrstatsTODO {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#'
#' @srrstatsTODO {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
#'
#' @srrstats {UL1.0} *Unsupervised Learning Software should explicitly document expected format (types or classes) for input data, including descriptions of types or classes which are not accepted; for example, specification that software accepts only numeric inputs in `vector` or `matrix` form, or that all inputs must be in `data.frame` form with both column and row names.*
#'
#' @srrstats {UL1.1} *Unsupervised Learning Software should provide distinct sub-routines to assert that all input data is of the expected form, and issue informative error messages when incompatible data are submitted.*
#'
#' @srrstatsTODO {UL1.2} *Unsupervised learning which uses row or column names to label output objects should assert that input data have non-default row or column names, and issue an informative message when these are not provided. (Such messages need not necessarily be provided by default, but should at least be optionally available.)*
#'
#' @srrstatsTODO {UL1.3} *Unsupervised Learning Software should transfer all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`, to corresponding aspects of return objects.*
#'
#' @srrstatsTODO {UL1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.*
#'
#' @srrstatsTODO {UL1.4} *Unsupervised Learning Software should explicitly document whether input data may include missing values.*
#'
#' @srrstatsTODO {UL1.5} *Functions in Unsupervised Learning Software which do not admit input data with missing values should provide informative error messages when data with missing values are submitted.*
#'
#' @srrstatsTODO {UL1.6} *Unsupervised Learning Software should document any assumptions made with regard to input data; for example assumptions about distributional forms or locations (such as that data are centred or on approximately equivalent distributional scales). Implications of violations of these assumptions should be both documented and tested, in particular:*
#'
#'
#' @srrstats {UL1.6b} *Examples or other documentation should not use `scale()` or equivalent transformations without explaining why scale is applied, and explicitly illustrating and contrasting the consequences of not applying such transformations.*
#'
#' @srrstatsTODO {UL2.3} *Unsupervised Learning Software should implement pre-processing routines to identify whether aspects of input data are perfectly collinear.*
#'
#' @srrstatsTODO {UL3.1} *Algorithms which apply sequential labels to input data (such as clustering or partitioning algorithms) should ensure that the sequence follows decreasing group sizes (so labels of "1", "a", or "A" describe the largest group, "2", "b", or "B" the second largest, and so on.)*
#'
#' @srrstatsTODO {UL3.2} *Dimensionality reduction or equivalent algorithms which label dimensions should ensure that that sequences of labels follows decreasing "importance" (for example, eigenvalues or variance contributions).*
#'
#' @srrstats {UL3.4} *Where applicable, Unsupervised Learning Software should implement routines to predict the properties (such as numerical ordinates, or cluster memberships) of additional new data without re-running the entire algorithm.*
#'
#' @srrstats {UL4.0} *Unsupervised Learning Software should return some form of "model" object, generally through using or modifying existing class structures for model objects, or creating a new class of model objects.*
#'
#' @srrstatsTODO {UL4.2} *The return object from Unsupervised Learning Software should include, or otherwise enable immediate extraction of, all parameters used to control the algorithm used.*
#'
#' @srrstats {UL4.2} *Model objects returned by Unsupervised Learning Software should implement or appropriately extend a default `print` method which provides an on-screen summary of model (input) parameters and methods used to generate results. The `print` method may also summarise statistical aspects of the output data or results.*
#'
#' @srrstatsTODO {UL4.2a} *The default `print` method should always ensure only a restricted number of rows of any result matrices or equivalent are printed to the screen.*
#'
#' @srrstats {UL7.0} *Inappropriate types of input data are rejected with expected error messages.*
#'
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {UL3.3} *Unsupervised Learning Software for which input data does not generally include labels (such as `array`-like data with no row names) should provide an additional parameter to enable cases to be labelled.*
#'
#' @srrstatsNA {UL4.1} *Unsupervised Learning Software may enable an ability to generate a model object without actually fitting values. This may be useful for controlling batch processing of computationally intensive fitting algorithms.*
#'
#' @srrstatsNA {UL1.6a} *Software which responds qualitatively differently to input data which has components on markedly different scales should explicitly document such differences, and implications of submitting such data.*
#'
#' @srrstatsNA {G1.4} *Software should include all code necessary to reproduce results which form the basis of performance claims made in associated publications.* This is not the appropriate place for that code. The paper has it's own standalone repo.
#'
#' @srrstatsNA {G1.5} *Software should include code necessary to compare performance claims with alternative implementations in other R packages.* I do not make any
#'
#' @srrstatsNA {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#'
#' @srrstatsNA {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#'
#' @srrstatsNA {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#'
#' @srrstatsNA {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#'
#' @srrstatsNA {G2.4d} *explicit conversion to factor via `as.factor()`*
#'
#' @srrstatsNA {G2.4e} *explicit conversion from factor via `as...()` functions*
#'
#' @srrstatsNA {G2.6} *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.*
#'
#' @srrstatsNA {G2.7} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
#'
#' @srrstatsNA {G2.8} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).*
#'
#' @srrstatsNA {G2.9} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
#'
#' @srrstatsNA {G2.10} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#'
#' @srrstatsNA {G2.11} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#'
#' @srrstatsNA {G2.11a} *error on missing data*
#'
#' @srrstatsNA {G2.11b} *ignore missing data with default warnings or messages issued*
#'
#' @srrstatsNA {G2.11c} *replace missing data with appropriately imputed values*
#'
#' @srrstatsNA {G2.12} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
#'
#' @srrstatsNA {G4.0} *Statistical Software which enables outputs to be written to local files should parse parameters specifying file names to ensure appropriate file suffices are automatically generated where not provided.*
#'
#' @srrstatsNA {G5.0} *Where applicable or practicable, tests should use standard data sets with known properties (for example, the [NIST Standard Reference Datasets](https://www.itl.nist.gov/div898/strd/), or data sets provided by other widely-used R packages).*
#'
#' @srrstatsNA {G5.1} *Data sets created within, and used to test, a package should be exported (or otherwise made generally available) so that users can confirm tests and run examples.*
#'
#' @srrstatsNA {G5.10} *Extended tests should included and run under a common framework with other tests but be switched on by flags such as as a `<MYPKG>_EXTENDED_TESTS=1` environment variable.*
#'
#' @srrstatsNA {G5.11} *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
#'
#' @srrstatsNA {G5.11a} *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
#'
#' @srrstatsNA {G5.12} *Any conditions necessary to run extended tests such as platform requirements, memory, expected runtime, and artefacts produced that may need manual inspection, should be described in developer documentation such as a `CONTRIBUTING.md` or `tests/README.md` file.*
#'
#' @srrstatsNA {UL2.0} *Routines likely to give unreliable or irreproducible results in response to violations of assumptions regarding input data (see UL1.6) should implement pre-processing steps to diagnose potential violations, and issue appropriately informative messages, and/or include parameters to enable suitable transformations to be applied (such as the `center` and `scale.` parameters of the [`stats::prcomp()`](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/prcomp.html) function).*
#'
#' @srrstatsNA {UL2.1} *Unsupervised Learning Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*
#'
#' @srrstatsNA {UL2.2} *For Unsupervised Learning Software which accepts missing values in input data, functions should implement explicit parameters controlling the processing of missing values, ideally distinguishing `NA` or `NaN` values from `Inf` values (for example, through use of `na.omit()` and related functions from the `stats` package).*
#'
#' @srrstatsNA {UL1.6a} *Software which responds qualitatively differently to input data which has components on markedly different scales should explicitly document such differences, and implications of submitting such data.*
#'
#' @srrstatsNA {UL3.3} *Unsupervised Learning Software for which input data does not generally include labels (such as `array`-like data with no row names) should provide an additional parameter to enable cases to be labelled.*
#'
#' @srrstatsNA {UL3.5} *Objects returned from Unsupervised Learning Software which labels, categorise, or partitions data into discrete groups should include, or provide immediate access to, quantitative information on intra-group variances or equivalent, as well as on inter-group relationships where applicable.*
#'
#' @srrstatsNA {UL4.1} *Unsupervised Learning Software may enable an ability to generate a model object without actually fitting values. This may be useful for controlling batch processing of computationally intensive fitting algorithms.*
#'
#' @srrstatsNA {UL4.3} *Unsupervised Learning Software should also implement `summary` methods for model objects which should summarise the primary statistics used in generating the model (such as numbers of observations, parameters of methods applied). The `summary` method may also provide summary statistics from the resultant model.* No. I don't like this.
#'
#' @srrstatsNA {UL6.0} *Objects returned by Unsupervised Learning Software should have default `plot` methods, either through explicit implementation, extension of methods for existing model objects, through ensuring default methods work appropriately, or through explicit reference to helper packages such as [`factoextra`](https://github.com/kassambara/factoextra) and associated functions.*
#'
#' @srrstatsNA {UL6.1} *Where the default `plot` method is **NOT** a generic `plot` method dispatched on the class of return objects (that is, through a `plot.<myclass>` function), that method dispatch should nevertheless exist in order to explicitly direct users to the appropriate function.*
#'
#' @srrstatsNA {UL6.2} *Where default plot methods include labelling components of return objects (such as cluster labels), routines should ensure that labels are automatically placed to ensure readability, and/or that appropriate diagnostic messages are issued where readability is likely to be compromised (for example, through attempting to place too many labels).*
#'
#' @srrstatsTODO {UL7.1} *Tests should demonstrate that violations of assumed input properties yield unreliable or invalid outputs, and should clarify how such unreliability or invalidity is manifest through the properties of returned objects.* LOL this will be really hard.
#'
#' @srrstatsNA {UL7.2} *Demonstrate that labels placed on output data follow decreasing group sizes (**UL3.1**)*
#'
#' @srrstatsNA {UL7.3} *Demonstrate that labels on input data are propagated to, or may be recovered from, output data (see **UL3.3**).*
#'
#' @srrstatsNA {UL7.4} *Demonstrate that submission of new data to a previously fitted model can generate results more efficiently than initial model fitting.*
#'
#' @srrstatsNA {UL7.5} *Batch processing routines should be explicitly tested, commonly via extended tests (see **G4.10**--**G4.12**).*
#'
#' @srrstatsNA {UL7.5a} *Tests of batch processing routines should demonstrate that equivalent results are obtained from direct (non-batch) processing.*
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @noRd
NULL
