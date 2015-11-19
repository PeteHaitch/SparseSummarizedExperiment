### -------------------------------------------------------------------------
### densify
###

# TODO: roxygen splits i,j by a newline in the .Rd file, which is ugly.
#' Densify sparse assays.
#'
#' A generic function to densify (expand) the elements of a SparseAssays object.
#'
#' @param x A \link{SparseAssays} object.
#' @param i,j Numeric or character vectors indicating which sparse assay
#'        (\code{i}) and samples (\code{j}) to extract and densify. At least
#'        one of \code{i} or \code{j} must be specified; see \sQuote{Details}.
#' @param ... Additional arguments, for use in specific methods.
#'
#' @details \strong{WARNING}: Since it is generally undesirable to
#'          simultaneously densify all sparse assays and samples, you must
#'          specify at least one of \code{i} and \code{j}. If you \emph{really}
#'          wish to simultaneously densify all sparse assays and samples, then
#'          use \code{densify(x, seq_len(length(x)), seq_len(ncol(x)))}. If
#'          \code{i} (resp. \code{j}) is missing then effectively
#'          \code{i = seq_len(length(x))} (resp. \code{j = seq_len(ncol(x))}).
#'
#' @return A \code{\link[S4Vectors]{SimpleList}} of length =
#'         \code{length(i)}, each containing a
#'         \code{\link[S4Vectors]{SimpleList}} of length = \code{length{j}},
#'         each containing a \code{matrix} of the densified data for that
#'         sample in that sparse assay.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @seealso
#' \itemize{
#'  \item \link{SimpleListSparseAssays} for examples.
#' }
#'
#' @aliases densify,SimpleListSparseAssays,character,character-method
#'          densify,SimpleListSparseAssays,character,missing-method
#'          densify,SimpleListSparseAssays,character,numeric-method
#'          densify,SimpleListSparseAssays,missing,character-method
#'          densify,SimpleListSparseAssays,missing,numeric-method
#'          densify,SimpleListSparseAssays,numeric,character-method
#'          densify,SimpleListSparseAssays,numeric,missing-method
#'          densify,SimpleListSparseAssays,numeric,numeric-method
#'
#' @examples
#' # See ?SimpleListSparseAssays
#'
#' @include SparseSummarizedExperiment.R S4Vectors-pkg.R GenomicRanges-pkg.R
#' SummarizedExperiment-pkg.R
#'
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("densify",
           function(x, i, j, ...)
             standardGeneric("densify"),
           signature = c("x", "i", "j")
)

### -------------------------------------------------------------------------
### sparsify
###

#' Sparsify a matrix-like object.
#'
#' A generic function to sparsify a matrix-like object, e.g.,
#' \code{\link[base]{matrix}}, \code{\link[base]{data.frame}}, or
#' \code{\link[data.table]{data.table}}.
#'
#' @param x A matrix-like object, e.g., \code{\link[base]{matrix}},
#'          \code{\link[base]{data.frame}}, or
#'          \code{\link[data.table]{data.table}}.
#' @param return_class The class of the returned object.
#' @param ... Additional arguments, for use in specific methods.
#'
#' @details While the concepts of "sparsify" and "densify" are inverse
#' operations, in general, \code{identical(densify(sparsify(x)), x)} is
#' \strong{\code{FALSE}}. This is because the return value of \code{sparsify()}
#' is generally not compatible with the signature of a
#' \code{\link{densify}()} method.
#'
#' @return An object with the concepts of a \code{key} and a \code{value}. The
#' specific return value is left to individual methods. For example, the
#' \code{sparsify,SimpleList,matrix-method} returns a
#' \link[S4Vectors]{SimpleList} object with an \code{integer} \code{key}
#' element and a \code{matrix} \code{value} element.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @seealso
#' \itemize{
#'  \item \link{SimpleListSparseAssays}
#' }
#'
#' @aliases sparsify,data.frame,character-method
#'          sparsify,data.table,character-method
#'          sparsify,matrix,character-method
#'
#' @examples
#' # See ?SimpleListSparseAssays
#'
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparsify",
           function(x, return_class, ...)
             standardGeneric("sparsify"),
           signature = c("x", "return_class")
)

### -------------------------------------------------------------------------
### sparseAssays
###
### NOTE: Documented in R/SparseSummarizedExperiment-class.R

# NOTE: Following assays(), sparseAssays() will not strip the dimnames if
#       withDimnames = FALSE but will simply fail to add them.
# NOTE: If the user wants sparseAssays as a ShallowSimpleListAssays object then
#       they should run as(sparseAssays(x), "ShallowSimpleListAssays"). The
#       returned object will not have rownames regardless of the value of
#       withDimnames.
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssays",
           function(x, ..., withDimnames = TRUE)
             standardGeneric("sparseAssays"),
           signature = "x"
)

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssays<-",
           function(x, ..., withDimnames = TRUE, value)
             standardGeneric("sparseAssays<-"),
           signature = c("x", "value")
)

### -------------------------------------------------------------------------
### sparseAssay
###
### NOTE: Documented in R/SparseSummarizedExperiment-class.R

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssay",
           function(x, i, ...)
             standardGeneric("sparseAssay"),
           signature = c("x", "i")
)

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssay<-",
           function(x, i, ..., value)
             standardGeneric("sparseAssay<-"),
           signature = c("x", "i", "value")

)

### -------------------------------------------------------------------------
### sparseAssayNames
###
### NOTE: Documented in R/SparseSummarizedExperiment-class.R

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssayNames",
           function(x, ...)
             standardGeneric("sparseAssayNames")
)

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssayNames<-",
           function(x, ..., value)
             standardGeneric("sparseAssayNames<-")
)

### -------------------------------------------------------------------------
### SparseSummarizedExperiment
###
### NOTE: Documented in R/SparseSummarizedExperiment-class.R

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("SparseSummarizedExperiment",
           function(sparseAssays, ...)
             standardGeneric("SparseSummarizedExperiment")
)
