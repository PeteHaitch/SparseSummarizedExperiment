### -------------------------------------------------------------------------
### densify
###
### NOTE: Documented in R/SparseAssays-class.R

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

### -------------------------------------------------------------------------
### saapply
###
### NOTE: Documented in R/SparseAsaays-class.R

# TODO: Will need to keep API up-to-date with BiocParallel::bplapply
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("SAapply",
           function(X, FUN, densify = TRUE, ..., BPREDO = list(),
                    BPPARAM = bpparam())
             standardGeneric("SAapply"),
           signature = "X"
)
