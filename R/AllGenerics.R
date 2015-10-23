### -------------------------------------------------------------------------
### sparseAssays
###

# NOTE: Following assays(), sparseAssays() will not strip the dimnames if
#       withDimnames = FALSE but will simply fail to add them.
# NOTE: The expand = TRUE argument returns a SimpleList, not a SparseAssays
#       object.
# NOTE: If the user wants sparseAssays as a ShallowSimpleListAssays object then
#       they should run as(sparseAssays(x), "ShallowSimpleListAssays"). The
#       returned object will not have rownames regardless of the value of
#       withDimnames. Note also that expand must be FALSE; the coercion to a
#       ShallowSimpleListAssays object automatically expands the
#       sparseAssays.
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssays", function(x, ..., withDimnames = TRUE, expand = FALSE) {
  standardGeneric("sparseAssays")
}, signature = "x")

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssays<-", function(x, ..., withDimnames = TRUE, value) {
  standardGeneric("sparseAssays<-")
}, signature = c("x", "value"))

### -------------------------------------------------------------------------
### sparseAssay
###

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssay", function(x, i, ..., withDimnames = TRUE,
                                   expand = FALSE) {
  standardGeneric("sparseAssay")
})

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssay<-", function(x, i, ..., withDimnames = TRUE, value) {
  standardGeneric("sparseAssay<-")
})

### -------------------------------------------------------------------------
### sparseAssayNames
###

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssayNames", function(x, ...) {
  standardGeneric("sparseAssayNames")
})

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("sparseAssayNames<-", function(x, ..., value) {
  standardGeneric("sparseAssayNames<-")
})

### -------------------------------------------------------------------------
### SparseSummarizedExperiment
###

#' @importFrom methods setGeneric
#'
#' @export
setGeneric("SparseSummarizedExperiment",
           function(sparseAssays, ...) {
             standardGeneric("SparseSummarizedExperiment")
})
