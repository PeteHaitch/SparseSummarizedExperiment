### -------------------------------------------------------------------------
### sparseAssays
###

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
