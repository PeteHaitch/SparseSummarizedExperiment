### -------------------------------------------------------------------------
### sparseAssays
###

#' @export
setGeneric("sparseAssays", function(x, ..., withDimnames = TRUE, expand = FALSE) {
  standardGeneric("sparseAssays")
}, signature = "x")

#' @export
setGeneric("sparseAssays<-", function(x, ..., withDimnames = TRUE, value) {
  standardGeneric("sparseAssays<-")
}, signature = c("x", "value"))

### -------------------------------------------------------------------------
### sparseAssay
###

#' @export
setGeneric("sparseAssay", function(x, i, ..., withDimnames = TRUE,
                                   expand = FALSE) {
  standardGeneric("sparseAssay")
})

#' @export
setGeneric("sparseAssay<-", function(x, i, ..., withDimnames = TRUE, value) {
  standardGeneric("sparseAssay<-")
})

### -------------------------------------------------------------------------
### sparseAssayNames
###

#' @export
setGeneric("sparseAssayNames", function(x, ...) {
  standardGeneric("sparseAssayNames")
})

#' @export
setGeneric("sparseAssayNames<-", function(x, ..., value) {
  standardGeneric("sparseAssayNames<-")
})

### -------------------------------------------------------------------------
### SparseSummarizedExperiment
###

setGeneric("SparseSummarizedExperiment",
           function(sparseAssays, ...) {
             standardGeneric("SparseSummarizedExperiment")
})
