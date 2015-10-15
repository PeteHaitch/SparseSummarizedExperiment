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

#' @export
setGeneric("SparseSummarizedExperiment",
           function(sparseAssays, ...) {
             standardGeneric("SparseSummarizedExperiment")
})

### -------------------------------------------------------------------------
### combine2
###

# TODO: Remove if/when a similar method is added to SummarizedExperiment
#       (see https://stat.ethz.ch/pipermail/bioc-devel/2015-October/008128.html)
#' @export
setGeneric("combine2", function(x, y, ...) {
  standardGeneric("combine2")
})
