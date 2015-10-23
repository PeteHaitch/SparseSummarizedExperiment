### =========================================================================
### SparseSummarizedExperiment objects
### -------------------------------------------------------------------------
###
### NOTE: The class hierarchy is as follows:
###       SummarizedExperiment0
###       ├── RangedSummarizedExperiment
###       │   ├── RangedSparseSummarizedExperiment
###       ├── SparseSummarizedExperiment

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseSummarizedExperiment class
###

#' SparseSummarizedExperiment objects
#'
#' @rdname SparseSummarizedExperiment
#'
#' @include SparseAssays-class.R
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment0
#' @importFrom methods setClass
#'
#' @export
setClass("SparseSummarizedExperiment",
         contains = "SummarizedExperiment0",
         representation = list(
           sparseAssays = "SparseAssays"
         ),
         prototype = list(
           sparseAssays = SparseAssays()
         )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

#' @importFrom S4Vectors setValidity2
setValidity2("SparseSummarizedExperiment",
             .valid.SSE)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# See R/RangedSparseSummarizedExperiment.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO: SparseSummarizedExperiment -> SummarizedExperiment0
# TODO: SparseSummarizedExperiment -> RangedSparseSummarizedExperiment

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssays", "SparseSummarizedExperiment",
          function(x, ..., withDimnames = TRUE, expand = FALSE) {
            .sparseAssays.SSE(x, ..., withDimnames = withDimnames, expand = expand)
          }
)

# TODO
#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays", "SparseSummarizedExperiment",
                 function(x, ..., value) {
                   stop("Not yet implemented")
                   .sparseAssaysReplace.SSE(x, ..., value)
                 }
)

## convenience for common use case

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "missing"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {
            .sparseAssay.SSE.missing(x = x, withDimnames = withDimnames,
                                 expand = expand)
          }
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "numeric"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {
            .sparseAssay.SSE.numeric(x, i, ..., withDimnames = withDimnames,
                                 expand = expand)
          }
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "character"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {
            .sparseAssay.SSE.character(x, i, ..., withDimnames = withDimnames,
                                   expand = expand)
          }
)

# TODO
# See assay<-,SummarizedExperiment-method, which has multiple methods; do
# I need something like this?
# What are valid signatures for this method?
# Can I call out to the `[[`,SparseAssays-method?
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay", "SparseSummarizedExperiment",
                 function(x, ..., value) {
                   stop("Not yet implemented")
                   .sparseAssayReplace.SSE(x, ..., value)
                 }
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssayNames", "SparseSummarizedExperiment",
          function(x, ...) {
            .sparseAssayNames.SSE(x)
          }
)

# TODO
#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssayNames", "SparseSummarizedExperiment",
                 function(x, ..., value) {
                   stop("Not yet implemented")
                   .sparseAssayNamesReplace.SSE(x)
                 }
)

# NOTE: The cannonical location for dim, dimnames. dimnames should be checked
#       for consistency (if non-null) and stripped from sparseAssays on
#       construction, or added from assays if dimnames are NULL in
#       <SparseSummarizedExperiment> but not sparseAssays. dimnames need to be
#       added on to sparse assays when sparseAssays() or sparseAssay() are
#       invoked.
# NOTE: dimnames and dimnames<- methods are inherited from
#       RangedSummarizedExperiment.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

#' @rdname SparseSummarizedExperiment
#'
#' @export
setMethod("[", "SparseSummarizedExperiment",
          function(x, i, j, ..., drop = TRUE) {
            .subsetSingleBracket.SSE(x, i, j, ..., drop = drop)
          }
)

#' @rdname SparseSummarizedExperiment
#'
#' @export
setReplaceMethod("[",
                 c("SparseSummarizedExperiment", "ANY", "ANY",
                   "SparseSummarizedExperiment"),
                 function(x, i, j, ..., value) {
                   .replaceSingleBracket.SSE(x, i, j, ..., value = value)
                 }
)

# NOTE: extractROWS() and replaceROWS() methods inherited from
#       SummarizedExperiment0 objects.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quick colData access.
###

# NOTE: There methods are inherited from SummarizedExperiment0 objects.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display.
###

#' @rdname SparseSummarizedExperiment
#'
#' @export
setMethod("show", "SparseSummarizedExperiment",
          function(object) {
            .show.SSE(object)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
###

# NOTE: Appropriate for objects with distinct features and identical samples.
#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("rbind", "SparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .rbind.SSE(args)
          }
)

# NOTE: Appropriate for objects with identical features and distinct samples.
#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("cbind", "SparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .cbind.SSE(args)
          }
)

#' @rdname SparseSummarizedExperiment
#'
#' @export
setMethod("combine",
          c("SparseSummarizedExperiment", "SparseSummarizedExperiment"),
          function(x, y, ...) {
            .combine.SSE(x, y, ...)
          }
)
