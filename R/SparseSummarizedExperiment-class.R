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
#' @include SparseAssays-class.R SSE-helpers.R
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
          .sparseAssays.SSE
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays",
                 c("SparseSummarizedExperiment", "SparseAssays"),
                 .sparseAssaysReplace.SSE
)

## convenience for common use case

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "missing"),
          .sparseAssay.SSE.missing
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "numeric"),
          .sparseAssay.SSE.numeric
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "character"),
          .sparseAssay.SSE.character
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "missing", "SimpleList"),
                 .sparseAssayReplace.SSE.missing
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "numeric", "SimpleList"),
                 .sparseAssayReplace.SSE.numeric
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "character", "SimpleList"),
                 .sparseAssayReplace.SSE.character
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssayNames", "SparseSummarizedExperiment",
          .sparseAssayNames.SSE
)

#' @rdname SparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssayNames",
                 c("SparseSummarizedExperiment", "character"),
                 .sparseAssayNamesReplace.SSE
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
          .subsetSingleBracket.SSE
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
          .show.SSE
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
          .combine.SSE
)
