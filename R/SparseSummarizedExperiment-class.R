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

#' @include SimpleListSparseAssays-class.R SSE-helpers.R
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
setValidity2("SparseSummarizedExperiment", .valid.SSE)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# See R/RangedSparseSummarizedExperiment.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# NOTE: Don't define this as an explicit as() method because it will break the
#       implicit coercion of SSE to SE;  this implicit/inherited coercoin drops
#       the sparseAssays slot whereas makeSEFromSSE() preserves it by expanding
#       the sparse assays and adding them to the assays slot. The
#       implicit/inherited coercion of SSE to SE is currently relied upon by
#       several functions in this package (most non-user facing).
#' @export
# setAs("SparseSummarizedExperiment", "SummarizedExperiment0",
#       .SSE.to.SE(from)
# )
makeSEFromSSE <- function(SSE, ...) {
  .SSE.to.SE(SSE)
}

#' @importClassesFrom S4Vectors SimpleList
#' @importFrom IRanges PartitioningByEnd
#' @importFrom GenomicRanges GRanges
#' @importFrom methods as setAs
#' @importFrom methods setAs
#'
#'
#' @export
setAs("SparseSummarizedExperiment", "RangedSparseSummarizedExperiment",
      function(from) {

        partitioning <- PartitioningByEnd(integer(length(from)),
                                          names = names(from))
        rowRanges <- relist(GRanges(), partitioning)
        SparseSummarizedExperiment(sparseAssays = from@sparseAssays,
                                   rowRanges = rowRanges,
                                   colData = from@colData,
                                   assays = as(from@assays, "SimpleList",
                                               strict = FALSE),
                                   metadata = from@metadata)
      }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssays", "SparseSummarizedExperiment",
          .sparseAssays.SSE
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays",
                 c("SparseSummarizedExperiment", "SparseAssays"),
                 .sparseAssaysReplace.SSE
)

## convenience for common use case

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "missing"),
          .sparseAssay.SSE.missing
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "numeric"),
          .sparseAssay.SSE.numeric
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "character"),
          .sparseAssay.SSE.character
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "missing", "SimpleList"),
                 .sparseAssayReplace.SSE.missing
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "numeric", "SimpleList"),
                 .sparseAssayReplace.SSE.numeric
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "character", "SimpleList"),
                 .sparseAssayReplace.SSE.character
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssayNames", "SparseSummarizedExperiment",
          .sparseAssayNames.SSE
)

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

#' @export
setMethod("[", "SparseSummarizedExperiment",
          .subsetSingleBracket.SSE
)

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

#' @export
setMethod("show", "SparseSummarizedExperiment",
          .show.SSE
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
###

# NOTE: Appropriate for objects with distinct features and identical samples.
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
#' @importFrom methods setMethod
#'
#' @export
setMethod("cbind", "SparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .cbind.SSE(args)
          }
)

#' @export
setMethod("combine",
          c("SparseSummarizedExperiment", "SparseSummarizedExperiment"),
          .combine.SSE
)
