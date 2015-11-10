### =========================================================================
### RangedSparseSummarizedExperiment objects
### -------------------------------------------------------------------------
###
### NOTE: The class hierarchy is as follows:
###       SummarizedExperiment0
###       ├── RangedSummarizedExperiment
###       │   ├── RangedSparseSummarizedExperiment
###       ├── SparseSummarizedExperiment

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangedSparseSummarizedExperiment class
###

#' @include SparseSummarizedExperiment-class.R
#'
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom methods setClass
#'
#' @export
setClass("RangedSparseSummarizedExperiment",
         contains = "RangedSummarizedExperiment",
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

# NOTE: Can re-use validity functions for SparseSummarizedExperiment objects.
#       If it is later determined that special checks of the rowRanges slot
#       are necessary then these can be added here (they shouldn't be since the
#       RangedSummarizedExperiment validity method should capture such issues).
.valid.RangedSparseSummarizedExperiment <- function(x) {
  c(.valid.SSE.sparseAssays_class(x),
    .valid.SSE.sparseAssays_dim(x))
}

# #' @include SparseSummarizedExperiment-class.R
#' @importFrom S4Vectors setValidity2
setValidity2("RangedSparseSummarizedExperiment",
             .valid.RangedSparseSummarizedExperiment)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

get_rownames_from_sparse_assays <- function(sparse_assays) {
  if (length(sparse_assays) == 0L) {
    return(NULL)
  }
  names(sparse_assays[[1L]][[1L]][["map"]])
}

#' @importClassesFrom SummarizedExperiment SummarizedExperiment0
#'                                         RangedSummarizedExperiment
#' @importFrom GenomicRanges GRangesList
#' @importFrom methods new setMethod
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importMethodsFrom S4Vectors endoapply
#'
#' @export
setMethod("SparseSummarizedExperiment", "SparseAssays",
          function(sparseAssays,
                   assays = SimpleList(),
                   rowRanges = GRangesList(),
                   colData = DataFrame(),
                   metadata = list()) {

            # Get sample names from sparseAssays if no colData supplied
            if (missing(colData) && length(sparseAssays) != 0L) {
              nms <- names(sparseAssays[[1]])
              if (is.null(nms) && length(sparseAssays[[1]]) != 0L)
                stop("'SparseSummarizedExperiment' sparse assay names must ", "
                     not be NULL")
              colData <- DataFrame(row.names = nms)
            }

            # NOTE: The cannonical location for dimnames, including sample
            # names, is in the colData. Therefore, to simplify things, we strip
            # them from the sparseAssays object.
            sparseAssays <- endoapply(sparseAssays, unname)

            # Construct the SummarizedExperiment
            if (missing(rowRanges)) {
              se <- SummarizedExperiment(assays = assays,
                                         colData = colData,
                                         metadata = metadata)
            } else {
              se <- SummarizedExperiment(assays = assays,
                                         rowRanges = rowRanges,
                                         colData = colData,
                                         metadata = metadata)
            }

            # Check that dimensions of assays and sparseAssays are compatible
            if (nrow(se) && nrow(sparseAssays)) {
              if (!identical(dim(se), dim(sparseAssays))) {
                stop("dimensions of 'assays' (", nrow(se), " x ", ncol(se),
                     ") and 'sparseAssays' (", nrow(sparseAssays), " x ",
                     ncol(sparseAssays), ") are not compatible")
              }
            } else {
              if (ncol(se) != ncol(sparseAssays)) {
                stop("ncol of 'assays' (", ncol(se), ") and 'sparseAssays' (",
                     ncol(sparseAssays), ") are not compatible")
              }
            }

            if (missing(rowRanges)) {
              # Need to update elementMetadata slot to have the valid dimensions
              se@elementMetadata <- DataFrame()
              se@elementMetadata@nrows <- nrow(sparseAssays)
            }

            if (missing(rowRanges)) {
              val <- new("SparseSummarizedExperiment",
                         se,
                         sparseAssays = sparseAssays)
            } else {
              val <- new("RangedSparseSummarizedExperiment",
                         se,
                         sparseAssays = sparseAssays)
            }

            # NOTE: rownames are taken from sparseAssays.
            # WARNING: This will override rownames present in assays argument
            #          and used by the SummarizedExperiment constructor.
            rownames(val) <- get_rownames_from_sparse_assays(sparseAssays)
            val
          }
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("SparseSummarizedExperiment", "missing",
          function(sparseAssays, ...) {
            SparseSummarizedExperiment(SparseAssays(), ...)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# NOTE: See R/SparseSummarizedExperiment.R for definition and rational behind
#       makeSEFromSSE() and why there is no as() method.
# setAs("RangedSparseSummarizedExperiment", "RangedSummarizedExperiment",
#       .SSE.to.SE(from)
# )

#' @importClassesFrom S4Vectors SimpleList
#' @importFrom IRanges PartitioningByEnd
#' @importFrom methods as setAs
#' @importFrom methods setAs
#'
#' @export
setAs("RangedSparseSummarizedExperiment", "SparseSummarizedExperiment",
      function(from) {

        new("SparseSummarizedExperiment",
            sparseAssays = from@sparseAssays,
            as(from, "SummarizedExperiment0")
        )
      }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssays", "RangedSparseSummarizedExperiment",
            .sparseAssays.SSE
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays",
                 c("RangedSparseSummarizedExperiment", "SparseAssays"),
                 .sparseAssaysReplace.SSE
)

## convenience for common use case

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "missing"),
          .sparseAssay.SSE.missing
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "numeric"),
          .sparseAssay.SSE.numeric
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "character"),
          .sparseAssay.SSE.character
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("RangedSparseSummarizedExperiment", "missing", "SimpleList"),
                 .sparseAssayReplace.SSE.missing
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("RangedSparseSummarizedExperiment", "numeric", "SimpleList"),
                 .sparseAssayReplace.SSE.numeric
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("RangedSparseSummarizedExperiment", "character", "SimpleList"),
                 .sparseAssayReplace.SSE.character
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssayNames", "RangedSparseSummarizedExperiment",
          .sparseAssayNames.SSE
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssayNames",
                 c("RangedSparseSummarizedExperiment", "character"),
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

#' @importFrom methods setMethod
#'
#' @export
setMethod("[", "RangedSparseSummarizedExperiment",
          .subsetSingleBracket.SSE
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("[",
                 c("RangedSparseSummarizedExperiment", "ANY", "ANY",
                   "RangedSparseSummarizedExperiment"),
                 .replaceSingleBracket.SSE
)

# NOTE: extractROWS() and replaceROWS() methods inherited from
#       SummarizedExperiment0 objects.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quick colData access.
###

# NOTE: There methods are inherited from SummarizedExperiment0 objects.
# TODO: Ensure the corresponding generics are imported by the NAMESPACE

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display.
###

#' @importFrom methods setMethod
#'
#' @export
setMethod("show", "RangedSparseSummarizedExperiment",
          .show.SSE
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
###

# NOTE: Appropriate for objects with distinct ranges and identical samples.
#' @importFrom methods setMethod
#'
#' @export
setMethod("rbind", "RangedSparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .rbind.SSE(args)
          }
)

# NOTE: Appropriate for objects with identical ranges and distinct samples.
#' @importFrom methods setMethod
#'
#' @export
setMethod("cbind", "RangedSparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .cbind.SSE(args)
          }
)


#' @importFrom methods setMethod
#'
#' @export
setMethod("combine",
          c("RangedSparseSummarizedExperiment", "RangedSparseSummarizedExperiment"),
          .combine.SSE
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs
###

# NOTE: granges,RangedSummarizedExperiment doesn't honour its contract to
#       return a GRanges object, e.g., the rowRanges slot could be a
#       GenomicTuples::GTuples object. A better definition might be
#       granges <- function(x) granges(rowRanges(x)).
