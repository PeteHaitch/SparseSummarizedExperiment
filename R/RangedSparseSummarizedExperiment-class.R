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

#' RangedSparseSummarizedExperiment objects
#'
#' @rdname RangedSparseSummarizedExperiment
#'
#' @include SparseAssays-class.R
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

#' @include SparseSummarizedExperiment-class.R
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

# TODO: Think about meaningful ways in which a RangedSparseSummarizedExperiment
# might be constructed. Note how SummarizedExperiment has several different
# constructors - a low-level function matched with several higher-level methods.

#' @rdname RangedSparseSummarizedExperiment
#'
#' @examples
#' sa <- SparseAssays(sparse_assays =
#'                      SimpleList(a1 =
#'                                   SimpleList(
#'                                     s1 = SimpleList(map =
#'                                                       as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5)),
#'                                                     data =
#'                                                       matrix(1:10, ncol = 2)),
#'                                     s2 = SimpleList(map =
#'                                                       as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA)),
#'                                                     data =
#'                                                       matrix(8:1, ncol = 2))),
#'                                 a2 =
#'                                   SimpleList(
#'                                     s1 = SimpleList(map =
#'                                                       as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1)),
#'                                                     data = matrix(1:2, ncol = 1)),
#'                                     s2 = SimpleList(map =
#'                                                       as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA)),
#'                                                     data = matrix(4:3, ncol = 1)))
#'
#'                      )
#' )
#' sse <- SparseSummarizedExperiment(sparseAssays = sa)
#' rr <- GRanges('chr1', IRanges(1:10, 2:11))
#' rsse <- SparseSummarizedExperiment(sparseAssays = sa,
#'                                    rowRanges = rr)
#' x <- rsse
#' y <- new("RangedSparseSummarizedExperiment")
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment0
#'                                         RangedSummarizedExperiment
#' @importFrom GenomicRanges GRangesList
#' @importFrom methods setMethod
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

#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("SparseSummarizedExperiment", "missing",
          function(sparseAssays, ...) {
            SparseSummarizedExperiment(SparseAssays(), ...)
          }
)

# TODO: SparseSummarizedExperiment,SimpleList-method and/or
#       SparseSummarizedExperiment,list-method? i.e. first try to coerce
#       SimpleList/list to a SparseAssays object and then call
#       SparseSummarizedExperiment,SparseAssays-method. When would this be
#       useful?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO: RangedSparseSummarizedExperiment -> RangedSummarizedExperiment

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssays", "RangedSparseSummarizedExperiment",
          function(x, ..., withDimnames = TRUE, expand = FALSE) {
            .sparseAssays.SSE(x, ..., withDimnames = withDimnames, expand = expand)
          }
)

# TODO
#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays", "RangedSparseSummarizedExperiment",
                 function(x, ..., value) {
                   stop("Not yet implemented")
                   .sparseAssaysReplace.SSE(x, ..., value)
                 }
)

## convenience for common use case

#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "missing"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {
            .sparseAssay.SSE.missing(x = x, withDimnames = withDimnames,
                                 expand = expand)
          }
)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "numeric"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {
            .sparseAssay.SSE.numeric(x, i, ..., withDimnames = withDimnames,
                                 expand = expand)
          }
)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "character"),
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
setReplaceMethod("sparseAssay", "RangedSparseSummarizedExperiment",
                 function(x, ..., value) {
                   stop("Not yet implemented")
                   .sparseAssayReplace.SSE(x, ..., value)
                 }
)


#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssayNames", "RangedSparseSummarizedExperiment",
          function(x, ...) {
            .sparseAssayNames.SSE(x)
          }
)

# TODO
#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssayNames", "RangedSparseSummarizedExperiment",
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

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("[", "RangedSparseSummarizedExperiment",
          function(x, i, j, ..., drop = TRUE) {
            .subsetSingleBracket.SSE(x, i, j, ..., drop = drop)
          }
)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setReplaceMethod("[",
                 c("RangedSparseSummarizedExperiment", "ANY", "ANY",
                   "RangedSparseSummarizedExperiment"),
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

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("show", "RangedSparseSummarizedExperiment",
          function(object) {
            .show.SSE(object)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
###

# NOTE: Appropriate for objects with distinct ranges and identical samples.
#' @rdname RangedSparseSummarizedExperiment
#'
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
#' @rdname RangedSparseSummarizedExperiment
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("cbind", "RangedSparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .cbind.SSE(args)
          }
)


#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("combine",
          c("RangedSparseSummarizedExperiment", "RangedSparseSummarizedExperiment"),
          function(x, y, ...) {
            .combine.SSE(x, y, ...)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs

# NOTE: granges,RangedSummarizedExperiment doesn't honour its contract to
#       return a GRanges object, e.g., the rowRanges slot could be a
#       GenomicTuples::GTuples object. A better definition might be
#       granges <- function(x) granges(rowRanges(x)).
