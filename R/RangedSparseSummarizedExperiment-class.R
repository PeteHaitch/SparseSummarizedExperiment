### =========================================================================
### RangedSparseSummarizedExperiment objects
### -------------------------------------------------------------------------
###
### NOTE: The class hierarchy is as follows:
###       SummarizedExperiment0
###       ├── RangedSummarizedExperiment
###       │   ├── RangedSparseSummarizedExperiment
###       ├── SparseSummarizedExperiment
###       │   ├── RangedSparseSummarizedExperiment
###       i.e. RangedSparseSummarizedExperiment is a subclass of both
###       SparseSummarizedExperiment and RangedSummarizedExperiment, although
###       SparseSummarizedExperiment has precedence.
###
### NOTE: Most methods are inherited from SparseSummarizedExperiment except
###       for those related to rowRanges, which are inherited from
###       RangedSummarizedExperiment.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangedSparseSummarizedExperiment class
###

# TODO: Document once SparseSummarizedExperiment object is documented.
#' RangedSparseSummarizedExperiment objects
#'
#' @include SparseSummarizedExperiment-class.R
#'
#' @aliases RangedSparseSummarizedExperiment
#'
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom methods setClass
#'
#' @export
setClass("RangedSparseSummarizedExperiment",
         contains = c("SparseSummarizedExperiment",
                      "RangedSummarizedExperiment")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# TODO: Need to generalise so as to work with all concrete subclasses of
#       SparseAssays objects.
# TODO: Should this be an internal method (i.e. with a "." prefix).
#' Get rownames of a SparseAssays object
#'
#' @param sparse_assays A SimpleListSparseAssays object.
get_rownames_from_sparse_assays <- function(sparse_assays) {
  if (length(sparse_assays) == 0L) {
    return(NULL)
  }
  names(sparse_assays[[1L]][[1L]][["key"]])
}

#' @param sparseAssays A \link{SparseAssays} object. All dimension names (if
#'        present) must be consistent across elements and with the row names of
#'        \code{rowRanges} and \code{colData}.
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#'
#' @rdname RangedSparseSummarizedExperiment-class
#'
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

#' @rdname RangedSparseSummarizedExperiment-class
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("SparseSummarizedExperiment", "missing",
          function(sparseAssays, ...) {
            SparseSummarizedExperiment(SparseAssays(), ...)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs
###

# NOTE: granges,RangedSummarizedExperiment doesn't honour its contract to
#       return a GRanges object, e.g., the rowRanges slot could be a
#       GenomicTuples::GTuples object. A better definition might be
#       granges <- function(x) granges(rowRanges(x)).
