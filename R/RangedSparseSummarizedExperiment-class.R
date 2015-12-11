### =========================================================================
### RangedSparseSummarizedExperiment objects
### -------------------------------------------------------------------------
###

#' @include SparseSummarizedExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangedSparseSummarizedExperiment class
###

# TODO: Document combine (noting it doesn't work well for GRangesList rowRanges
#       slot).

#' RangedSparseSummarizedExperiment objects
#'
#' @description The RangedSparseSummarizedExperiment class extends the
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} class by adding the
#' \code{sparseAssays} slot, which contains a \link{SparseAssays} object.
#'
#' RangedSummarizedSparseExperiment is a subclass of both
#' \link{SparseSummarizedExperiment} and
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}, with the former
#' having precedence. As such, all the methods documented in
#' \code{?}\link{SparseSummarizedExperiment} and
#' \code{?}\link[SummarizedExperiment]{RangedSummarizedExperiment} also work on
#' a RangedSparseSummarizedExperiment object. See
#' \link{SparseSummarizedExperiment} for details.
#'
#' @usage
#' ## Constructor
#'
#' SparseSummarizedExperiment(sparseAssays, ...)
#'
#' @details See \link{SparseAssays} and
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}.
#'
#' @section Constructor:
#' RangedSparseSummarizedExperiment instances are constructed using the
#' \code{SparseSummarizedExperiment} function with arguments outlined above.
#'
#' @section Accessors: See \link{SparseAssays} and
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}.
#'
#' @section GRanges compatibility (rowRanges access):
#' See \link[SummarizedExperiment]{RangedSummarizedExperiment}.
#'
#' @section Subsetting:
#' See \link[SummarizedExperiment]{RangedSummarizedExperiment}.
#'
#' @section Extension:
#' RangedSparseSummarizedExperiment is implemented as an S4 class, and can be
#' extended in the usual way, using
#' \code{contains = "RangedSparseSummarizedExperiment"} in the new class
#' definition.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @seealso
#' \itemize{
#'  \item \link{SparseSummarizedExperiment} objects.
#'  \item \link[SummarizedExperiment]{SummarizedExperiment} objects in the
#'    \pkg{SummarizedExperiment} package.
#'  \item \link{SparseAssays} and \link{SimpleListSparseAssays} objects.
#' }
#'
#' @aliases RangedSparseSummarizedExperiment
#'          SparseSummarizedExperiment
#'
#' @examples
#' sl1 <- SimpleList(
#' s1 = SimpleList(key = as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5)),
#'                 value = matrix(1:10, ncol = 2)),
#' s2 = SimpleList(key = as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA)),
#'                 value = matrix(8:1, ncol = 2)))
#'
#' sl2 <- SimpleList(
#'   s1 = SimpleList(key = as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1)),
#'                   value = matrix(1:2, ncol = 1)),
#'   s2 = SimpleList(key = as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA)),
#'                   value = matrix(4:3, ncol = 1)))
#' sa <- SparseAssays(SimpleList(sa1 = sl1, sa2 = sl2))
#' colData <- DataFrame(Genotype = c("WT", "KO"),
#'                      row.names = c("s1", "s2"))
#' sse <- SparseSummarizedExperiment(sparseAssays = sa,
#'                                   colData = colData)
#' rowRanges <- GRanges(rep(c("chr1", "chr2"), c(3, 7)),
#'                      IRanges(c(1, 5, 11, 15, 21, 100, 200, 300, 5000, 5010),
#'                              width = 50),
#'                      strand = rep(c("+", "-"), times = 5),
#'                      feature_id = paste0("f", 1:10))
#' rsse <- SparseSummarizedExperiment(sparseAssays = sa,
#'                                    rowRanges = rowRanges,
#'                                    colData = colData)
#' rsse
#' dim(rsse)
#' dimnames(rsse)
#' sparseAssayNames(rsse)
#' sparseAssay(rsse)
#' # densify the first sparse assay.
#' # In general its a bad idea to use densify = TRUE, but these data are small
#' # enough not to worry.
#' densify(sparseAssay(rsse), 1, 1:2)[[1]]
#' # TODO: Implement saapply
#' #sparseAssays(rsse) <- saapply(assays(rse), asinh)
#' sparseAssay(rsse)
#' # densify the first sparse assay
#' densify(sparseAssay(rsse), 1, 1:2)[[1]]
#'
#' rowRanges(rsse)
#' rowData(rsse)  # same as 'mcols(rowRanges(rsse))'
#'
#' rsse[, rsse$Genotype == "WT"]
#'
#' ## cbind() combines objects with the same features of interest
#' ## but different samples:
#' rsse1 <- rsse
#' rsse2 <- rsse1[, 1]
#' colnames(rsse2) <- "s3"
#' cmb1 <- cbind(rsse1, rsse2)
#' dim(cmb1)
#' dimnames(cmb1)
#'
#' ## rbind() combines objects with the same samples but different
#' ## features of interest:
#' rsse1 <- rsse
#' rsse2 <- rsse1[1:5, ]
#' rownames(rsse2) <- letters[1:nrow(rsse2)]
#' cmb2 <- rbind(rsse1, rsse2)
#' dim(cmb2)
#' dimnames(cmb2)
#'
#' ## combine() combines objects with potentially different genomic ranges of
#' ## interest and different samples, by finding matching genomic ranges:
#' rsse1 <- rsse[1:5, ]
#' names(rsse1) <- letters[1:5]
#' rsse2 <- rsse[3:8, 2]
#' names(rsse2) <- letters[3:8]
#' cmb3 <- combine(rsse1, rsse2)
#' dim(cmb3)
#' dimnames(cmb3)
#'
#' ## Coercion to/from SparseSummarizedExperiment:
#' sse <- as(rsse, "SparseSummarizedExperiment")
#' sse
#'
#' as(sse, "RangedSparseSummarizedExperiment")

#' ## Coercion to/from RangedSummarizedExperiment
#' ## Using as() drops the sparseAssays slot
#' rse <- as(rsse, "RangedSummarizedExperiment")
#' as(rse, "RangedSparseSummarizedExperiment")
#' ## But using makeSEFromSSE() preserves the sparseAssays slot by densifying and
#' ## storing it in the assays slot.
#' rse2 <- makeSEFromSSE(rsse)
#' assays(rse2)
#' ## However, converting back does not re-sparsify the sparse assays
#' rsse2 <- as(rse2, "RangedSparseSummarizedExperiment")
#' sparseAssays(rsse2)
#'
#' ## Setting rowRanges on a SparseSummarizedExperiment object turns it into a
#' ## RangedSparseSummarizedExperiment object:
#' sse2 <- sse
#' rowRanges(sse) <- rowRanges
#' sse  # RangedSparseSummarizedExperiment
#'
#' ## Sanity checks:
#' stopifnot(identical(assays(sse), assays(rsse)))
#' stopifnot(identical(dim(sse), dim(rsse)))
#' stopifnot(identical(dimnames(sse), dimnames(rsse)))
#' stopifnot(identical(rowData(sse), rowData(rsse)))
#' stopifnot(identical(colData(sse), colData(rsse)))
#'
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
#' Get rownames of a SparseAssays object
#'
#' @param sparse_assays A SimpleListSparseAssays object.
#'
#' @keywords internal
.get_rownames_from_sparse_assays <- function(sparse_assays) {
  if (length(sparse_assays) == 0L) {
    return(NULL)
  }
  names(sparse_assays[[1L]][[1L]][["key"]])
}

# TODO: Need to generalise so as to work with all concrete subclasses of
#       SparseAssays objects.
#' Set rownames of a SparseAssays object
#'
#' @param sparse_assays A SimpleListSparseAssays object.
#'
#' @keywords internal
.set_rownames_from_sparse_assays <- function(sparse_assays, rn) {
  if (length(sparse_assays) == 0L) {
    return(sparse_assays)
  }
  names(sparse_assays[[1L]][[1L]][["key"]]) <- rn
  sparse_assays
}

# # TODO: Comment on dimnames (see ?SummarizedExperiment)
#' @param sparseAssays A \link{SparseAssays} object.
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#'
#' @rdname RangedSparseSummarizedExperiment-class
#' @importFrom GenomicRanges GRangesList
#' @importFrom methods is new setMethod
#' @importFrom S4Vectors DataFrame SimpleList
#' @importMethodsFrom S4Vectors endoapply
#'
#' @export
setMethod("SparseSummarizedExperiment", "SparseAssays",
          function(sparseAssays,
                   assays = SimpleList(),
                   rowData = NULL,
                   rowRanges = GRangesList(),
                   colData = DataFrame(),
                   metadata = list()) {

            # Ensure assays is a SimpleList object.
            # NOTE: This code is based on that for the various
            #       SummarizedExperiment constructor methods.
            if (is.list(assays) && !is.matrix(assays)) {
              assays <- do.call(SimpleList, assays)
            } else if (is.matrix(assays)) {
              if (is.list(assays)) {
                # NOTE: special case -- matrix of lists
                assays <- list(assays)
              }
              assays <- SimpleList(assays)
            }

            # Get colnames, set in canonical location, and strip from elsewhere.
            # NOTE: These may be specified in one of three locations:
            #       1. rownames of colData,
            #       2. in the sparseAssays object
            #       3. in the assays object
            #       However, 1 is the canonical location in a
            #       SparseSummarizedExperiment object. Therefore, we check
            #       each location, check the compatibility of these names,
            #       and, if compatible, store these in 1 and strip from 2 and 3.
            cn1 <- rownames(colData)
            if (length(sparseAssays)) {
              cn2 <- unique(unlist(lapply(sparseAssays, names),
                                   use.names = FALSE))
            } else {
              cn2 <- NULL
            }
            if (length(assays)) {
              cn3 <- unique(unlist(lapply(assays, colnames), use.names = FALSE))
            } else {
              cn3 <- NULL
            }
            if (!identical(cn1, cn2) && !is.null(cn1) && !is.null(cn2)) {
              stop("'rownames(colData)' not identical to names of samples in ",
                   "'sparseAssays'.")
            }
            if (!identical(cn1, cn3) && !is.null(cn1) && !is.null(cn3)) {
              stop("'rownames(colData)' not identical to colnames of 'assays'.")
            }
            cn <- unique(c(cn1, cn2, cn3))
            if (is.null(cn) && (length(assays) || length(sparseAssays))) {
              stop("'SparseSummarizedExperiment' colnames must not be NULL.",
                   "\n  Please specify using 'row.names' of 'colData'.")
            }
            if (!missing(colData)) {
              rownames(colData) <- cn
            } else {
              colData <- DataFrame(row.names = cn)
            }
            sparseAssays <- endoapply(sparseAssays, unname)
            assays <- endoapply(assays, function(assay) {
              colnames(assay) <- NULL
            })

            # Get rownames, set in canonical location, and strip from elsewhere.
            # NOTE: These may be specified in one of three locations:
            #       1. names of rowRanges (resp. rownames of rowData) for RSE
            #          (resp. non-ranged SE).
            #       2. in the sparseAssays object
            #       3. in the assays object
            #       However, 1 is the canonical location in a
            #       SparseSummarizedExperiment object. Therefore, we check
            #       each location, check the compatibility of these names,
            #       and, if compatible, store these in 1 and strip from 2 and 3.
            if (!missing(rowRanges)) {
              rn1 <- names(rowRanges)
            } else {
              rn1 <- rownames(rowData)
            }
            # TODO: .get_rownames_from_sparse_assays() only works for
            #       SimpleListSparseAssays objects. Might be useful to define
            #       a dimnames() getter for each concrete subclass of
            #       SparseAssays.
            rn2 <- .get_rownames_from_sparse_assays(sparseAssays)
            rn3 <- SummarizedExperiment:::get_rownames_from_assays(assays)
            if (!identical(rn1, rn2) && !is.null(rn1) && !is.null(rn2)) {
              stop("'names(rowRanges)' or 'rownames(rowData)' not identical ",
                   "to ranges or feature names in 'sparseAssays'.")
            }
            if (!identical(rn1, rn3) && !is.null(rn1) && !is.null(rn3)) {
              stop("'names(rowRanges)' or 'rownames(rowData)' not identical ",
                   "to ranges or feature names in 'assays'.")
            }
            rn <- unique(c(rn1, rn2, rn3))
            if (!missing(rowRanges)) {
              names(rowRanges) <- rn
            } else {
              if (is.null(rowData)) {
                rowData <- new("DataFrame", nrows = nrow(sparseAssays))
              }
              rownames(rowData) <- rn
            }
            # TODO: .set_rownames_from_sparse_assays() only works for
            #       SimpleListSparseAssays objects. Might be useful to define
            #       a dimnames() setter for each concrete subclass of
            #       SparseAssays.
            sparseAssays <- .set_rownames_from_sparse_assays(sparseAssays, rn)
            assays <- endoapply(assays, function(assay) {
              rownames(assay) <- rn
            })

            # NOTE: We don't strip colnames from 'value' elements of the
            #       SparseAssays object in the sparseAssays slot because these
            #       are not stored elsewhere (these might be used to name the
            #       measurement, e.g., MM MU UM UU for methylation patterns)

            # Construct the SummarizedExperiment
            if (missing(rowRanges)) {
              se <- SummarizedExperiment(assays = assays,
                                         rowData = rowData,
                                         colData = colData,
                                         metadata = metadata)
            } else {
              se <- SummarizedExperiment(assays = assays,
                                         rowData = rowData,
                                         rowRanges = rowRanges,
                                         colData = colData,
                                         metadata = metadata)
            }

            # Check that dimensions of SummarizedExperiment and sparseAssays
            # are compatible.
            if (nrow(se) > 0 && (nrow(se) != nrow(sparseAssays))) {
              if (length(assays)) {
                stop("'nrow(sparseAssays)' != 'nrow(assays)' (",
                     nrow(sparseAssays), " != ", nrow(se), ").")
              } else {
                if (is(se, "RangedSummarizedExperiment")) {
                  stop("'nrow(sparseAssays)' != 'length(rowRanges)' (",
                       nrow(sparseAssays), " != ", nrow(se), ").")
                } else {
                  stop("'nrow(sparseAssays)' != 'nrow(rowData)' (",
                       nrow(sparseAssays), " != ", nrow(se), ").")
                }
              }
            } else if (ncol(se) != ncol(sparseAssays)) {
              if (length(assays)) {
                stop("'ncol(sparseAssays)' != 'ncol(assays)' (",
                     ncol(sparseAssays), " != ", ncol(se), ").")
              } else {
                stop("'ncol(sparseAssays)' != 'nrow(colData)' (",
                     ncol(sparseAssays), " != ", nrow(colData), ").")
              }
            }

            # TODO: Test if rowData is supplied with names
            # if (missing(rowRanges)) {
            #   # Need to update elementMetadata slot to have the valid dimensions
            #   se@elementMetadata <- rowData
            # }
            if (missing(rowRanges)) {
              val <- new("SparseSummarizedExperiment",
                         se,
                         sparseAssays = sparseAssays)
            } else {
              val <- new("RangedSparseSummarizedExperiment",
                         se,
                         sparseAssays = sparseAssays)
            }
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
