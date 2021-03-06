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

            # Construct minimal colData, necessary for ensuring
            # ncol(se) == ncol(sparseAssays)
            if (missing(colData) && length(assays) == 0L) {
              colData <- new("DataFrame", nrows = ncol(sparseAssays))
            }
            # Construct minimal rowData, necessary for ensuring
            # nrow(se) == nrow(sparseAssays)
            if (missing(rowData) && missing(rowRanges) &&
                length(assays) == 0L) {
              rowData <- new("DataFrame", nrows = nrow(sparseAssays))
            }

            # Construct the SE
            if ((missing(rowRanges))) {
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
            if (nrow(se) > 0 && nrow(sparseAssays) != 0 &&
                (nrow(se) != nrow(sparseAssays))) {
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

            # Check that dimnames of SummarizedExperiment and sparseAssays
            # are compatible.
            # TODO: What if this function is passed an invalid sparseAssays
            #       object, e.g, bad dimnames?
            ans_colnames <- colnames(se)
            colnames <- colnames(sparseAssays)
            if (is.null(ans_colnames)) {
              colnames(se) <- colnames(sparseAssays)
            } else {
              test <- is.null(colnames) || identical(colnames, ans_colnames)
              if (!test) {
                stop("sparseAssays colnames() must be NULL or equal colData ",
                     "rownames() / assays colnames().")
              }
            }
            ans_rownames <- rownames(se)
            rownames <- rownames(sparseAssays)
            if (is.null(ans_rownames)) {
              rownames(se) <- rownames(sparseAssays)
            } else {
              test <- is.null(rownames) || identical(rownames, ans_rownames)
              if (!test) {
                txt <- "sparseAssays rownames() must be NULL or equal rowData
                    rownames() / rowRanges names()"
                stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
              }
            }

            # Construct the SSE
            if (missing(rowRanges)) {
              new("SparseSummarizedExperiment",
                  se,
                  sparseAssays = sparseAssays)
            } else {
              new("RangedSparseSummarizedExperiment",
                  se,
                  sparseAssays = sparseAssays)
            }
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
