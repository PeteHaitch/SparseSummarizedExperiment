### =========================================================================
### Functionality that probably belongs in the SummarizedExperiment package
### -------------------------------------------------------------------------
###
### The functionality in this file is concerned with the development of a
### combine,SummarizedExperiment0,SummarizedExperiment0-method. This method
### essentially applies combine() to each slot of the SummarizedExperiment0
### object. Some of these slots are objects with class definition in the
### SummarizedExperiment package. Therefore, the functionality in this file
### probably belongs in the SummarizedExperiment package.

# NOTE: Not exported
.combine.NAMES <- function(x, y, ...) {
  shared_names <- intersect(x, y)
  c(x, setdiff(y, shared_names))
}

#' @rdname SummarizedExperiment
#'
#' @export
setMethod("combine", c("SummarizedExperiment0", "SummarizedExperiment0"),
          function(x, y, ...) {

            if (any(dim(y) == 0L)) {
              return(x)
            } else if (any(dim(x) == 0L)) {
              return(y)
            }

            # Give a helpful error message if the user tries to combine a
            # RangedSummarizedExperiment object to an non-ranged
            # SummarizedExperiment0.
            # NOTE: Can't simply check if object is SumamrizedExperiment0
            #       object since all RangedSummarizedExperiments are
            #       SummarizedExperiment0 objects but not all
            #       SummarizedExperiment0 objects are
            #       RangedSummarizedExperiment objects.
            if ((is(x, "RangedSummarizedExperiment") &&
                 !is(y, "RangedSummarizedExperiment")) ||
                (!is(x, "RangedSummarizedExperiment") &&
                 is(y, "RangedSummarizedExperiment"))) {
              stop("Cannot combine '", class(x), "' and '", class(y),
                   "' objects because only one of these has 'rowRanges'")
            }

            # Check colnames are set
            x_cn <- colnames(x)
            y_cn <- colnames(y)
            if (is.null(x_cn) || is.null(y_cn)) {
              stop("Cannot combine '", class(x), "' with NULL 'colnames'")
            }

            # Combine slots
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: This method doesn't work if rowRanges(x) or rowRanges(y)
              #       are GRangesList-derived objects since the
              #       combine,GRangesList,GRangesList-method currently requires
              #       that these have non-NULL names.
              if (is(rowRanges(x), "GRangesList") ||
                  is(rowRanges(y), "GRangesList")) {
                stop("Cannot combine '", class(x), "' objects with ",
                     "'GRangesList'-based 'rowRanges'")
              }
              # NOTE: Can't handle internal duplicate ranges in x or y.
              if (any(duplicated(rowRanges(x))) ||
                  any(duplicated(rowRanges(y)))) {
                stop("Cannot combine '", class(x), "' with internal ",
                     "duplicate 'rowRanges'")
              }
              # NOTE: mcols(x) and mcols(y)
              #       [i.e. mcols(rowRanges(x)) and mcols(rowRanges(y))]
              #       are separately combined. This could be simplified
              #       if combine,GRanges,GRanges-method had an
              #       ignore.mcols argument.
              x_em <- mcols(x, use.names = TRUE)
              mcols(x) <- NULL
              y_em <- mcols(y, use.names = TRUE)
              mcols(y) <- NULL

              # NOTE: names(rowRanges) are set to NULL
              rowRanges <- combine(rowRanges(x), rowRanges(y))
              # Set rownames based on findOverlaps() of rowRanges()
              # NOTE: We want ranges that are identical, hence
              #       type = "equal". Also, we want to identify identical
              #       zero-width ranges, hence minoverlap = 0
              x_ol <- findOverlaps(rowRanges(x), rowRanges,
                                   type = "equal", minoverlap = 0L)
              rownames(x) <- subjectHits(x_ol)
              rownames(x_em) <- subjectHits(x_ol)
              y_ol <- findOverlaps(rowRanges(y), rowRanges,
                                   type = "equal", minoverlap = 0L)
              rownames(y) <- subjectHits(y_ol)
              rownames(y_em) <- subjectHits(y_ol)
              mcols(rowRanges) <- combine(x_em, y_em)
            } else {
              if (anyDuplicated(names(x)) || anyDuplicated(names(y))) {
                stop("Cannot combine '", class(x), "' with internal ",
                     "duplicate 'names'")
              }
              NAMES <- .combine.NAMES(names(x), names(y))
            }
            colData <- combine(colData(x), colData(y))
            assays <- Assays(combine(assays(x, withDimnames = TRUE),
                                     assays(y, withDimnames = TRUE)))
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: elementMetadata slot of a RangedSummarizedExperiment
              #       object must be a zero-column DataFrame with nrow equal to
              #       the length of the RangedSummarizedExperiment objects at
              #       all times.
              elementMetadata <- DataFrame()
              elementMetadata@nrows <- length(rowRanges)
            } else {
              # NOTE: Using mcols() rather than slot(x, "elementMetadata") so
              #       that the combined objects have rownames on which to combine.
              elementMetadata <- combine(mcols(x, use.names = TRUE),
                                         mcols(y, use.names = TRUE))
              # NOTE: Once combined, drop rownames of elementMetadata since these
              #       are given by rownames(z) when z is a SummarizedExperiment0
              #       object.
              rownames(elementMetadata) <- NULL
            }
            metadata <- c(metadata(x), metadata(y))

            # Construct the combined SE
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: No need to update NAMES slot since it must be NULL in a
              #       valid RangedSummarizedExperiment object.
              BiocGenerics:::replaceSlots(x,
                                          rowRanges = rowRanges,
                                          colData = colData,
                                          assays = assays,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            } else {
              BiocGenerics:::replaceSlots(x,
                                          colData = colData,
                                          assays = assays,
                                          NAMES = NAMES,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            }
          }
)
