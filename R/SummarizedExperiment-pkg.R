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
#' @keywords internal
.combine.NAMES <- function(x, y, ...) {
  shared_names <- intersect(x, y)
  c(x, setdiff(y, shared_names))
}

# TODO: Avoid unnecessary (and possibly costly) object validation where
#       possible to safely do so.
#' Combining SummarizedExperiment0/RangedSummarizedExperiment objects
#'
#' Combine multiple \link[SummarizedExperiment]{SummarizedExperiment0} or
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} objects using a
#' union strategy.
#'
#' @details \link[SummarizedExperiment]{SummarizedExperiment0} objects are
#' combined based on the \code{names} of \code{x}, \code{y}, and \code{...},
#' while \link[SummarizedExperiment]{RangedSummarizedExperiment} are combined
#' based on finding matching genomic ranges. \strong{WARNING}: Does not
#' currently work if the \code{rowRanges} slot of \code{x}, \code{y}, or
#' \code{...} is a \link[GenomicRanges]{GRangesList} objets.
#'
#' @section Note for Developers:
#' Any class that extends the
#' \link[SummarizedExperiment]{SummarizedExperiment0} or
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} class by adding
#' additional slots will need to be careful when defining a combine() method
#' for the new class if it wants to call
#' \code{combine,SummarizedExperiment0-method} by inheritance through
#' \code{callNextMethod()}. Specifically, the \code{combine()} method will need
#' to first update these additional slots, update the corresponding slots in
#' \code{x} (thus likely making \code{x} an invalid object), and then calling
#' \code{callNextMethod()}. An alternative would be to set \code{check = FALSE}
#' when replacing the slots in the
#' \code{combine,SummarizedExperiment0,SummarizedExperiment0-method}, but this
#' would require that the validity of each slot was checked in some other way
#' to guard against generally returning of unvalidated
#' \link[SummarizedExperiment]{SummarizedExperiment0} or
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} objects.
#'
#' @param x A \link[S4Vectors]{DataFrame} object.
#' @param y A \link[S4Vectors]{DataFrame} object.
#' @param ... One or more \link[S4Vectors]{DataFrame} objects.
#'
#' @return A \link[S4Vectors]{DataFrame} object.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @examples
#' nrows <- 200; ncols <- 6
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
#'                      row.names=LETTERS[1:6])
#' se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
#'                             colData=colData)
#' names(se0) <- paste0("f", 1:200)
#' x <- se0[1:150, 1:4]
#' y <- se0[30:170, 2:6]
#' combine(x, y)
#' rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#'                      strand=sample(c("+", "-"), 200, TRUE),
#'                      feature_id=sprintf("ID%03d", 1:200))
#' rse <- unname(se0)
#' rowRanges(rse) <- rowRanges
#' x <- se0[1:150, 1:4]
#' y <- se0[30:170, 2:6]
#' combine(x, y)
#'
#' @importFrom methods is setMethod
#' @importFrom S4Vectors DataFrame
#' @importMethodsFrom IRanges findOverlaps
#' @importMethodsFrom S4Vectors mcols "mcols<-" metadata subjectHits
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

            if (!is(x, "RangedSummarizedExperiment") &&
                (is.null(names(x)) || is.null(names(y)))) {
              stop("Cannot combine '", class(x), "' objects with NULL 'names'")
            }

            # Check colnames are set
            x_cn <- colnames(x)
            y_cn <- colnames(y)
            if (is.null(x_cn) || is.null(y_cn)) {
              stop("Cannot combine '", class(x), "' objects with NULL 'colnames'")
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
                stop("Cannot combine '", class(x), "' objects with internal ",
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
                stop("Cannot combine '", class(x), "' objects with internal ",
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
