### =========================================================================
### Functionality that probably belongs in the SummarizedExperiment package
### -------------------------------------------------------------------------
###
### The functionality in this file is concerned with the development of a
### combine,SummarizedExperiment,SummarizedExperiment-method. This method
### essentially applies combine() to each slot of the SummarizedExperiment
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
#' Combining SummarizedExperiment/RangedSummarizedExperiment objects
#'
#' Combine multiple \link[SummarizedExperiment]{SummarizedExperiment} or
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} objects using a
#' union strategy.
#'
#' @details \link[SummarizedExperiment]{SummarizedExperiment} objects are
#' combined based on the \code{names} of \code{x}, \code{y}, and \code{...},
#' while \link[SummarizedExperiment]{RangedSummarizedExperiment} are combined
#' based on finding matching \code{rowRanges()} of \code{x}, \code{y}, and
#' \code{...}. \strong{WARNING}: Does not currently work if the
#' \code{rowRanges} slot of \code{x}, \code{y}, or \code{...} is a
#' \link[GenomicRanges]{GRangesList} objects.
#'
#' @section Note for Developers:
#' Any class that extends the
#' \link[SummarizedExperiment]{SummarizedExperiment} or
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} class by adding
#' additional slots will need to be careful when defining a combine() method
#' for the new class if it wants to call
#' \code{combine,SummarizedExperiment-method} by inheritance through
#' \code{callNextMethod()}. Specifically, the \code{combine()} method will need
#' to first update these additional slots, update the corresponding slots in
#' \code{x} (thus likely making \code{x} an invalid object), and then calling
#' \code{callNextMethod()}. An alternative would be to set \code{check = FALSE}
#' when replacing the slots in the
#' \code{combine,SummarizedExperiment,SummarizedExperiment-method}, but this
#' would require that the validity of each slot was checked in some other way
#' to guard against generally returning of unvalidated
#' \link[SummarizedExperiment]{SummarizedExperiment} or
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
#' colData <- DataFrame(Treatment = rep(c("ChIP", "Input"), 3),
#'                      row.names = LETTERS[1:6])
#' se <- SummarizedExperiment(assays = SimpleList(counts = counts),
#'                             colData = colData)
#' names(se) <- paste0("f", 1:200)
#' x <- se[1:150, 1:4]
#' y <- se[30:170, 2:6]
#' combine(x, y)
#' rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                      IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
#'                      strand = sample(c("+", "-"), 200, TRUE),
#'                      feature_id = sprintf("ID%03d", 1:200))
#' rse <- se
#' rownames(rse) <- NULL
#' rowRanges(rse) <- rowRanges
#' x <- rse[1:150, 1:4]
#' y <- rse[30:170, 2:6]
#' combine(x, y)
#'
#' @importFrom methods is setMethod
#' @importFrom S4Vectors DataFrame
#' @importMethodsFrom IRanges findOverlaps
#' @importMethodsFrom S4Vectors metadata subjectHits
#'
#' @export
setMethod("combine", c("SummarizedExperiment", "SummarizedExperiment"),
          function(x, y, ...) {

            if (all(dim(y) == 0L)) {
              return(x)
            } else if (all(dim(x) == 0L)) {
              return(y)
            }

            # Give a helpful error message if the user tries to combine a
            # RangedSummarizedExperiment object to a non-ranged
            # SummarizedExperiment.
            if ((is(x, "RangedSummarizedExperiment") &&
                 !is(y, "RangedSummarizedExperiment")) ||
                (!is(x, "RangedSummarizedExperiment") &&
                 is(y, "RangedSummarizedExperiment"))) {
              stop("Cannot combine '", class(x), "' and '", class(y),
                   "' objects because only one of these has a 'rowRanges' ",
                   "slot.")
            }

            # Check colnames are set (required for combining assays).
            x_cn <- colnames(x)
            y_cn <- colnames(y)
            if (is.null(x_cn) || is.null(y_cn)) {
              stop("Cannot combine '", class(x), "' objects with NULL ",
                   "'colnames()'")
            }

            # Cannot combine non-ranged SEs with NULL names (required for
            # combining assays).
            if (!is(x, "RangedSummarizedExperiment") &&
                (is.null(names(x)) || is.null(names(y)))) {
              stop("Cannot combine '", class(x), "' objects with NULL ",
                   "'names()'")
            }

            # NOTE: Currently, cannot combine SummarizedExperiments with
            # duplicate rows; this would require first combining the assays at
            # duplicate rows (e.g., summing the relevant rows). While this
            # might be doable, it seems in general to be a bad idea to have
            # duplicate rows in a SummarizedExperiment object.
            if (is(x, "RangedSummarizedExperiment")) {
              if (any(duplicated(x))) {
                stop("'combine(x = \"", class(x), "\", y = \"", class(y),
                     "\")'\n  'any(duplicated(x))' must be FALSE")
              }
              if (any(duplicated(y))) {
                stop("'combine(x = \"", class(x), "\", y = \"", class(y),
                     "\")'\n  'any(duplicated(y))' must be FALSE")
              }
            } else {
              if (anyDuplicated(names(x))) {
                stop("'combine(x = \"", class(x), "\", y = \"", class(y),
                     "\")'\n  'anyDuplicated(x)' must be 0 (FALSE)")
              }
              if (anyDuplicated(names(y))) {
                stop("'combine(x = \"", class(x), "\", y = \"", class(y),
                     "\")'\n  'anyDuplicated(y)' must be 0 (FALSE)")
              }
            }

            # Combine slots
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: combine,SummarizedExperiment,SummarizedExperiment-method
              #       doesn't currently work if rowRanges(x) or rowRanges(y)
              #       are GRangesList-derived objects because there is
              #       currently no combine,GRangesList,GRangesList-method.
              if (is(rowRanges(x), "GRangesList") ||
                  is(rowRanges(y), "GRangesList")) {
                stop("Cannot combine '", class(x), "' objects with ",
                     "'GRangesList'-based 'rowRanges'")
              }
              rowRanges <- combine(rowRanges(x), rowRanges(y))
            } else {
              rowData <- combine(rowData(x, use.names = TRUE),
                                 rowData(y, use.names = TRUE))
            }
            colData <- combine(colData(x), colData(y))
            x_assays <- assays(x, withDimnames = TRUE)
            y_assays <- assays(y, withDimnames = TRUE)

            x_assays_rn <-
              SummarizedExperiment:::.get_rownames_from_assays(x_assays)
            y_assays_rn <-
              SummarizedExperiment:::.get_rownames_from_assays(y_assays)

            assays_rm_rn <- FALSE
            if (is.null(x_assays_rn) || is.null(y_assays_rn)) {
              if (is(x, "RangedSummarizedExperiment")) {
                # NOTE: assay rownames are constructed by finding matches
                #       between the ranges of x (resp. y) and z.
                x_rn <- match(rowRanges(x), rowRanges)
                y_rn <- match(rowRanges(y), rowRanges)
                addRownames <- function(assay, rn) {
                  rownames(assay) <- rn
                  assay
                }
                x_assays <- endoapply(x_assays, addRownames, x_rn)
                y_assays <- endoapply(y_assays, addRownames, y_rn)
              }
              assays_rm_rn <- TRUE
            }
            assays <- combine(x_assays, y_assays)
            if (assays_rm_rn) {
              assays <- endoapply(assays, function(x) {
                rownames(x) <- NULL
                x
              })
            }
            metadata <- c(metadata(x), metadata(y))

            # Construct the combined SE
            if (is(x, "RangedSummarizedExperiment")) {
              SummarizedExperiment(assays = assays,
                                   rowRanges = rowRanges,
                                   colData = colData,
                                   metadata = metadata)
            } else {
              SummarizedExperiment(assays = assays,
                                   rowData = rowData,
                                   colData = colData,
                                   metadata = metadata)
            }
          }
)
