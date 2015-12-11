### =========================================================================
### Functionality that probably belongs in the GenomicRanges package
### -------------------------------------------------------------------------
###
### The functionality in this file is concerned with the development of a
### combine,SummarizedExperiment,SummarizedExperiment-method. This method
### essentially applies combine() to each slot of the SummarizedExperiment
### object. Some of these slots are objects with class definition in the
### GenomicRanges package. Therefore, the functionality in this file probably
### belongs in the GenomicRanges package.

# TODO: A ignore.mcols/use.mcols argument might be useful
#' Combining GRanges objects
#'
#' Combine multiple \link[GenomicRanges]{GenomicRanges} objects.
#'
#' @details Combines two or more \link[GenomicRanges]{GenomicRanges} objects
#' so that the resulting \link[GenomicRanges]{GenomicRanges} object contains
#' all elements of the original objects. Elements in the returned value are
#' unique, that is, an element represented in multiple arguments is represented
#' only only in the result. To perform this operation, \code{combine()} makes
#' sure that data in shared elements, \strong{including metadata columns
#' accessible with \code{\link[S4Vectors]{mcols}()}}, are identical in all the
#' \link[GenomicRanges]{GenomicRanges} objects. Data differences in shared
#' elements usually cause an error.
#'
#' \strong{WARNING}: \code{names()} of the \link[GenomicRanges]{GenomicRanges}
#' inputs are removed prior to combining, hence the returned value has
#' \code{names()} set to \code{NULL}.
#'
#' @param x A \link[GenomicRanges]{GenomicRanges} object.
#' @param y A \link[GenomicRanges]{GenomicRanges} object.
#' @param ... One or more \link[GenomicRanges]{GenomicRanges} objects.
#'
#' @return A \link[GenomicRanges]{GenomicRanges} object.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @examples
#' x <- GRanges(seqnames = c("chr1", "chr2"),
#'                ranges = IRanges(c(1, 11), c(11, 20)),
#'                seqinfo = Seqinfo(c("chr1", "chr2", "chr3")))
#' y <- GRanges(seqnames = c("chr2", "chr3"),
#'                ranges = IRanges(c(11, 21), c(20, 30)),
#'                seqinfo = Seqinfo(c("chr1", "chr2", "chr3")))
#' # x[2] (resp. y[1]) is only included once in the combined object
#' combine(x, y)
#'
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors findMatches mcols "mcols<-" queryHits
#'
#' @export
setMethod("combine", c("GenomicRanges", "GenomicRanges"),
          function(x, y, ...) {

            # NOTE: For each of x and y, need to check that any duplicate
            #       elements have identical mcols.
            x <- unname(x)
            x_is_dup <- duplicated(x)
            x_unique <- x[!x_is_dup]
            x_dup <- x[x_is_dup]
            x_ol <- findMatches(x_unique, x_dup)
            if (!identical(mcols(x_unique[queryHits(x_ol)]),
                           mcols(x_dup[subjectHits(x_ol)]))) {
              stop("'combine(x = \"", class(x), "\", y = \"", class(y),
                   "\")'\n  x contains duplicate ranges whose 'mcols()' ",
                   "differ.")
            }
            y <- unname(y)
            y_is_dup <- duplicated(y)
            y_unique <- y[!y_is_dup]
            y_dup <- y[y_is_dup]
            y_ol <- findMatches(y_unique, y_dup)
            if (!identical(mcols(y_unique[queryHits(y_ol)]),
                           mcols(y_dup[subjectHits(y_ol)]))) {
              stop("'combine(y = \"", class(y), "\", y = \"", class(y),
                   "\")'\n  y contains duplicate ranges whose 'mcols()' ",
                   "differ.")
            }

            # NOTE: unique() ignores mcols, so have to process these separately
            z <- unique(c(x_unique, y_unique, ignore.mcols = TRUE))
            x_unique_mcols <- mcols(x_unique, use.names = FALSE)
            y_unique_mcols <- mcols(y_unique, use.names = FALSE)
            # Create rownames for mcols based on z
            rownames(x_unique_mcols) <- match(x_unique, z)
            rownames(y_unique_mcols) <- match(y_unique, z)
            tryCatch({
              mcols(z) <- combine(x_unique_mcols, y_unique_mcols)
            }, error = function(err) {
              stop("\n'combine(x = \"", class(x), "\", y = \"", class(y),
                   "\")'\n  'mcols(x)' and 'mcols(y)' are not compatible.")
            })
            z
          }
)

# TODO: Would like a combine,GRangesList,GRangesList-method. However, I am not
#       familiar enough with GRangesList objects to come up with a useful
#       definition of what it means for these objects to be combine()-d.
