### =========================================================================
### Functionality that probably belongs in the GenomicRanges package
### -------------------------------------------------------------------------
###
### The functionality in this file is concerned with the development of a
### combine,SummarizedExperiment0,SummarizedExperiment0-method. This method
### essentially applies combine() to each slot of the SummarizedExperiment0
### object. Some of these slots are objects with class definition in the
### GenomicRanges package. Therefore, the functionality in this file probably
### belongs in the GenomicRanges package.

#' Combining GRanges objects
#'
#' Combine multiple \link[GenomicRanges]{GRanges} objects using a union strategy.
#'
#' @details Unlike other \code{combine} methods (e.g.,
#' \code{\link[BiocGenerics]{combine,matrix,matrix-method}}), this does not
#' make use of the \code{dimnames}/\code{names} of \code{x}, \code{y}, and
#' \code{...}. In fact, the \code{names} are ignored entirely and the method
#' effectively calls \code{unique(c(x, y, ...))}.
#'
#' @param x A \link[GenomicRanges]{GRanges} object.
#' @param y A \link[GenomicRanges]{GRanges} object.
#' @param ... One or more \link[GenomicRanges]{GRanges} objects.
#'
#' @return A \link[GenomicRanges]{GRanges} object.
#'
#' @author Peter Hickey, \url{peter.hickey@gmail.com}
#'
#' @examples
#' x <- GRanges(seqnames = c("chr1", "chr2"),
#'                ranges = IRanges(c(1, 11), c(11, 20)),
#'                seqinfo = Seqinfo(c("chr1", "chr2", "chr3")))
#' y <- GRanges(seqnames = c("chr2", "chr3"),
#'                ranges = IRanges(c(11, 21), c(20, 30)),
#'                seqinfo = Seqinfo(c("chr1", "chr2", "chr3")))
#' combine(x, y)
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("combine", c("GRanges", "GRanges"),
          function(x, y, ...) {

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            unique(c(x, y))
          }
)

# NOTE: Errors if any of the GRangesList objects have NULL names().
# TODO: How to handle unnamed GRangesList objects?
#' Combining GRangesList objects
#'
#' Combine multiple \link[GenomicRanges]{GRangesList} objects using a union
#' strategy.
#'
#' @details Objects are combined based on the \code{names} of \code{x},
#' \code{y}, and \code{...}.
#'
#' @param x A \link[GenomicRanges]{GRangesList} object.
#' @param y A \link[GenomicRanges]{GRangesList} object.
#' @param ... One or more \link[GenomicRanges]{GRangesList} objects.
#'
#' @return A \link[GenomicRanges]{GRanges} object.
#'
#' @author Peter Hickey, \url{peter.hickey@gmail.com}
#'
#' @examples
#' x <- GRangesList(
#'   "a" = GRanges(seqnames = "chr1",
#'                 ranges = IRanges(1, 10),
#'                 seqinfo = Seqinfo(c("chr1", "chr2", "chr3"))),
#'   "b" = GRanges(seqnames = "chr2",
#'                 ranges = IRanges(11, 20),
#'                 seqinfo = Seqinfo(c("chr1", "chr2", "chr3"))))
#' # NOTE: y$b is slightly different to x$b
#' y <- GRangesList(
#'   "b" = GRanges(seqnames = "chr2",
#'                 ranges = IRanges(12, 21),
#'                 seqinfo = Seqinfo(c("chr1", "chr2", "chr3"))),
#'   "c" = GRanges(seqnames = "chr3",
#'                 ranges = IRanges(21, 30),
#'                 seqinfo = Seqinfo(c("chr1", "chr2", "chr3"))))
#' # NOTE: the 'b' element of the combined object include both x$b and y$b
#' combine(x, y)
#'
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors mendoapply
#'
#' @export
setMethod("combine", c("GRangesList", "GRangesList"),
          function(x, y, ...) {

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            if (is.null(names(x)) || is.null(names(y))) {
              stop("'names' of 'x' and 'y' must be non-NULL when combining '",
                   class(x), "' objects")
            }

            shared_elements <- intersect(names(x), names(y))
            x[shared_elements] <- mendoapply(combine, x[shared_elements],
                                             y[shared_elements])
            c(x, y[setdiff(names(y), shared_elements)])
          }
)
