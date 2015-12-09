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

# TODO: Should this be a method defined for GRanges or GenomicRanges
# TODO: Update docs - if both x and y are named then match on names, otherwise
#       use match()
#' Combining GRanges objects
#'
#' Combine multiple \link[GenomicRanges]{GenomicRanges} objects.
#'
#' @details If and only if all of \code{x}, \code{y}, and \code{...} have
#' non-\code{NULL} \code{names()}, then objects are combined using these names.
#' Otherwise, names are first created based on finding identical ranges (using
#' \code{\link[S4Vectors]{findMatches}()}) and then the objects are combined
#' using these names. All of \code{x}, \code{y}, and \code{...} must also have
#' compatible \code{\link[S4Vectors]{mcols}()}.
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
#' combine(x, y)
#'
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors mcols "mcols<-"
#'
#' @export
setMethod("combine", c("GenomicRanges", "GenomicRanges"),
          function(x, y, ...) {

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            x_names <- names(x)
            y_names <- names(y)
            # NOTE: If either x or y has NULL names, use a findMatches()-based
            #       strategy to create the names.
            if (is.null(x_names) || is.null(y_names)) {
              if (!is.null(x_names)) {
                warning("'combine(x = \"", class(x), "\", y = \"", class(y),
                     "\")'\n  using 'findMatches()'-based naming strategy ",
                     "since 'names()' of some objects are NULL.")
              }
              # NOTE: ignore.mcols = TRUE since unique() ignores mcols.
              z <- unique(c(x, y, ignore.mcols = TRUE))
              x_names <- match(x, z)
              y_names <- match(y, z)
            }
            mxy <- match(x_names, y_names)
            myx <- match(y_names, x_names)
            x_shared <- which(!is.na(mxy))
            x_unique <- which(is.na(mxy))
            y_shared <- which(!is.na(myx))
            y_unique <- which(is.na(myx))
            if (!identical(x[x_shared], y[y_shared])) {
              # NOTE: names() should always be compatible if created using
              #       findMatches()-based strategy, but no guarantees about
              #       mcols().
              stop("'combine(x = \"", class(x), "\", y = \"", class(y),
                   "\")'\n  'names()' or 'mcols()' are not compatible.")
            }
            tryCatch({
              c(x, y[y_unique])
            }, error = function(err) {
              stop("\n'combine(x = \"", class(x), "\", y = \"", class(y),
                   "\")'\n  'mcols(x)' and 'mcols(y)' are not compatible.")
            })
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
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
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
            if (length(shared_elements)) {
              x[shared_elements] <- mendoapply(combine, x[shared_elements],
                                             y[shared_elements])
            }
            c(x, y[setdiff(names(y), shared_elements)])
          }
)
