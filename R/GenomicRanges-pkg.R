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

#' Code for GenomicRanges package
#'
#' @rdname GenomicRanges-pkg
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
#' @rdname GenomicRanges-pkg
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
