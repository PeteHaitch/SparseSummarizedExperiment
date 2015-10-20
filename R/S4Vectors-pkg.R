### =========================================================================
### Functionality that probably belongs in the S4Vectors package
### -------------------------------------------------------------------------
###
### The functionality in this file is concerned with the development of a
### combine,SummarizedExperiment0,SummarizedExperiment0-method. This method
### essentially applies combine() to each slot of the SummarizedExperiment0
### object. Some of these slots are objects with class definition in the
### S4Vectors package. Therefore, the functionality in this file probably
### belongs in the S4Vectors package.

#' @rdname S4Vectors
#'
#' @export
setMethod("combine", c("DataFrame", "DataFrame"),
          function(x, y, ...) {

            if (all(dim(x) == 0L) && all(dim(y) == 0L)) {
              return(x)
            } else if (all(dim(x) == 0L)) {
              return(y)
            } else if (all(dim(y) == 0L)) {
              return(x)
            }

            if (is.null(row.names(x)) || is.null(row.names(y))) {
              stop("'row.names' of 'x' and 'y' must be non-NULL")
            }
            as(combine(as.data.frame(x), as.data.frame(y)), "DataFrame")
          }
)

# TODO: How should elements unique to x or y be handled?
#' @rdname S4Vectors
#'
#' @export
setMethod("combine", c("SimpleList", "SimpleList"),
          function(x, y, ...) {

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            if (length(names(x)) != length(names(y))) {
              stop("'", class(x), "' objects have different number of elements:",
                   "\n\t", paste(names(x), collapse = " "), "\n\t",
                   paste(names(y), collapse = " "))
            }

            if (!all(names(x) == names(y))) {
              stop("'", class(x), "' objects have different element names:\n\t",
                   paste(names(x), collapse = " "), "\n\t",
                   paste(names(y), collapse = " "))
            }

            mendoapply(combine, x, y)
          }
)
