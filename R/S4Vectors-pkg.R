### =========================================================================
### Functionality that probably belongs in the S4Vectors package
### -------------------------------------------------------------------------
###
### The functionality in this file is concerned with the development of a
### combine,SummarizedExperiment,SummarizedExperiment-method. This method
### essentially applies combine() to each slot of the SummarizedExperiment
### object. Some of these slots are objects with class definition in the
### S4Vectors package. Therefore, the functionality in this file probably
### belongs in the S4Vectors package.

#' Combining DataFrame objects
#'
#' Combine multiple \link[S4Vectors]{DataFrame} objects using a union strategy.
#'
#' @details Objects are combined based on the \code{row.names} of \code{x},
#' \code{y}, and \code{...}.
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
#' x <- DataFrame(A = 1:10, B = letters[1:10], row.names = LETTERS[1:10])
#' y <- DataFrame(B = letters[6:10], C = 10:6, row.names = LETTERS[6:10])
#' combine(x, y)
#'
#' @importFrom methods as setMethod
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

#' Combining SimpleList objects
#'
#' Combine multiple \link[S4Vectors]{SimpleList} objects using a union strategy.
#'
#' @details Objects are combined \emph{element-wise} (i.e. using
#' \code{\link[S4Vectors]{mendoapply}}) based on the \code{names} of \code{x},
#' \code{y}, and \code{...}. This means that \code{x}, \code{y}, and all
#' elements of \code{...} must be \link[S4Vectors]{SimpleList} objects with the
#' same length and identical names, and that all elements of each
#' \link[S4Vectors]{SimpleList} must themselves have well-defined
#' \code{combine} methods.
#'
#' @param x A \link[S4Vectors]{SimpleList} object.
#' @param y A \link[S4Vectors]{SimpleList} object.
#' @param ... One or more \link[S4Vectors]{SimpleList} objects.
#'
#' @return A \link[S4Vectors]{SimpleList} object.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @examples
#' x <- SimpleList(a = DataFrame(A = 1:5, row.names = LETTERS[1:5]),
#'                 b = matrix(1:4, ncol = 2, byrow = TRUE,
#'                            dimnames = list(letters[1:2], LETTERS[1:2])))
#' y <- SimpleList(a = DataFrame(A = 3:7, row.names = LETTERS[3:7]),
#'                 b = matrix(1:8, ncol = 2, byrow = TRUE,
#'                            dimnames = list(letters[1:4], LETTERS[1:2])))
#' combine(x, y)
#'
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors mendoapply
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
