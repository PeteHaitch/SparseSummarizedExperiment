### =========================================================================
### sparsify a matrix, data.frame, or data.table
### -------------------------------------------------------------------------

#' @include AllGenerics.R
NULL

# Convert a matrix, data.frame, or data.table into a 'key' and 'value'
# elements. Basically, convert 'x' to a data.table object, use all columns as
# the keys, identify the unique rows of the data.table and map each row of
# 'x' to these unique rows.
#
# WARNING: The relative row-order of 'x' is not preserved in the returned 'data'.
#          However, the following should be TRUE when 'x' is
#          matrix:
#          identical(.densify.SimpleListSparseAssays.sample(sparsify(x)), x)
# NOTE: Returned object is stripped of dimnames
#' @param x A data.table.
#' @param value_class The class of the value element in the returned SimpleList.
#'
#' @return A SimpleList with an integer key element and a value element with
#'         class given by value_class.
#'
#' @importFrom data.table := .GRP .I key setkey setkeyv
#' @importFrom S4Vectors SimpleList
#'
.sparsify.SimpleList <- function(x,
                                 value_class =
                                   c("matrix", "data.frame", "data.table")) {


  # NOTE: Will otherwise get errors if data have zero rows.
  if (nrow(x)) {

    # Add an index for the original row number
    if (any(c(".myI", ".myKey")  %in% colnames(x))) {
      stop("'x' must not have a column named '.myI' or '.myKey'")
    }
    x[, .myI := .I]

    # Set the map (kind of like hashing the rows of the data.table since we use
    # all columns)
    my_key <- grep(".myI", colnames(x), value = TRUE, invert = TRUE)
    setkeyv(x, cols = my_key)

    # Create the map and data
    x[, .myKey := .GRP, by = key(x)]
    key <- setkey(x[, list(.myI, .myKey)], .myI)[, .myKey]
    value <- unique(x)[, c(".myI", ".myKey") := NULL]
  } else {
    value <- x
    key <- integer(0)
  }

  value <- unname(as.matrix(value))

  # Return the result
  SimpleList(key = key,
             value = value)
}
# To avoid WARNINGs about "Undefined global functions or variables" in
# R CMD check caused by the .sparsify.SimpleList() function.
#' @importFrom utils globalVariables
globalVariables(c(".myI", ".myKey"))

#' @importFrom data.table as.data.table
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparsify", c("matrix", "character"),
          function(x, return_class, ...) {
            # Convert input to data.table
            x <- as.data.table(x, keep.rownames = FALSE)
            .sparsify.SimpleList(x = x, ...)
}
)

#' @importFrom data.table setDT
#' @importFrom methods setMethod
#'
#' @export
setMethod("sparsify", c("data.frame", "character"),
          function(x, return_class, ...) {
            # Convert input to data.table by reference
            setDT(x)
            .sparsify.SimpleList(x = x, ...)
          }
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparsify", c("data.table", "character"),
          function(x, return_class, ...) {
            # Input is already a data.table
            .sparsify.SimpleList(x = x, ...)
          }
)