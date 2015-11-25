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
# NOTE: The sparsify generic and associated methods are currently internal
#       and are not exported. This is because I don't want to encourage users
#       to manually sparsify their data. Rather, I want them to use the
#       SparseAssays() constructor, along with the combine() method.
#
# WARNING: The relative row-order of 'x' is not preserved in the returned 'data'.
#          However, the following should be TRUE when 'x' is
#          matrix:
#          identical(.densify.SimpleListSparseAssays.sample(sparsify(x)), x)
# NOTE: Returned object is stripped of dimnames
#' @param x A data.table.
#'
#' @return A SimpleList with an integer key element and a matrix value element.
#'
#' @importFrom data.table := .GRP .I key setkey setkeyv
#' @importFrom S4Vectors SimpleList
#'
#' @keywords internal
.sparsify.SimpleList <- function(x) {


  # NOTE: Will otherwise get errors if data have zero rows.
  if (nrow(x)) {

    # Add an index for the original row number
    if (any(c(".myI", ".myKey")  %in% colnames(x))) {
      stop("'x' must not have a column named '.myI' or '.myKey'")
    }
    x[, .myI := .I]

    # Set the key (kind of like hashing the rows of the data.table since we use
    # all columns)
    my_key <- grep(".myI", colnames(x), value = TRUE, invert = TRUE)
    my_key <- grep("rn", my_key, value = TRUE, invert = TRUE)
    setkeyv(x, cols = my_key)

    # Create the key and value
    x[, .myKey := .GRP, by = key(x)]
    key <- setkey(x[, list(.myI, .myKey)], .myI)[, .myKey]
    value <- unique(x)[, c(".myI", ".myKey") := NULL]
  } else {
    value <- x
    key <- integer(0)
  }

  # Fix the dimnames
  if ("rn" %in% colnames(value)) {
    rn <- value[, rn]
    value <- as.matrix(value[, rn := NULL])
    names(key) <- rn[key]
  } else {
    value <- as.matrix(value)
  }
  # NOTE: Need to NULL-ify rownames differently depending on colnames,
  #       otherwise some downstream identical() checks can fail.
  if (identical(colnames(value), paste0("V", seq_len(ncol(value))))) {
    dimnames(value) <- NULL
  } else {
    rownames(value) <- NULL
  }

  # Handle the NA-row
  NA_idx <- which(!complete.cases(value))
  if (length(NA_idx)) {
    # Take care of NA rows
    stopifnot(length(NA_idx) == 1L)
    # Update value element by dropping NA row
    value <- value[-NA_idx, , drop = FALSE]
    # Update key element to replace index by NA for NA rows
    # TODO (longterm): Probably more efficient ways to do this
    if (!is.null(names(key))) {
      names(key)[key == NA_idx] <- NA
    }
    key[key == NA_idx] <- NA
    key[!is.na(key) & key > NA_idx] <- key[!is.na(key) & key > NA_idx] - 1L
  }

  # Return the result
  SimpleList(key = key, value = value)
}
# To avoid WARNINGs about "Undefined global functions or variables" in
# R CMD check caused by the .sparsify.SimpleList() function.
#' @importFrom utils globalVariables
globalVariables(c(".myI", ".myKey"))

#' @importFrom data.table as.data.table
#' @importFrom methods setMethod
#'
#' @keywords internal
setMethod("sparsify", c("matrix", "character"),
          function(x, return_class, ...) {
            if (!identical(return_class, "SimpleList")) {
              stop("'return_class' must be 'SimpleList', no others implemented")
            }
            # Convert input to data.table
            if ("rn" %in% colnames(x)) {
              stop("'x' must not have a column named 'rn'")
            }
            x <- as.data.table(x, keep.rownames = !is.null(rownames(x)))
            .sparsify.SimpleList(x = x)
}
)

#' @importFrom data.table setDT
#' @importFrom methods setMethod
#'
#' @keywords internal
setMethod("sparsify", c("data.frame", "character"),
          function(x, return_class, ...) {
            if (!identical(return_class, "SimpleList")) {
              stop("'return_class' must be 'SimpleList', no others implemented")
            }
            # Convert input to data.table by reference
            setDT(x)
            .sparsify.SimpleList(x = x)
          }
)

#' @importFrom methods setMethod
#'
#' @keywords internal
setMethod("sparsify", c("data.table", "character"),
          function(x, return_class, ...) {
            if (!identical(return_class, "SimpleList")) {
              stop("'return_class' must be 'SimpleList', no others implemented")
            }
            # Input is already a data.table
            .sparsify.SimpleList(x = x)
          }
)
