### =========================================================================
### SimpleListJointSparseAssays objects
### -------------------------------------------------------------------------
###
### This is just an idea at this stage. Within each sparse assay, there would
### be a 'key' and a 'value' element. The 'key' element would be a SimpleList or
### IntegerList object with one element per sample. The 'value' element would
### be a matrix containing all observed values across samples.
###
### The SimpleListJointSparseAssays class has a nested-list structure. This
### hierarchy is illustrated below for an example with two sparse assays and
### three samples.
###
### SimpleListJointSparseAssays
### ├── "sparse_assay_1"
### │   ├── "key"
###     │   ├── "sample1"
###     │   ├── "sample2"
###     │   ├── "sample3"
### │   ├── "value"
### ├── "sparse_assay_2"
### │   ├── "key"
###     │   ├── "sample1"
###     │   ├── "sample2"
###     │   ├── "sample3"
### │   ├── "value"
###
### Each key is an integer vector and all key elements must have identical
### length.
### Each value element is a matrix object (an assay with a 1-dimensional
### measurement is stored as a 1-column matrix). Each value element may have a
### different number of rows but the maximum value of the corresponding
### key elements must be less than or equal to the number of rows of that value
### element. A row of the value element may be pointed to multiple times by the
### key elements within the sparse assay.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SimpleListJointSparseAssays class
###

#' @importFrom methods setClass
#
#' @export
setClass("JointSparseAssays",
         contains = "SparseAssays"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimpleListJointSparseAssays <- function(x) {

  # Check that sparse assay level data has an element named 'key', an element
  # named 'value', and nothing else.
  element_names <- lapply(x, function(sparse_assay) {
    names(sparse_assay)
  })
  if (any(vapply(element_names, function(en) {
    !identical(en, c("key", "value")) && !identical(en, c("value", "key"))
  }, logical(length(1L))))) {
    return(paste0("Each sparse assay of a '", class(x), "' must have one ",
                  "element named 'map', an element named 'data', and nothing ",
                  " else."))
  }

  # Check all key elements are integer vectors.
  is_integer_key <- vapply(x, function(sparse_assay) {
    vapply(sparse_assay[["key"]], function(sample) {
      is(sample, "integer")
    }, logical(length(sparse_assay)))
  }, logical(length(x)))
  if (!isTRUE(all(is_integer_key))) {
    return(paste0("All 'key' elements of a '", class(x), "' object must be ",
                  "integer vectors."))
  }

  # Check all value elements are numeric matrix objects.
  is_numeric_matrix_value <- lapply(x, function(sparse_assay) {
    is(class(sparse_assay[["value"]]), "matrix") &&
      is.numeric(sparse_assay[["value"]])
  })
  if (!isTRUE(all(is_numeric_matrix_value))) {
    return(paste0("All 'value' elements of a '", class(x), "' object must be ",
                  "numeric matrix objects."))
  }

  # Check each sparse assay has the same number of samples.
  n_samples <- vapply(x, function(sparse_assay) {
    length(sparse_assay[["key"]])
  }, integer(length(x)))
  if (any(n_samples != n_samples[1L])) {
    return(paste0("All sparse assays of a '", class(x), "'must have an ",
                  "identical number of samples."))
  }

  # Check sample names are identical across sparse assays.
  sample_names <- lapply(x, function(sparse_assay) {
    names(sparse_assay[["key"]])
  })
  if (any(vapply(sample_names, function(sn, sn1) {
    !identical(sn, sn1)
  }, logical(1L), sn1 = sample_names[[1L]]))) {
    return(paste0("All sparse assays of a '", class(x), "'must have ",
                  "identical sample names."))
  }

  # Check all key elements have the same length
  key_length <- lapply(x, function(sparse_assay) {
    lapply(sparse_assay[["key"]], length)
  })
  ley_length <- unlist(key_length, use.names = FALSE)
  if (length(unique(unlist(key_length))) != 1L) {
    return(paste0("All 'key' elements of a '", class(x), "' object must have ",
                  "identical length."))
  }

  # Check that the maximum value in each key element is less than or equal to
  # the number of rows in each corresponding value element.
  key_max <- lapply(x, function(sparse_assay) {
    max(unlist(lapply(sparse_assay[["key"]], max, na.rm = TRUE),
               use.names = FALSE),
        na.rm = TRUE)
  })
  key_max <- unlist(key_max, use.names = FALSE)
  value_nrow <- lapply(x, function(sparse_assay) {
    nrow(sparse_assay[["value"]])
  })
  value_nrow <- unlist(value_nrow, use.names = FALSE)
  if (any(key_max > value_nrow)) {
    return(paste0("Maximum value in each 'key' element must be less than or ",
                  "equal to the number of rows in each corresponding 'value' ",
                  "element of a '", class(x), "' object."))
  }

  NULL
}

S4Vectors::setValidity2("JointSparseAssays", .valid.JointSparseAssays)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

SparseAssays(subclass = "SimpleListJointSparseAssays")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###
### NOTE: The following are defined via inheritance to the SparseAssays-method:
###       length, NROW, names, names<-, [[, [[<-
###       The following are specifically defined for SimpleListSparseAssays
###       objects: dim, [, [<-, rbind, cbind, combine

### dim

# NOTE: dim is defined by nrow = length of key and ncol = number of samples
setMethod("dim", "JointSparseAssays",
          function(x) {
            if (length(x) == 0L) {
              return(c(0L, 0L))
            } else {
              c(length(x[[1L]][["key"]][[1L]]),
                length(x[[1L]]))
            }
          }
)


### [

# TODO

# NOTE: This is going to be very fiddly.

# NOTE: This is old code. It might be better to start afresh than to try to
#       salvage this
# NOTE: Subsetting a SparseAssays object requires mapping the i to a new
# coordinate specified by the map element.
.extract_SimpleListJointSparseAssays_subset <- function(x, i, j) {

  if (!missing(i) && !missing(j)) {

    # normalize i
    i <- S4Vectors::normalizeSingleBracketSubscript(i, x, as.NSBS = FALSE)

    fun <- function(joint_sparse_assay) {
      stop("Not yet implemented")
    }
  } else if (!missing(i)) {

    # normalize i
    i <- S4Vectors::normalizeSingleBracketSubscript(i, x, as.NSBS = FALSE)

    fun <- function(joint_sparse_assay) {

      # Map i
      ii <- lapply(joint_sparse_assay[["map"]], function(map) {
        map[i]
      })

      # Extract and sparsify the data
      iii <- na.omit(unlist(ii, use.names = FALSE))
      data <- joint_sparse_assay[["data"]][iii, ]
      sparsified <- .sparsify(data)

      # Update the map
      # Should have all(lengths(map) == length(i))
      # PROBLEM: The data.table version of .sparsify() will reorder the
      #          rows of the input. I need to figure out how to update the map;
      #          the way of doing things doesn't work.
      stop("Not yet implemented")

      stopifnot(all(lengths(map) == length(i)))

      SimpleList(map = map, data = data)

    }
  } else if (!missing(j)) {

    fun <- function(joint_sparse_assay) {

      # Map i (drop those map elements no longer required)
      stop("Not yet implemented")

    }

  }
  endoapply(x, fun)

}

setMethod("[", "JointSparseAssays",
          function(x, i, j, ..., drop = FALSE) {
            if (drop) {
              warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
            }
            stop("Not yet implemented")
            .extract_SimpleListJointSparseAssays_subset(x, i, j)
          }
)

### [<-

# TODO

# NOTE: This is going to be very fiddly.

### rbind/cbind

# TODO

### combine

# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO: .densify.SimpleListJointSparseAssays.sample(),
#       .densify.SimpleListJointSparseAssays,
#       densify,SimpleListJointSparseAssays-method.
#
# TODO: as(SimpleListJointSparseAssays, SimpleListSparseAssays) and
#       as(SimpleListSparseAssays, SimpleListJointSparseAssays)?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TODOs
###

# TODO: Check import/export tags
