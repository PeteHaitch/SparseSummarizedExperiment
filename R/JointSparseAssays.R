### =========================================================================
### JointSparseAssays objects
### -------------------------------------------------------------------------
###
### This is just an idea at this stage. Within each sparse assay, there would
### be a 'map' and a 'data' element. The 'map' element would be a SimpleList
### object with one element per sample. The 'data' element would be a matrix
### containing all observed values across samples.
### The JointSparseAssays API consists of:
###   (a) The JointSparseAssays() constructor function.
###   (b) length, names, names<-, [[, [[<-, dim, [, [<-, rbind, cbind
###
### The JointSparseAssays class has a nested-list structure. This hierarchy is
### illustrated below for an example with two sparse assays and three samples.
###
### SparseAssays
### ├── "sparse_assay_1"
### │   ├── "map"
###     │   ├── "sample1"
###     │   ├── "sample2"
###     │   ├── "sample3"
### │   ├── "data"
### ├── "sparse_assay_2"
### │   ├── "map"
###     │   ├── "sample1"
###     │   ├── "sample2"
###     │   ├── "sample3"
### │   ├── "data"
###
### Each map is an integer vector and all map elements must have identical
### length.
### Each data element is a matrix object (an assay with a 1-dimensional
### variable is stored as a 1-column matrix). Each data element may have a
### different number of rows but the maximum value of the corresponding
### map elements must be less than or equal to the number of rows of that data
### element. A row of the data element may be pointed to multiple times by the
### map elements within the sparse assay.
###
### NOTE: It may be useful to make JointSparseAssays a VIRTUAL class à la the
### Assays class in the SummarizedExperiment package. I have opted for the
### simpler implementation, at least while this package is in the experimental
### stage.
###
### NOTE: JointSparseAssays only payoff when you get more than one measurement
### per-feature, per-sample, e.g., a vector-valued measurement. The payoff is
### greater if there are lots of features with the same measurement within a
### sparse assay and/or lots of NAs per-sample. Compared to a SparseAssays
### object, the payoff increases if most samples have the same measurements
### per sparse assay.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### JointSparseAssays class
###

# JointSparseAssays objects
#
# @rdname JointSparseAssays
#
# @export
setClass("JointSparseAssays",
         contains = "SimpleList"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.JointSparseAssays <- function(x) {

  # The prototype is basically a zero-length SimpleList
  if (length(x) == 0L) {
    return(NULL)
  } else {

    # Check that sparse assay level data has an element named 'map', an element
    # named 'data', and nothing else.
    element_names <- lapply(x, function(sparse_assay) {
      names(sparse_assay)
    })
    if (any(vapply(element_names, function(en) {
      !identical(en, c("map", "data")) && !identical(en, c("data", "map"))
    }, logical(length(1))))) {
      return(paste0("Each sparse assay must have one element named 'map', an ",
                    "element named 'data', and nothing else."))
    }

    # Check all map elements are integer vectors.
    map_class <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay[["map"]], class)
    })
    map_class <- unlist(map_class, use.names = FALSE)
    if (any(map_class != "integer")) {
      return("All 'map' elements must be integer vectors.")
    }

    # Check all data elements are matrix objects.
    data_class <- lapply(x, function(sparse_assay) {
      class(sparse_assay[["data"]])
    })
    data_class <- unlist(data_class, use.names = FALSE)
    if (any(data_class != "matrix")) {
      return("All 'data' elements must be matrix objects.")
    }

    # TODO: Check all data elements have integer or double storage.mode. See
    # functions that call storage.mode for where this is implicitly assumed.

    # Check sample names are identical across sparse assays.
    sample_names <- lapply(x, function(sparse_assay) {
      names(sparse_assay[["map"]])
    })
    if (any(vapply(sample_names, function(sn, sn1) {
      !identical(sn, sn1)
    }, logical(1), sn1 = sample_names[[1]]))) {
      return("All sparse assays must have identical sample names.")
    }

    # Check all map elements have the same length
    map_length <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay[["map"]], length)
    })
    map_length <- unlist(map_length, use.names = FALSE)
    if (length(unique(map_length)) != 1L) {
      return("All 'map' elements must have identical length.")
    }

    # Check that the maximum value in each map element is less than or equal to
    # the number of rows in each corresponding data element.
    map_max <- lapply(x, function(sparse_assay) {
      max(unlist(lapply(sparse_assay[["map"]], max, na.rm = TRUE),
                 use.names = FALSE),
          na.rm = TRUE)
    })
    map_max <- unlist(map_max, use.names = FALSE)
    data_nrow <- lapply(x, function(sparse_assay) {
      nrow(sparse_assay[["data"]])
    })
    data_nrow <- unlist(data_nrow, use.names = FALSE)
    if (any(map_max > data_nrow)) {
      return(paste0("Maximum value in each 'map' element must be less than or ",
                    "equal to the number of rows in each corresponding 'data' ",
                    "element."))

    }

  }

  # TODO: Check each sparse assay has the same number of samples

  NULL
}

S4Vectors::setValidity2("JointSparseAssays", .valid.JointSparseAssays)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# NOTE: This constructor isn't much more helpful than calling
# new("JointSparseAssays", x), where x is an appropriately structured
# SimpleList object.

# SparseAssays
#
# @rdname SparseAssays
#
# @examples
# jsa <- JointSparseAssays(sparse_assays = SimpleList(
#   a1 = SimpleList(
#     map = SimpleList(
#       s1 = as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5)),
#       s2 = as.integer(c(NA, NA, 6, 7, NA, NA, 8, 9, NA, NA))
#     ),
#     data = rbind(matrix(1:10, ncol = 2),
#                  matrix(8:1, ncol = 2))),
#   a2 = SimpleList(
#     map = SimpleList(
#       s1 = as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1)),
#       s2 = as.integer(c(3, 3, 3, 4, NA, NA, NA, NA, NA, NA))
#   ),
#   data = rbind(matrix(1:2, ncol = 1),
#                matrix(4:3, ncol = 1)))
# ))
# jsa1 <- JointSparseAssays(sparse_assays = SimpleList(
#   a1 = SimpleList(
#     map = SimpleList(
#       s1 = as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5))),
#     data = matrix(1:10, ncol = 2)),
#   a2 = SimpleList(
#     map = SimpleList(
#       s1 = as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1))),
#   data = matrix(1:2, ncol = 1))
# ))
# jsa2 <- JointSparseAssays(sparse_assays = SimpleList(
#   a1 = SimpleList(
#     map = SimpleList(
#       s2 = as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA))
#     ),
#     data = matrix(8:1, ncol = 2)),
#   a2 = SimpleList(
#     map = SimpleList(
#       s2 = as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA))
#   ),
#   data = matrix(4:3, ncol = 1))
# ))
# @export
JointSparseAssays <- function(sparse_assays = SimpleList()) {
  ans <- new("JointSparseAssays", sparse_assays)
  validObject(ans)
  ans
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###
### The following are defined via inheritance to the SimpleList-method:
###   - length
###   - names
###   - names<-
###   - [[
###
### The following are specifically defined for JointSparseAssays objects:
###   - [[<-
###   - dim
###   - NROW
###   - [
###   - [<-
###

### [[<-

# NOTE: Can't defer to [[<-,SimpleList-method because it doesn't validate
# the modified object.
# @rdname JointSparseAssays
#
# @export
setReplaceMethod("[[", "JointSparseAssays",
                 function(x, i, j, ..., value) {
                   joint_sparse_assays <- S4Vectors::setListElement(x, i, value)
                   validObject(joint_sparse_assays)
                   joint_sparse_assays
                 }
)

### dim

# NOTE: dim is defined by nrow = length of map and ncol = number of samples
# @rdname JointSparseAssays
#
# @export
setMethod("dim", "JointSparseAssays",
          function(x) {
            if (length(x) == 0L) {
              return(c(0L, 0L))
            } else {
              c(length(x[[1]][["map"]][[1]]),
                length(x[[1]]))
            }
          }
)

### NROW

# NOTE: NROW,JointSparseAssays-method would otherwise defer to
# NROW,Vector-method which calls length(). This matters because otherwise
# NROW(x) != nrow(x) when x is a JointSparseAssays object.
# @rdname JointSparseAssays
#
# @export
setMethod("NROW", "JointSparseAssays",
          function(x) {
            nrow(x)
          }
)

### [

# TODO

# NOTE: This is going to be very fiddly.

# NOTE: Subsetting a SparseAssays object requires mapping the i to a new
# coordinate specified by the map element.
.extract_JointSparseAssays_subset <- function(x, i, j) {

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


# @rdname JointSparseAssays
#
# @export
setMethod("[", "JointSparseAssays",
          function(x, i, j, ..., drop = FALSE) {
            if (drop) {
              warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
            }
            stop("Not yet implemented")
            .extract_JointSparseAssays_subset(x, i, j)
          }
)

### [<-

# TODO

# NOTE: This is going to be very fiddly.

### rbind/cbind

# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO
