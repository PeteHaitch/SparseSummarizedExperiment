### =========================================================================
### SparseAssays objects
### -------------------------------------------------------------------------
###
###
### The SparseAssays API consists of:
###   (a) The SparseAssays() constructor function.
###   (b) length, names, names<-, [[, [[<-, dim, [, [<-, rbind, cbind
###
### The SparseAssays class has a nested-list structure. This hierarchy is
### illustrated below for an example with two sparse assays and three samples.
###
### SparseAssays
### ├── "sparse_assay_1"
### │   ├── "sample_1"
###     │   ├── map_1_1
###     │   ├── data_1_1
### │   ├── "sample_2"
###     │   ├── map_1_2
###     │   ├── data_1_2
### │   ├── "sample3"
###     │   ├── map_1_3
###     │   ├── data_1_3
### ├── "sparse_assay_2"
### │   ├── "sample_1"
###     │   ├── map_2_1
###     │   ├── data_2_1
### │   ├── "sample_2"
###     │   ├── map_2_2
###     │   ├── data_2_2
### │   ├── "sample_3"
###     │   ├── map_2_3
###     │   ├── data_2_3
###
###
### Each map is an integer vector and all map_*_* elements must have identical
### length.
### Each data element is a matrix object (an assay with a 1-dimensional
### variable is stored as a 1-column matrix). Each data_*_* element may have a
### different number of rows but the maximum number of rows must be less than
### or equal to the length of the map_*_* elements. A row of the data_*_*
### element may be pointed to multiple times by the map_*_* element within the
### same sample and asssay.
###
### NOTE: It may be useful to make SparseAssays a VIRTUAL class à la the
### Assays class in the SummarizedExperiment package. I have opted for the
### simpler implementation, at least while this package is in the experimental
### stage.
###
### NOTE: SparseAssays only payoff when you get more than one measurement
### per-feature, per-sample. The payoff is greater if there are lots of
### features with the same measurement within a sample and/or lots of NAs
### per-sample.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseAssays class
###

#' SparseAssays objects
#'
#' @rdname SparseAssays
#'
#' @export
setClass("SparseAssays",
         contains = "SimpleList"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SparseAssays <- function(x) {

  # The prototype is basically a zero-length SimpleList
  if (length(x) == 0L) {
    return(NULL)
  } else {

    # Check that sample level data has an element named 'map', an element named
    # 'data', and nothing else.
    element_names <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay, function(sample) {
        names(sample)
      })
    })
    element_names <- unlist(element_names, recursive = FALSE, use.names = FALSE)
    if (any(vapply(element_names, function(en) {
      !identical(en, c("map", "data")) && !identical(en, c("data", "map"))
    }, logical(length(1))))) {
      return(paste0("All sample-level data within each sparse assay must ",
                    "have one element named 'map', an element named 'data', ",
                    "and nothing else."))
    }

    # Check all map elements are integer vectors.
    map_class <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay, function(sample) {
        class(sample[["map"]])
      })
    })
    map_class <- unlist(map_class, use.names = FALSE)
    if (any(map_class != "integer")) {
      return("All 'map' elements must be integer vectors.")
    }

    # Check all data elements are matrix objects.
    data_class <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay, function(sample) {
        class(sample[["data"]])
      })
    })
    data_class <- unlist(data_class, use.names = FALSE)
    if (any(data_class != "matrix")) {
      return("All 'data' elements must be matrix objects.")
    }

    # TODO: Check all data elements have integer or double storage.mode. See
    # functions that call storage.mode for where this is implicitly assumed.

    # Check sample names are identical across sparse assays.
    sample_names <- lapply(x, function(sparse_assay) {
      names(sparse_assay)
    })
    if (any(vapply(sample_names, function(sn, sn1) {
      !identical(sn, sn1)
    }, logical(1), sn1 = sample_names[[1]]))) {
      return("All sparse assays must have identical sample names.")
    }

    # Check all map elements have the same length
    map_length <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay, function(sample) {
        length(sample[["map"]])
      })
    })
    if (length(unique(unlist(map_length))) != 1L) {
      return("All 'map' elements must have identical length.")
    }

    # Check all data elements within each sparse assay have the same number of
    # columns.
    # NOTE: More generally, if data were an n-dimensional array rather than a
    # 2-dimensional matrix, would require that all dimensions except the number
    # of rows were identical within a sparse assay.
    data_ncol <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay, function(sample) {
        ncol(sample[["data"]])
      })
    })
    if (any(vapply(data_ncol, function(sparse_assay_ncol) {
      length(unique(unlist(sparse_assay_ncol))) != 1L
    }, logical(1)))) {
      return(paste0("All 'data' elements within each sparse assay must have ",
                    "the same number of columns."))
    }

    # Check that the maximum value in each map element is less than or equal to
    # the number of rows in each corresponding data element.
    map_max <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay, function(sample) {
        # Suppress warning about taking max of an empty vector
        suppressWarnings(max(sample[["map"]], na.rm = TRUE))
      })
    })
    map_max <- unlist(map_max, use.names = FALSE)
    data_nrow <- lapply(x, function(sparse_assay) {
      lapply(sparse_assay, function(sample) {
        nrow(sample[["data"]])
      })
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

S4Vectors::setValidity2("SparseAssays", .valid.SparseAssays)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# NOTE: This constructor isn't much more helpful than calling
# new("SparseAssays", x), where x is an appropriately structured SimpleList
# object.
#
#' SparseAssays
#'
#' @rdname SparseAssays
#'
#' @examples
#' sa <- SparseAssays(sparse_assays =
#'                      SimpleList(a1 =
#'                                   SimpleList(
#'                                     s1 = SimpleList(map =
#'                                                       as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5)),
#'                                                     data =
#'                                                       matrix(1:10, ncol = 2)),
#'                                     s2 = SimpleList(map =
#'                                                       as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA)),
#'                                                     data =
#'                                                       matrix(8:1, ncol = 2))),
#'                                 a2 =
#'                                   SimpleList(
#'                                     s1 = SimpleList(map =
#'                                                       as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1)),
#'                                                     data = matrix(1:2, ncol = 1)),
#'                                     s2 = SimpleList(map =
#'                                                       as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA)),
#'                                                     data = matrix(4:3, ncol = 1)))
#'
#'                      )
#' )
#'
#'
#' value <- SparseAssays(sparse_assays =
#'                         SimpleList(a1 =
#'                                      SimpleList(
#'                                        s1 = SimpleList(map =
#'                                                          as.integer(c(NA, 1, NA, 2, 3)),
#'                                                        data =
#'                                                          matrix(c(3, 14, 15, 8, 19, 20), ncol = 2)),
#'                                        s2 = SimpleList(map =
#'                                                          as.integer(c(NA, 1, 2, NA, NA)),
#'                                                        data =
#'                                                          matrix(c(6, 15, 2, 11), ncol = 2))),
#'                                    a2 =
#'                                      SimpleList(
#'                                        s1 = SimpleList(map =
#'                                                          as.integer(c(NA, 1, NA, NA, 1)),
#'                                                        data = matrix(11, ncol = 1)),
#'                                        s2 = SimpleList(map =
#'                                                          as.integer(c(NA, NA, NA, NA, NA)),
#'                                                        data = matrix(, ncol = 1)))
#'
#'                         )
#' )
#'
#' @export
SparseAssays <- function(sparse_assays = SimpleList()) {
  ans <- new("SparseAssays", sparse_assays)
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
### The following are specifically defined for SparseAssays objects:
###   - [[<-
###   - dim
###   - NROW
###   - [
###   - [<-
###


### [[<-

# NOTE: Can't defer to [[<-,SimpleList-method because it doesn't validate
# the modified object.
#' @rdname SparseAssays
#'
#' @export
setReplaceMethod("[[", "SparseAssays",
                 function(x, i, j, ..., value) {
                   sparse_assays <- S4Vectors::setListElement(x, i, value)
                   validObject(sparse_assays)
                   sparse_assays
                 }
)

### dim

# NOTE: dim is defined by nrow = length of map and ncol = number of samples
#' @rdname SparseAssays
#'
#' @export
setMethod("dim", "SparseAssays",
          function(x) {
            if (length(x) == 0L) {
              return(c(0L, 0L))
            } else {
              c(length(x[[1]][[1]][["map"]]),
                length(x[[1]]))
            }
          }
)

### NROW

# NOTE: NROW,SparseAssays-method would otherwise defer to NROW,Vector-method
# which calls length(). This matters because otherwise NROW(x) != nrow(x) when
# x is a SparseAssays object.
#' @rdname SparseAssays
#'
#' @export
setMethod("NROW", "SparseAssays",
          function(x) {
            nrow(x)
          }
)

### [

# NOTE: Subsetting a SparseAssays object requires mapping the i to a new
# coordinate specified by the map element.
.extract_SparseAssays_subset <- function(x, i, j) {

  if (!missing(i) && !missing(j)) {

    # normalize i
    i <- S4Vectors::normalizeSingleBracketSubscript(i, x, as.NSBS = FALSE)

    fun <- function(sparse_assay) {
      endoapply(sparse_assay[j], function(sample) {
        # Map i
        ii <- na.omit(sample[["map"]][i])
        # Extract using mapped i
        data <- sample[["data"]][ii, , drop = FALSE]
        # Sparsify the data
        sparsified <- .sparsify(data)
        # Update the map
        # Should have length(map) == length(i)
        if (!is.null(attr(ii, "na.action"))) {
          map <- rep(NA_integer_, length(i))
          map[-attr(ii, "na.action")] <- sparsified[["map"]]
        } else {
          map <- sparsified[["map"]]
        }
        stopifnot(length(map) == length(i))

        SimpleList(map = map,
                   data = sparsified[["data"]])
      })
    }
  } else if (!missing(i)) {

    # normalize i
    i <- S4Vectors::normalizeSingleBracketSubscript(i, x, as.NSBS = FALSE)

    fun <- function(sparse_assay) {
      endoapply(sparse_assay, function(sample) {
        # Map i
        ii <- na.omit(sample[["map"]][i])
        # Extract using mapped i
        data <- sample[["data"]][ii, , drop = FALSE]
        # Sparsify the data
        sparsified <- .sparsify(data)
        # Update the map
        # Should have length(map) == length(i)
        if (!is.null(attr(ii, "na.action"))) {
          map <- rep(NA_integer_, length(i))
          map[-attr(ii, "na.action")] <- sparsified[["map"]]
        } else {
          map <- sparsified[["map"]]
        }
        stopifnot(length(map) == length(i))

        SimpleList(map = map,
                   data = sparsified[["data"]])
      })
    }
  } else if (!missing(j)) {
    fun <- function(sparse_assay) {
      sparse_assay[j]
    }
  }
  endoapply(x, fun)
}

#' @rdname SparseAssays
#'
#' @export
setMethod("[", "SparseAssays",
          function(x, i, j, ..., drop = FALSE) {
            if (drop) {
              warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
            }
            .extract_SparseAssays_subset(x, i, j)
          }
)

### [<-

# IDEA (when !missing(i)):
# (1) Create new data element by
#     rbind(sample[["data"]], v_sample[["data"]])
# (2) Create new map element
#   (a) Add nrow(sample[["data"]]) to all elements of v_sample[["map"]]
#       to account for rbind() operation.
#   (b) sample[["map"]][i] <- v_sample[["map"]].
# (3) Subset data to only contain the required rows by
#     data[na.omit(unique(map)), , drop = FALSE]
# (4) Sparsify the data
# (5) Re-map the non-NA values of map by the "sparsified" map.
#
# NOTE: Steps 3-5 are "sparsifying" the data. The "expansion" of
# (map, data) at (2) and (5) should be identical, even though the
# individual elements may not be identical.
#
# IDEA (when missing(i)):
# (1) Simply replace the j-th sample(s) (map, data)-pair by that given in value.
.replace_SparseAssays_subset <- function(x, i, j, value) {

  if (!missing(i) && !missing(j)) {
    # Sanity check j
    if (length(j) != ncol(value)) {
      stop("length(j) != ncol(value)")
    }

    fun <- function(sparse_assay, v_sparse_assay) {
      sparse_assay[j] <- mendoapply(function(sample, v_sample) {

        # (1)
        data <- rbind(sample[["data"]], v_sample[["data"]])
        # (2)
        # NOTE: NAs are correctly propogated since NA + number = NA.
        vsm_updated <- v_sample[["map"]] + nrow(sample[["data"]])
        map <- sample[["map"]]
        map[i] <- vsm_updated
        # (3) and (4)
        sparsified <- .sparsify(data[na.omit(unique(map)), , drop = FALSE])
        # (5)
        new_lvls <- sparsified[["map"]]
        old_lvls <- na.omit(unique(map))
        map[!is.na(map)] <- new_lvls[match(map[!is.na(map)], old_lvls)]
        SimpleList(map = map,
                   data = sparsified[["data"]])
      }, sample = sparse_assay[j],
      v_sample = v_sparse_assay) # No need to subset v_sparse_assay by j
      sparse_assay
    }
  } else if (!missing(i)) {
    fun <- function(sparse_assay, v_sparse_assay) {
      sparse_assay <- mendoapply(function(sample, v_sample) {

        # (1)
        data <- rbind(sample[["data"]], v_sample[["data"]])
        # (2)
        # NOTE: NAs are correctly propogated since NA + number = NA.
        vsm_updated <- v_sample[["map"]] + nrow(sample[["data"]])
        map <- sample[["map"]]
        map[i] <- vsm_updated
        # (3) and (4)
        sparsified <- .sparsify(data[na.omit(unique(map)), , drop = FALSE])
        # (5)
        new_lvls <- sparsified[["map"]]
        old_lvls <- na.omit(unique(map))
        map[!is.na(map)] <- new_lvls[match(map[!is.na(map)], old_lvls)]
        SimpleList(map = map,
                   data = sparsified[["data"]])
      }, sample = sparse_assay, v_sample = v_sparse_assay)
      sparse_assay
    }
  } else if (!missing(j)) {
    # Sanity check j
    if (length(j) > ncol(value)) {
      stop("j > ncol(value)")
    }
    # Sanity check number of replacement rows
    # NOTE: Only necessary if missing(i)
    if (nrow(x) != nrow(value)) {
      stop("Cannot replace on j if nrow(x) != nrow(value)")
    }

    fun <- function(sparse_assay, v_sparse_assay) {
      sparse_assay[j] <- v_sparse_assay
      sparse_assay
    }
  }

  # Normalize value
  value <- S4Vectors::normalizeSingleBracketReplacementValue(value, x, i)

  # Loop over each sparse assay and do replacement
  val <- mendoapply(fun, x, value)

  # NOTE: Sanity check (shouldn't be necessary and may kill performance, but
  # until I have good unit tests in place this stays).
  validObject(val)

  val
}

#' @rdname SparseAssays
#'
#' @examples
#' sa_ <- sa
#' # Seems to work
#' sa[1, ] <- value[5, ]
#'
#' # Seems to work
#' sa <- sa_
#' sa[1, 1] <- value[5, 1]
#'
#' # Rightfully fails validObject(val)
#' sa <- sa_
#' \dontrun{sa[, 1] <- value[, 1]}
#'
#' # Fails but not for obvious reason
#' \dontrun{sa[1, 1] <- value[5, 1:2]}
#'
#' # Seems to work (basically does an rbind)
#' sa[11, ] <- value[1, ]
#'
#' # Rightfully fails validObject(val)
#' \dontrun{sa[11, 1] <- value[1, 1]}
#' @export
setReplaceMethod("[", "SparseAssays",
                 function(x, i, j, ..., value) {
                   .replace_SparseAssays_subset(x, i, j, value)
                 }
)

### rbind/cbind

.bind_SparseAssays <- function(lst, bind) {
  # This is copied from SummarizedExperiment:::.bind_Assays(). I'm not sure
  # how this could be passed a zero-length list, but I keep it here until I
  # know it's safe to remove.
  if (length(lst) == 0L) {
    return(SparseAssays())
  }

  # If the list has only a single element then just return that element.
  if (length(lst) == 1L) {
    return(lst[[1]])
  }

  lens <- sapply(lst, length)
  len1 <- lens[1L]
  if (any(lens != len1)) {
    stop("elements in sparse assays must have the same length")
  }
  if (len1 == 0L) {
    return(SparseAssays())
  }

  # Check that samples names are unique for cbind and are identical for rbind
  sample_names <- lapply(lst, function(e) {
    lapply(e, names)
  })
  # Don't need check the first against itself
  sample_names_identical <- vapply(sample_names[-1], function(sn, sn1) {
    identical(sn, sn1)
  }, FUN.VALUE = logical(1L), sn1 = sample_names[[1]])
  if (identical(bind, cbind)) {
    if (any(sample_names_identical)) {
      stop(paste0("Sample names (if present) must be unique when calling ",
                  "cbind() on 'SparseAssays'"))
    }
  } else {
    if (!all(sample_names_identical)) {
      stop(paste0("Sample names (if present) must be identical when calling ",
                  "rbind() on 'SparseAssays'"))
    }
  }

  # Check all elements of lst have the same sparse assay names
  sparse_assay_names <- lapply(lst, names)
  if (any(vapply(sparse_assay_names, function(san, san1) {
    !identical(san, san1)
  }, logical(1), san1 = sparse_assay_names[[1]]))) {
    stop("All 'SparseAssay' objects must have the same sparse assay names.")
  }
  sparse_assay_names <- sparse_assay_names[[1]]

  if (identical(bind, rbind)) {
    # If rbind-ing, need to check that all data elements within each sparse
    # assay have the same number of columns.
    same_ncol <- lapply(sparse_assay_names, function(san) {
      l_sparse_assay <- lapply(lst, "[[", san)
      ncol <- lapply(l_sparse_assay, function(sparse_assay) {
        lapply(sparse_assay, function(sample) {
          ncol(sample[["data"]])
        })
      })
      ncol <- unlist(ncol, use.names = FALSE)
      all(ncol == ncol[1])
    })
    same_ncol <- unlist(same_ncol, use.names = FALSE)
    ncol <- unlist(ncol, use.names = FALSE)
    if (any(!same_ncol)) {
      stop(paste0("Can only rbind 'SparseAssays' objects where the data ",
                  "elements within each sparse assay have the same number of ",
                  "columns."))
    }

    # rbind,SparseAssays-method uses the SparseAssays,`[<-`-method to
    # recursively add the next SparseAssays object to the end of the already
    # rbind-ed SparseAssays objects.
    # It's not the most efficient way, but it works and avoids repeating
    # much of the code used by SparseAssays,`[<-`-method.
    # This assumes that length(lst) > 1, which it should be given above checks
    # on length(lst).
    val <- lst[[1]]
    for (idx in seq.int(from = 2L, to = length(lst), by = 1L)) {
      val_nrow <- nrow(val)
      i <- seq.int(from = val_nrow + 1,
                   to = val_nrow + nrow(lst[[idx]]),
                   by = 1L)
      val[i, ] <- lst[[idx]]
    }
  } else {
    # If cbind()-ing, need to check that all map elements have the same length.
    same_length <- lapply(sparse_assay_names, function(san) {
      l_sparse_assay <- lapply(lst, "[[", san)
      length <- lapply(l_sparse_assay, function(sparse_assay) {
        lapply(sparse_assay, function(sample) {
          length(sample[["map"]])
        })
      })
      length <- unlist(length, use.names = FALSE)
      all(length == length[1])
    })
    same_length <- unlist(same_length, use.names = FALSE)
    length <- unlist(length, use.names = FALSE)
    if (any(!same_length)) {
      stop(paste0("Can only cbind 'SparseAssays' objects where the map ",
                  "elements within each sparse assay have the same length."))
    }
    val <- lapply(sparse_assay_names, function(san) {
      l_sparse_assay <- lapply(lst, "[[", san)
      do.call("c", l_sparse_assay)
    })
    val <- as(setNames(SimpleList(val), sparse_assay_names), "SparseAssays")
  }

  val
}

# NOTE: Can't defer to rbind,SimpleList-method because it in turn defers to
# rbind,ANY-method, which fails because a SparseAssays object cannot be coerced
# to a vector (and even if it could, the resulting operation probably wouldn't
# make sense).
#' @rdname SparseAssays
#'
#' @export
setMethod("rbind", "SparseAssays",
          function(..., deparse.level = 1) {
            .bind_SparseAssays(unname(list(...)), rbind)
          }
)

# NOTE: Can't defer to cbind,SimpleList-method because it in turn defers to
# cbind,ANY-method, which returns a matrix with elements being list objects
# (and doesn't make much sense for SparseAssays objects).
# WARNING: Not really a cbind. It adds new elements to the 'sparse_assay'-level
# SimpleList.
#' @rdname SparseAssays
#'
#' @export
setMethod("cbind", "SparseAssays",
          function(..., deparse.level = 1) {
            .bind_SparseAssays(unname(list(...)), cbind)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

.expand.SparseAssays.sample <- function(sample, na.rm = FALSE) {
  if (na.rm) {
    stop("'na.rm = TRUE' not yet implemented.")
  }
  sample[["data"]][sample[["map"]], , drop = FALSE]
}

# NOTE: Not a method for the time being; also cannot be called expand() because
# of IRanges::expand().
.expand <- function(x, na.rm = FALSE, assays = FALSE) {

  if (na.rm) {
    stop("'na.rm = TRUE' not yet implemented.")
  }

  if (assays) {
    l <- lapply(x, function(sparse_assay) {
      # A kludge to guess whether the data are integer or numeric. If multiple
      # data storage modes are found then assume numeric.
      data_storage_mode <- lapply(sparse_assay, function(sample) {
        storage.mode(sample[["data"]])
      })
      data_storage_mode <- unlist(data_storage_mode)
      if (all(data_storage_mode == "integer")) {
        val <- array(NA_integer_,
                     dim = c(nrow(x),
                             ncol(x),
                             ncol(sparse_assay[[1]][["data"]])),
                     dimnames = list(NULL, names(sparse_assay), NULL))
      } else {
        val <- array(NA_real_,
                     dim = c(nrow(x),
                             ncol(x),
                             ncol(sparse_assay[[1]][["data"]])),
                     dimnames = list(NULL, names(sparse_assay), NULL))
      }

      # TODO (longterm): Investigate a Rcpp version
      # Fill val with the "expanded" data
      for (sample in seq_along(sparse_assay)) {
        val[ , sample, ] <- sparse_assay[[sample]][["data"]][
          sparse_assay[[sample]][["map"]], , drop = FALSE]
      }
      val
    })
    SummarizedExperiment::Assays(l)

  } else {
    lapply(x, function(sparse_assay) {
      lapply(sparse_assay, .expand.SparseAssays.sample)
    })
  }
}

#' @rdname SparseAssays
#'
#' @name as
#'
#' @export
setAs("SparseAssays", "Assays",
      function(from) {
        .expand(from, na.rm = FALSE, assays = TRUE)
      }
)

# TODO: Move to utils.R or similar
.rowHash <- function(m) {
  # NOTE: This code uses the parts of the digest::digest() function that are
  # necessary for my needs. This approach is approximately 4-times faster than
  # using apply(X = m, MARGIN = 1L, digest).
  # TODO: Check how to properly call a NativeSymbolInfo from another package.
  # See http://r-pkgs.had.co.nz/src.html#clang. Perhaps the definition of
  # digest should be done via .onLoad()
  digest <- getNativeSymbolInfo(name = "digest", PACKAGE = "digest")
  apply(X = m, MARGIN = 1L, FUN = function(object) {
    object <- serialize(object, connection = NULL, ascii = FALSE)

    .Call(digest, object, 1L, -1L, 14L, 0L, 0L)
  })
}

# Convert a matrix, data.frame, or data.table into a 'map' and 'data'
# elements. Basically, convert 'x' to a data.table object, use all columns as
# the keys, identify the unique rows of the data.table and map each row of
# 'x' to these unique rows.
#
# NOTE: The following should be TRUE when 'x' is matrix:
#       identical(.expand.SparseAssays.sample(.sparsify(x)), x)
# NOTE: Returned object is stripped of dimnames
.sparsify <- function(x, data_class = c("matrix", "data.frame", "data.table")) {

  # Convert input to data.table
  if (is(x, "data.frame")) {
    # Modifiy by reference
    data.table::setDT(x)
  } else if (is(x, "matrix")) {
    x <- data.table::as.data.table(x, keep.rownames = FALSE)
  } else if (is(x, "data.table")) {
    # Nothing to do
  } else {
    stop("'x' must be a 'matrix', 'data.frame', or 'data.table' object")
  }

  data_class <- match.arg(data_class)

  # Add an index for the original row number
  if (any(c(".myI", ".myMap")  %in% colnames(x))) {
    stop("'x' must not have a column named '.myI' or '.myMap'")
  }
  x[, .myI := .I]

  # Set the key (kind of like hashing the rows of the data.table since we use all columns)
  my_key <- grep(".myI", colnames(x), value = TRUE, invert = TRUE)
  data.table::setkeyv(x, cols = my_key)

  # Create the map and data
  x[, .myMap := .GRP, by = key(x)]
  map <- data.table::setkey(x[, .(.myI, .myMap)], .myI)[, .myMap]
  data <- unique(x)[, c(".myI", ".myMap") := NULL]
  if (identical(data_class, "matrix")) {
    data <- unname(as.matrix(data))
  } else if (identical(data_class, "data.frame")) {
    data <- unname(as.data.frame(data))
  }

  # Return the result
  SimpleList(map = map,
             data = data)
}

# Convert a matrix into 'map' and 'data' elements.
# Basically hash each row of the matrix, check for duplicates amongst the
# hashes, and, if there are any, find the matches using the hashes.
# NOTE: The following should be TRUE:
# identical(.expand.SparseAssays.sample(.sparsify(m)), m)
# WARNING: This is **much** slower than .sparsify(), especially as nrow(m)
#          grows. It will ultimately be removed from the package.
.sparsify_old <- function(m, data_class = class(m)) {

  if (!is.matrix(m) && !is.data.frame(m)) {
    stop("'m' must be a matrix or data.frame.")
  }

  if (data_class != "matrix" && data_class != "data.frame") {
    stop("'data_class' must be 'matrix' or 'data.frame'")
  }

  hash <- .rowHash(m)

  if (base::anyDuplicated.default(hash)) {
    # map <- selfmatch(hash)
    map <- match(hash, unique(hash))
    data <- m[!base::duplicated(hash), , drop = FALSE]
    # Remove row.names (for when m is a data.frame)
    row.names(data) <- NULL
    if (!identical(class(data), data_class)) {
      if (data_class == "data.frame") {
        # NOTE: as(data, "data.frame") doesn't work
        data <- as.data.frame.matrix(data)
      } else {
      data <- as(data, data_class)
      }
    }
  } else {
    map <- seq_len(nrow(m))
    data <- m
  }
  SimpleList(map = map, data = data)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs
###

# NOTE: rev,Assays-method doesn't work! Nor does rev,SparseAssays-method.
# NOTE: x[] errors if x is an Assays object. Should be a no-op (I think). Also
# errors if x is a SparseAssays object.
