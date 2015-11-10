### =========================================================================
### SparseAssays objects
### -------------------------------------------------------------------------
###
### The SparseAssays class has a nested-list structure. This hierarchy is
### illustrated below for an example with two sparse assays and three samples.
###
### SimpleListSparseAssays
### ├── "sparse_assay_1"
### │   ├── "sample_1"
###     │   ├── key_1_1
###     │   ├── value_1_1
### │   ├── "sample_2"
###     │   ├── key_1_2
###     │   ├── value_1_2
### │   ├── "sample3"
###     │   ├── key_1_3
###     │   ├── value_1_3
### ├── "sparse_assay_2"
### │   ├── "sample_1"
###     │   ├── key_2_1
###     │   ├── value_2_1
### │   ├── "sample_2"
###     │   ├── key_2_2
###     │   ├── value_2_2
### │   ├── "sample_3"
###     │   ├── key_2_3
###     │   ├── value_2_3
###
### Each key is an integer vector and all key_*_* elements must have identical
### length.
### Each value element is a matrix object (an assay with a 1-dimensional
### measurement is stored as a 1-column matrix). Each value_*_* element may
### have a different number of rows but the maximum number of rows must be less
### than or equal to the length of the key_*_* elements. A row of the value_*_*
### element may be pointed to multiple times by the key_*_* element within the
### same sample and asssay.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseAssays class
###

# NOTE: Lossless back and forth coercion from/to SimpleList are automatically
#       taken care of by automatic methods defined by the methods package.
#' @include SparseAssays-class.R
#'
#' @importFrom methods setClass
#'
#' @export
setClass("SimpleListSparseAssays",
         contains = c("SparseAssays", "SimpleList")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimpleListSparseAssays <- function(x) {

  sparse_assays <- as(x, "SimpleList", strict = FALSE)

  if (!is(sparse_assays, "SimpleList")) {
    return("'sparseAssays' must be a SimpleList object")
  }
  if (length(sparse_assays) == 0L) {
    return(NULL)
  }

  # Check that sample level data has an element named 'key', an element named
  # 'value', and nothing else.
  element_names <- lapply(sparse_assays, function(sparse_assay) {
    lapply(sparse_assay, function(sample) {
      names(sample)
    })
  })
  element_names <- unlist(element_names, recursive = FALSE, use.names = FALSE)
  if (any(vapply(element_names, function(en) {
    !identical(en, c("key", "value")) && !identical(en, c("value", "key"))
  }, logical(1L)))) {
    return(paste0("All sample-level data within each sparse assay of a '",
                  class(x), "' object must have one element named 'key', an ",
                  "element named 'value' and nothing else."))
  }

  # Check all key elements are integer vectors.
  is_integer_key <- vapply(sparse_assays, function(sparse_assay) {
    vapply(sparse_assay, function(sample) {
      is(sample[["key"]], "integer")
    }, logical(1L))
  }, logical(length(sparse_assays[[1L]])))
  if (!isTRUE(all(is_integer_key))) {
    return(paste0("All 'key' elements of a '", class(x), "' object must be ",
                  "integer vectors."))
  }

  # Check all data elements are numeric matrix objects.
  # NOTE: This could probably be relaxed to allow non-numeric storage.mode but
  # will require changes to functions that call storage.mode() for where it is
  # implicitly assumed that the value is either 'integer' or 'double'.
  is_numeric_matrix_value <- vapply(sparse_assays, function(sparse_assay) {
    vapply(sparse_assay, function(sample) {
      # The check on nrow is because if the matrix has zero rows then it will
      # be coerced to a logical matrix even if it was originally numeric (e.g.,
      # following subsetting with the [-method).
      is(sample[["value"]], "matrix") &&
        (is.numeric(sample[["value"]]) || nrow(sample[["value"]]) == 0L)
    }, logical(1L))
  }, logical(length(sparse_assays[[1L]])))
  if (!isTRUE(all(is_numeric_matrix_value))) {
    return(paste0("All 'value' elements of a '", class(x), "' object must be ",
                  "numeric matrix objects."))
  }

  # Check each sparse assay has the same number of samples.
  n_samples <- vapply(sparse_assays, length, integer(1L))
  if (any(n_samples != n_samples[1])) {
    return(paste0("All sparse assays of a '", class(x), "'must have an ",
                  "identical number of samples."))
  }

  # Check sample names are identical across sparse assays.
  sample_names <- lapply(sparse_assays, function(sparse_assay) {
    names(sparse_assay)
  })
  if (any(vapply(sample_names, function(sn, sn1) {
    !identical(sn, sn1)
  }, logical(1L), sn1 = sample_names[[1L]]))) {
    return(paste0("All sparse assays of a '", class(x), "'must have ",
                  "identical sample names."))
  }

  # Check all key elements have the same length.
  key_length <- lapply(sparse_assays, function(sparse_assay) {
    lapply(sparse_assay, function(sample) {
      length(sample[["key"]])
    })
  })
  if (length(unique(unlist(key_length))) != 1L) {
    return(paste0("All 'key' elements of a '", class(x), "' object must have ",
                  "identical length."))
  }

  # Check all value elements within each sparse assay have the same number of
  # columns.
  # NOTE: More generally, if data were an n-dimensional array rather than a
  # 2-dimensional matrix, would require that all dimensions except the number
  # of rows were identical within a sparse assay.
  value_ncol <- lapply(sparse_assays, function(sparse_assay) {
    lapply(sparse_assay, function(sample) {
      ncol(sample[["value"]])
    })
  })
  if (any(vapply(value_ncol, function(sparse_assay_ncol) {
    length(unique(unlist(sparse_assay_ncol))) != 1L
  }, logical(1L)))) {
    return(paste0("All 'data' elements within each sparse assay of a '",
                  class(x), "' object must have identical ncol."))
  }

  # Check that the maximum value in each key element is less than or equal to
  # the number of rows in each corresponding value element.
  key_max <- lapply(sparse_assays, function(sparse_assay) {
    lapply(sparse_assay, function(sample) {
      # Suppress warning about taking max of an empty vector
      suppressWarnings(max(sample[["key"]], na.rm = TRUE))
    })
  })
  key_max <- unlist(key_max, use.names = FALSE)
  value_nrow <- lapply(sparse_assays, function(sparse_assay) {
    lapply(sparse_assay, function(sample) {
      nrow(sample[["value"]])
    })
  })
  value_nrow <- unlist(value_nrow, use.names = FALSE)
  if (any(key_max > value_nrow)) {
    return(paste0("Maximum value in each 'key' element must be less than or ",
                  "equal to the number of rows in each corresponding 'value' ",
                  "element of a '", class(x), "' object."))
  }

  NULL
}


#' @importFrom S4Vectors setValidity2
setValidity2("SimpleListSparseAssays", .valid.SimpleListSparseAssays)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###
### NOTE: The following are defined via inheritance to the SparseAssays-method:
###       length, NROW, names, names<-, [[, [[<-
###       The following are specifically defined for SimpleListSparseAssays
###       objects: dim, [, [<-, rbind, cbind, combine, densify

### dim

# NOTE: dim is defined by nrow = length of key, and ncol = number of samples
#' @importFrom methods setMethod
#'
#' @export
setMethod("dim", "SimpleListSparseAssays",
          function(x) {
            if (length(x) == 0L) {
              return(c(0L, 0L))
            } else {
              c(length(x[[1L]][[1L]][["key"]]),
                length(x[[1L]]))
            }
          }
)

### [

# NOTE: Subsetting a SimpleListSparseAssays object requires mapping the i to a
# new coordinate specified by the key element.
#' @importFrom S4Vectors normalizeSingleBracketSubscript SimpleList
#' @importFrom stats na.omit
#' @importMethodsFrom S4Vectors endoapply
.extract_SimpleListSparseAssays_subset <- function(x, i, j) {

  if (!missing(i) && !missing(j)) {

    # normalize i
    i <- normalizeSingleBracketSubscript(i, x, as.NSBS = FALSE)

    fun <- function(sparse_assay) {
      endoapply(sparse_assay[j], function(sample) {
        # Map i using key
        ii <- na.omit(sample[["key"]][i])
        # Extract using mapped i
        data <- sample[["value"]][ii, , drop = FALSE]
        # Sparsify the data
        sparsified <- .sparsify(data)
        # Update the key
        # Should have length(key) == length(i)
        if (!is.null(attr(ii, "na.action"))) {
          key <- rep(NA_integer_, length(i))
          key[-attr(ii, "na.action")] <- sparsified[["key"]]
        } else {
          key <- sparsified[["key"]]
        }
        stopifnot(length(key) == length(i))

        SimpleList(key = key,
                   value = sparsified[["value"]])
      })
    }
  } else if (!missing(i)) {

    # normalize i
    i <- normalizeSingleBracketSubscript(i, x, as.NSBS = FALSE)

    fun <- function(sparse_assay) {
      endoapply(sparse_assay, function(sample) {
        # Map i using key
        ii <- na.omit(sample[["key"]][i])
        # Extract using mapped i
        data <- sample[["value"]][ii, , drop = FALSE]
        # Sparsify the data
        sparsified <- .sparsify(data)
        # Update the key
        # Should have length(key) == length(i)
        if (!is.null(attr(ii, "na.action"))) {
          key <- rep(NA_integer_, length(i))
          key[-attr(ii, "na.action")] <- sparsified[["key"]]
        } else {
          key <- sparsified[["key"]]
        }
        stopifnot(length(key) == length(i))

        SimpleList(key = key,
                   value = sparsified[["value"]])
      })
    }
  } else if (!missing(j)) {
    fun <- function(sparse_assay) {
      sparse_assay[j]
    }
  }
  endoapply(x, fun)
}

#' @importFrom methods setMethod
#'
#' @export
setMethod("[", "SimpleListSparseAssays",
          function(x, i, j, ..., drop = FALSE) {
            if (drop) {
              warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
            }
            .extract_SimpleListSparseAssays_subset(x, i, j)
          }
)

### [<-

# IDEA (when !missing(i)):
# (1) Create new value element by
#     rbind(sample[["value"]], v_sample[["value"]])
# (2) Create new key element
#   (a) Add nrow(sample[["value"]]) to all elements of v_sample[["key"]]
#       to account for rbind() operation.
#   (b) sample[["key"]][i] <- v_sample[["key"]].
# (3) Subset value to only contain the required rows by
#     value[na.omit(unique(key)), , drop = FALSE]
# (4) Sparsify the value
# (5) Re-map the non-NA values of key by the "sparsified" key.
#
# NOTE: Steps 3-5 are "sparsifying" the value. The "expansion" of
# (key, value) at (2) and (5) should be identical, even though the
# individual elements may not be identical.
#
# IDEA (when missing(i)):
# (1) Simply replace the j-th sample(s) (key, value)-pair by that given in value.
#' @importFrom methods validObject
#' @importFrom S4Vectors SimpleList
#' @importFrom stats na.omit
#' @importMethodsFrom S4Vectors mendoapply
#'                              normalizeSingleBracketReplacementValue
.replace_SimpleListSparseAssays_subset <- function(x, i, j, value) {

  if (!missing(i) && !missing(j)) {
    # Sanity check j
    if (length(j) != ncol(value)) {
      stop("length(j) != ncol(value)")
    }

    fun <- function(sparse_assay, v_sparse_assay) {
      sparse_assay[j] <- mendoapply(function(sample, v_sample) {

        # (1)
        # NOTE: call this matrix 'tmp_value' so as not to cause confusion with
        #       'value' (which is a SimpleListSparseAssays object).
        tmp_value <- rbind(sample[["value"]], v_sample[["value"]])
        # (2)
        # NOTE: NAs are correctly propogated since NA + number = NA.
        vsm_updated <- v_sample[["key"]] + nrow(sample[["value"]])
        key <- sample[["key"]]
        key[i] <- vsm_updated
        # (3) and (4)
        sparsified <- .sparsify(tmp_value[na.omit(unique(key)), , drop = FALSE])
        # (5)
        new_lvls <- sparsified[["key"]]
        old_lvls <- na.omit(unique(key))
        key[!is.na(key)] <- new_lvls[match(key[!is.na(key)], old_lvls)]
        SimpleList(key = key,
                   value = sparsified[["value"]])
      }, sample = sparse_assay[j],
      v_sample = v_sparse_assay) # No need to subset v_sparse_assay by j
      sparse_assay
    }
  } else if (!missing(i)) {
    fun <- function(sparse_assay, v_sparse_assay) {
      sparse_assay <- mendoapply(function(sample, v_sample) {

        # (1)
        # NOTE: call this matrix 'tmp_value' so as not to cause confusion with
        #       'value' (which is a SimpleListSparseAssays object).
        tmp_value <- rbind(sample[["value"]], v_sample[["value"]])
        # (2)
        # NOTE: NAs are correctly propogated since NA + number = NA.
        vsm_updated <- v_sample[["key"]] + nrow(sample[["value"]])
        key <- sample[["key"]]
        key[i] <- vsm_updated
        # (3) and (4)
        sparsified <- .sparsify(tmp_value[na.omit(unique(key)), , drop = FALSE])
        # (5)
        new_lvls <- sparsified[["key"]]
        old_lvls <- na.omit(unique(key))
        key[!is.na(key)] <- new_lvls[match(key[!is.na(key)], old_lvls)]
        SimpleList(key = key,
                   value = sparsified[["value"]])
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
  value <- normalizeSingleBracketReplacementValue(value, x, i)

  # Loop over each sparse assay and do replacement
  val <- mendoapply(fun, x, value)

  # NOTE: Sanity check (shouldn't be necessary and may kill performance, but
  # until I have good unit tests in place this stays).
  validObject(val)

  val
}

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("[", "SimpleListSparseAssays",
                 function(x, i, j, ..., value) {
                   .replace_SimpleListSparseAssays_subset(x, i, j, value)
                 }
)

### rbind/cbind

#' @importFrom methods as
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
.bind_SimpleListSparseAssays <- function(lst, bind) {

  # NOTE: This is copied from SummarizedExperiment:::.bind_Assays(). I'm not
  # sure how this could be passed a zero-length list, but I keep it here until
  # I know it's safe to remove.
  if (length(lst) == 0L) {
    return(SparseAssays())
  }

  # If the list has only a single element then just return that element.
  if (length(lst) == 1L) {
    return(lst[[1L]])
  }

  lens <- sapply(lst, length)
  len1 <- lens[1L]
  if (any(lens != len1)) {
    # TODO: I'm not sure that this error message if very informative/accurate.
    stop("elements in sparse assays must have the same length")
  }
  if (len1 == 0L) {
    return(SparseAssays())
  }

  # Check that samples names are unique for cbind and are identical for rbind
  sample_names <- lapply(lst, function(e) {
    lapply(e, names)
  })
  # Don't need check the first element against itself
  sample_names_identical <- vapply(sample_names[-1L], function(sn, sn1) {
    identical(sn, sn1)
  }, FUN.VALUE = logical(1L), sn1 = sample_names[[1L]])
  if (identical(bind, cbind)) {
    if (any(sample_names_identical)) {
      stop("Sample names (if present) must be unique when calling 'cbind()' ",
           "on '", class(lst[[1L]]), "'.")
    }
  } else {
    if (!all(sample_names_identical)) {
      stop("Sample names (if present) must be identical when calling ",
           "'rbind()' on '", class(lst[[1L]]), "'")
    }
  }

  # Check all elements of lst have the same sparse assay names
  sparse_assay_names <- lapply(lst, names)
  if (any(vapply(sparse_assay_names, function(san, san1) {
    !identical(san, san1)
  }, logical(1L), san1 = sparse_assay_names[[1L]]))) {
    stop("All '", class(lst[[1L]]), "' objects must have the same sparse ",
         "assay names.")
  }
  sparse_assay_names <- sparse_assay_names[[1L]]

  if (identical(bind, rbind)) {
    # If rbind-ing, need to check that all data elements within each sparse
    # assay have the same number of columns.
    same_ncol <- lapply(sparse_assay_names, function(san) {
      l_sparse_assay <- lapply(lst, "[[", san)
      ncol <- lapply(l_sparse_assay, function(sparse_assay) {
        lapply(sparse_assay, function(sample) {
          ncol(sample[["value"]])
        })
      })
      ncol <- unlist(ncol, use.names = FALSE)
      all(ncol == ncol[1L])
    })
    same_ncol <- unlist(same_ncol, use.names = FALSE)
    if (!isTRUE(any(same_ncol))) {
      stop("Can only rbind '", class(lst[[1L]]), "' objects where the 'value' ",
           "elements within each sparse assay have the same number of columns.")
    }

    # NOTE: rbind,SimpleListSparseAssays-method uses the
    #       SparseAssays,`[<-`-method to recursively add the next
    #       SimpleListSparseAssays object to the end of the already rbind-ed
    #       SimpleListSparseAssays objects. It's not the most efficient method,
    #       but it works and avoids repeating much of the code used by
    #       SparseAssays,`[<-`-method.
    # This assumes that length(lst) > 1, which it should be given above checks
    # on length(lst).
    val <- lst[[1L]]
    for (idx in seq.int(from = 2L, to = length(lst), by = 1L)) {
      val_nrow <- nrow(val)
      i <- seq.int(from = val_nrow + 1,
                   to = val_nrow + nrow(lst[[idx]]),
                   by = 1L)
      val[i, ] <- lst[[idx]]
    }
  } else {
    # If cbind()-ing, need to check that all key elements have the same length.
    same_length <- lapply(sparse_assay_names, function(san) {
      l_sparse_assay <- lapply(lst, "[[", san)
      length <- lapply(l_sparse_assay, function(sparse_assay) {
        lapply(sparse_assay, function(sample) {
          length(sample[["key"]])
        })
      })
      length <- unlist(length, use.names = FALSE)
      all(length == length[1])
    })
    same_length <- unlist(same_length, use.names = FALSE)
    length <- unlist(length, use.names = FALSE)
    if (any(!same_length)) {
      stop("Can only cbind '", class(lst[[1L]]), "' objects where the 'key' ",
           "elements within each sparse assay have the same length.")
    }
    val <- lapply(sparse_assay_names, function(san) {
      l_sparse_assay <- lapply(lst, "[[", san)
      do.call("c", l_sparse_assay)
    })
    val <- as(setNames(SimpleList(val), sparse_assay_names), class(lst[[1L]]))
  }

  val
}

# NOTE: Can't defer to rbind,SimpleList-method because it in turn defers to
# rbind,ANY-method, which fails because a SimpleListSparseAssays object cannot
# be coerced to a vector (and even if it could, the resulting operation
# probably wouldn't make sense).
#' @importFrom methods setMethod
#'
#' @export
setMethod("rbind", "SimpleListSparseAssays",
          function(..., deparse.level = 1) {
            .bind_SimpleListSparseAssays(unname(list(...)), rbind)
          }
)

# NOTE: Can't defer to cbind,SimpleList-method because it in turn defers to
# cbind,ANY-method, which returns a matrix with elements being list objects
# (and doesn't make much sense for SimpleListSparseAssays objects).
# WARNING: Not really a cbind. It adds new elements to the 'sparse_assay'-level
#          SimpleList.
#' @importFrom methods setMethod
#'
#' @export
setMethod("cbind", "SimpleListSparseAssays",
          function(..., deparse.level = 1) {
            .bind_SimpleListSparseAssays(unname(list(...)), cbind)
          }
)

### combine

# Combine sample-level key and value elements
# TODO (longterm): Combine without an expand-then-sparsify operation.
#' @importFrom S4Vectors SimpleList
#' @importFrom stats complete.cases
.combine_sample_level.SimpleListSparseAssays <- function(x, y) {

  xe <- .densify.SimpleListSparseAssays.sample(x)
  dimnames(xe) <- list(names(x[["key"]]), paste0("V", seq_len(ncol(xe))))
  ye <- .densify.SimpleListSparseAssays.sample(y)
  dimnames(ye) <- list(names(y[["key"]]), paste0("V", seq_len(ncol(ye))))
  # NOTE: If the key is all NAs then the 'expanded' data are logical NA,
  #       which will cause problems when we try to combine this with a
  #       non-logical matrix.
  if (storage.mode(xe) == "logical") {
    storage.mode(xe) <- storage.mode(ye)
  }
  if (storage.mode(ye) == "logical") {
    storage.mode(ye) <- storage.mode(xe)
  }

  sparsified <- .sparsify(combine(xe, ye))

  key <- sparsified[["key"]]
  value <- sparsified[["value"]]
  NA_idx <- which(!complete.cases(value))
  if (length(NA_idx)) {
    # Take care of NA rows
    stopifnot(length(NA_idx) == 1L)
    # Update value element by dropping NA row
    value <- value[-NA_idx, , drop = FALSE]
    # Update key element to replace index by NA for NA rows
    # TODO (longterm): Probably more efficient ways to do this
    key[key == NA_idx] <- NA
    key[!is.na(key) & key > NA_idx] <- key[!is.na(key) & key > NA_idx] - 1L
  }

  SimpleList(key = key, value = value)
}

# NOTE: Can't defer to combine,SimpleList,SimpleList-method because
#       SparseAssays objects may have a different number of elements at the
#       second level (corresponding to different samples).
# NOTE: Requires that both x and y have the identical number of sparse assays
#       with identical names.
#' @importFrom methods setMethod
#' @importFrom S4Vectors SimpleList
#' @importMethodsFrom S4Vectors endoapply mendoapply
#'
#' @export
setMethod("combine", c("SimpleListSparseAssays", "SimpleListSparseAssays"),
          function(x, y, ...) {

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            if (is.null(dimnames(x)) || is.null(dimnames(y))) {
              stop("'x' and 'y' must have non-NULL dimnames")
            }

            mendoapply(function(x_sa, y_sa) {

              if (is.null(names(x_sa)) || is.null(names(y_sa))) {
                stop("sample names must be non-NULL when combining '", class(x),
                     "' objects")
              }

              # Identify shared samples
              shared_samples <- intersect(names(x_sa), names(y_sa))
              rownames <- unique(unlist(c(lapply(x_sa, function(sample) {
                names(sample[["key"]])}),
                lapply(y_sa, function(sample) {
                  names(sample[["key"]])})
              ), use.names = FALSE))

              # Update shared samples
              x_sa[shared_samples] <- mendoapply(
                .combine_sample_level.SimpleListSparseAssays,
                x_sa[shared_samples],
                y_sa[shared_samples])
              # Update samples unique to x
              x_u <- setdiff(names(x_sa), shared_samples)
              fake <- endoapply(x_sa[x_u], function(e) {
                missing_rownames <- setdiff(rownames, names(e[["key"]]))
                key <- rep(NA_integer_, length(missing_rownames))
                names(key) <- missing_rownames
                value <- matrix(ncol = ncol(e[["value"]]),
                                dimnames =
                                  list(NULL,
                                       paste0("V",
                                              seq_len(ncol(e[["value"]])))))
                storage.mode(value) <- storage.mode(e[["value"]])
                SimpleList(key = key, value = value)
              })
              x_sa[x_u] <- mendoapply(
                .combine_sample_level.SimpleListSparseAssays,
                x_sa[x_u],
                fake)

              # Update samples unique to y
              y_u <- setdiff(names(y_sa), shared_samples)
              fake <- endoapply(y_sa[y_u], function(e) {
                missing_rownames <- setdiff(rownames, names(e[["key"]]))
                key <- rep(NA_integer_, length(missing_rownames))
                names(key) <- missing_rownames
                value <- matrix(ncol = ncol(e[["value"]]),
                                dimnames =
                                  list(NULL,
                                       paste0("V",
                                              seq_len(ncol(e[["value"]])))))
                storage.mode(value) <- storage.mode(e[["value"]])
                SimpleList(key = key, value = value)
              })
              # NOTE: Append to x_sa
              x_sa <- c(x_sa,
                        mendoapply(.combine_sample_level.SimpleListSparseAssays,
                                   y_sa[y_u], fake))

              x_sa
            }, x, y)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# NOTE: Useful if wanting to densify/expand one sample's worth of data at a
#       time.
# TODO: Might be useful to export so that developers' can make use of it. If
#       exported, then need to come up with a better name.
.densify.SimpleListSparseAssays.sample <- function(sample) {
  sample[["value"]][sample[["key"]], , drop = FALSE]
}

# NOTE: Not a method for the time being; also cannot be called expand() because
# of IRanges::expand().
#' @importClassesFrom GenomicRanges ShallowSimpleListAssays
#' @importFrom SummarizedExperiment Assays
.densify.SimpleListSparseAssays <- function(x, ShallowSimpleListAssays = FALSE) {


  if (ShallowSimpleListAssays) {
    l <- lapply(x, function(sparse_assay) {
      # A kludge to guess whether the data are integer or numeric. If multiple
      # data storage modes are found then assume numeric.
      data_storage_mode <- lapply(sparse_assay, function(sample) {
        storage.mode(sample[["value"]])
      })
      data_storage_mode <- unlist(data_storage_mode)
      if (all(data_storage_mode == "integer")) {
        val <- array(NA_integer_,
                     dim = c(nrow(x),
                             ncol(x),
                             ncol(sparse_assay[[1L]][["value"]])),
                     dimnames = list(NULL, names(sparse_assay), NULL))
      } else {
        val <- array(NA_real_,
                     dim = c(nrow(x),
                             ncol(x),
                             ncol(sparse_assay[[1L]][["value"]])),
                     dimnames = list(NULL, names(sparse_assay), NULL))
      }

      # TODO (longterm): Investigate a Rcpp version
      # Fill val with the "expanded" data
      for (sample in seq_along(sparse_assay)) {
        val[ , sample, ] <- sparse_assay[[sample]][["value"]][
          sparse_assay[[sample]][["key"]], , drop = FALSE]
      }
      val
    })
    return(Assays(l))

  } else {
    lapply(x, function(sparse_assay) {
      lapply(sparse_assay, .densify.SimpleListSparseAssays.sample)
    })
  }
}

#' @importFrom methods setMethod
#'
#' @export
setMethod("densify", "SimpleListSparseAssays",
          .densify.SimpleListSparseAssays)

#' @importFrom methods setAs
#'
#' @export
setAs("SimpleListSparseAssays", "ShallowSimpleListAssays",
      function(from) {
        .densify.SimpleListSparseAssays(from, ShallowSimpleListAssays = TRUE)
      }
)

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
#' @importFrom data.table := .GRP .I as.data.table key setDT setkey setkeyv
#' @importFrom S4Vectors SimpleList
#'
.sparsify <- function(x, data_class = c("matrix", "data.frame", "data.table")) {

  # Convert input to data.table
  if (is(x, "data.frame")) {
    # Modifiy by reference
    setDT(x)
  } else if (is(x, "matrix")) {
    x <- as.data.table(x, keep.rownames = FALSE)
  } else if (is(x, "data.table")) {
    # Nothing to do
  } else {
    stop("'x' must be a 'matrix', 'data.frame', or 'data.table' object")
  }

  data_class <- match.arg(data_class)

  # NOTE: Will get errors if data have zero rows.
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


  if (identical(data_class, "matrix")) {
    value <- unname(as.matrix(value))
  } else if (identical(data_class, "data.frame")) {
    value <- unname(as.data.frame(value))
  }

  # Return the result
  SimpleList(key = key,
             value = value)
}
# To avoid WARNINGs about "Undefined global functions or variables" in
# R CMD check caused by the .sparsify() function.
#' @importFrom utils globalVariables
globalVariables(c(".myI", ".myKey"))

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparsify", "matrix",
          function(x, data_class = "matrix") {
            .sparsify(x = x, data_class = data_class)
          }
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparsify", "data.frame",
          function(x, data_class = "data.frame") {
            .sparsify(x = x, data_class = data_class)
          }
)

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparsify", "data.table",
          function(x, data_class = "data.table") {
            .sparsify(x = x, data_class = data_class)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TODOs
###

# TODO: show,SimpleListSparseAssays prints truncated class name.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs
###

# NOTE: rev,Assays-method doesn't work! Nor does rev,SparseAssays-method.
# NOTE: x[] errors if x is an Assays object. Should be a no-op (I think). Also
# errors if x is a SparseAssays object.
