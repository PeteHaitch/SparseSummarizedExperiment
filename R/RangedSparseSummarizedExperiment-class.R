### =========================================================================
### RangedSparseSummarizedExperiment objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangedSparseSummarizedExperiment class
###

#' RangedSparseSummarizedExperiment objects
#'
#' @rdname RangedSparseSummarizedExperiment
#'
#' @include SparseAssays-class.R
#'
#' @export
setClass("RangedSparseSummarizedExperiment",
         contains = "RangedSummarizedExperiment",
         representation = list(
           sparseAssays = "SparseAssays"
         ),
         prototype = list(
           sparseAssays = SparseAssays()
         )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.RangedSparseSummarizedExperiment.sparseAssays_class <- function(x) {

  if (!is(x@sparseAssays, "SparseAssays")) {
    return("'sparseAssays' slot must contain a 'SparseAssays' object.")
  }
  NULL
}

.valid.RangedSparseSummarizedExperiment.sparseAssays_nrow <- function(x) {

  if (length(x@sparseAssays) == 0L) {
    return(NULL)
  }

  if (nrow(x@sparseAssays) != length(x)) {
    return("'sparseAssays' nrow differs from 'mcols' nrow")
  }
  NULL
}

.valid.RangedSparseSummarizedExperiment.sparseAssays_ncol <- function(x) {
  if (length(x@sparseAssays) == 0L) {
    return(NULL)
  }

  if (ncol(x@sparseAssays) != nrow(colData(x))) {
    return("'sparseAssays' ncol differs from 'colData' nrow")
  }
  NULL
}

.valid.RangedSparseSummarizedExperiment.sparseAssays_dim <- function(x) {
  c(.valid.RangedSparseSummarizedExperiment.sparseAssays_nrow(x),
    .valid.RangedSparseSummarizedExperiment.sparseAssays_ncol(x))
}

.valid.RangedSparseSummarizedExperiment <- function(x) {
  c(.valid.RangedSparseSummarizedExperiment.sparseAssays_class(x),
    .valid.RangedSparseSummarizedExperiment.sparseAssays_dim(x))
}

setValidity2("RangedSparseSummarizedExperiment",
             .valid.RangedSparseSummarizedExperiment)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# TODO: Think about meaningful ways in which a RangedSparseSummarizedExperiment
# might be constructed. Note how SummarizedExperiment has several different
# constructors - a low-level function matched with several higher-level methods.

#' @rdname RangedSparseSummarizedExperiment
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
#' rse <- SummarizedExperiment(rowRanges = GRanges('chr1', IRanges(1:10, 2:11)),
#'                             colData = DataFrame(row.names = c("s1", "s2")))
#' rsse <- SparseSummarizedExperiment(sparseAssays = sa,
#'                                    rowRanges = rowRanges(rse),
#'                                    colData = colData(rse))
#' x <- rsse
#' y <- new("RangedSparseSummarizedExperiment")
#'
#' @export
setMethod("SparseSummarizedExperiment", "SparseAssays",
          function(sparseAssays,
                   assays = SimpleList(),
                   rowRanges = GRangesList(),
                   colData = DataFrame(),
                   metadata = list()) {

            # Get sample names from sparseAssays if no colData supplied
            if (missing(colData) && length(sparseAssays) != 0L) {
              nms <- names(sparseAssays[[1]])
              if (is.null(nms) && length(sparseAssays[[1]]) != 0L)
                stop(paste0("'SparseSummarizedExperiment' sparse assay names ",
                            "must not be NULL"))
              colData <- S4Vectors::DataFrame(row.names = nms)
            }

            # NOTE: The cannonical location for dimnames, including sample
            # names, is in the colData. Therefore, to simplify things, we strip
            # them from the sparseAssays object.
            sparseAssays <- endoapply(sparseAssays, unname)

            # Construct the SummarizedExperiment
            se <- SummarizedExperiment::SummarizedExperiment(assays,
                                                             rowRanges,
                                                             colData,
                                                             metadata)
            if (missing(rowRanges)) {
              stop("'SparseSummarizedExperiment0' class not yet implemented")
              new("SparseSummarizedExperiment0",
                  se,
                  sparseAssays = sparseAssays)
            } else {
              new("RangedSparseSummarizedExperiment",
                  se,
                  sparseAssays = sparseAssays)

            }
          }
)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("SparseSummarizedExperiment", "missing",
          function(sparseAssays, ...) {
            SparseSummarizedExperiment(SparseAssays(), ...)
          }
)

# TODO: SparseSummarizedExperiment,SimpleList-method and/or
#       SparseSummarizedExperiment,list-method? i.e. first try to coerce
#       SimpleList/list to a SparseAssays object and then call
#       SparseSummarizedExperiment,SparseAssays-method. When would this be
#       useful?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO: RangedSparseSummarizedExperiment -> RangedSummarizedExperiment

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

# NOTE: Following assays(), sparseAssays() will not strip the dimnames if
#       withDimnames = FALSE but will simply fail to add them.
# NOTE: The expand = TRUE argument returns a SimpleList, not a SparseAssays
#       object.
# NOTE: If the user wants sparseAssays as a ShallowSimpleListAssays object then
#       they should run as(sparseAssays(x), "ShallowSimpleListAssays"). The
#       returned object will not have rownames regardless of the value of
#       withDimnames. Note also that expand must be FALSE; the coercion to a
#       ShallowSimpleListAssays object automatically expands the
#       sparseAssays.
#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("sparseAssays", "RangedSparseSummarizedExperiment",
          function(x, ..., withDimnames = TRUE, expand = FALSE) {

            val <- x@sparseAssays

            if (expand) {
              val <- SimpleList(.expand(val))
              if (withDimnames) {
                val <- endoapply(val, function(e) {
                  names(e) <- colnames(x)
                  endoapply(e, "rownames<-", rownames(x))
                })
              }
            } else {
              if (withDimnames) {
                val <- endoapply(val, function(sparse_assay) {
                  sparse_assay <- endoapply(sparse_assay, function(sample) {
                    names(sample[["map"]]) <- names(x)
                    sample
                  })
                  names(sparse_assay) <- colnames(x)
                  sparse_assay
                })
              }
            }

            val
          }
)

# TODO
#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setReplaceMethod("sparseAssays", "RangedSparseSummarizedExperiment",
                 function(x, ..., value) {
                   stop("Not yet implemented")
                 }
)

## convenience for common use case

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "missing"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {
            # Don't want to expand all the sparseAssays, just the one being
            # extracted, so don't expand just yet.
            sparse_assays <- sparseAssays(x, ..., withDimnames = withDimnames,
                                          expand = FALSE)
            if (length(sparse_assays) == 0L)
              stop("'sparseAssay(<", class(x), ">, i=\"missing\", ...) ",
                   "length(sparseAssays(<", class(x), ">)) is 0'")
            val <- sparse_assays[[1]]

            if (expand) {
              val <- stats::setNames(SparseAssays(SimpleList(val)),
                                     sparseAssayNames(x)[1])
              val <- as(val, "ShallowSimpleListAssays")[[1]]
              if (withDimnames) {
                dimnames(val) <- dimnames(x)
              }
            }
            val
          }
)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "numeric"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {

            tryCatch({
              # Don't want to expand all the sparseAssays, just the one being
              # extracted, so don't expand just yet.
              val <- sparseAssays(x, ..., withDimnames = withDimnames,
                                  expand = FALSE)[[i]]
            }, error = function(err) {
              stop("'sparseAssay(<", class(x), ">, i=\"numeric\", ...)' ",
                   "invalid subscript 'i'\n", conditionMessage(err))
            })

            if (expand) {
              val <- stats::setNames(SparseAssays(SimpleList(val)),
                                     sparseAssayNames(x)[i])
              # extract first element, not i-th element, because this only has
              # length 1.
              val <- as(val, "ShallowSimpleListAssays")[[1]]
              if (withDimnames) {
                dimnames(val) <- dimnames(x)
              }
            }
            val
          }
)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("sparseAssay", c("RangedSparseSummarizedExperiment", "character"),
          function(x, i, ..., withDimnames = TRUE, expand = FALSE) {

            msg <- paste0("'sparseAssay(<", class(x), ">, i=\"character\",",
                          "...)' invalid subscript 'i'")
            val <- tryCatch({
              # Don't want to expand all the sparseAssays, just the one being
              # extracted, so don't expand just yet.
              sparseAssays(x, ..., withDimnames = withDimnames,
                           expand = FALSE)[[i]]
            }, error = function(err) {
              stop(msg, "\n", conditionMessage(err))
            })
            if (is.null(val)) {
              stop(msg, "\n'i' not in names(sparseAssays(<", class(x), ">))")
            }

            if (expand) {
              val <- stats::setNames(SparseAssays(SimpleList(val)),
                                     sparseAssayNames(x)[i])
              # extract first element, not i-th element, because this only has
              # length 1.
              val <- as(val, "ShallowSimpleListAssays")[[1]]
              if (withDimnames) {
                dimnames(val) <- dimnames(x)
              }
            }
            val
          }
)

# TODO: setReplaceMethod("sparseAssay", "RangedSummarizedExperiment").
# See assay<-,RangedSummarizedExperiment-method, which has multiple methods; do
# I need something like this?
# What are valid signatures for this method?
# Can I call out to the `[[`,SparseAssays-method?

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("sparseAssayNames", "RangedSparseSummarizedExperiment",
          function(x, ...) {
            names(sparseAssays(x))
          }
)

# TODO
#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setReplaceMethod("sparseAssayNames", "RangedSparseSummarizedExperiment",
                 function(x, ..., value) {
                   stop("Not yet implemented")
                 }
)

# NOTE: The cannonical location for dim, dimnames [here 'dimnames' almost
#       always means 'colnames' since it doesn't really make sense for a
#       SparseAssays object to have rownames (unlike an Assays object)].
#       dimnames should be checked for consistency (if non-null) and stripped
#       from sparseAssays on construction, or added from assays if dimnames
#       are NULL in <SparseSummarizedExperiment> but not sparseAssays.
#       dimnames need to be added on to sparse assays when sparseAssays() or
#       sparseAssay() are invoked.
# NOTE: dimnames and dimnames<- methods are inherited from
#       RangedSummarizedExperiment.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

# TODO: Why aren't I replacing the metadata slot (I don't think the [,SE-method
#       does this either)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("[", "RangedSparseSummarizedExperiment",
          function(x, i, j, ..., drop = TRUE) {

            if (length(drop) != 1L || (!missing(drop) && drop)) {
              warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
            }

            # Nothing to do if both i and j are missing
            if (missing(i) && missing(j)) {
              return(x)
            }

            # Subset the sparseAssays slot
            if (!missing(i) && !missing(j)) {
              ans_sparseAssays <- x@sparseAssays[i, j, drop = FALSE]
            } else if (!missing(i)) {
              ans_sparseAssays <- x@sparseAssays[i, , drop = FALSE]
            } else if (!missing(j)) {
              ans_sparseAssays <- x@sparseAssays[, j, drop = FALSE]
            }

            # Subset the rest of the object via callNextMethod()
            ans_rse <- callNextMethod()

            # Replace slots
            BiocGenerics:::replaceSlots(x, ...,
                                        sparseAssays = ans_sparseAssays,
                                        elementMetadata = ans_rse@elementMetadata,
                                        rowRanges = ans_rse@rowRanges,
                                        colData = ans_rse@colData,
                                        assays = ans_rse@assays,
                                        check = FALSE)
          }
)

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setReplaceMethod("[",
                 c("RangedSparseSummarizedExperiment", "ANY", "ANY",
                   "RangedSparseSummarizedExperiment"),
                 function(x, i, j, ..., value) {

                   # Nothing to do if both i and j are missing
                   if (missing(i) && missing(j)) {
                     return(x)
                   }

                   # Replace the sparseAssays slot
                   if (!missing(i) && !missing(j)) {
                     # NOTE: The use of local() is copied from `[<-`,SE-method.
                     ans_sparseAssays <- local({
                       sa <- x@sparseAssays
                       sa[i, j] <- value@sparseAssays
                       sa
                     })
                   } else if (!missing(i)) {
                     # NOTE: The use of local() is copied from `[<-`,SE-method.
                     ans_sparseAssays <- local({
                       sa <- x@sparseAssays
                       sa[i, ] <- value@sparseAssays
                       sa
                     })
                   } else if (!missing(j)) {
                     # NOTE: The use of local() is copied from `[<-`,SE-method.
                     ans_sparseAssays <- local({
                       sa <- x@sparseAssays
                       sa[, j] <- value@sparseAssays
                       sa
                     })
                   }

                   # Replace the rest of the object
                   # TODO: callNextMethod() doesn't work; why?
                   # ans_rse <- callNextMethod()

                   if (!missing(i) && !missing(j)) {
                     ans_rse <- as(x, "RangedSummarizedExperiment")
                     ans_rse[i, j] <- as(value, "RangedSummarizedExperiment")
                   } else if (!missing(i)) {
                     ans_rse <- as(x, "RangedSummarizedExperiment")
                     ans_rse[i, ] <- as(value, "RangedSummarizedExperiment")
                   } else if (!missing(j)) {
                     ans_rse <- as(x, "RangedSummarizedExperiment")
                     ans_rse[, j] <- as(value, "RangedSummarizedExperiment")
                   }

                   # Replace slots
                   val <- BiocGenerics:::replaceSlots(x, ...,
                                                      sparseAssays = ans_sparseAssays,
                                                      elementMetadata = ans_rse@elementMetadata,
                                                      rowRanges = ans_rse@rowRanges,
                                                      colData = ans_rse@colData,
                                                      assays = ans_rse@assays,
                                                      metadata = ans_rse@metadata,
                                                      check = FALSE)
                   msg <- .valid.RangedSparseSummarizedExperiment.sparseAssays_dim(val)

                   if (!is.null(msg)) {
                     msg
                   }
                   val
                 }
)

# NOTE: extractROWS() and replaceROWS() methods inherited from
#       SummarizedExperiment0 objects.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quick colData access.
###

# NOTE: There methods are inherited from SummarizedExperiment0 objects.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display.
###

# NOTE: Based on show,SummarizedExperiment0-method
#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("show", "RangedSparseSummarizedExperiment",
          function(object) {
            selectSome <- S4Vectors:::selectSome
            scat <- function(fmt, vals = character(), exdent = 2, ...) {
              vals <- ifelse(nzchar(vals), vals, "''")
              lbls <- paste(S4Vectors:::selectSome(vals), collapse = " ")
              txt <- sprintf(fmt, length(vals), lbls)
              cat(strwrap(txt, exdent = exdent, ...), sep = "\n")
            }

            cat("class:", class(object), "\n")
            cat("dim:", dim(object), "\n")

            # metadata()
            expt <- names(metadata(object))
            if (is.null(expt)) {
              expt <- character(length(S4Vectors::metadata(object)))
            }
            scat("metadata(%d): %s\n", expt)

            # sparseAssays()
            nms <- sparseAssayNames(object)
            if (is.null(nms)) {
              nms <- character(length(sparseAssays(object, withDimnames = FALSE)))
            }
            scat("sparseAssays(%d): %s\n", nms)

            # assays()
            nms <- assayNames(object)
            if (is.null(nms)) {
              nms <- character(length(assays(object, withDimnames = FALSE)))
            }
            scat("assays(%d): %s\n", nms)

            # rownames()
            dimnames <- dimnames(object)
            dlen <- sapply(dimnames, length)
            if (dlen[[1]]) {
              scat("rownames(%d): %s\n", dimnames[[1]])
            } else {
              scat("rownames: NULL\n")
            }

            # mcols()
            mcolnames <- names(mcols(object))
            fmt <- "metadata column names(%d): %s\n"
            if (is(object, "RangedSummarizedExperiment")) {
              fmt <- paste("rowRanges", fmt)
            }
            scat(fmt, mcolnames)

            # colnames()
            if (dlen[[2]]) {
              scat("colnames(%d): %s\n", dimnames[[2]])
            } else {
              cat("colnames: NULL\n")
            }

            # colData()
            scat("colData names(%d): %s\n", names(colData(object)))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine.
###

# NOTE: Appropriate for objects with different ranges and same samples.
#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("rbind", "RangedSparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .rbind.RangedSparseSummarizedExperiment(args)
          }
)

.rbind.RangedSparseSummarizedExperiment <- function(args) {

  # rbind sparseAssays slot
  sparseAssays <- do.call(rbind, lapply(args, sparseAssays))

  # rbind the rest of the object.
  # NOTE: Can't use callNextMethod() because I'm using a local function
  rse <- do.call(rbind, lapply(args, as, "RangedSummarizedExperiment"))

  BiocGenerics:::replaceSlots(args[[1L]],
                              sparseAssays = sparseAssays,
                              elementMetadata = rse@elementMetadata,
                              rowRanges = rse@rowRanges,
                              colData = rse@colData,
                              assays = rse@assays,
                              metadata = rse@metadata,
                              check = FALSE)
}

# NOTE: Appropriate for objects with same ranges and different samples.
#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("cbind", "RangedSparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .cbind.RangedSparseSummarizedExperiment(args)
          }
)

.cbind.RangedSparseSummarizedExperiment <- function(args) {

  # cbind sparseAssays slot
  # NOTE: cbind,sparseAssays-method isn't strictly a cbind (see
  #       R/SparseAssays-class.R)
  sparseAssays <- do.call(cbind, lapply(args, sparseAssays))

  # cbind the rest of the object.
  # NOTE: Can't use callNextMethod() because I'm using a local function
  rse <- do.call(cbind, lapply(args, as, "RangedSummarizedExperiment"))

  BiocGenerics:::replaceSlots(args[[1L]],
                              sparseAssays = sparseAssays,
                              elementMetadata = rse@elementMetadata,
                              rowRanges = rse@rowRanges,
                              colData = rse@colData,
                              assays = rse@assays,
                              metadata = rse@metadata,
                              check = FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs

# NOTE: granges,RangedSummarizedExperiment doesn't honour its contract to
#       return a GRanges object, e.g., the rowRanges slot could be a
#       GenomicTuples::GTuples object. A better definition might be
#       granges <- function(x) granges(rowRanges(x)).
