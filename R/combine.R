### -------------------------------------------------------------------------
### combine
###

# TODO: Can we add additional arguments to specific combine() methods, e.g.,
#       ignore.mcols or use.mcols to combine,GRanges-method,
#       combine,GRangesList-method, combine,SummarizedExperiment0-method?
#       Related, should this argument be named ignore.mcols, like c(), or
#       use.mcols, like granges()?
# NOTE: All methods are based on matching names or row.names, as appropriate.
#       Consequently, these must be non-NULL in order to work. This follows
#       the behaviour of combine,matrix-method and combine,data.frame-method.
#       Currently, we error if rownames/names are NULL, but in certain
#       cirumstances we could construct meaningful rownames/names prior to
#       calling combine(), which would allow combine() to work in a more
#       general setting.
#' @rdname NULL
#'
#' @export
setMethod("combine", c("DataFrame", "DataFrame"),
  function(x, y, ...) {
    if (is.null(row.names(x)) || is.null(row.names(y))) {
      stop("'row.names' of 'x' and 'y' must be non-NULL")
      # TODO: Could we combine DataFrame objects by using merge()?
    }
    # Coerce to data.frame, use combine,data.frame-method, coerce back to
    # DataFrame.
    # TODO: This coercion incurs a copy, right? Should I basically copy the code
    # from combine,data.frame-method in order to avoid this copy?
    as(combine(as.data.frame(x), as.data.frame(y)), "DataFrame")
  }
)

#' @rdname NULL
#'
#' @export
setMethod("combine", c("GRanges", "GRanges"),
          function(x, y, ...) {
            if (is.null(names(x)) || is.null(names(y))) {
              stop("'names' of 'x' and 'y' must be non-NULL")
              # TODO: Could combine GRanges objects and unique-ify instead of
              #       erroring.
              unique(c(x, y, ignore.mcols = FALSE))
            }

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            shared_ranges <- intersect(names(x), names(y))
            # TODO: If ignore.mcols can never be TRUE remove this conditional
            # If ignore.mcols = TRUE then don't check equality of mcols
            # if (!ignore.mcols) {
            #   ok <- all.equal(x[shared_ranges], y[shared_ranges])
            # } else {
            #   xx <- x[shared_ranges]
            #   mcols(xx) <- NULL
            #   yy <- y[shared_ranges]
            #   mcols(yy) <- NULL
            #   ok <- all.equal(xx, yy)
            # }
            ok <- all.equal(x[shared_ranges], y[shared_ranges])

            if (!isTRUE(ok)) {
              stop("GRanges shared elements differ: ", ok)
            }

            c(x, y[setdiff(names(y), shared_ranges)], ignore.mcols = FALSE)
          }
)

# TODO: Test combine,GRangesList-method by itself and in conjunction with
#       combine,SummarizedExperiment0 method.
#' @rdname NULL
#'
#' @export
setMethod("combine", c("GRangesList", "GRangesList"),
          function(x, y, ..., ignore.mcols = FALSE) {
            if (is.null(names(x)) || is.null(names(y))) {
              stop("'names' of 'x' and 'y' must be non-NULL")
              # TODO: Could combine GRangesList objects and unique-ify instead
              #       of erroring.
              # NOTE: unique,GRangesList-method doesn't work do what I want. It
              #       checks for duplicates *within* each element of the
              #       GRangesList. What I want is to identify duplicate
              #       elements, hence this construction.
              # WARNING: This is potentially slow since it uses unique for
              #          lists (see ?unique) and a clunky set of coercions
              GRangesList(unique(as(c(x, y, ignore.mcols = FALSE), "list")))
            }

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            shared_elements <- intersect(names(x), names(y))
            # TODO: If ignore.mcols can never be TRUE remove this conditional
            # If ignore.mcols = TRUE then don't check equality of mcols
            # if (!ignore.mcols) {
            #   ok <- all.equal(x[shared_elements], y[shared_elements])
            # } else {
            #   xx <- x[shared_elements]
            #   mcols(xx) <- NULL
            #   yy <- y[shared_elements]
            #   mcols(yy) <- NULL
            #   ok <- all.equal(xx, yy)
            # }
            ok <- all.equal(x[shared_ranges], y[shared_ranges])

            if (!isTRUE(ok)) {
              stop("'GRangesList' shared elements differ: ", ok)
            }

            c(x, y[setdiff(names(y), shared_elements)], ignore.mcols = FALSE)
          }
)



# NOTE: Not a method since NAMES is just a character vector.
.combine.NAMES <- function(x, y, ...) {
  shared_names <- intersect(x, y)
  c(x, setdiff(y, shared_names))
}

#' @rdname RangedSparseSummarizedExperiment
#'
#' @export
setMethod("combine", c("SummarizedExperiment0", "SummarizedExperiment0"),
          function(x, y, ...) {

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            # Give a helpful error message if the user tries to combine a
            # RangedSummarizedExperiment object to an non-ranged
            # SummarizedExperiment.
            # NOTE: Can't simply check if object is SumamrizedExperiment0
            #       object since all RangedSummarizedExperiments are
            #       SummarizedExperiment0 objects but not all
            #       SummarizedExperiment0 objects are
            #       RangedSummarizedExperiment objects.
            if ((is(x, "RangedSummarizedExperiment") &&
                 !is(y, "RangedSummarizedExperiment")) ||
                (!is(x, "RangedSummarizedExperiment") &&
                 is(y, "RangedSummarizedExperiment"))) {
              stop("Cannot combine '", class(x), "' and '", class(y),
                   "' objects because only one of these has 'rowRanges'.")
            }

            x_dim <- dimnames(x)
            y_dim <- dimnames(y)

            if (is.null(x_dim) || is.null(y_dim) ||
                any(sapply(x_dim, is.null)) || any(sapply(y_dim, is.null))) {
              stop("'", class(x), "' must have dimnames for 'combine'")
              # TODO: Could add rownames based on match()/findOverlaps(), but
              #       what about if colnames are missing?
            }

            # Combine NAMES/rowRanges
            if (is(x, "RangedSummarizedExperiment")) {
              rowRanges <- combine(rowRanges(x), rowRanges(y))
            } else {
              NAMES <- .combine.NAMES(names(x), names(y))
            }

            # Combine colData
            colData <- combine(colData(x), colData(y))

            # Check that each object has the same assays
            if (!identical(assayNames(x), assayNames(y))) {
              stop("All '", x, "' objects must have identical 'assayNames'")
            }

            # Combine assays
            # NOTE: This applies combine() to each element of the assays slot
            #       in each object (e.g., combine,matrix-method if all elements
            #       are matrix objects).
            # TODO: Are non-matrix objects allowed as elements in an Assays
            #       object? If so, will need a combine() method for each of
            #       the allowed classes.
            # TODO: Use this 2-step process or should I define a
            #       combine,ShallowSimpleListAssays-method and use
            #       combine(slot(x, "assays"), slot(y, "assays"))? The trouble
            #       with the latter approach is that the rownames aren't added
            #       when using slot(x, "assays") whereas
            #       assays(x, withDimanmes = TRUE) does this  (albeit requiring
            #       a copy).
            assays <- mendoapply(combine, assays(x, withDimnames = TRUE),
                                 assays(y, withDimnames = TRUE))
            assays <- Assays(assays)

            # Combine elementMetadata
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: elementMetadata slot of a RangedSummarizedExperiment
              #       object must be a zero-column DataFrame at all times.
              elementMetadata <- DataFrame()
              elementMetadata@nrows <- length(rowRanges)
            } else {
              # NOTE: Using mcols() rather than slot(x, "elementMetadata") so
              #       rownames are added to the returned DataFrame.
              elementMetadata <- combine(mcols(x, use.names = TRUE),
                                         mcols(y, use.names = TRUE))
              # NOTE: Drop rownames of elementMetadata since these are given by
              # rownames(z) when z is a SummarizedExperiment0 object.
              rownames(elementMetadata) <- NULL
            }

            # Combine metadata
            # NOTE: metadata from all objects are combined into a list with no
            #       name checking.
            metadata <- c(metadata(x), metadata(y))

            # Construct the combined SE
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: No need to update NAMES slot since it must be NULL in a
              #       valid RangedSummarizedExperiment object.
              BiocGenerics:::replaceSlots(x,
                                          rowRanges = rowRanges,
                                          colData = colData,
                                          assays = assays,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            } else {
              BiocGenerics:::replaceSlots(x,
                                          colData = colData,
                                          assays = assays,
                                          NAMES = NAMES,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            }
          }
)

# TODO: Remove if/when a similar method is added to SummarizedExperiment
#       (see https://stat.ethz.ch/pipermail/bioc-devel/2015-October/008128.html)
# combineSE <- function(x, y, ..., nomatch = NA, use.mcols = FALSE) {
#
#   args <- unname(list(x, y, ...))
#
#   # Check that each object has unique colnames
#   colnames <- unlist(lapply(args, colnames))
#   if (anyDuplicated(colnames)) {
#     stop("Cannot combine ", class(args[[1]]), " objects with ", "duplicate ",
#          "'colnames'")
#   }
#
#   # Check that each object has the same assays
#   an <- lapply(args, assayNames)
#   if (any(sapply(an, function(x, y) any(is.na(match(x, y))),
#                  y = an[[1]]))) {
#     stop("All ", class(args[[1]]), " objects must have identical 'assayNames'")
#   }
#
#   # Combine rowRanges or NAMES slot
#   if (is(args[[1L]], "RangedSummarizedExperiment")) {
#     # NOTE: rowRanges(x) includes elementMetadata(x) since
#     #       elementMetadata(x) = elementMetadata(rowRanges(x)).
#     #       This is the reason for the use.mcols if-else block.
#     if (use.mcols) {
#       all_rowRanges <- do.call(c, lapply(args, rowRanges))
#     } else {
#       all_rowRanges <- do.call(c, lapply(args, function(x) {
#         rr <- rowRanges(x)
#         elementMetadata(rr) <- NULL
#         rr
#       }))
#     }
#     rowRanges <- unique(all_rowRanges)
#     nr <- length(rowRanges)
#   } else {
#     if (any(vapply(args, function(x) is.null(x@NAMES), logical(1)))) {
#       stop("Cannot combine ", class(args[[1]]), " objects with ",
#            "'NAMES' set to NULL")
#     }
#     all_NAMES <- do.call(c, lapply(args, function(x) x@NAMES))
#     NAMES <- unique(all_NAMES)
#     nr <- length(NAMES)
#   }
#
#   # Combine colData
#   colData <- do.call(rbind, lapply(args, colData))
#
#   # Combine elementMetadata
#   if (use.mcols) {
#     if (is(args[[1L]], "RangedSummarizedExperiment")) {
#       elementMetadata <- mcols(rowRanges)
#     } else {
#       # IDEA: Create DataFrame with all_NAMES/all_rowRanges in one
#       #       column and elementMetadata in others, unique-ify, and
#       #       check that the number of unique rows equals nr.
#       # WARNING: This will be slow for large SummarizedExperiment0
#       #          objects
#       elementMetadata <- unique(
#         cbind(DataFrame(all_NAMES),
#               do.call(rbind, lapply(args, elementMetadata))))
#       # Now drop the all_NAMES column (assumes it is the first column)
#       elementMetadata <- elementMetadata[, -c(1L), drop = FALSE]
#       # Sanity check
#       if (nrow(elementMetadata) != nr) {
#         stop("'elementMetadata' must match across ",class(args[[1]]), " objects")
#       }
#     }
#   } else {
#     elementMetadata <- DataFrame()
#     elementMetadata@nrows <- nr
#   }
#
#   # Create assays of the correct dimension (fill with 'nomatch')
#   # First, create the empty combined assay using the appropriate
#   # storage.mode (guessed from the storage.mode of the assay in the
#   # first sample).
#   nomatch <- lapply(seq_along(an[[1]]), function(i) {
#     storage.mode(nomatch) <- storage.mode(assay(args[[1]],
#                                                 i,
#                                                 withDimnames = FALSE))
#     nomatch
#   })
#   assays <- lapply(nomatch, function(nm) {
#     matrix(nm, nrow = nr, ncol = length(colnames))
#   })
#   names(assays) <- an[[1]]
#
#   # NOTE: I suspect that there are faster and more efficient ways to
#   # combine the assays, perhaps at the C-level.
#   if (is(args[[1L]], "RangedSummarizedExperiment")) {
#     for (j in seq_along(args)) {
#       ol <- findOverlaps(args[[j]], rowRanges, type = "equal")
#       for (i in seq_along(assays)) {
#         assays[[i]][subjectHits(ol),
#                     match(colnames(args[[j]]), colnames)] <-
#           assay(args[[j]], i, withDimnames = FALSE)
#       }
#     }
#   } else {
#     for (j in seq_along(args)) {
#       ol <- match(args[[j]]@NAMES, NAMES)
#       for (i in seq_along(assays)) {
#         assays[[i]][ol, match(colnames(args[[j]]), colnames)] <-
#           assay(args[[j]], i, withDimnames = FALSE)
#       }
#     }
#   }
#   assays <- Assays(assays)
#
#   # Combine metadata
#   metadata <- do.call(c, lapply(args, metadata))
#
#   if (is(args[[1L]], "RangedSummarizedExperiment")) {
#     # No need to replace elementMetadata slot since it is part of
#     # rowRanges.
#     BiocGenerics:::replaceSlots(args[[1L]],
#                                 rowRanges = rowRanges,
#                                 colData = colData,
#                                 assays = assays,
#                                 metadata = metadata)
#   } else {
#     BiocGenerics:::replaceSlots(args[[1L]],
#                                 NAMES = NAMES,
#                                 colData = colData,
#                                 assays = assays,
#                                 metadata = metadata,
#                                 elementMetadata = elementMetadata)
#   }
# }


# #' Combining SummarizedExperiment objects
# #'
# #' @rdname combine
# #'
# #' @export
# setMethod("combine", c("SummarizedExperiment0", "SummarizedExperiment0"),
#           function(x, y, ..., nomatch = NA, use.mcols = FALSE) {
#             combineSE(x, y, ..., nomatch = nomatch, use.mcols = use.mcols)
#           }
# )
