### =========================================================================
### SSE helper functions
### -------------------------------------------------------------------------
###
### These are helper functions for methods defined for
### SparseSummarizedExperiment and RangedSparseSummarizedExperiment objects,
### collectively abbrevaited as SSE objects.
###
### None of these functions are exported but their functionality is made
### available via the appropriate S4 method.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

#' @importFrom methods is
#' @keywords internal
.valid.SSE.sparseAssays_class <- function(x) {

  if (!is(x@sparseAssays, "SparseAssays")) {
    return("'sparseAssays' slot must contain a 'SparseAssays' object.")
  }
  NULL
}

#' @keywords internal
.valid.SSE.sparseAssays_nrow <- function(x) {

  if (length(x@sparseAssays) == 0L) {
    return(NULL)
  }

  if (nrow(x@sparseAssays) != length(x)) {
    return("'sparseAssays' nrow differs from 'mcols' nrow")
  }
  NULL
}

#' @keywords internal
.valid.SSE.sparseAssays_ncol <- function(x) {
  if (length(x@sparseAssays) == 0L) {
    return(NULL)
  }

  if (ncol(x@sparseAssays) != nrow(colData(x))) {
    return("'sparseAssays' ncol differs from 'colData' nrow")
  }
  NULL
}

#' @keywords internal
.valid.SSE.sparseAssays_dim <- function(x) {
  c(.valid.SSE.sparseAssays_nrow(x),
    .valid.SSE.sparseAssays_ncol(x))
}

#' @keywords internal
.valid.SSE <- function(x) {
  c(.valid.SSE.sparseAssays_class(x),
    .valid.SSE.sparseAssays_dim(x))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

#' @keywords internal
#'
#' @importClassesFrom GenomicRanges ShallowSimpleListAssays
#' @importFrom methods as is
.SSE.to.SE <- function(from) {

  extra_assays <- as(sparseAssays(from), "ShallowSimpleListAssays")
  assays <- Assays(c(assays(from),
                     as(extra_assays, "SimpleList", strict = FALSE)))
  if (is(from, "RangedSparseSummarizedExperiment")) {
    from <- as(from, "RangedSummarizedExperiment")
  } else {
    from <- as(from, "SummarizedExperiment0")
  }
  BiocGenerics:::replaceSlots(from,
                              assays = assays)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

#' @importFrom S4Vectors SimpleList
#' @importMethodsFrom S4Vectors endoapply
#' @keywords internal
.sparseAssays.SSE <- function(x, ..., withDimnames = TRUE, expand = FALSE) {
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

#' @keywords internal
.sparseAssaysReplace.SSE <- function(x, ..., withDimnames = TRUE, value) {

  # NOTE: withDimnames arg allows
  # names(sparseAssays(se, withDimnames = FALSE)) <- value

  ok <- vapply(value, function(sa, x_dimnames) {
    # TODO: Replace with dimnames() when there is a
    # dimnames,SparseAssay[[1]]-method (i.e., one that acts on an element
    # of a SparseAssays object).
    sa_dimnames <- list(names(sa[[1]][["map"]]),
                        names(sa))
    (is.null(sa_dimnames[[1L]]) ||
      identical(sa_dimnames[[1L]], x_dimnames[[1L]]) &&
      (is.null(sa_dimnames[[2L]]) ||
         identical(sa_dimnames[[2L]], x_dimnames[[2L]])))
  }, logical(1L), x_dimnames = dimnames(x))

  if (!all(ok)) {
    stop("current and replacement 'dimnames' differ")
  }
  # NOTE: .SummarizedExperiment.assays.replace uses check = FALSE due to
  #       some unusual behaviour by packages that depend on the
  #       SummarizedExperiment package.
  x <- BiocGenerics:::replaceSlots(x, sparseAssays = value, check = TRUE)

  x
}

## convenience for common use case

#' @importClassesFrom GenomicRanges ShallowSimpleListAssays
#' @importFrom methods as
#' @importFrom stats setNames
#' @keywords internal
.sparseAssay.SSE.missing <- function(x, ..., withDimnames = TRUE, expand = FALSE) {
  # Don't want to expand all the sparseAssays, just the one being
  # extracted, so don't expand just yet.
  sparse_assays <- sparseAssays(x, ..., withDimnames = withDimnames,
                                expand = FALSE)
  if (length(sparse_assays) == 0L)
    stop("'sparseAssay(<", class(x), ">, i=\"missing\", ...) ",
         "length(sparseAssays(<", class(x), ">)) is 0'")
  val <- sparse_assays[[1]]

  if (expand) {
    val <- setNames(SparseAssays(SimpleList(val)),
                    sparseAssayNames(x)[1])
    val <- as(val, "ShallowSimpleListAssays")[[1]]
    if (withDimnames) {
      dimnames(val) <- dimnames(x)
    }
  }
  val
}

#' @importClassesFrom GenomicRanges ShallowSimpleListAssays
#' @importFrom methods as
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @keywords internal
.sparseAssay.SSE.numeric <- function(x, i, ..., withDimnames = TRUE, expand = FALSE) {
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
    val <- setNames(SparseAssays(SimpleList(val)),
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

#' @importClassesFrom GenomicRanges ShallowSimpleListAssays
#' @importFrom methods as
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @keywords internal
.sparseAssay.SSE.character <- function(x, i, ..., withDimnames = TRUE, expand = FALSE) {

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
    val <- setNames(SparseAssays(SimpleList(val)),
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

#' @keywords internal
.sparseAssayReplace.SSE.missing <- function(x, i, ..., value) {

  if (length(sparseAssays(x, withDimnames = FALSE)) == 0L) {
    stop("'sparseAssay(<", class(x), ">) <- value' ", "length(sparseAssays(<",
         class(x), ">)) is 0")
  }
  sparseAssays(x)[[1]] <- value
  x
}

#' @keywords internal
.sparseAssayReplace.SSE.numeric <- function(x, i, ..., value) {

  sparseAssays(x, ...)[[i]] <- value
  x
}

#' @keywords internal
.sparseAssayReplace.SSE.character <- function(x, i, ..., value) {

  sparseAssays(x, ...)[[i]] <- value
  x
}

#' @keywords internal
.sparseAssayNames.SSE <- function(x, ...) {
  names(sparseAssays(x, withDimnames = FALSE))
}

#' @keywords internal
.sparseAssayNamesReplace.SSE <- function(x, ..., value) {
  names(sparseAssays(x, withDimnames = FALSE)) <- value
  x
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

#' @importFrom methods is setMethod
#' @keywords internal
.subsetSingleBracket.SSE <- function(x, i, j, ..., drop = TRUE) {

  if (length(drop) != 1L || (!missing(drop) && drop)) {
    warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
  }

  # Nothing to do if both i and j are missing
  if (missing(i) && missing(j)) {
    return(x)
  }

  # Subset the sparseAssays slot
  # NOTE: Don't use the sparseAssays() accessor since can modify
  #       the returned object under its default settings (e.g.,
  #       withDimnames = TRUE).
  if (!missing(i) && !missing(j)) {
    ans_sparseAssays <- x@sparseAssays[i, j, drop = FALSE]
  } else if (!missing(i)) {
    ans_sparseAssays <- x@sparseAssays[i, , drop = FALSE]
  } else if (!missing(j)) {
    ans_sparseAssays <- x@sparseAssays[, j, drop = FALSE]
  }

  # NOTE: Can't use callNextMethod() because I'm using a .local function and
  #       not directly inside a method definition.
  if (is(x, "RangedSparseSummarizedExperiment")) {
    as_class <- "RangedSummarizedExperiment"
  } else {
    as_class <- "SummarizedExperiment0"
  }
  # NOTE: drop is ignored by both `[`,SummarizedExperiment0,ANY,ANY-method,
  #       so no need to pass it down.
  if (!missing(i) && missing(j)) {
    ans_se <- as(x, as_class)[i, j]
  } else if (!missing(i)) {
    ans_se <- as(x, as_class)[i, ]
  } else if (!missing(j)) {
    ans_se <- as(x, as_class)[, j]
  }

  # Replace slots
  # NOTE: No need to replace the metadata slot since it isn't subset by
  #       "[".
  if (is(x, "RangedSparseSummarizedExperiment")) {
    BiocGenerics:::replaceSlots(x, ...,
                                sparseAssays = ans_sparseAssays,
                                elementMetadata = ans_se@elementMetadata,
                                rowRanges = ans_se@rowRanges,
                                colData = ans_se@colData,
                                assays = ans_se@assays,
                                check = FALSE)
  } else {
    BiocGenerics:::replaceSlots(x, ...,
                                sparseAssays = ans_sparseAssays,
                                elementMetadata = ans_se@elementMetadata,
                                NAMES = ans_se@NAMES,
                                colData = ans_se@colData,
                                assays = ans_se@assays,
                                check = FALSE)
  }
}

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'                                         SummarizedExperiment0
#' @importFrom methods as is setReplaceMethod
#' @keywords internal
.replaceSingleBracket.SSE <- function(x, i, j, ..., value) {

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
  # NOTE: Can't use callNextMethod() because I'm using a .local function and
  #       not directly inside a method definition.
  if (is(x, "RangedSparseSummarizedExperiment")) {
    as_class <- "RangedSummarizedExperiment"
  } else {
    as_class <- "SummarizedExperiment0"
  }

  if (!missing(i) && !missing(j)) {
    ans_se <- as(x, as_class)
    ans_se[i, j] <- as(value, as_class)
  } else if (!missing(i)) {
    ans_se <- as(x, as_class)
    ans_se[i, ] <- as(value, as_class)
  } else if (!missing(j)) {
    ans_se <- as(x, as_class)
    ans_se[, j] <- as(value, as_class)
  }

  # Replace slots
  if (is(x, "RangedSparseSummarizedExperiment")) {
    val <- BiocGenerics:::replaceSlots(x, ...,
                                       sparseAssays = ans_sparseAssays,
                                       elementMetadata = ans_se@elementMetadata,
                                       rowRanges = ans_se@rowRanges,
                                       colData = ans_se@colData,
                                       assays = ans_se@assays,
                                       metadata = ans_se@metadata,
                                       check = FALSE)
  } else {
    val <- BiocGenerics:::replaceSlots(x, ...,
                                       sparseAssays = ans_sparseAssays,
                                       elementMetadata = ans_se@elementMetadata,
                                       NAMES = ans_se@NAMES,
                                       colData = ans_se@colData,
                                       assays = ans_se@assays,
                                       metadata = ans_se@metadata,
                                       check = FALSE)
  }
  msg <- .valid.SSE.sparseAssays_dim(val)

  if (!is.null(msg)) {
    msg
  }
  val
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display.
###

# NOTE: Based on show,SummarizedExperiment0-method
#' @importFrom methods is setMethod
#' @importMethodsFrom S4Vectors mcols metadata
#' @importMethodsFrom SummarizedExperiment assayNames assays colData
#' @keywords internal
.show.SSE <- function(object) {
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
    expt <- character(length(metadata(object)))
  }
  scat("metadata(%d): %s\n", expt)

  # sparseAssays()
  nms <- sparseAssayNames(object)
  if (is.null(nms)) {
    nms <- character(length(sparseAssays(object,
                                         withDimnames = FALSE)))
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
###

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'                                         SummarizedExperiment0
#' @keywords internal
.rbind.SSE <- function(args) {

  # rbind sparseAssays slot
  sparseAssays <- do.call(rbind, lapply(args, sparseAssays))

  # rbind the rest of the object.
  # NOTE: Can't use callNextMethod() because I'm using a .local function and
  #       not directly inside a method definition.
  if (is(args[[1L]], "RangedSparseSummarizedExperiment")) {
    as_class <- "RangedSummarizedExperiment"
  } else {
    as_class <- "SummarizedExperiment0"
  }
  se <- do.call(rbind, lapply(args, as, as_class))

  if (is(args[[1L]], "RangedSparseSummarizedExperiment")) {
    BiocGenerics:::replaceSlots(args[[1L]],
                                sparseAssays = sparseAssays,
                                elementMetadata = se@elementMetadata,
                                rowRanges = se@rowRanges,
                                colData = se@colData,
                                assays = se@assays,
                                metadata = se@metadata,
                                check = FALSE)
  } else {
    BiocGenerics:::replaceSlots(args[[1L]],
                                sparseAssays = sparseAssays,
                                elementMetadata = se@elementMetadata,
                                NAMES = se@NAMES,
                                colData = se@colData,
                                assays = se@assays,
                                metadata = se@metadata,
                                check = FALSE)
  }
}

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @keywords internal
.cbind.SSE <- function(args) {

  # cbind sparseAssays slot
  # NOTE: cbind,sparseAssays-method isn't strictly a cbind (see
  #       R/SparseAssays-class.R)
  sparseAssays <- do.call(cbind, lapply(args, sparseAssays))

  # cbind the rest of the object.
  # NOTE: Can't use callNextMethod() because I'm using a .local function and
  #       not directly inside a method definition.
  if (is(args[[1L]], "RangedSparseSummarizedExperiment")) {
    as_class <- "RangedSummarizedExperiment"
  } else {
    as_class <- "SummarizedExperiment0"
  }
  se <- do.call(cbind, lapply(args, as, as_class))

  if (is(args[[1L]], "RangedSparseSummarizedExperiment")) {
    BiocGenerics:::replaceSlots(args[[1L]],
                                sparseAssays = sparseAssays,
                                elementMetadata = se@elementMetadata,
                                rowRanges = se@rowRanges,
                                colData = se@colData,
                                assays = se@assays,
                                metadata = se@metadata,
                                check = FALSE)
  } else {
    BiocGenerics:::replaceSlots(args[[1L]],
                                sparseAssays = sparseAssays,
                                elementMetadata = se@elementMetadata,
                                NAMES = se@NAMES,
                                colData = se@colData,
                                assays = se@assays,
                                metadata = se@metadata,
                                check = FALSE)
  }
}

# TODO: There's quite a bit of room for optimising this, e.g., there's a lot of
#       coercion and validity checking that likely adds a fair bit of overhead.
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'                                         SummarizedExperiment0
#' @importMethodsFrom IRanges findOverlaps
#' @importFrom methods as is setMethod
#' @importMethodsFrom S4Vectors endoapply subjectHits
#' @importMethodsFrom SummarizedExperiment rowRanges
#' @keywords internal
.combine.SSE <- function(x, y, ...) {
  if (any(dim(y) == 0L)) {
    return(x)
  } else if (any(dim(x) == 0L)) {
    return(y)
  }

  # Update the part of the object that are derived from
  # SummarizedExperiment0/RangedSummarizedExperiment.
  if (is(x, "RangedSparseSummarizedExperiment")) {
    se <- combine(as(x, "RangedSummarizedExperiment"),
                  as(y, "RangedSummarizedExperiment"))
  } else {
    se <- combine(as(x, "SummarizedExperiment0"),
                  as(y, "SummarizedExperiment0"))
  }

  # Update the sparseAssays slot
  x_sa <- sparseAssays(x, withDimnames = TRUE)
  y_sa <- sparseAssays(y, withDimnames = TRUE)
  if (is(x, "RangedSparseSummarizedExperiment")) {
    x_ol <- findOverlaps(rowRanges(x), rowRanges(se),
                         type = "equal", minoverlap = 0L)
    y_ol <- findOverlaps(rowRanges(y), rowRanges(se),
                         type = "equal", minoverlap = 0L)
    # A kludge to update the "rownames" of the sparseAssays objects
    # so that they are combined using the findOverlaps()-based
    # rownames.
    x_sa <- endoapply(x_sa, function(sparse_assay) {
      endoapply(sparse_assay, function(sample) {
        names(sample[["map"]]) <- subjectHits(x_ol)
        sample
      })
    })
    y_sa <- endoapply(y_sa, function(sparse_assay) {
      endoapply(sparse_assay, function(sample) {
        names(sample[["map"]]) <- subjectHits(y_ol)
        sample
      })
    })
  }
  sparseAssays <- combine(x_sa, y_sa)

  # Construct the combined SSE
  if (is(x, "RangedSparseSummarizedExperiment")) {
    BiocGenerics:::replaceSlots(x,
                                sparseAssays = sparseAssays,
                                rowRanges = se@rowRanges,
                                colData = se@colData,
                                assays = se@assays,
                                NAMES = se@NAMES,
                                elementMetadata = se@elementMetadata,
                                metadata = se@metadata)
  } else {
    BiocGenerics:::replaceSlots(x,
                                sparseAssays = sparseAssays,
                                colData = se@colData,
                                assays = se@assays,
                                NAMES = se@NAMES,
                                elementMetadata = se@elementMetadata,
                                metadata = se@metadata)
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs
###

# TODO: The `assay<-()` replacement methods for SummarizedExperiment0 don't
#       set withDimnames = FALSE when checking length of assays, which
#       likely slows things down somewhat since it incurs a copy.
