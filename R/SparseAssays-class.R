### =========================================================================
### SparseAssays objects
### -------------------------------------------------------------------------
###
### The SparseAssays API consists of:
###   (a) The SparseAssays() constructor function.
###   (b) Lossless back and forth coercion from/to SimpleList. The coercion
###       method from SimpleList doesn't need (and should not) validate the
###       returned object.
###   (c) length, NROW, names, names<-, [[, [[<-
###   (d) dim, [, [<-, rbind, cbind, combine, densify
###
### A SparseAssays concrete subclass needs to implement (b) (required) plus
### the methods in (d) (required). The methods in (c) are inherited from the
### SimpleList class. Each element of a SparseAssays object is referred to as a
### "sparse assay" (lowercase).
###
### The SparseSummarizedExperiment package currently implements the
### SimpleListSparseAssays concrete subclass and includes specs for a
### SimpleListJointSparseAssays concrete subclass, although this is not yet
### implemented (TODO: Update as needed).
###
### IMPORTANT: Methods that return a modified SparseAssays object (a.k.a.
### endomorphisms), that is, [ as well as replacement methods names<-, [[<-,
### and [<-, must respect the copy-on-change contract. With objects that
### don't make use of references internally, the developer doesn't need to
### take any special action for that because it's automatically taken care of
### by R itself. However, for objects that do make use of references internally
### (e.g. environments, external pointers, pointer to a file on disk, etc...),
### the developer needs to be careful to implement endomorphisms with
### copy-on-change semantics. This can be achieved by performaing a full (deep)
### copy of the object before modifying it instead of trying to modify it
### in-place. Note that the full (deep) copy is not always necessary in order
### to achieve copy-on-change semantics: it's enough (and often preferrable for
### performance reasons) to copy only the parts of the objects that need to
### be modified.
###
###
### NOTE: SparseAssays only payoff compared to SummarizedExperiment::Assays
### when you get more than one measurement per-feature, per-sample. The payoff
### is greater if there are lots of features with the same measurement within a
### sample and/or lots of NAs per-sample.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseAssays class
###

#' @include AllGenerics.R
#'
#' @importFrom methods setClass
#'
#' @export
setClass("SparseAssays")

### Validity

# TODO: Note what should should be validated by the validity methods of a
#       concrete subclass.

#' @importFrom methods as is
.valid.SparseAssays <- function(x) {

  sparse_assays <- as(x, "SimpleList", strict = FALSE)

  if (!is(sparse_assays, "SimpleList")) {
    return("'sparseAssays' must be a SimpleList object")
  }

  NULL
}

#' @importFrom S4Vectors setValidity2
setValidity2("SparseAssays", .valid.SparseAssays)

### Constructor

#' @importClassesFrom S4Vectors SimpleList
#' @importFrom methods is
#' @importFrom S4Vectors SimpleList
.normarg.sparse_assays <- function(sparse_assays) {
  if (!is(sparse_assays, "SimpleList")) {
    if (is.list(assays)) {
      sparse_assays <- SimpleList(sparse_assays)
    } else {
      stop("'sparse_assays' must be a SimpleList or list")
    }
  }
  sparse_assays
}

#' @importFrom methods as validObject
#' @importFrom S4Vectors SimpleList
#'
#' @export
SparseAssays <- function(sparse_assays = SimpleList(), subclass) {
  if (missing(subclass)) {
    subclass <- "SimpleListSparseAssays"
  }
  # TODO: Some check that subclass is a valid concrete subclass of the
  #       virtual SparseAssays class.
  sparse_assays <- .normarg.sparse_assays(sparse_assays)

  ans <- as(sparse_assays, subclass)
  validObject(ans)
  ans
}

### Accessors

#' @importFrom methods selectMethod
.SL_get_length <- selectMethod("length", "SimpleList")
#' @importFrom methods as setMethod
#'
#' @export
setMethod("length", "SparseAssays",
          function(x) {
            sparse_assays <- as(x, "SimpleList", strict = FALSE)
            .SL_get_length(sparse_assays)
          }
)

# NOTE: NROW,SparseAssays-method would otherwise defer to NROW,Vector-method
#       which calls length(). This matters because otherwise NROW(x) != nrow(x)
#       when x is a SparseAssays object.
#' @importFrom methods setMethod
setMethod("NROW", "SparseAssays",
          function(x) {
            nrow(x)
          }
)

#' @importFrom methods selectMethod
.SL_get_names <- selectMethod("names", "SimpleList")
#' @importFrom methods as setMethod
#'
#' @export
setMethod("names", "Assays",
          function(x) {
            sparse_assays <- as(x, "SimpleList", strict = FALSE)
            .SL_get_names(sparse_assays)
          }
)

#' @importFrom methods setMethod
.SL_set_names <- selectMethod("names<-", "SimpleList")
#' @importFrom methods as setMethod
#'
#' @export
setReplaceMethod("names", "Assays",
                 function(x, value) {
                   sparse_assays <- as(x, "SimpleList", strict = FALSE)
                   sparse_assays <- .SL_set_names(sparse_assays, value)
                   as(sparse_assays, class(x))
                 }
)

#' @importFrom methods as setMethod
#' @importMethodsFrom S4Vectors getListElement
#'
#' @export
setMethod("[[", "Assays",
          function(x, i, j, ...) {
            sparse_assays <- as(x, "SimpleList", strict = FALSE)
            getListElement(sparse_assays, i)
          }
)

#' @importFrom methods as setMethod validObject
#' @importMethodsFrom S4Vectors setListElement
#'
#' @export
setReplaceMethod("[[", "Assays",
                 function(x, i, j, ..., value) {
                   sparse_assays <- as(x, "SimpleList", strict = FALSE)
                   sparse_assays <- setListElement(sparse_assays, i, value)
                   ans <- as(sparse_assays, class(x))
                   validObject(ans)
                   ans
                 }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TODOs
###

# TODO: Add some higher level accessors for working with data from a
#       SparseAssays object, e.g., to access keys, values, densified/expanded
#       data, saaply (SparseAssay apply?), etc.
