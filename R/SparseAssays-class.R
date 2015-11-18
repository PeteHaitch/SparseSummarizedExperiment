### =========================================================================
### SparseAssays objects
### -------------------------------------------------------------------------

#' @include AllGenerics.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseAssays class
###

#' SparseAssays objects
#'
#' @description The SparseAssays virtual class and its methods provide a
#' formal abstraction of the sparseAssays slot of
#' \link{SparseSummarizedExperiment} and
#' \link{RangedSparseSummarizedExperiment} objects.
#'
#' \link{SimpleListSparseAssays} and SimpleListJointSparseAssays (\strong{not
#' yet implemented}) are concrete subclasses of SparseAssays with the former
#' being currently the default implementation of SparseAssays objects. Other
#' implementations (e.g. disk-based, environment-based) could easily be added.
#'
#' Note that these classes are not meant to be used directly by the end-user
#' and the material in this man page is aimed at package developers.
#'
#' @details SparseAssays objects have a list-like semantics with elements
#' containing key and value elements.
#'
#' The SparseAssays API consists of:
#' \itemize{
#'  \item (a) The \code{SparseAssays()} constructor function.
#'  \item (b) Lossless back and forth coercion from/to
#'  \code{\link[S4Vectors]{SimpleList}}. The coercion
#'  method from \code{\link[S4Vectors]{SimpleList}} doesn't need (and should
#'  not) validate the returned object.
#'  \item (c) \code{\link{length}}, \code{\link{NROW}}, \code{\link{names}},
#'  \code{\link{names<-}}, \code{\link{[[}}, \code{\link{[[<-}}
#'  \item (d) \code{\link{dim}}, \code{\link{[}}, \code{\link{[<-}},
#'  \code{\link{rbind}}, \code{\link{cbind}}, \code{\link{combine}},
#'  \code{\link{densify}}
#' }
#'
#' A SparseAssays concrete subclass needs to implement (b) (required) plus
#' the methods in (d) (required). The methods in (c) are inherited from the
#' \code{\link[S4Vectors]{SimpleList}} class. Each element of a SparseAssays
#' object is referred to as a "sparse assay" (lowercase).
#'
#' \strong{IMPORTANT}: Methods that return a modified SparseAssays object
#' (a.k.a. endomorphisms), that is, \code{[} as well as replacement methods
#' \code{names<-}, \code{[[<-}, and \code{[<-}, must respect the
#' \emph{copy-on-change contract}.
#' With objects that don't make use of references internally, the developer
#' doesn't need to take any special action for that because it's automatically
#' taken care of by R itself. However, for objects that do make use of
#' references internally (e.g. environments, external pointers, pointer to a
#' file on disk, etc...), the developer needs to be careful to implement
#' endomorphisms with copy-on-change semantics. This can easily be achieved by
#' performaing a full (deep) copy of the object before modifying it instead of
#' trying to modify it in-place. Note that the full (deep) copy is not always
#' necessary in order to achieve copy-on-change semantics: it's enough (and
#' often preferrable for performance reasons) to copy only the parts of the
#' objects that need to be modified.
#'
#' SparseAssays has currently 1 implementation formalized by concrete subclass
#' \link{SimpleListSparseAssays}. There are written specs for a second
#' formalization,  SimpleListJointSparseAssays, although this is not yet
#' implemented.
#'
#' The sparseAssays slot of a \link{SparseSummarizedExperiment} object contains
#' an instance of \link{SimpleListSparseAssays}.
#'
#' \strong{NOTE}: SparseAssays only payoff compared to
#' \code{SummarizedExperiment::\link[SummarizedExperiment]{Assays}} when you
#' get more than one measurement per-feature, per-sample. The payoff is greater
#' when there are lots of features with the same measurement (normally within a
#' sample, although SimpleListJointSparseAssays should allow this constraint to
#' be removed) and/or lots of NAs per-sample.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @seealso
#' \itemize{
#'  \item \link{SparseSummarizedExperiment} objects.
#'  \item \link[S4Vectors]{SimpleList} objects in the \pkg{S4Vectors} package.
#' }
#'
#' @examples
#' # See ?SimpleListSparseAssays
#'
#' @importFrom methods setClass
#'
#' @export
setClass("SparseAssays")

### Validity

# TODO: Note what should should be validated by the validity methods of a
#       concrete subclass.

#' @param x A SparseAssays object.
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

#' Normalise argument passed to SparseAssays constructor.
#'
#' @param sparse_assays A SimpleList or list.
#'
#' @return A SimpleList or an error if x is not a SimpleList or list instance.
#'
#' @importClassesFrom S4Vectors SimpleList
#' @importFrom methods is
#' @importFrom S4Vectors SimpleList
.normarg.sparse_assays <- function(sparse_assays) {
  if (!is(sparse_assays, "SimpleList")) {
    if (is.list(sparse_assays)) {
      sparse_assays <- SimpleList(sparse_assays)
    } else {
      stop("'sparse_assays' must be a SimpleList or list")
    }
  }
  sparse_assays
}

#' @param sparse_assays A SimpleList or list of sparse assays.
#' @param subclass The concrete subclass to be instantiated. The default is
#'        \link{SimpleListSparseAssays}.
#'
#' @return For \code{SparseAssays()}, an instance of a concrete SparseAssays
#' subclass given by \code{subclass}.
#'
#' @rdname SparseAssays-class
#'
#' @importFrom methods as validObject
#' @importFrom S4Vectors SimpleList
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
#' @param x A SparseAssays object.
#'
#' @rdname SparseAssays-class
#'
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
#' @inheritParams length,SparseAssays-method
#'
#' @rdname SparseAssays-class
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod("NROW", "SparseAssays",
          function(x) {
            nrow(x)
          }
)

#' @importFrom methods selectMethod
.SL_get_names <- selectMethod("names", "SimpleList")
#' @inheritParams length,SparseAssays-method
#'
#' @rdname SparseAssays-class
#'
#' @importFrom methods as setMethod
#'
#' @export
setMethod("names", "SparseAssays",
          function(x) {
            sparse_assays <- as(x, "SimpleList", strict = FALSE)
            .SL_get_names(sparse_assays)
          }
)

#' @importFrom methods setMethod
.SL_set_names <- selectMethod("names<-", "SimpleList")
#' @inheritParams length,SparseAssays-method
#' @param value An object of a class specified in the S4 method signature or as
#'        outlined in \sQuote{Details}.
#'
#' @rdname SparseAssays-class
#'
#' @importFrom methods as setMethod
#'
#' @export
setReplaceMethod("names", "SparseAssays",
                 function(x, value) {
                   sparse_assays <- as(x, "SimpleList", strict = FALSE)
                   sparse_assays <- .SL_set_names(sparse_assays, value)
                   as(sparse_assays, class(x))
                 }
)

#' @inheritParams length,SparseAssays-method
#' @param i Numeric or character vector of length 1 indicating which sparse
#'        assay to select.
#' @param j Not used by \code{[[,SparseAssays,ANY,ANY-method} or
#'        \code{[[<-,SparseAssays,ANY,ANY-method}.
#' @param ... Additional arguments, for use in specific methods.
#'
#' @rdname SparseAssays-class
#'
#' @importFrom methods as setMethod
#' @importMethodsFrom S4Vectors getListElement
#'
#' @export
setMethod("[[", "SparseAssays",
          function(x, i, j, ...) {
            sparse_assays <- as(x, "SimpleList", strict = FALSE)
            getListElement(sparse_assays, i)
          }
)

# TODO: Should really inherit params from [[,SparseAssays,ANY,ANY-method, but I
#       can't get this to work.
#' @inheritParams length,SparseAssays-method
#' @inheritParams names<-,SparseAssays-method
#'
#' @rdname SparseAssays-class
#'
#' @importFrom methods as setMethod validObject
#' @importMethodsFrom S4Vectors setListElement
#'
#' @export
setReplaceMethod("[[", "SparseAssays",
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
