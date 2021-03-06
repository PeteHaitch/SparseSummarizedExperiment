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
#' @param sparse_assays A SimpleList or list that can be used to construct a
#'        SparseAssays instance; see \sQuote{Examples}.
#' @param subclass The concrete subclass to be instantiated. The default is
#'        \link{SimpleListSparseAssays}.
#' @param x A SparseAssays object.
# TODO: Should x[[i, j]] return an error/warning since j has no effect?
# TODO: roxygen splits i,j by a newline in the .Rd file, which is ugly.
#' @param i,j For \code{[[, [[<-}, \code{i} is a numeric or character vector of
#'        length 1 indicating which sparse assay to select. \code{j} is not
#'        used.
#'
#'        For \code{densify}, \code{i} and \code{j} are numeric or character
#'        vectors indicating which sparse assays (\code{i}) and samples
#'        (\code{j}) to extract and densify. At least one of \code{i} or
#'        \code{j} must be provided; see \sQuote{Densify}.
#' @param value An object of a class specified in the S4 method signature or as
#'        outlined in \sQuote{Details}.
#' @param withRownames A \code{logical(1)}, indicating whether rownames should
#'        be applied to densified sparse assay elements. Setting
#'        \code{withRownames = FALSE} increases the speed and memory efficiency
#'        with which sparse assays are extracted. Note that colnames are always
#'        added.
#'
#'        For \code{SAaaply}, \code{withRownames} has no effect if
#'        \code{densify = FALSE}.
#' @param X A SparseAssays object.
#' @param FUN The function to be applied to each element of \code{X}: see
#' \sQuote{Applying a function to a SparseAssays object (\code{SAapply})}. In
#' the case of functions like \code{+}, \code{\%*\%}, the function name must be
#' backquoted or quoted.
#' @param densify A \code{logical(1)}, indicating whether the sparse data need
#'        to be densified prior to applying \code{FUN}.
#' @param sparsify A \code{logical(1)}, indicating whether the result should
#'        be sparsified following the application of \code{FUN}. By default,
#'        \code{sparsify = !densify}, that is, sparse data will remain sparse
#'        and densified data will remain densified.
#' @param ... Optional arguments to \code{FUN} or additional arguments, for use
#' in specific methods.
#' @param BPREDO,BPPARAM See \code{?\link[BiocParallel]{bplapply}}.
#'
#' @usage
#' ## Constructor
#'
#' SparseAssays(sparse_assays = SimpleList(), subclass)
#'
#' ## Accessors
#'
#' \S4method{length}{SparseAssays}(x)
#'
#' \S4method{NROW}{SparseAssays}(x)
#'
#' \S4method{names}{SparseAssays}(x)
#'
#' \S4method{names}{SparseAssays}(x) <- value
#'
#' \S4method{[[}{SparseAssays,ANY,ANY}(x, i, j, ...)
#'
#' \S4method{[[}{SparseAssays,ANY,ANY}(x, i, j, ...) <- value
#'
#' ## Densify a SparseAssays object
#'
#' densify(x, i, j, ..., withRownames = TRUE)
#'
#' ## Apply a function to a SparseAssays object
#'
#' SAapply(X, FUN, densify = TRUE, sparsify = !densify,
#'         withRownames = TRUE, ...,
#'         BPREDO = list(), BPPARAM = bpparam())
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
#'  \code{\link{names<-}}, \code{\link{[[}}, \code{\link{[[<-}}.
#'  \item (d) \code{\link{dim}}, \code{\link{dimnames}}, \code{\link{[}},
#'  \code{\link{[<-}}, \code{\link{rbind}}, \code{\link{cbind}},
#'  \code{\link{combine}}, \code{\link{densify}}, \code{\link{SAapply}}.
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
#' @section Dimensions:
#' The dimensions of a SparseAssays object are defined by nrow = length of
#' features (usually the length of the key), and ncol = number of samples.
#'
#' @section Subsetting:
#' Subsetting with \code{[} uses \code{i} to subset rows/features in each
#' sparse assay and \code{j} to subset samples in each sparse assay.
#' \strong{NOTE}: Use \code{[[} with \code{i} to select the \code{i}-th sparse
#' assay.
#'
#' @section Combining:
#' SparseAssays objects can be combined in three different ways.
#' \enumerate{
#'  \item \code{rbind} Suitable for when each object has the same samples.
#'  \item \code{cbind} Suitable for when each object has unique samples.
#'  \item \code{combine} Suitable in either case, \strong{however}, requires
#'  that \code{dimnames} are set on each object and that all objects have an
#'  identical number of sparse assays with identical names.
#' }
#'
#' @section Densify:
#' SparseAssays objects can be \emph{densified} (expanded) using the
#' \code{densify()} method. For each sample, the densified data for a single
#' sparse assay is returned as a matrix. Therefore, the \code{densify} generic
#' returns a \link[S4Vectors]{SimpleList} of length = \code{length(i)}, each
#' containing a \code{\link[S4Vectors]{SimpleList}} of length =
#' \code{length{j}}, each containing a \code{matrix} of the densified data for
#' that sample in that sparse assay.
#'
#' \strong{WARNING}: It is generally advisable to not simulatenously densify
#' all sparse assays in all samples since the entire point of using
#' SparseAssays is to use a more memory-efficient storage of the data.
#' Therefore, users must provide at least one of \code{i} (to select sparse
#' assays) and \code{j} (to select samples). If you \emph{really} wish to
#' simultaneously densify all sparse assays and samples, then use
#' \code{densify(x, seq_along(x), seq_len(ncol(x)))}. If \code{i} (resp.
#' \code{j}) is missing then effectively \code{i = seq_along(x)} (resp.
#' \code{j = seq_len(ncol(x))}).
#'
#' @section Coercion:
#' SparseAssays objects can be coerced into a
#' \link[SummarizedExperiment]{ShallowSimpleListAssays} object (from the
#' \pkg{SummarizedExperiment} package); this will also densify the object.
#' This can be done using \code{as(x, "ShallowSimpleListAssays")}, where
#' \code{x} is a SparseAssays object. \emph{WARNING}: The resulting
#' \link[SummarizedExperiment]{ShallowSimpleListAssays} object will typically
#' require much more memory than the equivalent SparseAssays object.
#'
#' @section Applying a function to a SparseAssays object (\code{SAapply}):
#' A common use case is to apply a function to a SparseAssays object. For
#' example, we might wish to compute the column-wise mean(s) for each sample
#' in a sparse assay. \code{SAapply} is designed to do this in an efficient
#' manner with an interface that is modelled on the \code{\link{lapply}}
#' functional in base R.
#'
#' \code{SAapply} takes a SparseAssays object (\code{X}) and
#' applies a single function (\code{FUN}) to each sample in each sparse assay.
#' It is worth emphasising that this means that the same function is applied to
#' all samples and sparse assays in \code{X} (use \code{sparseAssay()} with the
#' \code{i} argument to extract specific sparse assays).
#'
#' While it is desirable to apply \code{FUN} to the data in its sparse form,
#' this is not always possible and the data may need to be densified prior to
#' \code{FUN} being applied. The \code{SAapply} method simplifies this process
#' in two ways:
#'
#' \enumerate{
#'  \item \code{SAapply} allows the user to pass a function, \code{FUN}, that
#'    works on sparse or dense data. The \code{densify} argument specifies
#'    whether the data need to be densified prior to \code{FUN} being applied.
#'  \item If the data need to be densified, then \code{SAaaply} does this in a
#'    memory-efficient manner. For example, it will serially densify each
#'    sample in each sparse assay and apply \code{FUN} before moving onto the
#'    next sample's data (this is appropriately generalised if the user
#'    specifies a non-serial \code{\link[BiocParallel]{BiocParallelParam}}
#'    backend via the \code{BPPARAM} argument).
#' }
#'
#' Parallelisation is implemented via the \pkg{BiocParallel} package. Please
#' consult its documentation for further details on parallelisation options,
#' in particular the \code{?\link[BiocParallel]{BiocParallelParam}} help page.
#'
#' Finally, the \code{sparsify} argument determines the class of the return
#' value of \code{SAapply()}. If \code{sparsify = FALSE}, the return value is a
#' nested \link[base]{list} where the first level is the sparse assay and the
#' second level is the sample-level data as a dense \link[base]{matrix}. If
#' \code{sparsify = TRUE}, the return value is a SparseAssays object
#' with the same concrete subclass as \code{X}. By default,
#' \code{sparsify = !densify}, that is, sparse data will remain sparse and
#' densified data will remain densified. The use of \code{densify = TRUE}
#' allows the output of \code{SAapply()} to be used as the \code{value} in a
#' call to \code{sparseAssays(x) <- value}; see \sQuote{Examples}.
#'
#' \strong{NOTE}: The generic is called \code{SAapply} rather than
#' \code{saapply} to reduce the confusion/typo-rate with \code{\link{sapply}}.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @seealso
#' \itemize{
#'  \item \link{SimpleListSparseAssays} objects, the current default concrete
#'    subclass of the SparseAssays virtual class.
#'  \item \link{SparseSummarizedExperiment} objects, which use a SparseAssays
#'    object in the \code{sparseAssays} slot.
#' }
#'
#' @aliases SparseAssays
#'          SAapply
#'          densify
#'          NROW,SparseAssays-method
#'          SAapply,SimpleListSparseAssays-method
#'          [[,SparseAssays,ANY,ANY-method
#'          [[<-,SparseAssays,ANY,ANY-method
#'          densify,SimpleListSparseAssays,character,character-method
#'          densify,SimpleListSparseAssays,character,missing-method
#'          densify,SimpleListSparseAssays,character,numeric-method
#'          densify,SimpleListSparseAssays,missing,character-method
#'          densify,SimpleListSparseAssays,missing,missing-method
#'          densify,SimpleListSparseAssays,missing,numeric-method
#'          densify,SimpleListSparseAssays,numeric,character-method
#'          densify,SimpleListSparseAssays,numeric,missing-method
#'          densify,SimpleListSparseAssays,numeric,numeric-method
#'          length,SparseAssays-method
#'          names,SparseAssays-method
#'          names<-,SparseAssays-method
#'
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
  if (length(sparse_assays) == 0L) {
    return(NULL)
  }

  # NOTE: Unlike SummarizedExperiment:::.valid.Assays(), we can't check dims.
  #       because the dim() method is subclass-dependent.

  NULL
}

#' @importFrom S4Vectors setValidity2
setValidity2("SparseAssays", .valid.SparseAssays)

### Constructor

#' Normalise argument passed to SparseAssays constructor.
#'
#' If a matrix is supplied as 'sparse_assays', this function assumes a very
#' specific structure of the sparsified data. Namely, it assumes that the
#' NA-row, if present, is the first row of the value element. This should be
#' TRUE because data.table::setkey() is consistent with
#' base::sort(x, na.last = FALSE), however, we include an explict check and
#' return an error if this assumption is violated.
#'
#' @param sparse_assays A SimpleList, list, or matrix.
#'
#' @return A SimpleList or an error if x is not a SimpleList or list instance.
#'
#' @keywords internal
#'
#' @importClassesFrom S4Vectors SimpleList
#' @importFrom methods is
#' @importFrom S4Vectors SimpleList
.normarg.sparse_assays <- function(sparse_assays) {
  if (!is(sparse_assays, "SimpleList")) {
    if (is.list(sparse_assays)) {
      sparse_assays <- SimpleList(sparse_assays)
    } else if (is.matrix(sparse_assays)) {
      if (identical(sparse_assays, matrix())) {
        return(SimpleList())
      }
      sparsified <- sparsify(sparse_assays, "SimpleList")
      value <- sparsified[["value"]]
      key <- sparsified[["key"]]
      names(key) <- rownames(sparse_assays)
      colnames(value) <- colnames(sparse_assays)
      # TODO: any(is.na(x)) is slow: (1) it creates a matrix of the same
      #       dimension as x; (2) it checks each element of this matrix.
      #       There are surely faster ways to do this. One solution would be to
      #       use Rcpp (or .Call() to avoid the Rcpp dependency).
      if (isTRUE(all(is.na(value[1L, ])))) {
        # Need to drop the first row of the value, which is the NA-row, and
        # finesse the key.
        key <- key - 1L
        key[key == 0] <- NA_integer_
        value <- value[-1L, ]
      } else if (any(is.na(value))) {
          stop("Could not sparsify data. Please file an issue at ",
               "https://github.com/PeteHaitch/SparseSummarizedExperiment")
      }
      # TODO: This horrible nested structure is specific to
      #       SimpleListSparseAssays.
      sparse_assays <- SimpleList(
        SimpleList(SimpleList(key = key, value = value)))
    } else {
      stop("'sparse_assays' must be a SimpleList, list, or matrix")
    }
  }
  sparse_assays
}

#' @importFrom methods as validObject
#' @importFrom S4Vectors SimpleList
#' @export
SparseAssays <- function(sparse_assays = SimpleList(), subclass) {
  if (missing(subclass)) {
    subclass <- "SimpleListSparseAssays"
  }
  # Check whether subclass is a valid subclass of SparseAssays.
  # NOTE: Kind of convoluted, but extends() doesn't seemt to work.
  error_msg <- paste0("'subclass' error: '", subclass,
                      "' does not extend 'SparseAssays'")
  is_valid_subclass <- try(inherits(new(subclass), "SparseAssays"),
      silent = TRUE)
  if (is(is_valid_subclass, "try-error") || !isTRUE(is_valid_subclass)) {
    stop(error_msg)
  }

  # Normalise the arguments
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
#'
#' @export
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
setMethod("names", "SparseAssays",
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
setReplaceMethod("names", "SparseAssays",
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
setMethod("[[", "SparseAssays",
          function(x, i, j, ...) {
            sparse_assays <- as(x, "SimpleList", strict = FALSE)
            getListElement(sparse_assays, i)
          }
)

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

# TODO: Would it be useful to have higher-level accessors for working with data
#       from a SparseAssays object? E.g., to access keys or values in such a way
#       that works with all concrete subclasses of SparseAssays.
