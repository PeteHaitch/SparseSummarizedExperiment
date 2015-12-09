### =========================================================================
### SparseSummarizedExperiment objects
### -------------------------------------------------------------------------
###

#' @include SimpleListSparseAssays-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseSummarizedExperiment class
###

#' SparseSummarizedExperiment objects
#'
#' @description The SparseSummarizedExperiment class extends the
#' \link[SummarizedExperiment]{SummarizedExperiment} class by adding the
#' \code{sparseAssays} slot, which contains a \link{SparseAssays} object.
#'
#' Note that \link[SummarizedExperiment]{SummarizedExperiment} is the parent
#' of the SparseSummarizedExperiment class which means that all methods
#' documented in \code{?}\link[SummarizedExperiment]{SummarizedExperiment}
#' also work on a SparseSummarizedExperiment. Note also that
#' SparseSummarizedExperiment is a parent of the
#' \link{RangedSparseSummarizedExperiment} class which means that all the
#' methods documented below also work on a
#' \link{RangedSparseSummarizedExperiment}. See
#' \sQuote{Implementation and Extension} for details.
#'
#' @usage
#' ## Constructor
#'
#' # See ?RangedSparseSummarizedExperiment for the constructor function.
#'
#' ## Accessors
#'
#' sparseAssayNames(x, ...)
#' sparseAssayNames(x, ...) <- value
#' sparseAssays(x, ..., withDimnames = TRUE)
#' sparseAssays(x, ..., withDimnames = TRUE) <- value
#' sparseAssay(x, i, ...)
#' sparseAssay(x, i, ...) <- value
#'
#' ## Subsetting
#'
#' \S4method{[}{SparseSummarizedExperiment}(x, i, j, ...)
#'
#' ## Combining
#'
#' \S4method{combine}{SparseSummarizedExperiment,SparseSummarizedExperiment}(x, y, ...)
#'
#' @param x,y A SparseSummarizedExperiment object.
#' @param ... For \code{sparseAssay}, \code{...} may contain
#'        \code{withDimnames}, which is forwarded to \code{sparseAssays}.
#'
#'        For \code{cbind}, \code{rbind}, and \code{combine}, \code{...}
#'        contains SparseSummarizedExperiment objects to be combined.
#'
#'        For other accessors, ignored.
#' @param i,j For \code{sparseAssay}, \code{sparseAssay<-}, \code{i} is a
#'        numeric or character scalar; see \sQuote{Accessors} for additional
#'        constraints.
#'
#'        For \code{[,SparseSummarizedExperiment},
#'        \code{[,SparseSummarizedExperiment<-}, \code{i}, \code{j} are
#'        subscripts that can act to subset the rows (features) and columns
#'        (samples) of \code{x}, that is the sparse assay elements of
#'        \code{sparseAssays} and the \code{matrix} elements of \code{assays}.
#'
#'        For \code{[[,SparseSummarizedExperiment},
#'        \code{[[,SparseSummarizedExperiment<-}, \code{i}, is a scalar index
#'        (e.g., \code{character(1)} or \code{integer(1)}) into a column of
#'        \code{colData}.
#' @param withDimnames A \code{logical(1)}, indicating whether dimnames should
#'        be applied to extracted sparse assay elements. Setting
#'        \code{withDimnames = FALSE} increases the speed and memory efficient
#'        with which sparse assays are extracted. \code{withDimnames = TRUE}
#'        in the setter \code{sparseAssays<-} allows efficient complex
#'        assignments (e.g., updating names of sparse assays,
#'        \code{names(sparseAssays(x, withDimnames = FALSE)) <- ...} is more
#'        efficient that \code{names(sparseAssays(x)) <- ...}); it does not
#'        influence actual assignment of dimnames to sparse assays.
#'        \strong{NOTE}: For this particular example, it is simpler and just as
#'        efficient to use \code{sparseAssayNames(x) <- ...}.
#' @param value An object of a class specified in the S4 method signature or as
#'        outlined in \sQuote{Details}.
#'
#' @details These details assume familiarity with the
#' \link[SummarizedExperiment]{SummarizedExperiment} class; please first read
#' this linked documentation.
#'
#' The SparseSummarizedExperiment class is meant for \emph{sparse}
#' numeric data derived from a sequencing experiment. These data are stored as
#' a \link{SparseAssays} object in the \code{sparseAssays} slot of the
#' SparseSummarizedExperiment. In this instance, \emph{sparse} means data where
#' there are multiple measurements per-feature, per-sample and where
#' measurements with the same value (including missing values) are frequently
#' observed. \strong{NOTE}: SparseSummarizedExperiment objects only payoff
#' compared to \code{\link[SummarizedExperiment]{SummarizedExperiment}} when
#' this condition is satisfied.
#'
#' A SparseSummarizedExperiment object can also store non-sparse data by
#' storing these data in the \code{assays} slot, as would be done in a
#' \link[SummarizedExperiment]{SummarizedExperiment} object.
#'
#' The \emph{sparse data} are accessed by using the \code{sparseAssays}
#' funcion, described below. This returns a \link{SimpleList} object.
#'
#' For an example of where SparseSummarizedExperiment objects are useful,
#' please see the MethPat class in the \pkg{MethylationTuples} package
#' (currently GitHub-only,
#' \email{https://github.com/PeteHaitch/MethylationTuples/}).
#'
#' @section Constructor:
#' SparseSummarizedExperiment instances are constructed using the
#' \code{SparseSummarizedExperiment} function documented in
#' \code{?}\link{RangedSparseSummarizedExperiment}.
#'
#' @section Accessors:
#' All the accessors documented in
#' \code{?}\link[SummarizedExperiment]{SummarizedExperiment} are also
#' applicable to SparseSummarizedExperiment objects. In addition, when \code{x}
#' is a SparseSummarizedExperiment objects, the following accessors are
#' applicable.
#'
#' \describe{
#' # TODO: Check equivalence claim
#'  \item{\code{sparseAssays(x)}, \code{sparseAssays(x) <- value}:}{Get or set
#'    the sparse assays. Unlike \code{\link[SummarizedExperiment]{assays}(x)},
#'    \code{sparseAssays(x)} does not coerce the returned object to a
#'    \link[S4Vectors]{SimpleList} object but preserves it as the concrete
#'    \link{SparseAssays} subclass. \code{value} is a \link{SparseAssays}
#'    object with the same dimensions as \code{x} or a \link{SimpleList} object
#'    (which will be coerced to a \code{SparseAssays} object and must then have
#'    the same dimensions as \code{x}).}
#'
#' # TODO: Check equivalence claim
#' # TODO: Could/should I allow value to be a SimpleList
#'  \item{\code{sparseAssay(x, i)}, \code{sparseAssay(x, i) <- value}:}{Get or
#'    set the \code{i}th (default first) sparse assay elements. Unlike
#'    \code{\link[SummarizedExperiment]{assay}(x, i)},
#'    \code{sparseAssay(x, i)} allows vector \code{i} and preserves the
#'    returned object as the concrete  \link{SparseAssays} subclass.
#'    \code{value} must be a \link{SparseAssays} object (with the same concrete
#'    subclass) of the same dimension as \code{x}, and with dimension names
#'    \code{NULL} or consistent with those of \code{x} and \code{length} equal
#'    to \code{length(i)}.}
#'
#'  \item{\code{sparseAssayNames(x)}, \code{sparseAssayNames(x) <- value}:}{Get
#'    or set the names of \code{sparseAssays()} elements.}
#' }
#'
#' @section Subsetting:
#' Subsetting behaviour is inherited from methods defined for
#' SummarizedExperiment methods; see
#' \code{?}\link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @section Combining:
#' SparseSummarizedExperiment objects can be combined in three different ways.
#' \enumerate{
#'  \item \code{rbind} Suitable for when each object has the same samples.
#'  \item \code{cbind} Suitable for when each object has unique samples.
#'  \item \code{combine} Suitable in either case, \strong{however}, requires
#'  that \code{dimnames} are set on each object and that all objects have an
#'  identical number of sparse assays with identical names.
#' }
#'
#' \code{cbind()} and \code{rbind()} behaviour is inherited from methods
#' defined for SummarizedExperiment methods; see
#' \code{?}\link[SummarizedExperiment]{SummarizedExperiment}. The
#' \code{sparseAssays} slot is appropriately handled in a \code{cbind()} or
#' \code{rbind()}; see \code{\link{cbind,SimpleListSparseAssays-method}} and
#' \code{\link{rbind,SimpleListSparseAssays-method}} for details.
#
#  # TODO: Update if this functionality is moved to the SummarizedExperiment pkg
#' Additionally, the \pkg{SparseSummarizedExperiment} defines
#' \code{\link[BiocGenerics]{combine}} methods for both
#' \link[SummarizedExperiment]{SummarizedExperiment} and
#' SparseSummarizedExperiment objects. The \code{sparseAssays} slot is
#' appropriately handled in a \code{combine()}; see
#' \code{\link{combine,SimpleListSparseAssays,SimpleListSparseAssays-method}}
#' for details.
#'
#' @section Coercion:
#' Coercion from a SparseSummarizedExperiment (resp.
#' \link{RangedSparseSummarizedExperiment}) to a
#' \link[SummarizedExperiment]{SummarizedExperiment} (resp.
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}) can be done in one
#' of two ways. The first method uses implicit coercion, e.g., if \code{x} is a
#' SparseSummarizedExperiment object then
#' \code{as(x, "SparseSummarizedExperiment")} coerces it to a
#' \link[SummarizedExperiment]{SummarizedExperiment} \strong{but drops the
#' \code{sparseAssays} slot}. The second method uses an
#' explicit coercion to coerce the \link{SparseAssays} object in
#' \code{sparseAssays} slot into a \link{Assays} object and adds it to the
#' \code{assays} slot of the resulting object, \code{makeSEFromSSE(x)}.
#'
#' @section Implementation and Extension:
#' This section contains advanced material meant for package developers.
#'
#' The SparseSummarizedExperiment/RangedSparseSummarizedExperiment class
#' hierarchy is as follows:
#' \preformatted{
#' SummarizedExperiment
#' |-- RangedSummarizedExperiment
#' |   |-- RangedSparseSummarizedExperiment
#' |-- SparseSummarizedExperiment
#' |   |-- RangedSparseSummarizedExperiment
#' }
#'
#' That is, the \link{RangedSparseSummarizedExperiment} is a subclass of both
#' SparseSummarizedExperiment and
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}, although
#' SparseSummarizedExperiment takes precedence.
#'
#' SparseSummarizedExperiment is implemented as an S4 class, and can be
#' extended in the usual way, using
#' \code{contains = "SparseSummarizedExperiment"} in the new class definition.
#' Similarly, the RangedSparseSummarizedExperiment can be extended using
#' \code{contains = "RangedSparseSummarizedExperiment"} in the new class
#' definition.
#'
#' In addition, the representation of the \code{sparseAssays} slot of
#' SparseSummarizedExperiment is as a virtual class, \link{SparseAssays}. This
#' allows derived classes (\code{contains = "SparseAssays"}) to easily
#' implement alternative requirement for the sparse assays, e.g., backed by
#' file-based storage like NetCDF or the \pkg{ff} package, while re-using
#' the existing SparseSummarizedExperiment class without modification. See
#' \link{SparseAssays} for more information.
#'
#  # TODO: Update docs following changes to sparseAssay(), e.g., densify, or
#          addition of saapply().
#' The current \code{sparseAssays} slot is implemented as a
#' \link{SimpleListSparseAssays} object.
#'
#' It is generally advisable to work with
#' the sparse representation of the data wherever possible, but there are
#' times when the \emph{densified} version of the sparse data are required.
#' This can be achieved using the \code{densify = TRUE} argument in
#' \code{sparseAssays()} and \code{sparseAssay()}. Note, however, that it is
#' generally unadvisable to simultaneously densify all sparse assays and
#' samples; see \code{\link{densify}}.
#'
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}
#'
#' @seealso
#' \itemize{
#'  \item \link{RangedSparseSummarizedExperiment} objects.
#'  \item \link[SummarizedExperiment]{SummarizedExperiment} objects in the
#'    \pkg{SummarizedExperiment} package.
#'  \item \link{SparseAssays} and \link{SimpleListSparseAssays} objects.
#' }
#'
#' @aliases makeSEFromSSE
#'          sparseAssays
#'          sparseAssays,SparseSummarizedExperiment-method
#'          sparseAssays<-
#'          sparseAssays<-,SparseSummarizedExperiment,list-method
#'          sparseAssays<-,SparseSummarizedExperiment,SimpleList-method
#'          sparseAssays<-,SparseSummarizedExperiment,SparseAssays-method
#'          sparseAssay
#'          sparseAssay,SparseSummarizedExperiment,character-method
#'          sparseAssay,SparseSummarizedExperiment,missing-method
#'          sparseAssay,SparseSummarizedExperiment,numeric-method
#'          sparseAssay<-
#'          sparseAssay<-,SparseSummarizedExperiment,character,SimpleList-method
#'          sparseAssay<-,SparseSummarizedExperiment,missing,SimpleList-method
#'          sparseAssay<-,SparseSummarizedExperiment,numeric,SimpleList-method
#'          sparseAssayNames
#'          sparseAssayNames,SparseSummarizedExperiment-method
#'          sparseAssayNames<-
#'          sparseAssayNames<-,SparseSummarizedExperiment,character-method
#'          [,SparseSummarizedExperiment-method
#'          [,SparseSummarizedExperiment,ANY-method
#'      [<-,SparseSummarizedExperiment,ANY,ANY,SparseSummarizedExperiment-method
#'          show,SparseSummarizedExperiment-method
#'          rbind,SparseSummarizedExperiment-method
#'          cbind,SparseSummarizedExperiment-method
#'          combine,SparseSummarizedExperiment,SparseSummarizedExperiment-method
#' @examples
#' sl1 <- SimpleList(
#' s1 = SimpleList(key = as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5)),
#'                 value = matrix(1:10, ncol = 2)),
#' s2 = SimpleList(key = as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA)),
#'                 value = matrix(8:1, ncol = 2)))
#'
#' sl2 <- SimpleList(
#'   s1 = SimpleList(key = as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1)),
#'                   value = matrix(1:2, ncol = 1)),
#'   s2 = SimpleList(key = as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA)),
#'                   value = matrix(4:3, ncol = 1)))
#' sa <- SparseAssays(SimpleList(sa1 = sl1, sa2 = sl2))
#'
#' colData <- DataFrame(Genotype = c("WT", "KO"),
#'                      row.names = c("s1", "s2"))
#' sse <- SparseSummarizedExperiment(sparseAssays = sa,
#'                                   colData = colData)
#' sse
#' dim(sse)
#' dimnames(sse)
#' sparseAssay(sse)
#' # densify the first sparse assay.
#' # In general its a bad idea to use densify = TRUE, but these data are small
#' # enough not to worry.
#' # TODO: Should I use sparseAssay() or sparseAssays() in the example; check
#' #       out SummarizedExperiment examples.
#' densify(sparseAssay(sse), 1, 1:2)[[1]]
#' SAapply(sparseAssays(sse), function(x) x^2)
#' \dontrun{
#'  # Need sparsify = TRUE to use the replace method
#'  sparseAssays(sse) <- SAapply(sparseAssays(sse), function(x) x^2)
#' }
#' sparseAssays(sse) <- SAapply(sparseAssays(sse), function(x) x^2,
#'                              sparsify = TRUE)
#' densify(sparseAssays(sse), 1:2, 1:2)
#'
#' sparseAssay(sse)
#' # densify the first sparse assay
#' densify(sparseAssay(sse), 1, 1:2)[[1]]
#'
#' sse[, sse$Genotype == "WT"]
#'
#' ## cbind() combines objects with the same features of interest
#' ## but different samples:
#' sse1 <- sse
#' sse2 <- sse1[, 1]
#' colnames(sse2) <- "s3"
#' cmb1 <- cbind(sse1, sse2)
#' dim(cmb1)
#' dimnames(cmb1)
#'
#' ## rbind() combines objects with the same samples but different
#' ## features of interest:
#' sse1 <- sse
#' sse2 <- sse1[1:5, ]
#' rownames(sse2) <- letters[1:nrow(sse2)]
#' cmb2 <- rbind(sse1, sse2)
#' dim(cmb2)
#' dimnames(cmb2)
#'
#' ## combine() combines objects with potentially different features of interest
#' ## and different samples, by matching on names:
#' sse1 <- sse[1:5, ]
#' names(sse1) <- letters[1:5]
#' sse2 <- sse[3:8, 2]
#' names(sse2) <- letters[3:8]
#' cmb3 <- combine(sse1, sse2)
#' dim(cmb3)
#' dimnames(cmb3)
#'
#' @importFrom methods setClass
#'
#' @export
setClass("SparseSummarizedExperiment",
         contains = "SummarizedExperiment",
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

  sparseAssays_nrow <- nrow(x@sparseAssays)
  rowData_nrow <- length(x)
  if (sparseAssays_nrow != rowData_nrow) {
    txt <- sprintf(
      "\n  nb of rows in 'sparseAssays' (%d) must equal nb of rows in 'rowData' (%d)",
      sparseAssays_nrow, rowData_nrow)
    return(txt)
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

#' @importFrom S4Vectors setValidity2
setValidity2("SparseSummarizedExperiment", .valid.SSE)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# See R/RangedSparseSummarizedExperiment.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# NOTE: Don't define this as an explicit as() method because it will break the
#       implicit coercion of SSE to SE; this implicit/inherited coercion drops
#       the sparseAssays slot whereas makeSEFromSSE() preserves it by expanding
#       the sparse assays and adding them to the assays slot. The
#       implicit/inherited coercion of SSE to SE is currently relied upon by
#       several functions in this package (most non-user facing).

#' @param x A SparseSummarizedExperiment object.
#'
#' @keywords internal
#'
#' @importFrom methods as is
.SSE.to.SE <- function(x) {

  # NOTE: Use withDimnames = TRUE so that sample names are retrieved, but then
  #       strip from extra_assays like the SummarizedExperiment() constructor.
  extra_assays <- sparseAssays(x, withDimnames = TRUE)
  extra_assays <- endoapply(extra_assays, function(sparse_assay) {
    names(sparse_assay) <- NULL
    endoapply(sparse_assay, function(sample) {
      names(sample[["key"]]) <- NULL
      sample
    })
  })
  extra_assays <- as(extra_assays, "ShallowSimpleListAssays")
  assays <- Assays(c(assays(x, withDimnames = FALSE),
                     as(extra_assays, "SimpleList", strict = FALSE)))
  if (is(x, "RangedSparseSummarizedExperiment")) {
    x <- as(x, "RangedSummarizedExperiment")
  } else {
    x <- as(x, "SummarizedExperiment")
  }
  BiocGenerics:::replaceSlots(x,
                              assays = assays)
}

#' @export
makeSEFromSSE <- function(x) {
  .SSE.to.SE(x)
}

.from_SSE_to_RSSE <- function(from) {
  se <- as(from, "SummarizedExperiment")
  rse <- as(se, "RangedSummarizedExperiment")
  new("RangedSparseSummarizedExperiment",
      rse,
      sparseAssays = from@sparseAssays)
}

setAs("SparseSummarizedExperiment", "RangedSparseSummarizedExperiment",
      .from_SSE_to_RSSE
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

#' @importFrom methods as
#' @importFrom S4Vectors SimpleList
#' @importMethodsFrom S4Vectors endoapply
#' @keywords internal
.sparseAssays.SSE <- function(x, ..., withDimnames = TRUE) {

  sparse_assays <- x@sparseAssays

  if (withDimnames) {
    sparse_assays <- endoapply(sparse_assays, function(sparse_assay) {
      sparse_assay <- endoapply(sparse_assay, function(sample) {
        names(sample[["key"]]) <- names(x)
        sample
      })
      names(sparse_assay) <- colnames(x)
      sparse_assay
    })
  }

  sparse_assays
}

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssays", "SparseSummarizedExperiment",
          .sparseAssays.SSE
)

#' @keywords internal
.sparseAssaysReplace.SSE <- function(x, ..., withDimnames = TRUE, value) {

  # NOTE: withDimnames arg allows
  # names(sparseAssays(se, withDimnames = FALSE)) <- value

  ok <- vapply(value, function(sa, x_dimnames) {
    # TODO: Replace with dimnames() if/when there is a
    # dimnames,SparseAssay[[1]]-method (i.e., one that acts on an element
    # of a SparseAssays object).
    sa_dimnames <- list(names(sa[[1]][["key"]]),
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

# TODO: Unsure whether this should be SparseAssays or SimpleListSparseAssays
#       in signature.
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays",
                 c("SparseSummarizedExperiment", "SparseAssays"),
                 .sparseAssaysReplace.SSE
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays",
                 c("SparseSummarizedExperiment", "SimpleList"),
                 function(x, ..., withDimnames, value) {
                   value <- SparseAssays(value)
                   .sparseAssaysReplace.SSE(x,
                                            ...,
                                            withDimnames = withDimnames,
                                            value = value)
                 }
)

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssays",
                 c("SparseSummarizedExperiment", "list"),
                 function(x, ..., withDimnames, value) {
                   value <- SparseAssays(value)
                   .sparseAssaysReplace.SSE(x,
                                            ...,
                                            withDimnames = withDimnames,
                                            value = value)
                 }
)

## convenience for common use case

#' @importFrom methods as
#' @importFrom stats setNames
#' @keywords internal
.sparseAssay.SSE.missing <- function(x, i, ...) {

  # Don't want to densify all the sparseAssays, just the one being
  # extracted, so don't densify just yet.
  # if (!missing(j)) {
  #   sparse_assays <- sparseAssays(x, j, densify = FALSE, ...)
  # } else {
  sparse_assays <- sparseAssays(x, ...)
  # }
  if (length(sparse_assays) == 0L)
    stop("'sparseAssay(<", class(x), ">, i=\"missing\", ...) ",
         "length(sparseAssays(<", class(x), ">)) is 0'")
  subclass <- class(x@sparseAssays)
  # NOTE: Need strict = TRUE otherwise [,SimpleListSparseAssays-method is
  #       called instead of [,SimpleList-method
  as(as(sparse_assays, "SimpleList", strict = TRUE)[1L], subclass)
}

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "missing"),
          .sparseAssay.SSE.missing
)

#' @importFrom methods as
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @keywords internal
.sparseAssay.SSE.numeric <- function(x, i, ...) {
  tryCatch({
    subclass <- class(x@sparseAssays)
    # NOTE: Need strict = TRUE otherwise [,SimpleListSparseAssays-method is
    #       called instead of [,SimpleList-method
    as(as(sparseAssays(x, ...), "SimpleList", strict = TRUE)[i], subclass)
  }, error = function(err) {
    stop("'sparseAssay(<", class(x), ">, i=\"numeric\", ...)' ",
         "invalid subscript 'i'\n", conditionMessage(err))
  })
}

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "numeric"),
          .sparseAssay.SSE.numeric
)

#' @importFrom methods as
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @keywords internal
.sparseAssay.SSE.character <- function(x, i, ...) {

  msg <- paste0("'sparseAssay(<", class(x), ">, i=\"character\",",
                "...)' invalid subscript 'i'")
  val <- tryCatch({
    subclass <- class(x@sparseAssays)
    # NOTE: Need strict = TRUE otherwise [,SimpleListSparseAssays-method is
    #       called instead of [,SimpleList-method
    as(as(sparseAssays(x, ...), "SimpleList", strict = TRUE)[i], subclass)
  }, error = function(err) {
    stop(msg, "\n", conditionMessage(err))
  })
  # TODO: Is this strictly necessary (and does it even get called)?
  if (is.null(val)) {
    stop(msg, "\n'i' not in names(sparseAssays(<", class(x), ">))")
  }
  val
}

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssay", c("SparseSummarizedExperiment", "character"),
          .sparseAssay.SSE.character
)

#' @keywords internal
.sparseAssayReplace.SSE.missing <- function(x, i, ..., value) {

  # NOTE: Need strict = TRUE otherwise [,SimpleListSparseAssays-method is
  #       called instead of [,SimpleList-method
  if (length(sparseAssays(x, withDimnames = TRUE)) == 0L) {
    stop("'sparseAssay(<", class(x), ">) <- value' ", "length(sparseAssays(<",
         class(x), ">)) is 0")
  }
  sparseAssays(x)[[1]] <- value
  x
}

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "missing", "SimpleList"),
                 .sparseAssayReplace.SSE.missing
)

#' @keywords internal
.sparseAssayReplace.SSE.numeric <- function(x, i, ..., value) {

  sparseAssays(x, ...)[[i]] <- value
  x
}

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "numeric", "SimpleList"),
                 .sparseAssayReplace.SSE.numeric
)

#' @keywords internal
.sparseAssayReplace.SSE.character <- function(x, i, ..., value) {

  sparseAssays(x, ...)[[i]] <- value
  x
}

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssay",
                 c("SparseSummarizedExperiment", "character", "SimpleList"),
                 .sparseAssayReplace.SSE.character
)

#' @keywords internal
.sparseAssayNames.SSE <- function(x, ...) {
  names(sparseAssays(x, withDimnames = FALSE))
}

#' @importFrom methods setMethod
#'
#' @export
setMethod("sparseAssayNames", "SparseSummarizedExperiment",
          .sparseAssayNames.SSE
)

#' @keywords internal
.sparseAssayNamesReplace.SSE <- function(x, ..., value) {
  names(sparseAssays(x, withDimnames = FALSE)) <- value
  x
}

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("sparseAssayNames",
                 c("SparseSummarizedExperiment", "character"),
                 .sparseAssayNamesReplace.SSE
)

# NOTE: The cannonical location for dim, dimnames. dimnames should be checked
#       for consistency (if non-null) and stripped from sparseAssays on
#       construction, or added from assays if dimnames are NULL in
#       <SparseSummarizedExperiment> but not sparseAssays. dimnames need to be
#       added on to sparse assays when sparseAssays() or sparseAssay() are
#       invoked.
# NOTE: dimnames and dimnames<- methods are inherited from
#       RangedSummarizedExperiment.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

#' @importFrom methods as is
#' @keywords internal
.subsetSingleBracket.SSE <- function(x, i, j, ..., drop = FALSE) {

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
    as_class <- "SummarizedExperiment"
  }

  if (!missing(i) && !missing(j)) {
    ans_se <- as(x, as_class)[i, j, drop = drop]
  } else if (!missing(i)) {
    ans_se <- as(x, as_class)[i, , drop = drop]
  } else if (!missing(j)) {
    ans_se <- as(x, as_class)[, j, drop = drop]
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

#' @importFrom methods setMethod
#'
#' @export
setMethod("[", "SparseSummarizedExperiment",
          .subsetSingleBracket.SSE
)

#' @importFrom methods as is
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
    as_class <- "SummarizedExperiment"
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

#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("[",
                 c("SparseSummarizedExperiment", "ANY", "ANY",
                   "SparseSummarizedExperiment"),
                 function(x, i, j, ..., value) {
                   .replaceSingleBracket.SSE(x, i, j, ..., value = value)
                 }
)

# NOTE: extractROWS() and replaceROWS() methods inherited from
#       SummarizedExperiment objects.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quick colData access.
###

# NOTE: There methods are inherited from SummarizedExperiment objects.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display.
###

# NOTE: Based on show,SummarizedExperiment-method
#' @importFrom methods is
#' @importMethodsFrom S4Vectors mcols metadata
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

#' @importFrom methods setMethod
#'
#' @export
setMethod("show", "SparseSummarizedExperiment",
          .show.SSE
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
###

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
    as_class <- "SummarizedExperiment"
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


# NOTE: Appropriate for objects with distinct features and identical samples.
#' @importFrom methods setMethod
#'
#' @export
setMethod("rbind", "SparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .rbind.SSE(args)
          }
)

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
    as_class <- "SummarizedExperiment"
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

# NOTE: Appropriate for objects with identical features and distinct samples.
#' @importFrom methods setMethod
#'
#' @export
setMethod("cbind", "SparseSummarizedExperiment",
          function(..., deparse.level = 1) {
            args <- unname(list(...))
            .cbind.SSE(args)
          }
)

# TODO: There's quite a bit of room for optimising this, e.g., there's a lot of
#       coercion and validity checking that likely adds a fair bit of overhead.
#' @importMethodsFrom IRanges findOverlaps
#' @importFrom methods as is setMethod
#' @importMethodsFrom S4Vectors endoapply subjectHits
#' @keywords internal
.combine.SSE <- function(x, y, ...) {
  if (any(dim(y) == 0L)) {
    return(x)
  } else if (any(dim(x) == 0L)) {
    return(y)
  }

  # Update the part of the object that are derived from
  # SummarizedExperiment/RangedSummarizedExperiment.
  if (is(x, "RangedSparseSummarizedExperiment")) {
    se <- combine(as(x, "RangedSummarizedExperiment"),
                  as(y, "RangedSummarizedExperiment"))
  } else {
    se <- combine(as(x, "SummarizedExperiment"),
                  as(y, "SummarizedExperiment"))
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
        names(sample[["key"]]) <- subjectHits(x_ol)
        sample
      })
    })
    y_sa <- endoapply(y_sa, function(sparse_assay) {
      endoapply(sparse_assay, function(sample) {
        names(sample[["key"]]) <- subjectHits(y_ol)
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

#' @importFrom methods setMethod
#'
#' @export
setMethod("combine",
          c("SparseSummarizedExperiment", "SparseSummarizedExperiment"),
          .combine.SSE
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous NOTEs
###

# TODO: The `assay<-()` replacement methods for SummarizedExperiment don't
#       set withDimnames = FALSE when checking length of assays, which
#       likely slows things down somewhat since it incurs a copy.
