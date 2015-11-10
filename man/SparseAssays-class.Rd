\name{SparseAssays-class}
\docType{class}

\alias{class:SparseAssays}
\alias{SparseAssays-class}
\alias{SparseAssays}
% methods that defer to SimpleList
\alias{length,SparseAssays-method}
\alias{names,SparseAssays-method}
\alias{names<-,SparseAssays-method}
\alias{[[,SparseAssays,ANY,ANY-method}
\alias{[[<-,SparseAssays,ANY,ANY-method}
% NOTE: No explicit coercion implemented for SimpleList <-> SparseAssays
%       because SSparseAssays are already SimpleList via inheritance.
% TODO: Should I document it anyway?

% methods defined specifically for SparseAssays
\alias{NROW,SparseAssays-method}

% An implemented concrete subclass, SimpleListSparseAssays
\alias{class:SimpleListSparseAssays}
\alias{SimpleSparseListAssays-class}
\alias{SimpleSparseListAssays}
% method defined specifically for SimpleListSparseAssays
\alias{coerce,SimpleListSparseAssays,ShallowSimpleListAssays-method}
\alias{dim,SimpleListSparseAssays-method}
\alias{[,SimpleListSparseAssays-method}
\alias{[<-,SimpleListSparseAssays-method}
\alias{rbind,SimpleListSparseAssays-method}
\alias{cbind,SimpleListSparseAssays-method}
\alias{combine,SimpleListSparseAssays-method}
\alias{densify,SimpleListSparseAssays-method}
% TODO: Uncomment once implemented
% \alias{coerce,SimpleListSparseAssays, SimpleListJointSparseAssays}

% An unimplemented concrete subclass, SimpleListJointSparseAssays
\alias{class:SimpleListJointSparseAssays}
\alias{SimpleListJointSparseAssays-class}
\alias{SimpleListJointSparseAssays}
% TODO: Uncomment the below methods once implemented
% method defined specifically for SimpleListSparseAssays
%\alias{coerce,SimpleListJointSparseAssays,ShallowSimpleListAssays-method}
%\alias{dim,SimpleListJointSparseAssays-method}
%\alias{[,SimpleListJointSparseAssays-method}
%\alias{[<-,SimpleListJointSparseAssays-method}
%\alias{rbind,SimpleListJointSparseAssays-method}
%\alias{cbind,SimpleListJointSparseAssays-method}
%\alias{combine,SimpleListJointSparseAssays-method}
%\alias{densify,SimpleListJointSparseAssays-method}
%\alias{coerce,SimpleListJointSparseAssays, SimpleListSparseAssays}

\title{SparseAssays objects}

\description{
  The SparseAssays virtual class and its methods provide a formal abstraction
  of the sparseAssays slot of \link{SparseSummarizedExperiment} and
  \link{RangedSparseSummarizedExperiment} objects.

  SimpleListSparseAssays and SimpleListJointSparseAssays (\strong{not yet
  implemented}) are concrete subclasses of SparseAssays with the former being
  currently the default implementation of SparseAssays objects. Other
  implementations (e.g. disk-based, environment-based) could easily be added.

  Note that these classes are not meant to be used directly by the end-user
  and the material in this man page is aimed at package developers.
}

\details{
  SparseAssays objects have a list-like semantics with elements containing
  key and value elements.

  The SparseAssays API consists of:
  \itemize{
   \item (a) The SparseAssays() constructor function.
   \item (b) Lossless back and forth coercion from/to SimpleList. The coercion
             method from SimpleList doesn't need (and should not) validate the
             returned object.
   \item (c) length, NROW, names, names<-, [[, [[<-
   \item (d) dim, [, [<-, rbind, cbind, combine, densify
  }

  A SparseAssays concrete subclass needs to implement (b) (required) plus
  the methods in (d) (required). The methods in (c) are inherited from the
  \code{\link[S4Vectors]{SimpleList}} class. Each element of a SparseAssays
  object is referred to as a "sparse assay" (lowercase).

  IMPORTANT: Methods that return a modified SparseAssays object (a.k.a.
  endomorphisms), that is, \code{[} as well as replacement methods
  \code{names<-}, \code{[[<-}, and \code{[<-}, must respect the
  \emph{copy-on-change contract}.
  With objects that don't make use of references internally, the developer
  doesn't need to take any special action for that because it's automatically
  taken care of by R itself. However, for objects that do make use of
  references internally (e.g. environments, external pointers, pointer to a
  file on disk, etc...), the developer needs to be careful to implement
  endomorphisms with copy-on-change semantics.
  This can easily be achieved by performaing a full (deep) copy of the object
  before modifying it instead of trying to modify it in-place. Note that the
  full (deep) copy is not always necessary in order to achieve copy-on-change
  semantics: it's enough (and often preferrable for performance reasons) to
  copy only the parts of the objects that need to be modified.

  SparseAssays has currently 1 implementation formalized by concrete subclass
  SimpleListSparseAssays. There are written specs for a second formalization,
  SimpleListJointSparseAssays, although this is not yet implemented.


  3 implementations which are formalized by concrete
  subclasses SimpleListAssays, ShallowSimpleListAssays, and AssaysInEnv.
  ShallowSimpleListAssays is the default. AssaysInEnv is a \emph{broken}
  alternative to ShallowSimpleListAssays that does NOT respect the
  \emph{copy-on-change contract}. It is only provided for illustration
  purposes (see source file Assays-class.R for the details).

  The sparseAssays slot of a \link{SparseSummarizedExperiment} object contains
  an instance of SimpleListSparseAssays.
}

\author{Peter Hickey, \url{peter.hickey@gmail.com}}

\seealso{
  \itemize{
    \item \link{SparseSummarizedExperiment} objects.

    \item \link[S4Vectors]{SimpleList} objects in the \pkg{S4Vectors} package.
  }
}

\examples{

% UP TO HERE: Get old examples from docs

## ---------------------------------------------------------------------
## DIRECT MANIPULATION OF SparseAssays OBJECTS
## ---------------------------------------------------------------------
sl1 <- SimpleList(
  s1 = SimpleList(key = as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5)),
                  value = matrix(1:10, ncol = 2)),
  s2 = SimpleList(key = as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA)),
                  value = matrix(8:1, ncol = 2)))

sl2 <- SimpleList(
  s1 = SimpleList(key = as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1)),
                  value = matrix(1:2, ncol = 1)),
  s2 = SimpleList(key = as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA)),
                  value = matrix(4:3, ncol = 1)))
sa <- SparseAssays(SimpleList(sl1, sl2))
sa

as(sa, "SimpleList")

length(sa)
sa[[2]]
dim(sa)

sa2 <- sa[-4, 2]
sa2
length(sa2)
sa2[[2]]
dim(sa2)

names(sa)
names(sa) <- c("sa1", "sa2")
names(sa)
sa[["sa2"]]

rbind(sa, sa)
\dontrun{
# ERROR: cbind-ing requires unique sample names
cbind(sa, sa)
}
\dontrun{
# ERROR: missing dimnames (which can't because there is no
#        dimnames,SparseAssays-method.
combine(sa[1:7, 1], sa[3:10, 2])
}
}
