#' An experimental implementation of a container modeled on the
#' SummarizedExperiment class for storing sparse genomic data.
#'
#' \pkg{SparseSummarizedExperiment} defines general purpose containers for
#' storing sparse genomic data. It aims to provide functionality for sparse
#' genomic data analogous to that available in the SummarizedExperiment
#' Bioconductor package.
#'
#' @docType package
#' @name SparseSummarizedExperiment-package
#'
# NOTE: For simplicity, just import the entire BiocGenerics package
#' @import BiocGenerics
#'
NULL

# TODO: What's the correct way to import the C-routines from digest (if digest
#       is still required)?

### NOTE: The SparseSummarizedExperiment class and associated methods
###       are inspired by the of the SummarizedExperiment0 class and
###       associated methods in the SummarizedExperiment package. There is also
###       a RangedSparseSummarizedExperiment class and associated methods that
###       is inspired by the RangedSummarizedExperiment class and associated
###       methods in the SummarizedExperiment package. The
###       RangedSummarizedExperiment class extends the SummarizedExperiment0
###       class by addition of the rowRanges slot. However, the
###       RangedSparseSummarizedExperiment class does not directly extend the
###       SparseSummarizedExperiment class. While these classes and methods
###       could have been written in such a way, I decided it was easier for
###       the RangedSparseSummarizedExperiment class to extend the
###       RangedSummarizedExperiment class and the SparseSummarizedExperiment
###       class to extend the SummarizedExperiment0 class. This hierarchy is
###       illustrated below:
###
###       SummarizedExperiment0
###       ├── RangedSummarizedExperiment
###       │   ├── RangedSparseSummarizedExperiment
###       ├── SparseSummarizedExperiment
###
###       The main benefit of this hierarchy is that all rowRanges-based
###       methods automatically work via inheritance for the
###       RangedSparseSummarizedExperiment objects.
###
###       The rules of S4 inheritance means that methods defined for
###       SparseSummarizedExperiment objects won't automatically be
###       inherited for RangedSparseSummarizedExperiment objects. This is
###       somewhat unfortunate from a developers' perspective since it requires
###       some duplication of effort to define methods that work for both
###       classes. I have tried to re-use basic functions as part of S4-methods
###       for these objects where possible. The user shouldn't have to worry
###       about this at all.
