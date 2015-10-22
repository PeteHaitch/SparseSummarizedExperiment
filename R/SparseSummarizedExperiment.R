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
