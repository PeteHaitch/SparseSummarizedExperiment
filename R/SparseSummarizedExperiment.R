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
#' @import methods
#  May need to be more selecting when importing from data.table since it
#  defines some functions that may clobber Bioconductor ones, e.g., shift()
#' @import data.table
#' @import BiocGenerics
#' @importFrom S4Vectors SimpleList setValidity2 endoapply mendoapply
#'             normalizeSingleBracketSubscript
#'             normalizeSingleBracketReplacementValue DataFrame metadata
#' @import SummarizedExperiment
#  TODO: How to just "import" the C-routines of digest?
#' @importFrom GenomicRanges GRangesList
#' @import digest
#' @importFrom stats setNames
NULL
