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
#' @import BiocGenerics
#' @importFrom S4Vectors SimpleList setValidity2 endoapply mendoapply
#'             normalizeSingleBracketSubscript
#'             normalizeSingleBracketReplacementValue DataFrame metadata
#' @import SummarizedExperiment
#  TODO: How to just "import" the C-routines of digest?
#' @importFrom GenomicRanges GRangesList
#' @import digest
#' @importFrom stats setNames
#' @import data.table
NULL
