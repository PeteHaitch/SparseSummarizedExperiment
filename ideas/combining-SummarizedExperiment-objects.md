# Combining SummarizedExperiment objects
Peter Hickey  
20 October 2015  



# Motivation

I often find myself with multiple `SE` objects (I'm using `SE` as a shorthand for the `SummarizedExperiment0` and `RangedSummarizedExeriment` classes), each with potentially non-distinct samples and potentially non-overlapping features/ranges. Currently, it is difficult to combine these objects; `rbind()` can only combine objects with the same samples but distinct features/ranges and `cbind()` can only combine objects with the same features/ranges but distinct samples. I think it would be useful to have a "combine" method for `SE` objects that handles the most general situation where each object has potentially non-distinct samples and potentially non-overlapping features/ranges.

__NOTE__: All `RangedSummarizedExperiment` objects are `SummarizedExperiment0` objects but the converse is not true.

# Methods

My initial attempt was posted to [https://gist.github.com/PeteHaitch/8993b096cfa7ccd08c13/2c82a279ee7d56cec3e9dadb8024cbb9183da576](https://gist.github.com/PeteHaitch/8993b096cfa7ccd08c13/2c82a279ee7d56cec3e9dadb8024cbb9183da576) and cross-posted to the Bioconductor development mailing list ([https://stat.ethz.ch/pipermail/bioc-devel/2015-October/008128.html](https://stat.ethz.ch/pipermail/bioc-devel/2015-October/008128.html)). Following feedback from Martin Morgan ([https://stat.ethz.ch/pipermail/bioc-devel/2015-October/008144.html](https://stat.ethz.ch/pipermail/bioc-devel/2015-October/008144.html)), I revised to make the following major changes:

1. Provide a `combine()` method for each component of the `SE` object. This means I wrote: `combine,DataFrame,DataFrame-method` (for the `colData` and `elementMetadata` slots), `combine,GRanges,GRanges-method` and `combine,GRangesList,GRangesList-method` (for the `rowRanges` slot), and `combine,SimpleList,SimpleList-method` (for the `assays` slot). There is also a helper function `.combine.NAMES()` to combine the `NAMES` slot of a non-ranged `SE`.
2. All `combine()` methods are initially designed to "combine" based on the `names`/`rownames`/`dimnames` of the objects. This follows the existing `combine,data.frame,data.frame-method` and `combine,matrix,matrix-method` from the _BiocGenerics_ package. 

With respect to (2), in some cases I realised there were equivalent definitions that did not rely on these "names" that were faster, and other scenarios where it might be desirable not to rely on these "names". For example, in my experience, it is common for a `RangedSummarizedExperiment` to have `NULL` names, i.e. `identical(names(x), NULL)` is `TRUE` where `x` is a `RangedSummarizedExperiment` (note that `names(x)` calls `names(rowRanges(x))` and so it is the `rowRanges(x)` that commonly have `NULL` names). A consequence of this is that we cannot combine these `RangedSummarizedExperiment` objects based on "names" alone. However, there is an obvious way to construct meaningful `rownames` for these objects by a `findOverlaps()`-based naming scheme; this is implemented in `combine,SummarizedExperiment0,SummarizedExperiment0-method` for `RangedSummarizedExperiment` objects.

Comments on all aspects of the implementation are appreciated. I have highlighted particular concerns with __RFC__.

## Setup


```r
suppressPackageStartupMessages(library(SummarizedExperiment))
```

## `combine,DataFrame,DataFrame-method`

This coerces the `DataFrame` to a `data.frame`, calls the `combine,data.frame,data.frame-method`, and then coerces back to a `DataFrame`.

__RFC__: Do these coercions incur a copy operation(s)? Can these be avoided, e.g., should I basically copy the code from `combine,data.frame-method` to define `combine,DataFrame,DataFrame-method`?


```r
setMethod("combine", c("DataFrame", "DataFrame"),
  function(x, y, ...) {
    
    if (all(dim(x) == 0L) && all(dim(y) == 0L)) {
      return(x)
    } else if (all(dim(x) == 0L)) {
      return(y)
    } else if (all(dim(y) == 0L)) {
      return(x)
    }
    
    if (is.null(row.names(x)) || is.null(row.names(y))) {
      stop("'row.names' of 'x' and 'y' must be non-NULL")
    }
    as(combine(as.data.frame(x), as.data.frame(y)), "DataFrame")
  }
)
#> [1] "combine"
```

## `combine,GRanges,GRanges-method`

There are two strategies available (let `x` and `y` be `GRanges` objects in what follows):

1. Rely on `names(x)` and `names(y)` to identify "shared" elements.
2. Ignore `names(x)` and `names(y)` and simply identify unique elements of a combined `x` with `y`.

Option 1 is more restrictive (since it requires `names(x)` and `names(y)` to be non-`NULL`) but matches the behaviour of other `combine` methods that are defined in terms of the `names`/`rownames`/`dimnames` of the objects.

Option 2 will work on all `GRanges` objects even when `names(x)` or `names(y)` are `NULL` since `unique,GRanges-method` ignores `names`. The `names` of the "shared" elements are taken from `x` in the returned object. This option is approximately 2x faster.

I prefer option 2 since it is more general and faster.

__RFC__: Which option is preferrable?

__RFC__: Could/should we add a `ignore.mcols` argument to `combine,GRanges,GRanges-method` like that available in `c,GRanges,GRanges-method`? Also, should this argument be `ignore.mcols` or `use.mcols` (like that available in `granges,GRanges-method`)?


```r
setMethod("combine", c("GRanges", "GRanges"),
          function(x, y, ...) {
            
            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }
            
            if (is.null(names(x)) || is.null(names(y))) {
              stop("'names' of 'x' and 'y' must be non-NULL")
            }

            shared_elements <- intersect(names(x), names(y))
            ok <- all.equal(x[shared_elements], y[shared_elements])
            if (!isTRUE(ok)) {
              stop("GRanges shared elements differ: ", ok)
            }

            c(x, y[setdiff(names(y), shared_elements)], ignore.mcols = FALSE)
          }
)
#> [1] "combine"
```


```r
setMethod("combine", c("GRanges", "GRanges"),
          function(x, y, ...) {
            
            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }
            
            unique(c(x, y))
          }
)
#> [1] "combine"
```

## `combine,GRangesList,GRangesList-method`

__RFC__: What does it mean to combine two `GRangesList` objects? For example, consider two `GRangesList` objects `x` and `y`; what should the result of calling `combine(x, y)` be? Since these are `List`-derived objects, it seems to make sense to `mendoapply()` `combine` to each element of `x` and `y`. This works when `x` and `y` have identical element names (case A), but not when `x` or `y` contains elements not found in the other (case B).


```r
x <- GRangesList(gr1 = GRanges("chr1", IRanges(1, 10)))
x
#> GRangesList object of length 1:
#> $gr1 
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1   [1, 10]      *
#> 
#> -------
#> seqinfo: 1 sequence from an unspecified genome; no seqlengths

y <- GRangesList(gr1 = GRanges("chr1", IRanges(12, 16)))
y
#> GRangesList object of length 1:
#> $gr1 
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1  [12, 16]      *
#> 
#> -------
#> seqinfo: 1 sequence from an unspecified genome; no seqlengths

# Case A (works)
mendoapply(combine, x, y)
#> GRangesList object of length 1:
#> $gr1 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1  [ 1, 10]      *
#>   [2]     chr1  [12, 16]      *
#> 
#> -------
#> seqinfo: 1 sequence from an unspecified genome; no seqlengths

# Case B (expect error)
y[["gr2"]] <- GRanges("chr1", IRanges(66, 99))
mendoapply(combine, x, y)
#> Error in validObject(.Object): invalid class "GRangesList" object: 1: 'x@elementMetadata' is not parallel to 'x'
#> invalid class "GRangesList" object: 2: 'names(x)' must be NULL or have the length of 'x'
```

It seems like in case B that the user might want for `y$gr2` to be included in the returned value, so perhaps the following definition is suitable:


```r
setMethod("combine", c("GRangesList", "GRangesList"),
          function(x, y, ...) {
            
            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }           
            
            if (is.null(names(x)) || is.null(names(y))) {
              stop("'names' of 'x' and 'y' must be non-NULL")
            }
            
            shared_elements <- intersect(names(x), names(y))
            x[shared_elements] <- mendoapply(combine, x[shared_elements], 
                                             y[shared_elements])
            c(x, y[setdiff(names(y), shared_elements)])
          }
)
#> [1] "combine"
```

This will appropriately combine `GRangesList` objects in case B but requires that both objects have non-`NULL` `names`:


```r
# Works
combine(x, y)
#> GRangesList object of length 2:
#> $gr1 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1  [ 1, 10]      *
#>   [2]     chr1  [12, 16]      *
#> 
#> $gr2 
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames   ranges strand
#>   [1]     chr1 [66, 99]      *
#> 
#> -------
#> seqinfo: 1 sequence from an unspecified genome; no seqlengths

# Expect error
combine(unname(x), unname(y))
#> Error in combine(unname(x), unname(y)): 'names' of 'x' and 'y' must be non-NULL
```

Unlike `combine,GRanges,GRanges-method`, it seems difficult to extend `combine,GRangesList,GRangesList-method` to work when either `names(x)` or `names(y)` is `NULL`.

__RFC__: Any suggestions on how to extend `combine,GRangesList,GRangesList-method` to work when either `names(x)` or `names(y)` is `NULL`? Is this even desirable? Without this, we won't be able to use `combine()` on `DESeqDataSet` objects, for example.

## `combine,SimpleList,SimpleList-method`

This calls `mendoapply()` with `combine()` on `x` and `y`. This means it will only work if a suitable `combine()` method is defined for `class(x)` and `class(y)`. Unlike the proposed `combine,GRangesList,GRangesList-method`, `combine,SimpleList,SimpleList-method` requires that `x` and `y` have an identical number of elements with identical names.


```r
setMethod("combine", c("SimpleList", "SimpleList"),
          function(x, y, ...) {

            if (length(y) == 0L) {
              return(x)
            } else if (length(x) == 0L) {
              return(y)
            }

            if (length(names(x)) != length(names(y))) {
              stop("'", class(x), "' objects have different number of elements:",
                   "\n\t", paste(names(x), collapse = " "), "\n\t",
                   paste(names(y), collapse = " "))
            }

            if (!all(names(x) == names(y))) {
              stop("'", class(x), "' objects have different element names:\n\t",
                   paste(names(x), collapse = " "), "\n\t",
                   paste(names(y), collapse = " "))
            }

            mendoapply(combine, x, y)
          }
)
#> [1] "combine"
```

__RFC__: Should this be `combine,List,List-method` or might this make it too general?


## `combine,SummarizedExperiment0,SummarizedExperiment0-method`

This applies the `combine` to each slot of the `SummarizedExperiment0` objects. There is a helper function, `.combine.NAMES()` to combine the `NAMES` slot of non-ranged `SummarizedExperiment0` objects; this isn't a method since the `NAMES` slot is `NULL` or a `character()` object rather than a formal S4 class.


```r
.combine.NAMES <- function(x, y, ...) {
  shared_names <- intersect(x, y)
  c(x, setdiff(y, shared_names))
}
```

There are two strategies available (let `x` and `y` be `SummarizedExperiment0` objects in what follows):

1. Rely on `dimnames(x)` and `dimnames(y)` to identify "shared" features/ranges and samples.
2. If `x` and `y` are `RangedSummarizedExperiment`-based objects, then ignore `rownames(x)` and `rownames(y)` and simply identify unique ranges of a combined `x` with `y`. This still requires that both `x` and `y` have `colnames`. Currently, this will not work for a `RangedSummarizedExperiment`-based object with `GRangesList`-based `rowRanges` (this would require a version of `combine,GRangesList,GRangesList-method` that doesn't rely on `names`).

Each option has limitations.

Option 1 is restrictive since it requires `dimnames(x)` and `dimnames(y)` to be non-`NULL`, but matches the behaviour of other `combine` methods that are defined in terms of the `names`/`rownames`/`dimnames` of the objects.

Option 2 will work on wider class of `RangedSummarizedExperiments` objects since the `names` of such objects will commonly be `NULL` (it will not work for `RangedSummarizedExperiment` objects with `GRangesList`-derived `rowRanges` slots since the `combine` method only ).  This option is generally slower. The `names` of the "shared" elements are taken from `x` in the returned object. Option 2 will effectively use Option 1 for non-ranged `SummarizedExperiment0`-based objects.

I prefer option 2 since it is more general for a common use case.

A third option is a hybrid: use the first option if `names(x)` and `names(y)` are both non-`NULL` and the second option otherwise. This could lead to surprising behaviour, however, if the user expects the option 2 to be used (e.g., they are unaware that their objects have non-`NULL` `names`).

__RFC__: Which option is preferrable?


```r
setMethod("combine", c("SummarizedExperiment0", "SummarizedExperiment0"),
          function(x, y, ...) {

            if (any(dim(y) == 0L)) {
              return(x)
            } else if (any(dim(x) == 0L)) {
              return(y)
            }

            # Give a helpful error message if the user tries to combine a
            # RangedSummarizedExperiment object to an non-ranged
            # SummarizedExperiment0.
            # NOTE: Can't simply check if object is SumamrizedExperiment0
            #       object since all RangedSummarizedExperiments are
            #       SummarizedExperiment0 objects but not all
            #       SummarizedExperiment0 objects are
            #       RangedSummarizedExperiment objects.
            if ((is(x, "RangedSummarizedExperiment") &&
                 !is(y, "RangedSummarizedExperiment")) ||
                (!is(x, "RangedSummarizedExperiment") &&
                 is(y, "RangedSummarizedExperiment"))) {
              stop("Cannot combine '", class(x), "' and '", class(y),
                   "' objects because only one of these has 'rowRanges'.")
            }

            # Check dimnames are set
            x_dim <- dimnames(x)
            y_dim <- dimnames(y)
            if (is.null(x_dim) || is.null(y_dim) ||
                any(sapply(x_dim, is.null)) || any(sapply(y_dim, is.null))) {
              stop("'", class(x), "' must have dimnames for 'combine'")
            }

            # Combine slots
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: Can't handle internal duplicate ranges in x or y.
              if (any(duplicated(rowRanges(x))) || 
                  any(duplicated(rowRanges(y)))) {
                stop("Cannot combine '", class(x), "' with internal ", 
                     "duplicate 'rowRanges'")
              }
              rowRanges <- combine(rowRanges(x), rowRanges(y))
            } else {
              if (anyDuplicated(names(x)) || anyDuplicated(names(y))) {
                stop("Cannot combine '", class(x), "' with internal ", 
                     "duplicate 'names'")
              }
              NAMES <- .combine.NAMES(names(x), names(y))
            }
            colData <- combine(colData(x), colData(y))
            assays <- Assays(combine(assays(x, withDimnames = TRUE),
                                     assays(y, withDimnames = TRUE)))
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: elementMetadata slot of a RangedSummarizedExperiment
              #       object must be a zero-column DataFrame with nrow equal to 
              #       the length of the RangedSummarizedExperiment objects at 
              #       all times.
              elementMetadata <- DataFrame()
              elementMetadata@nrows <- length(rowRanges)
            } else {
              # NOTE: Using mcols() rather than slot(x, "elementMetadata") so
              #       that the combined objects have rownames on which to combine.
              elementMetadata <- combine(mcols(x, use.names = TRUE),
                                         mcols(y, use.names = TRUE))
              # NOTE: Once combined, drop rownames of elementMetadata since these 
              #       are given by rownames(z) when z is a SummarizedExperiment0 
              #       object.
              rownames(elementMetadata) <- NULL
            }
            metadata <- c(metadata(x), metadata(y))

            # Construct the combined SE
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: No need to update NAMES slot since it must be NULL in a
              #       valid RangedSummarizedExperiment object.
              BiocGenerics:::replaceSlots(x,
                                          rowRanges = rowRanges,
                                          colData = colData,
                                          assays = assays,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            } else {
              BiocGenerics:::replaceSlots(x,
                                          colData = colData,
                                          assays = assays,
                                          NAMES = NAMES,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            }
          }
)
#> [1] "combine"
```


```r
setMethod("combine", c("SummarizedExperiment0", "SummarizedExperiment0"),
          function(x, y, ...) {

            if (any(dim(y) == 0L)) {
              return(x)
            } else if (any(dim(x) == 0L)) {
              return(y)
            }
            
            # Give a helpful error message if the user tries to combine a
            # RangedSummarizedExperiment object to an non-ranged
            # SummarizedExperiment0.
            # NOTE: Can't simply check if object is SumamrizedExperiment0
            #       object since all RangedSummarizedExperiments are
            #       SummarizedExperiment0 objects but not all
            #       SummarizedExperiment0 objects are
            #       RangedSummarizedExperiment objects.
            if ((is(x, "RangedSummarizedExperiment") &&
                 !is(y, "RangedSummarizedExperiment")) ||
                (!is(x, "RangedSummarizedExperiment") &&
                 is(y, "RangedSummarizedExperiment"))) {
              stop("Cannot combine '", class(x), "' and '", class(y),
                   "' objects because only one of these has 'rowRanges'")
            }

            # Check colnames are set
            x_cn <- colnames(x)
            y_cn <- colnames(y)
            if (is.null(x_cn) || is.null(y_cn)) {
              stop("Cannot combine '", class(x), "' with NULL 'colnames'")
            }
            
            # Combine slots
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: This method doesn't work if rowRanges(x) or rowRanges(y) 
              #       are GRangesList-derived objects since the 
              #       combine,GRangesList,GRangesList-method currently requires 
              #       that these have non-NULL names.
              if (is(rowRanges(x), "GRangesList") || 
                  is(rowRanges(y), "GRangesList")) {
                stop("Cannot combine '", class(x), "' objects with ", 
                     "'GRangesList'-based 'rowRanges'")
              }
              # NOTE: Can't handle internal duplicate ranges in x or y.
              if (any(duplicated(rowRanges(x))) || 
                  any(duplicated(rowRanges(y)))) {
                stop("Cannot combine '", class(x), "' with internal ", 
                     "duplicate 'rowRanges'")
              }
              # NOTE: mcols(x) and mcols(y) 
              #       [i.e. mcols(rowRanges(x)) and mcols(rowRanges(y))] 
              #       are separately combined. This could be simplified 
              #       if combine,GRanges,GRanges-method had an 
              #       ignore.mcols argument.
              x_em <- mcols(x, use.names = TRUE)
              mcols(x) <- NULL
              y_em <- mcols(y, use.names = TRUE)
              mcols(y) <- NULL
              
              # NOTE: names(rowRanges) are set to NULL
              rowRanges <- combine(rowRanges(x), rowRanges(y))
              # Set rownames based on findOverlaps() of rowRanges()
              # NOTE: We want ranges that are identical, hence 
              #       type = "equal". Also, we want to identify identical 
              #       zero-width ranges, hence minoverlap = 0
              x_ol <- findOverlaps(rowRanges(x), rowRanges, 
                                   type = "equal", minoverlap = 0L)
              rownames(x) <- subjectHits(x_ol)
              rownames(x_em) <- subjectHits(x_ol)
              y_ol <- findOverlaps(rowRanges(y), rowRanges, 
                                   type = "equal", minoverlap = 0L)
              rownames(y) <- subjectHits(y_ol)
              rownames(y_em) <- subjectHits(y_ol)
              mcols(rowRanges) <- combine(x_em, y_em)
            } else {
              if (anyDuplicated(names(x)) || anyDuplicated(names(y))) {
                stop("Cannot combine '", class(x), "' with internal ", 
                     "duplicate 'names'")
              }
              NAMES <- .combine.NAMES(names(x), names(y))
            }
            colData <- combine(colData(x), colData(y))
            assays <- Assays(combine(assays(x, withDimnames = TRUE),
                                     assays(y, withDimnames = TRUE)))
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: elementMetadata slot of a RangedSummarizedExperiment
              #       object must be a zero-column DataFrame with nrow equal to 
              #       the length of the RangedSummarizedExperiment objects at 
              #       all times.
              elementMetadata <- DataFrame()
              elementMetadata@nrows <- length(rowRanges)
            } else {
              # NOTE: Using mcols() rather than slot(x, "elementMetadata") so
              #       that the combined objects have rownames on which to combine.
              elementMetadata <- combine(mcols(x, use.names = TRUE),
                                         mcols(y, use.names = TRUE))
              # NOTE: Once combined, drop rownames of elementMetadata since these 
              #       are given by rownames(z) when z is a SummarizedExperiment0 
              #       object.
              rownames(elementMetadata) <- NULL
            }
            metadata <- c(metadata(x), metadata(y))

            # Construct the combined SE
            if (is(x, "RangedSummarizedExperiment")) {
              # NOTE: No need to update NAMES slot since it must be NULL in a
              #       valid RangedSummarizedExperiment object.
              BiocGenerics:::replaceSlots(x,
                                          rowRanges = rowRanges,
                                          colData = colData,
                                          assays = assays,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            } else {
              BiocGenerics:::replaceSlots(x,
                                          colData = colData,
                                          assays = assays,
                                          NAMES = NAMES,
                                          elementMetadata = elementMetadata,
                                          metadata = metadata)
            }
          }
)
#> [1] "combine"
```

__RFC__: Should the check of whether the user is trying to combine a non-ranged `SummarizedExperiment0`-based object with a `RangedSummarizedExeriment`-based object be the more restrictive `if(class(x) != class(y)) stop()` (this is what `combine,eSet,eSet-method` uses)?

__RFC__: Are non-`matrix` objects allowed as elements in an `Assays` object? If so, will need a `combine()` method for each of these allowed classes.

__RFC__: `metadata` from all objects are combined into a list with no name checking (following `cbind,SummarizedExperiment0-method` and `rbind,SummarizedExperiment0-method`). Does this seem reasonable?

__RFC__: If `nrow(x)` (resp. `nrow(y)`) is zero (i.e. no features/ranges in that object) or `ncol(x)` (resp. `ncol(y)`) is zero (i.e. no samples in that object) then `y` (resp. `x`) is returned (or `x` if the condition occurs for both samples); does this seem reasonable? In particular, the zero-column object could have features/ranges that the user wants to add to the other object (or the zero-row object could have samples that the user wants to add to the other object). How to handle these scenarios?

__RFC__: I thought that `colnames` were a requirement of a `SummarizedExperiment0` object, but this is not true, at least according to the validity method (the `SummarizedExperiment()` constructor __does__ enforce non-`NULL` `colnames`). Is this a bug in the validity method?


```r
example("RangedSummarizedExperiment", echo = FALSE)
rse
#> class: RangedSummarizedExperiment 
#> dim: 200 6 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowRanges metadata column names(1): feature_id
#> colnames(6): A B ... E F
#> colData names(1): Treatment
validObject(rse)
#> [1] TRUE

colnames(rse) <- NULL
rse
#> class: RangedSummarizedExperiment 
#> dim: 200 6 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowRanges metadata column names(1): feature_id
#> colnames: NULL
#> colData names(1): Treatment
validObject(rse)
#> [1] TRUE
```

__RFC__: Option 2 uses `findOverlaps()` with `type = "equal"` (because we want identical ranges) and `minoverlap = 0` (because otherwise zero-width identical ranges are not found). Are the missed ranges when `minoverlap = 1` (the default) intended behaviour? For example:


```r
x <- GRanges("chr1", IRanges(10, width = 0))
x
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1   [10, 9]      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# No hits
findOverlaps(x, x, type = "equal")
#> Hits object with 0 hits and 0 metadata columns:
#>    queryHits subjectHits
#>    <integer>   <integer>
#>   -------
#>   queryLength: 1
#>   subjectLength: 1

# 1 hit
findOverlaps(x, x, type = "equal", minoverlap = 0L)
#> Hits object with 1 hit and 0 metadata columns:
#>       queryHits subjectHits
#>       <integer>   <integer>
#>   [1]         1           1
#>   -------
#>   queryLength: 1
#>   subjectLength: 1
```

__RFC__: Both Option 1 and 2 will return an error if either `x` or `y` contains "internal" duplicate features/ranges. I don't know how else you would handle this situation; suppose `x` contains an "internal" duplicate, it would first require an "internal combine" of `x` followed by a "combine"  with `y` (see example below).

## Discuss the signature of the `combine` generic

My understanding is that the `combine()` generic defined in _BiocGenerics_ does not allow for additional arguments (at least, additional arguments that are not objects to be "combined"). It might be useful to be able to pass `ignore.mcols` to `combine,GRanges,GRanges-method`, `combine,GRangesList,GRangesList`, and `combine,SummarizedExperiment0,SummarizedExperiment0-method`. I believe this would require the `combine()` generic to be redefined.

# Examples

I use data from the `SummarizedExperiment::SummarizedExeriment` and `SummarizedExperimentRangedSummarizedExperiment` man pages.


```r
# Need to set seed because example("RangedSummarizedExperiment") and 
# example("SummarizedExperiment") both call runif().
set.seed(667)
```

## `SummarizedExperiment0`


```r
# Get some example data
example("SummarizedExperiment0", echo = FALSE)

# NOTE: SummarizedExperiment0 objects must have non-NULL NAMES slot in order 
#       to combine
names(se0) <- as.character(1:200)

# Create data to combine
A <- se0[1:4, "A"]
B <- se0[3:6, "B"]
C <- se0[5:8, "C"]
BC <- se0[c(4:7, 9), c("B", "C")]
D <- se0[7:10, "D"]
E <- se0[9:12, "E"]
F <- se0[11:14, "F"]

# Sanity check: identical to cbind when given compatible arguments
all.equal(cbind(se0[, 1], se0[, 2], se0[, 3]),
          combine(se0[, 1], se0[, 2], se0[, 3]))
#> [1] TRUE

# Combining objects with non-overlapping features and non-distinct samples
z <- combine(A, B, C, BC)
z
#> class: SummarizedExperiment0 
#> dim: 9 3 
#> metadata(0):
#> assays(1): counts
#> rownames(9): 1 2 ... 8 9
#> metadata column names(0):
#> colnames(3): A B C
#> colData names(1): Treatment
assay(z)
#>          A        B        C
#> 1 8.274185       NA       NA
#> 2 9.165680       NA       NA
#> 3 9.654497 9.698581       NA
#> 4 9.554457 9.773346 9.774785
#> 5       NA 7.065869 9.516848
#> 6       NA 9.552036 9.338527
#> 7       NA 5.872714 9.788340
#> 8       NA       NA 9.424306
#> 9       NA 9.901994 9.297545

# Incomplete elementMetadata are elegantly handled
a <- A
mcols(a) <- DataFrame(J = seq_len(nrow(a)))
z <- combine(a, B)
assay(z)
#>          A        B
#> 1 8.274185       NA
#> 2 9.165680       NA
#> 3 9.654497 9.698581
#> 4 9.554457 9.773346
#> 5       NA 7.065869
#> 6       NA 9.552036
mcols(z)
#> DataFrame with 6 rows and 1 column
#>           J
#>   <integer>
#> 1         1
#> 2         2
#> 3         3
#> 4         4
#> 5        NA
#> 6        NA

# "Internal duplicates" are not allowed
x <- se0[1:2, 1]
# Create an "internal duplicate" feature name
names(x)[2] <- names(x)[1]
y <- se0[10, 2]
# Expect error
combine(x, y)
#> Error in combine(x, y): Cannot combine 'SummarizedExperiment0' with internal duplicate 'names'
```

## `RangedSummarizedExperiment`


```r
# Get some example data
example("RangedSummarizedExperiment", echo = FALSE)

# Create 
A <- rse[1:4, "A"]
B <- rse[3:6, "B"]
C <- rse[5:8, "C"]
BC <- rse[c(4:7, 9), c("B", "C")]
D <- rse[7:10, "D"]
E <- rse[9:12, "E"]
F <- rse[11:14, "F"]

# Sanity check: identical to cbind when given compatible arguments
all.equal(cbind(rse[, 1], rse[, 2], rse[, 3]),
          combine(rse[, 1], rse[, 2], rse[, 3]))
#> [1] TRUE

# Combining objects with non-overlapping features and non-distinct samples
z <- combine(A, B, C, BC)
z
#> class: RangedSummarizedExperiment 
#> dim: 9 3 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowRanges metadata column names(1): feature_id
#> colnames(3): A B C
#> colData names(1): Treatment
assay(z)
#>              A        B        C
#>  [1,] 9.283199       NA       NA
#>  [2,] 8.786685       NA       NA
#>  [3,] 9.412699 8.916581       NA
#>  [4,] 9.851115 9.261356 9.770541
#>  [5,]       NA 9.248532 7.105869
#>  [6,]       NA 9.225918 8.958373
#>  [7,]       NA 8.642004 8.422512
#>  [8,]       NA       NA 8.675099
#>  [9,]       NA 8.332671 9.868225

# Incomplete elementMetadata are elegantly handled
a <- A
colnames(a) <- "a"
mcols(a) <- DataFrame(J = seq_len(nrow(a)))
z <- combine(a, B)
assay(z)
#>             a        B
#> [1,] 9.283199       NA
#> [2,] 8.786685       NA
#> [3,] 9.412699 8.916581
#> [4,] 9.851115 9.261356
#> [5,]       NA 9.248532
#> [6,]       NA 9.225918
mcols(z)
#> DataFrame with 6 rows and 2 columns
#>           J  feature_id
#>   <integer> <character>
#> 1         1          NA
#> 2         2          NA
#> 3         3       ID003
#> 4         4       ID004
#> 5        NA       ID005
#> 6        NA       ID006

# "Internal duplicates"" are not allowed
x <- rse[1:2, 1]
# Create an "internal duplicate" ranges within x
rowRanges(x)[2] <- rowRanges(x)[1]
y <- rse[10, 2]
# Expect error
combine(x, y)
#> Error in combine(x, y): Cannot combine 'RangedSummarizedExperiment' with internal duplicate 'rowRanges'

# Some examples of what happens when rownames are available on only on x or y
A <- rse[1:4, "A"]
B <- rse[3:6, "B"]
z <- combine(A, B)
rownames(z)
#> NULL
assay(z)
#>             A        B
#> [1,] 9.283199       NA
#> [2,] 8.786685       NA
#> [3,] 9.412699 8.916581
#> [4,] 9.851115 9.261356
#> [5,]       NA 9.248532
#> [6,]       NA 9.225918

# "Shared" ranges have names taken from 'x'
names(A) <- letters[seq_along(A)]
names(B) <- LETTERS[seq_along(B)]
zz <- combine(A, B)
rownames(zz)
#> [1] "a" "b" "c" "d" "C" "D"
assay(zz)
#>          A        B
#> a 9.283199       NA
#> b 8.786685       NA
#> c 9.412699 8.916581
#> d 9.851115 9.261356
#> C       NA 9.248532
#> D       NA 9.225918

# elementMetadata are still combined even though rownames differ for "shared"
# ranges
mcols(A) <- DataFrame(J = seq_along(A))
zzz <- combine(A, B)
rownames(zzz)
#> [1] "a" "b" "c" "d" "C" "D"
assay(zzz)
#>          A        B
#> a 9.283199       NA
#> b 8.786685       NA
#> c 9.412699 8.916581
#> d 9.851115 9.261356
#> C       NA 9.248532
#> D       NA 9.225918

# Can end up with an object with "incomplete" rownames
rownames(A) <- NULL
zzzz <- combine(A, B)
rownames(zzzz)
#> [1] ""  ""  ""  ""  "C" "D"
assay(zzzz)
#>          A        B
#>   9.283199       NA
#>   8.786685       NA
#>   9.412699 8.916581
#>   9.851115 9.261356
#> C       NA 9.248532
#> D       NA 9.225918
```

# Summary

The `combine()` method for `SE` objects addresses the aim of being able to combine multiple `SE` objects when they have potentially different features/ranges and potentially different samples. Furthermore, I have also defined `combine,DataFrame,DataFrame-method`, `combine,GRanges,GRanges-method`, `combine,GRangesList,GRangesList-method`, and `combine,SimpleList,SimpleList-method`.

The chief limitation is an incomplete `combine,GRangesList,GRangesList-method`; we need to be able to appropriately "combine" `GRangesList` objects when `names` are `NULL`. Without this, the `combine,SummarizedExperiment0,SummarizedExperiment0-method` is also incomplete.

It is perhaps worth noting that `identical(combine(x, y), combine(y, x))` is generally `FALSE`.

## Futher work

- Complete support of `RangedSummarizedExperiment` objects where the `rowRanges` are `GRangesList`-derived objects, e.g., `DESeqDataSet` objects.
- Unit tests
- Documentation

# Session info


```r
devtools::session_info()
#> Session info --------------------------------------------------------------
#>  setting  value                                             
#>  version  R Under development (unstable) (2015-10-13 r69511)
#>  system   x86_64, darwin13.4.0                              
#>  ui       X11                                               
#>  language (EN)                                              
#>  collate  en_AU.UTF-8                                       
#>  tz       Australia/Melbourne                               
#>  date     2015-10-20
#> Packages ------------------------------------------------------------------
#>  package              * version date      
#>  Biobase              * 2.31.0  2015-10-14
#>  BiocGenerics         * 0.17.0  2015-10-14
#>  devtools               1.9.1   2015-09-11
#>  digest                 0.6.8   2014-12-31
#>  evaluate               0.8     2015-09-18
#>  formatR                1.2.1   2015-09-18
#>  GenomeInfoDb         * 1.7.0   2015-10-14
#>  GenomicRanges        * 1.21.32 2015-10-14
#>  htmltools              0.2.6   2014-09-08
#>  IRanges              * 2.4.0   2015-10-14
#>  knitr                  1.11    2015-08-14
#>  magrittr               1.5     2014-11-22
#>  memoise                0.2.1   2014-04-22
#>  rmarkdown              0.8.1   2015-10-10
#>  S4Vectors            * 0.9.1   2015-10-17
#>  stringi                0.5-5   2015-06-29
#>  stringr                1.0.0   2015-04-30
#>  SummarizedExperiment * 1.1.0   2015-10-14
#>  XVector                0.11.0  2015-10-14
#>  yaml                   2.1.13  2014-06-12
#>  zlibbioc               1.17.0  2015-10-14
#>  source                                      
#>  Bioconductor                                
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  Bioconductor                                
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  Github (Bioconductor-mirror/IRanges@b0b5fff)
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  Bioconductor                                
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  Bioconductor
```
