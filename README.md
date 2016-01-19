# R package: SparseSummarizedExperiment

[![Project Status: Wip - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/0.1.0/wip.svg)](http://www.repostatus.org/#wip)
[![Linux Build Status](https://travis-ci.org/PeteHaitch/SparseSummarizedExperiment.svg?branch=master)](https://travis-ci.org/PeteHaitch/SparseSummarizedExperiment)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/github/PeteHaitch/SparseSummarizedExperiment?svg=true)](https://ci.appveyor.com/project/PeteHaitch/SparseSummarizedExperiment)
[![Coverage Status](https://img.shields.io/codecov/c/github/PeteHaitch/SparseSummarizedExperiment/master.svg)](https://codecov.io/github/PeteHaitch/SparseSummarizedExperiment?branch=master)

---

This is an __experimental__ package that defines S4 classes and methods for 
sparse genomic data. The API (and, indeed, much of the infrastructure) is based 
on the [_SummarizedExperiment_ Bioconductor package](http://bioconductor.org/packages/SummarizedExperiment/).

Fair warning, while the package is now reasonably stable with most functionality documented and ever-improving test coverage, the final API is not locked in.

---

## Installation
 
_SparseSummarizedExperiment_ is in development and can only be installed using 
the development version of Bioconductor. Please first read 
[these instructions on installing the development version of Bioconductor](http://www.bioconductor.org/developers/how-to/useDevel/). 

```r
devtools::install_github("PeteHaitch/SparseSummarizedExperiment")
```

## License

The MIT License (MIT)

Copyright (c) 2015-2016 Peter Hickey

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
