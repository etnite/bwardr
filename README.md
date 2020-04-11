# bwardr

Author: [Brian Ward](https://brianpward.net/)  
Email: [brian@brianpward.net](mailto:brian@brianpward.net)  
Github: [etnite](https://github.com/etnite)  
License: [MIT](https://opensource.org/licenses/MIT)

[![Travis build status](https://travis-ci.org/etnite/bwardr.svg?branch=master)](https://travis-ci.org/etnite/bwardr)

This is a personal R package containing a collection of functions that I have written over time to help in my own research. I work in the field of wheat genetics, and therefore this package contains a mixture of functions relating to bioinformatics, genetics, and agronomic analyses. I have made it publicly available in case it might be of use to others, but please note that it lacks the polish of a CRAN-hosted package. Individual functions do have documentation, but there is no manual. Being that this is an R package, it is suited to moderately-sized genomic datasets. Those working with very large files will likely be better served using the usual purpose-built C/C++ or Java executables.

This package can be installed using the devtools package:

```r
library(devtools)
install_github("etnite/bwardr")
```

Or to install with vignettes:

```r
library(devtools)
install_github("etnite/bwardr", build_vignettes = TRUE)
```

## Functionality

Some demonstrations are forthcoming! However, the package includes (but is not limited to) the following functions:

* eigstrat() - an R port of the [Eigenstrat](https://github.com/DReichLab/EIG/tree/master/EIGENSTRAT) algorithm for delineating population structure using principal component analysis of a SNP matrix
* stand_str() - a function for standardizing strings to a consistent format, to facilitate successful intersections of sets
* format_qtlmap_geno() - a function for converting VCF or PLINK binary PED (BED) files storing data on a biparental population into files suitable for input to [R/qtl](http://www.rqtl.org/) 
* fdr_thresh() - finds an approximate p-value threshold corresponding to a specified false-discovery-rate threshold
* mc_train_val_multienv() - generate training/validation splits in multi-environment phenotypic datasets
