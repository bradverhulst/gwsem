# gwsem

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jpritikin/gwsem.svg?branch=master)](https://travis-ci.org/jpritikin/gwsem)
[![Codecov test coverage](https://codecov.io/gh/jpritikin/gwsem/branch/master/graph/badge.svg)](https://codecov.io/gh/jpritikin/gwsem?branch=master)
[![cran version](http://www.r-pkg.org/badges/version/gwsem)](https://cran.r-project.org/package=gwsem)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/gwsem)](http://cranlogs.r-pkg.org/badges/gwsem)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/gwsem)](http://cranlogs.r-pkg.org/badges/grand-total/gwsem)
<!-- badges: end -->

The goal of gwsem is to provide users with the opportunity to analyze the complex, interconnected array of risk factors, biomarkers, environmental antecedents, comorbid disorders, and other health outcomes on a genome-wide basis using structural equation modeling techniques.

## Installation

GW-SEM utilizes the optimization function of OpenMx.

You can install the released version of OpenMx from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("OpenMx")
```

If you want a new version of OpenMx, you can follow the instruction to build it from source [HERE](https://openmx.ssri.psu.edu).

GW-SEM is currently under development. To install gwsem you will need a few dependencies found below:

``` r
install.packages("roxygen2")
install.packages("qqman")
install.packages("rpf")
install.packages("RcppEigen")
install.packages("StanHeaders")
install.packages("BH")
install.packages("testthat")
install.packages("digest")
install.packages("Rcpp")
install.packages("lifecycle")
install.packages("devtools")
```

Then, you can clone the source code and do

```
git clone https://github.com/jpritikin/gwsem
cd gwsem 
./tools/rox
R CMD INSTALL .
Rscript tools/test.R
```

You cannot use **devtools** `install_github` because it does not run **roxygen2**.
