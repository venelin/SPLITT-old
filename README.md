
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis build status](https://travis-ci.org/venelin/SPLITT.svg?branch=master)](https://travis-ci.org/venelin/SPLITT) [![Coverage status](https://codecov.io/gh/venelin/SPLITT/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/SPLITT?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/SPLITT?color=blue)](http://cran.r-project.org/web/packages/SPLITT) [![Downloads](http://cranlogs.r-pkg.org/badges/SPLITT?color=blue)](http://cran.rstudio.com/package=SPLITT)

SPLITT: Serial and Parallel LIneage Traversal of Trees
======================================================

SPLITT is a one-header C++ library that makes it easy to program fast pre-order and post-order traversal algorithms over tree-like data. The library is written in C++11. If [OpenMP](https://www.openmp.org) is enabled during compilation, SPLITT can take advantage of a multi-core CPU and/or vectorized instructions to execute the tree traversal in parallel.

The SPLITT.h file
-----------------

The SPLITT library represents a C++ header file called "SPLITT.h". Technically, downloading [the latest version of the file from github](https://github.com/venelin/SPLITT/raw/master/src/SPLITT.h) and including it in one or more of your C++ files is all you need to start using it. However, before you can do anything useful with SPLITT you will need to get familiar with some of its classes, i.e. its programming interface.

The SPLITT R-package
--------------------

The R-package called with the same name (SPLITT) bundles the SPLITT.h header together with some helper functions for creating R-packages using SPLITT and this documentation. You will find it useful to install this package if you wish to create R-pacages that call SPLITT via Rcpp modules. This is explained in the [Calling SPLITT from an R-package guide](https://venelin.github.io/SPLITT/articles/SPLITTRcppModules.html).

User guides and resources
=========================

The [Get started guide](https://venelin.github.io/SPLITT/articles/SPLITT.html) provides a general introduction to the library and a roadmap throughout the documentation.

The research article "Parallel Likelihood Calculation for Phylogenetic Comparative Models: the SPLITT C++ Library" introduces the library and reports a large-scale performance benchmark. The example application to the PMM model depicted on Fig. 1 in this article is used as a coding example in the SPLITT online documentation. The article is currently in peer review for a publication and is available as a preprint from [BioRxiv](https://www.biorxiv.org/content/early/2018/10/29/235739).

-   The SPLITT documentation is available from the [SPLITT web page](https://venelin.github.io/SPLITT).
-   The SPLITT source code is located in the [SPLITT github repository](https://github.com/venelin/SPLITT).
-   The classes of the SPLITT API are documented in the [SPLITT Reference page](https://venelin.github.io/SPLITT/reference/SPLITT.html).
-   Feature requests, bugs, etc can be reported in the [SPLITT issues list](https://github.com/venelin/SPLITT/issues).

Citing SPLITT
=============

To give credit to the SPLITT library in a publication, please cite the following article:

Mitov, V., & Stadler, T. (2018). Parallel Likelihood Calculation for Phylogenetic Comparative Models: the SPLITT C++ Library. bioRxiv, 235739. <http://doi.org/10.1101/235739>

Used R-packages
===============

The SPLITT R-package uses the following 3rd party R-packages:

-   For calling C++ from R: Rcpp v1.0.0 (Eddelbuettel et al. 2018);
-   For tree processing in R: ape v5.2 (Paradis et al. 2018);
-   For unit-testing: testthat v2.0.1 (Wickham 2018);
-   For documentation and web-site generation: roxygen2 v6.1.1 (Wickham, Danenberg, and Eugster 2018), pkgdown v1.1.0.9000 (Wickham and Hesselberth 2018);

References
==========

Eddelbuettel, Dirk, Romain Francois, JJ Allaire, Kevin Ushey, Qiang Kou, Nathan Russell, Douglas Bates, and John Chambers. 2018. *Rcpp: Seamless R and C++ Integration*. <https://CRAN.R-project.org/package=Rcpp>.

Paradis, Emmanuel, Simon Blomberg, Ben Bolker, Joseph Brown, Julien Claude, Hoa Sien Cuong, Richard Desper, et al. 2018. *Ape: Analyses of Phylogenetics and Evolution*. <https://CRAN.R-project.org/package=ape>.

Wickham, Hadley. 2018. *Testthat: Unit Testing for R*. <https://CRAN.R-project.org/package=testthat>.

Wickham, Hadley, and Jay Hesselberth. 2018. *Pkgdown: Make Static Html Documentation for a Package*.

Wickham, Hadley, Peter Danenberg, and Manuel Eugster. 2018. *Roxygen2: In-Line Documentation for R*. <https://CRAN.R-project.org/package=roxygen2>.
