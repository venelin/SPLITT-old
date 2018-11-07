
<!-- README.md is generated from README.Rmd. Please edit that file -->
SPLITT: Serial and Parallel LIneage Traversal of Trees
======================================================

SPLITT is a one-header C++ library that makes it easy to program fast pre-order and post-order traversal algorithms over tree-like data. The library is written in C++11. If [OpenMP](https://www.openmp.org) is enabled during compilation, SPLITT can take advantage of a multi-core CPU and/or vectorized instructions to execute the tree traversal in parallel.

The SPLITT.h file
-----------------

The SPLITT library represents a C++ header file called "SPLITT.h". Technically, downloading [the latest version of the file from github](https://github.com/venelin/SPLITT/raw/master/src/SPLITT.h) and including it in one or more of your C++ files is all you need to start using it. However, before you can do anything useful with SPLITT you will need to get familiar with some of its classes, i.e. its programming interface.

The SPLITT R-package
--------------------

The R-package called with the same name (SPLITT) bundles the SPLITT.h header together with its documentation and examples. You don't need to install this package nor do you need to know R or to have an R environment on your system for using the SPLITT library.

User guides and resources
-------------------------

The user guides and technical reference for the library are available from the [SPLITT web-page](https://venelin.github.io/SPLITT/index.html).

-   The [Getting started](https://venelin.github.io/SPLITT/articles/SPLITT.html) guide shows how to compile and run a console application using SPLITT.
-   The [Writing a traversal specification](https://venelin.github.io/SPLITT/articles/SPLITTTraversalSpecification.html) guide shows how to define the application specific data-type and node-traversal operations.
-   The [Running a traversal](https://venelin.github.io/SPLITT/articles/SPLITTRunTraversal.html) guide shows how to create and run a traversal task.
-   The [View from above](https://venelin.github.io/SPLITT/articles/SPLITTClasses.html) guide overviews the SPLITT classes.
-   The [Calling SPLITT from an R-package](https://venelin.github.io/SPLITT/articles/SPLITTRcppModules.html) guide shows how to use Rcpp to wrap SPLITT traversal functionality into an R-package.

The research article "Parallel Likelihood Calculation for Phylogenetic Comparative Models: the SPLITT C++ Library" provides a general overview of SPLITT. The example application to the PMM model depicted on Fig. 1 in this article is used as a coding example in the SPLITT online documentation. The article is currently in peer review for a publication and is available as a preprint from [BioRxiv](https://www.biorxiv.org/content/early/2018/10/29/235739).

The SPLITT source code is located in the [SPLITT github repository](https://github.com/venelin/SPLITT).

Feature requests, bugs, etc can be reported in the [SPLITT issues list](https://github.com/venelin/SPLITT/issues).

Citing SPLITT
-------------

To give credit to the SPLITT library in a publication, please cite the following article:

Mitov, V., & Stadler, T. (2018). Parallel Likelihood Calculation for Phylogenetic Comparative Models: the SPLITT C++ Library. bioRxiv, 235739. <http://doi.org/10.1101/235739>

Used R-packages
---------------

The SPLITT R-package uses the following 3rd party R-packages:

-   For calling C++ from R: Rcpp v0.12.19 (Eddelbuettel et al. 2017), RcppArmadillo v0.9.100.5.0 (Eddelbuettel, Francois, and Bates 2016)
-   For tree processing in R: ape v5.2 (Paradis et al. 2016);
-   For reporting the test-benchmark: data.table v1.11.8 (Dowle and Srinivasan 2016);
-   For unit-testing: testthat v2.0.1 (Wickham 2016).

References
----------

Dowle, Matt, and Arun Srinivasan. 2016. *Data.table: Extension of ‘Data.frame‘*. <https://CRAN.R-project.org/package=data.table>.

Eddelbuettel, Dirk, Romain Francois, JJ Allaire, Kevin Ushey, Qiang Kou, Nathan Russell, Douglas Bates, and John Chambers. 2017. *Rcpp: Seamless R and C++ Integration*. <https://CRAN.R-project.org/package=Rcpp>.

Eddelbuettel, Dirk, Romain Francois, and Doug Bates. 2016. *RcppArmadillo: ’Rcpp’ Integration for the ’Armadillo’ Templated Linear Algebra Library*. <https://CRAN.R-project.org/package=RcppArmadillo>.

Paradis, Emmanuel, Simon Blomberg, Ben Bolker, Julien Claude, Hoa Sien Cuong, Richard Desper, Gilles Didier, et al. 2016. *Ape: Analyses of Phylogenetics and Evolution*. <https://CRAN.R-project.org/package=ape>.

Wickham, Hadley. 2016. *Testthat: Unit Testing for R*. <https://CRAN.R-project.org/package=testthat>.
