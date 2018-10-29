
<!-- README.md is generated from README.Rmd. Please edit that file -->
SPLITT: Serial and Parallel LIneage Traversal of Trees
======================================================

SPLITT is a one-header C++ library that makes it easy to program fast pre-order and post-order traversal algorithms over tree-like data structures. The library is written in C++11. If [OpenMP](https://www.openmp.org) is enabled during compilation, SPLITT can take advantage of a multi-core CPU and/or vectorized instructions to execute the tree traversal in parallel.

The SPLITT.h file
-----------------

The SPLITT library represents a C++ header file called "SPLITT.h". Technically, downloading [the latest version of the file from github](https://github.com/venelin/SPLITT/raw/master/src/SPLITT.h) and including it in one or more of your C++ files is all you need to start using it. However, before you can do anything useful with SPLITT you will need to get familiar with some of its classes, i.e. its programming interface.

The SPLITT R-package
--------------------

The R-package called with the same name (SPLITT) bundles the SPLITT.h header together with its documentation and examples. You don't need to install this package nor do you need to know R or to have an R environment on your system for using the SPLITT library. Nevertheless, if your final goal is to write a C++ extension that uses SPLITT as a fast back-end for an R-package of your own, the SPLITT package provides an example how you can do that.

Resources
---------

The user guides and technical reference for the library are available from the [SPLITT web-page](https://venelin.github.io/SPLITT/index.html).

The [Getting started](https://venelin.github.com/SPLITT/articles/SPLITT.html) guide shows how to write a console application using SPLITT.

The [SPLITT classes](https://venelin.github.com/SPLITT/articles/SPLITTClasses.html) guide overviews the SPLITT classes.

The [Calling SPLITT from an R-package](https://venelin.github.com/SPLITT/articles/SPLITTRcppModules.html) guide shows how to use Rcpp to wrap SPLITT traversal functionality into an R-package.

The research article "Parallel Likelihood Calculation for Phylogenetic Comparative Models: the SPLITT C++ Library" provides a general overview of SPLITT. The example application to the PMM model depicted on Fig. 1 is in this article is used as a coding example in the SPLITT online documentation. The article is currently undergoing peer review for a publication and is available as a preprint from [BioRxiv](https://www.biorxiv.org/content/early/2017/12/18/235739).

The SPLITT source code is located in the [SPLITT github repository](https://github.com/venelin/SPLITT).
