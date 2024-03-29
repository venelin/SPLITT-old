% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SPLITT::Tree::LengthType}
\alias{SPLITT::Tree::LengthType}
\title{Abstract type for the branch-lengths in a \code{\link[=SPLITT::Tree]{Tree}}.}
\description{
\code{typedef Length LengthType;}

A public typedef in class \code{\link[=SPLITT::Tree]{Tree}}. 
  A synonym for the template argument Length. Defines a type that can be 
  associated with a branch. Can be a basic type, e.g. \code{double}, but also
  a composite of several attributes on a branch, such as a \code{double} 
  length and an \code{int} color.
}
\seealso{
\code{\link[=SPLITT::Tree]{Tree}}

\link{SPLITT}
}
