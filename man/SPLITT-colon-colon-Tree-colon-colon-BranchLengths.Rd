% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SPLITT::Tree::BranchLengths}
\alias{SPLITT::Tree::BranchLengths}
\title{Get a reference to the internal vector of branch lengths.}
\value{
\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const&}; 
  a const reference to the internally stored vector of branch lengths, in the order of the end-node ids.
}
\description{
\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const& BranchLengths() const;}
}
\seealso{
\code{\link[=SPLITT::Tree]{Tree}}

\link{SPLITT}

Other public methods in SPLITT::Tree: \code{\link{SPLITT::Tree::FindChildren}},
  \code{\link{SPLITT::Tree::FindIdOfNode}},
  \code{\link{SPLITT::Tree::FindIdOfParent}},
  \code{\link{SPLITT::Tree::FindNodeWithId}},
  \code{\link{SPLITT::Tree::HasBranchLengths}},
  \code{\link{SPLITT::Tree::LengthOfBranch}},
  \code{\link{SPLITT::Tree::OrderNodesPosType}},
  \code{\link{SPLITT::Tree::OrderNodes}},
  \code{\link{SPLITT::Tree::SetBranchLengths}},
  \code{\link{SPLITT::Tree::SetLengthOfBranch}},
  \code{\link{SPLITT::Tree::Tree}},
  \code{\link{SPLITT::Tree::num_nodes}},
  \code{\link{SPLITT::Tree::num_tips}}
}
\concept{public methods in SPLITT::Tree}
