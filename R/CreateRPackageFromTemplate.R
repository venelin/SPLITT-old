# CreateRPackageFromTemplate.R
# SPLITT
# 
# Copyright 2018 Venelin Mitov
# 
# This file is part of SPLITT: a generic C++ library for Serial and Parallel
# Lineage Traversal of Trees.
# 
# SPLITT is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# 
# SPLITT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with SPLITT.  If not, see
# <http://www.gnu.org/licenses/>.
# 
# @author Venelin Mitov

#' Create an R-package based on a template
#' @param packageName a character string denoting the name of the package.
#' @param path a character string denoting a directory where to hold the package.
#' @param class a character string denoting the name of a TraversalSpecification 
#' class to be implemented in the package. 
#' @details This function needs Internet connection to clone the 
#' \link{https://github.com/venelin/PMMUsingSPLITT.git} repository in \code{path}. 
#' Then, it replaces the name of the package and the the C++ namespace 
#' \code{PMMUsingSPLITT} with \code{packageName} and the name of the class 
#' \code{AbcPMM} with \code{className}. All these arguments must
#' be valid C++ identifiers. The \code{packageName} parameter must be a valid 
#' name for an R-package. 
#' @return If the function runs correctly nothing. Otherwise an error is raised.
CreateRPackageFromTemplate <- function(
  packageName, path = ".", class = "MyTraversalTask", 
  overwrite = FALSE) {
  
  dp <- file.path(path, packageName)
  
  if(file.exists(file.path(path, packageName))) {
    if(overwrite) {
      unlink(dp, recursive = TRUE) 
    } else {
      stop(paste0("Directory ", dp,
                  " already exists. Set overwrite to TRUE to overwrite."))
    }
  }
  
  repo <- clone("https://github.com/venelin/PMMUsingSPLITT.git", dp)
  unlink(file.path(path, packageName, ".git"), recursive = TRUE)
  
  ## Replace packageName and class name everywhere
  for( f in list.files(dp, recursive = TRUE) ){
    fp <- file.path(dp, f)
    if(length(grep("favicon.ico", fp)) > 0) {
      # skip binary icon file.
      next
    }
    x <- readLines(fp)
    y <- gsub( "PMMUsingSPLITT", packageName, x )
    z <- gsub( "AbcPMM", class, y)
    cat(z, file=fp, sep="\n")
  }
  
  ## Rename files containing packageName and class in their name:
  for( f in list.files(dp, recursive = TRUE) ){
    fp <- file.path(dp, f)
    if(length(grep("favicon.ico", fp)) > 0) {
      # skip binary icon file.
      next
    }
    fp2 <- fp
    if(length(grep("PMMUsingSPLITT", fp)) > 0) {
      fp2 <- gsub("PMMUsingSPLITT", packageName, fp)
    } 
    
    if(length(grep("AbcPMM", fp)) > 0) {
      fp2 <- gsub("AbcPMM", class, fp2)
    }
    
    if(!identical(fp2, fp)) {
      file.rename(fp, fp2)
    }
  }
  
  repo <- init(file.path(path, packageName))
  
  add(repo, "*")
  status(repo)
}