## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE, out.width = "792px", out.height = "612px"--------------
knitr::include_graphics("figures/Uml1.png")

## ----create-references, echo=FALSE, include=FALSE, eval=TRUE-------------
treeProcessing <- c("ape")
data <- c("data.table")
poumm <- c("POUMM")
testing <- c("testthat")
boot <- c("boot")
 
packagesUsed <- c(treeProcessing, data, poumm, boot, testing)

printPackages <- function(packs) {
  res <- ""
  for(i in 1:length(packs)) {
    res <- paste0(res, paste0(packs[i], ' v', packageVersion(packs[i]), ' [@R-', packs[i], ']'))
    if(i < length(packs)) {
      res <- paste0(res, ', ')
    }
  }
  res
}

# Write bib information (this line is executed manually and the bib-file is edited manually after that)
#knitr::write_bib(packagesUsed, file = "./REFERENCES-R.bib")

