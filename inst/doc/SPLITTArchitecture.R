## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- out.width = "792px", out.height = "612px"--------------------------
knitr::include_graphics("figures/UmlDiagram8.pdf")

## ----fig1, dpi=70, eval=FALSE, include=TRUE, results="hide",fig.width=7.2, fig.height=7.2----
#  #newick <- '(((1:4,2:4.5)11:4.3,3:4)10:1.4,(4:10.25,((5:4,6:3.8)14:4.2,7:4.1)13:4.5,(8:6,9:6)15:6)12:2.2)0;'
#  newick <- '(((1:4,2:4.5)11:4.3,3:4)10:3.8,(((5:3,6:2.8)14:4.2,7:4.1)13:4.5,(8:5,9:5)15:6)12:4.2)0;'
#  tree <- read.tree(text=newick)
#  
#  par(mfrow=c(1,5))
#  par(mar=c(0,0,0,0))
#  
#  plotParallelPruningOrder(tree, TRUE)
#  
#  fig.a <- recordPlot()
#  
#  newick <- "(((((1:3.2,2:3.2)10:3.2,3:3.2)11:3.2,4:3.2)12:3.2,5:3.2)13:3.2,6:3.2)14;"
#  tree <- read.tree(text=newick)
#  
#  par(mfrow=c(1,6))
#  par(mar=c(0,0,0,0))
#  plotParallelPruningOrder(tree, TRUE)
#  
#  fig.b <- recordPlot()
#  
#  cowplot::plot_grid(fig.a, fig.b, nrow=2, labels = c("a", "b"))

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

