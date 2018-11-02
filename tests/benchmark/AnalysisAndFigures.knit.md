---
output: 
  bookdown::pdf_document2:
    keep_tex: true
    fig_caption: true
    latex_engine: pdflatex
    template: template-sysbio.tex
thanks: ""
geometry: margin=1in
fontfamily: mathpazo
fontsize: 11pt
graphics: true
tables: true
---




<!-- Analyse benchmark data -->

```r
load("../../tests/benchmark/Trees.RData")
library(ggtree)

pTypes <- c(p0.5="p = 0.5 (balanced)", 
            p0.1="p = 0.1", 
            p0.01="p = 0.01", 
            p.c99="p = 0.01/N (ladder)")

treeList <- lapply(
  trees[N==1000 & num_traits == 1, list(tree), keyby=pSymb][list(names(pTypes))]$tree, 
  function(tr) {
    tr$edge.length[] <- 1
    ladderize(tr)
  })

class(treeList) <- "multiPhylo"

ggtree(treeList, size = 0.1) + facet_wrap(~factor(sapply(.id, function(i) pTypes[i]), levels = pTypes), scales="free", ncol=4)
```




```r
ggplot(times1[expr != "gc" & 
               type != "binary.with.star" &
               order != "split" & 
               order != "visit-ranges" &
               !sapply(impl, startsWith, "abc") &
               #!sapply(as.character(expr), endsWith, "split") &
               !mode %in% c("auto", "hybrid")
             ]) +
  geom_line(aes(x=N, y=mean, col=Implementation, linetype=Mode, shape=Order), size=0.3) +
  geom_point(aes(x=N, y=mean, col=Implementation, linetype=Mode, shape=Order), size = 0.9) +
  scale_shape_manual(values = c(`postorder`=20, `queue-based`=2, `range-based`=4, `prune-ranges`=22)) + 
  scale_y_continuous(limits = c(0.02, 4000), 
                     trans=log10_trans(), 
                     breaks = c(seq(0.02, 0.1, by=.04), 
                                seq(0.2, 1, by=.4), 
                                seq(2, 10, by=4), 
                                seq(20, 100, by=40),
                                seq(200, 1000, by=400), 
                                seq(2000, 10000, by=4000)),
                     minor_breaks = c(seq(0.04, 0.1, by=.4), 
                                      seq(0.4, 1, by=.4), 
                                      1.5,
                                      seq(4, 10, by=4), 
                                      15,
                                      seq(40, 100, by=40),
                                      150,
                                      seq(400, 1000, by=400),
                                      1500,
                                      seq(4000, 10000, by=4000),
                                      15000),
                     labels = comma) + 
  scale_x_continuous(limits = c(50, 2.6e5),
                     trans=log10_trans(), 
                     breaks = 10^(2:5), 
                     minor_breaks = NULL,
                     labels = comma) + 
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=2, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Time [ms]") +
  facet_grid(.~pUnbalanced)
```

```
## Warning: Ignoring unknown aesthetics: shape
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

![Likelihood calculation times for R and C++ implementations of the pruning algorithm.](AnalysisAndFigures_files/figure-latex/fig-plot-times-1.pdf) 


```r
ggplot(times[
  Computer == "MacBookPro" &
    Package %in% c("Rphylopars (C++)", "PCMBaseCpp (C++)") & 
    Mode %in% c("serial", "parallel range")]) + 
  geom_line(
    aes(x=N, 
        y=time, 
        col=Package, 
        #shape=as.factor(k),
        linetype=Mode
        ), size=0.3) +
  geom_point(
    aes(x=N, 
        y=time, 
        col=Package, 
        #shape=as.factor(k),
        linetype=Mode
        ), size = 0.9) +
  scale_y_continuous(limits = c(0.6, 40000), 
                     trans=log10_trans(), 
                     breaks = c(seq(0.02, 0.1, by=.04), 
                                seq(0.2, 1, by=.4), 
                                seq(2, 10, by=4), 
                                seq(20, 100, by=40),
                                seq(200, 1000, by=400), 
                                seq(2000, 10000, by=4000),
                                seq(20000, 100000, by=40000)),
                     minor_breaks = c(seq(0.04, 0.1, by=.4), 
                                      seq(0.4, 1, by=.4), 
                                      1.5,
                                      seq(4, 10, by=4), 
                                      15,
                                      seq(40, 100, by=40),
                                      150,
                                      seq(400, 1000, by=400),
                                      1500,
                                      seq(4000, 10000, by=4000),
                                      15000,
                                      seq(40000, 100000, by=40000),
                                      150000),
                     labels = comma) + 
  scale_x_continuous(limits = c(50, 2.6e5),
                     trans=log10_trans(), 
                     breaks = 10^(2:5), 
                     minor_breaks = NULL,
                     labels = comma) + 
  scale_color_discrete(name = "Implementation") +
  scale_shape_discrete(name = "Number of traits") +
  scale_linetype_manual(values = c(serial = 1, `parallel range` = 2)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(ncol=2, order=1), 
         #shape=guide_legend(ncol=2, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  ylab("Time [ms]") +
  facet_grid(factor(k, levels = k, labels=(paste0("k=", k, " traits")))~pUnbalanced)
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

```
## Warning: Removed 9 rows containing missing values (geom_point).
```

![Multivariate likelihood calculation times on MacBookPro](AnalysisAndFigures_files/figure-latex/fig-plot-times-multivariate-MacBookPro-1.pdf) 



```r
ggplot(times[
  Computer == "Euler" &
  #package=="PCMBaseCpp (C++)" &
  package=="POUMM (C++)" &
               !(Mode %in% c("parallel queue", 
                             "hybrid range", 
                             "parallel range except"
                             )) &
               #!pUnbalanced %in% c("p = 0.01/N (ladder)") &
               k == 1]) +
  geom_line(aes(x=omp_max_threads, y=time, 
                linetype=Mode), size = 0.3) +
  geom_point(aes(x=omp_max_threads, y=time, 
                 linetype=Mode), size = 0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma,
                     limits = c(0, NA)) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 24),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Time [ms]") +
  facet_wrap(paste0("N = ",N)~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-times-euler-1-1.pdf) 


```r
pl <- ggplot(times[
  Computer == "Euler" &
  #package=="PCMBaseCpp (C++)" &
  package=="PCMBaseCpp (C++)" &
               !(Mode %in% c(#"parallel queue", 
                             "hybrid range", 
                             "parallel range except"
                             )) &
               #!pUnbalanced %in% c("p = 0.01/N (ladder)") &
               k == 1]) +
  geom_line(aes(x=omp_max_threads, y=time, 
                linetype=Mode), size = 0.3) +
  geom_point(aes(x=omp_max_threads, y=time, 
                 linetype=Mode), size = 0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma,
                     limits = c(0, NA)) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 24),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range` = 2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Time [ms]") +
  facet_wrap(paste0("N = ",N)~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```



```r
ggplot(times[Computer == "Euler" & 
               k==1 & Package=="POUMM (C++)" & 
               Mode %in% c("serial", "parallel range", "parallel queue"),
             list(Mode, speedup=c(1, time[1]/time[-1])), 
             keyby=list(pUnbalanced, Algorithm, N, k, omp_max_threads)]) +
  geom_abline(aes(slope = 1, intercept = 0), linetype=1, size=0.3, color="darkgrey") +
  geom_abline(aes(slope = 0.5, intercept = 0), linetype=1, size=0.3, color="red") +
  geom_line(aes(x=omp_max_threads, y=speedup, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=speedup, 
                 linetype=Mode), size=0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma,
                     limits = c(0, 12),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 24),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range`=2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Parallel speedup [x]") +
  facet_wrap(paste0("N = ", N)~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-speedup-euler-1-1.pdf) 


```r
ggplot(times[Computer == "Euler" & 
               k==1 & Package=="PCMBaseCpp (C++)" & 
               Mode %in% c("serial", "parallel range", "parallel queue"),
             list(Mode, speedup=c(1, time[1]/time[-1])), 
             keyby=list(pUnbalanced, Algorithm, N, k, omp_max_threads)]) +
  geom_abline(aes(slope = 1, intercept = 0), linetype=1, size=0.3, color="darkgrey") +
  geom_abline(aes(slope = 0.5, intercept = 0), linetype=1, size=0.3, color="red") +
  geom_line(aes(x=omp_max_threads, y=speedup, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=speedup, 
                 linetype=Mode), size=0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma,
                     limits = c(0, 12),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 24),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range`=2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Parallel speedup [x]") +
  facet_wrap(paste0("N = ", N)~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

![Parallel speedup, Euler multivariate 1 trait.](AnalysisAndFigures_files/figure-latex/fig-plot-speedup-euler-multi-1-1.pdf) 


```r
num_traits <- 4

ggplot(times[
  Computer == "Euler" &
  package %in% c("PCMBaseCpp (C++)", "Rphylopars (C++)") &
               !(Mode %in% c(#"parallel queue", 
                             "hybrid range", 
                             "parallel range except"
                             )) &
               #!pUnbalanced %in% c("p = 0.01/N (ladder)") &
               k == num_traits]) +
  geom_line(aes(x=omp_max_threads, y=time, 
  #geom_line(aes(x=nCores, y=time, 
                #col=Package, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=time, 
  #geom_point(aes(x=nCores, y=time, 
                 #col=Package, 
                 linetype=Mode), size = 0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma, 
                     limits = c(0, NA)) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 20),
                     breaks = c(1, 2, seq(4, 20, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range` = 2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Time [ms]") +
  facet_wrap(N~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

```
## Warning: Removed 3 rows containing missing values (geom_path).
```

```
## Warning: Removed 48 rows containing missing values (geom_point).
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-times-euler-4-1.pdf) 


```r
num_traits <- 4
ggplot(times[
  Computer == "Euler" &
    k==num_traits & 
    Package=="PCMBaseCpp (C++)" & 
    Mode %in% c("serial", "parallel range", "parallel queue"), 
  list(Mode, speedup=c(1, time[1]/time[-1])), 
  keyby=list(pUnbalanced, Algorithm, N, k, omp_max_threads)]) +
  geom_abline(aes(slope = 1, intercept = 0), linetype=1, size=0.3, color="darkgrey") +
  geom_abline(aes(slope = 0.5, intercept = 0), linetype=1, size=0.3, color="red") +
  geom_line(aes(x=omp_max_threads, y=speedup, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=speedup, 
                 linetype=Mode), size=0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma,
                     limits = c(0, 12),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 24),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range`=2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Parallel speedup [x]") +
  facet_wrap(paste0("N = ", N)~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-speedup-euler-4-1.pdf) 



```r
num_traits <- 8

ggplot(times[
  Computer == "Euler" & 
  package %in% c("PCMBaseCpp (C++)", "Rphylopars (C++)") &
               !(Mode %in% c(#"parallel queue", 
                             "hybrid range", 
                             "parallel range except"
                             )) &
               #!pUnbalanced %in% c("p = 0.01/N (ladder)") &
               k == num_traits]) +
  geom_line(aes(x=omp_max_threads, y=time, 
  #geom_line(aes(x=nCores, y=time, 
                #col=Package, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=time, 
  #geom_point(aes(x=nCores, y=time, 
                 #col=Package, 
                 linetype=Mode), size=0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma, 
                     limits = c(0, NA)) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 20),
                     breaks = c(1, 2, seq(4, 20, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range` = 2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Time [ms]") +
  facet_wrap(N~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

```
## Warning: Removed 3 rows containing missing values (geom_path).
```

```
## Warning: Removed 48 rows containing missing values (geom_point).
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-times-euler-8-1.pdf) 


```r
num_traits <- 8
ggplot(times[
  Computer == "Euler" & 
    k==num_traits & 
    Package=="PCMBaseCpp (C++)" & 
    Mode %in% c("serial", "parallel range", "parallel queue"), 
  list(Mode, speedup=c(1, time[1]/time[-1])), 
  keyby=list(pUnbalanced, Algorithm, N, k, omp_max_threads)]) +
  geom_abline(aes(slope = 1, intercept = 0), linetype=1, size=0.3, color="darkgrey") +
  geom_abline(aes(slope = 0.5, intercept = 0), linetype=1, size=0.3, color="red") +
  geom_line(aes(x=omp_max_threads, y=speedup, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=speedup, 
                 linetype=Mode), size=0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma,
                     limits = c(0, 12),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 24),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range`=2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Parallel speedup [x]") +
  facet_wrap(paste0("N = ", N)~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-speedup-euler-8-1.pdf) 


```r
num_traits <- 16

ggplot(times[
  Computer == "Euler" &
    package %in% c("PCMBaseCpp (C++)", "Rphylopars (C++)") &
    !(Mode %in% c(#"parallel queue", 
      "hybrid range", 
      "parallel range except"
    )) &
    #!pUnbalanced %in% c("p = 0.01/N (ladder)") &
    k == num_traits]) +
  geom_line(aes(x=omp_max_threads, y=time, 
                #geom_line(aes(x=nCores, y=time, 
                #col=Package, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=time, 
                 #geom_point(aes(x=nCores, y=time, 
                 #col=Package, 
                 linetype=Mode), size = 0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma, 
                     limits = c(0, NA)) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 20),
                     breaks = c(1, 2, seq(4, 20, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range` = 2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Time [ms]") +
  facet_wrap(N~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

```
## Warning: Removed 3 rows containing missing values (geom_path).
```

```
## Warning: Removed 48 rows containing missing values (geom_point).
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-times-euler-16-1.pdf) 


```r
num_traits <- 16
ggplot(times[
  Computer == "Euler" & 
    k==num_traits & 
    Package=="PCMBaseCpp (C++)" & 
    Mode %in% c("serial", "parallel range", "parallel queue"), 
  list(Mode, speedup=c(1, time[1]/time[-1])), 
  keyby=list(pUnbalanced, Algorithm, N, k, omp_max_threads)]) +
  geom_abline(aes(slope = 1, intercept = 0), linetype=1, size=0.3, color="darkgrey") +
  geom_abline(aes(slope = 0.5, intercept = 0), linetype=1, size=0.3, color="red") +
  geom_line(aes(x=omp_max_threads, y=speedup, 
                linetype=Mode), size=0.3) +
  geom_point(aes(x=omp_max_threads, y=speedup, 
                 linetype=Mode), size=0.6) +
  #scale_shape_manual(values = c(`Quadratic polynomial`=20, `3-point`=2)) + 
  scale_y_continuous(labels = comma,
                     limits = c(0, 12),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL) + 
  scale_x_continuous(name = "Number of CPU cores",
                     limits = c(1, 24),
                     breaks = c(1, 2, seq(4, 24, by=4)),
                     minor_breaks = NULL,
                     labels = comma) +
  scale_linetype_manual(values = c(serial = 1, `parallel range`=2, `parallel queue`=3)) +
  ggplot2::theme_grey() + 
  guides(linetype=guide_legend(order=1), 
         shape=guide_legend(ncol=1, order=2),
         col=guide_legend(ncol=2, order=3)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.grid = element_line(size = 0.4),
        strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        line = element_line(size=0.6)) +
  
  ylab("Parallel speedup [x]") +
  facet_wrap(paste0("N = ", N)~pUnbalanced, scales="free")
```

```
## Warning: Ignoring unknown aesthetics: linetype
```

![Likelihood calculation times for R and C++ implementations of multivariate OU caclulation.](AnalysisAndFigures_files/figure-latex/fig-plot-speedup-euler-16-1.pdf) 



```r
tableCreateCache <- times1[type=="binary" & 
                             !is.na(timeCreateCache) & 
                             !sapply(Implementation, function(.) startsWith(as.character(.), 'ln')) & 
                             !Implementation=="diversitree (R)", 
                           list(timeCreateCache=min(timeCreateCache)), 
                          keyby=list(N, pSymb, Implementation)]

table <- rbindlist(lapply(10^{2:5}, function(NN) 
  tableCreateCache[N==NN, list(N=NN, Implementation=Implementation[pSymb=="p0.5"], 
                                `p=0.5`=.SD[pSymb=="p0.5", timeCreateCache*1000],
                                `p=0.1`=.SD[pSymb=="p0.1", timeCreateCache*1000],
                                `p=0.01`=.SD[pSymb=="p0.01", timeCreateCache*1000],
                                `p=0.01/N`=.SD[pSymb=="p.c99", timeCreateCache*1000]
                                )]
  
))

require(xtable)
xt <- xtable(table)
options(digits=6)
print(xt)
```
