# cluster_analysis
Tools for expression based analyses of microorganisms
* clustering
* regulatory motif discovery
* sequence search
* exploratory analysis workflow capture

Coding conventions:
* env global in RData output has the main data items for a study
* globals use . as word separator, e.g. env$meme.data

Dependencies:
* MEME, 4.10.2 (important!)
* R

```R
> sessionInfo()
R version 3.1.2 (2014-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
[1] MASS_7.3-37       flexclust_1.3-4   modeltools_0.2-21 lattice_0.20-29  
[5] shinyBS_0.50      shiny_0.11.0.9001

loaded via a namespace (and not attached):
 [1] BiocGenerics_0.12.1  digest_0.6.8         GenomeInfoDb_1.2.4  
 [4] GenomicRanges_1.18.4 htmltools_0.2.6      httpuv_1.3.2        
 [7] IRanges_2.0.1        mime_0.2             parallel_3.1.2      
[10] R6_2.0.1             Rcpp_0.11.3          S4Vectors_0.4.0     
[13] xtable_1.7-4         XVector_0.6.0       
```
