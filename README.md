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
[1] parallel  stats4    grid      stats     graphics  grDevices utils
[8] datasets  methods   base

other attached packages:
[1] DESeq2_1.6.3            RcppArmadillo_0.4.600.0 Rcpp_0.11.3
[4] GenomicRanges_1.18.4    GenomeInfoDb_1.2.4      IRanges_2.0.1
[7] S4Vectors_0.4.0         BiocGenerics_0.12.1     pdist_1.2
[10] seqLogo_1.32.1          reshape2_1.4.1          gridExtra_0.9.1
[13] ggplot2_1.0.0           sROC_0.1-2              flexclust_1.3-4
[16] modeltools_0.2-21       lattice_0.20-29         shiny_0.10.2.2

loaded via a namespace (and not attached):
[1] acepack_1.3-3.3      annotate_1.44.0      AnnotationDbi_1.28.1
[4] base64enc_0.1-2      BatchJobs_1.5        BBmisc_1.8
[7] Biobase_2.26.0       BiocParallel_1.0.0   brew_1.0-6
[10] checkmate_1.5.1      cluster_1.15.3       codetools_0.2-9
[13] colorspace_1.2-4     DBI_0.3.1            digest_0.6.8
[16] fail_1.2             foreach_1.4.2        foreign_0.8-62
[19] Formula_1.1-2        genefilter_1.48.1    geneplotter_1.44.0  
[22] gtable_0.1.2         Hmisc_3.14-6         htmltools_0.2.6
[25] httpuv_1.3.2         iterators_1.0.7      latticeExtra_0.6-26
[28] locfit_1.5-9.1       MASS_7.3-37          mime_0.2
[31] munsell_0.4.2        nnet_7.3-8           plyr_1.8.1
[34] proto_0.3-10         R6_2.0.1             RColorBrewer_1.1-2  
[37] RJSONIO_1.3-0        rpart_4.1-8          RSQLite_1.0.0
[40] scales_0.2.4         sendmailR_1.2-1      splines_3.1.2
[43] stringr_0.6.2        survival_2.37-7      tools_3.1.2
[46] XML_3.98-1.1         xtable_1.7-4         XVector_0.6.0
```
