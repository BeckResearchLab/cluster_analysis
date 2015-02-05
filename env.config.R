# to be sourced

# metadata
env$study.title <- "5GB1"
# sample
env$file.sample.info <- "sample_info.xls"
env$file.sample.ordering <- "sample_ordering.xls"
env$file.sample.tracks <- "sample_tracks.xls"
# gene
env$file.genes.annotations <- "5G.annotations.xls"
env$file.genes.upstream.seqs <- "5G.upstream.tab"

# data 
env$file.rpkm <- "5G.rpkm.xls"
env$file.counts <- "5G.counts.xls"
env$reference.sample <- "FM23_TR3"

# settings for upstream data
env$upstream.start <- -1
env$upstream.end <- -301

# pathes
env$dir.root <- "/Users/dacb/work/5G/cluster_analysis"
env$dir.output <- "cluster_analysis.dir"
env$dir.motif.plots <- "motif_plots.dir"
env$dir.my.cluster <- "my_cluster.dir"
env$url.prefix <- "http://127.0.0.1:4202"

# clustering
env$flexclust.seed <- 8576
env$flexclust.k.values <- c(seq(10, 150, by <- 10))
env$flexclust.cores <- 8
env$flexclust.nrep <- 8

# motif prediction
env$path.to.meme <- "meme"
env$meme.base.args <- "-dna -maxsize 1000000 -evt 1e10 -minw 6 -maxw 25 -mod zoops -nostatus -text"
env$meme.nmotifs <- 4
env$file.upstream.fa <- "upstream.fa"
env$file.meme.bfile <- "5G.upstream.bfile"
# meme.jobs can be fed to gnu parallel or similar for dispatch
env$file.meme.jobs <- "meme.jobs"
env$file.meme.txt <- "meme.txt"

# environment cache pathes
env$file.env.cache <- "env.cache.RData"
env$file.env.samples.cache <- "env.samples.cache.RData"
env$file.env.genes.cache <- "env.genes.cache.RData"
env$file.env.cluster.ensemble.cache <- "env.cluster.ensemble.cache.RData"
env$file.env.meme.data.cache <- "env.meme.data.cache.RData"
env$file.env.meme.sites.cache <- "env.meme.sites.cache.RData"

# varia, for plotting, etc.
env$motif.colors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")

# database stuff
env$mysql.database <- "cluster_analysis"
env$mysql.log.table <- "log"
env$mysql.log.like.view <- "log_likes"
