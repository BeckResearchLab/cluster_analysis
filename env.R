source("env.samples.R")
source("env.genes.R")
source("env.cluster.ensemble.R")
source("env.meme.data.R")
source("env.meme.sites.R")

source("utilities.R")

env <- new.env(parent = emptyenv())

env$study.title <- "5GB1"
env$file.sample.info <- "sample_info.xls"
env$file.sample.ordering <- "sample_ordering.xls"
env$file.sample.tracks <- "sample_tracks.xls"
env$file.rpkm <- "5G.rpkm.xls"
env$file.counts <- "5G.counts.xls"
env$reference.sample <- "FM23_TR3"
env$file.genes.annotations <- "5G.annotations.xls"
env$file.genes.upstream.seqs <- "5G.upstream.tab"
env$upstream.start <- -1
env$upstream.end <- -301
env$dir.output <- "cluster_analysis.dir"
env$flexclust.seed <- 8576
env$flexclust.k.values <- c(seq(10, 150, by <- 10))
env$flexclust.cores <- 8
env$flexclust.nrep <- 8
env$file.upstream.fa <- "upstream.fa"
env$path.to.meme <- "meme"
env$meme.base.args <- "-dna -maxsize 1000000 -evt 1e10 -minw 6 -maxw 25 -mod zoops -nostatus -text"
env$meme.nmotifs <- 4
env$file.meme.bfile <- "5G.upstream.bfile"
env$file.meme.jobs <- "meme.jobs"
env$file.meme.txt <- "meme.txt"
env$motif.colors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")
env$file.env.cache <- "env.cache.RData"
env$file.env.sample.cache <- "env.sample.cache.RData"
env$file.env.genes.cache <- "env.genes.cache.RData"
env$file.env.cluster.ensemble.cache <- "env.cluster.ensemble.cache.RData"
env$file.env.meme.data.cache <- "env.meme.data.cache.RData"
env$file.env.meme.sites.cache <- "env.meme.sites.cache.RData"

# include.meme: should the meme.data and meme.sites be loaded?
env.assemble <- function(include.meme = T) {
	# reconstruct the environment by loading each piece or executing code
	# if exists load
	if (file.exists(env$file.env.sample.cache)) {
		cat("reading samples cache...\n")
		env$samples <- readRDS(env$file.env.sample.cache)
	} else {
		# or compute
		env$samples <- env.samples.setup(env$file.sample.info, env$file.sample.ordering, env$file.sample.tracks, 
			env$file.rpkm, env$file.counts, env$reference.sample)
		saveRDS(env$samples, env$file.env.sample.cache);
	}
	# repeat as above for gene data
	if (file.exists(env$file.env.genes.cache)) {
		cat("reading genes cache...\n")
		env$genes <- readRDS(env$file.env.genes.cache)
	} else {
		env$genes <- env.genes.setup(env$file.genes.annotations, env$file.genes.upstream.seqs)
		saveRDS(env$genes, env$file.env.genes.cache)
	}
	# repeat as above for clustering data
	if (file.exists(env$file.env.cluster.ensemble.cache)) {
		cat("reading cluster ensemble cache...\n")
		env$cluster.ensemble <- readRDS(env$file.env.cluster.ensemble.cache)
		} else {
		env$cluster.ensemble <- env.cluster.ensemble.setup(flexclust.seed, env$samples$log.ratio, env$flexclust.k.values, env$flexclust.nrep, env$flexclust.cores)
	saveRDS(env$cluster.ensemble, env$file.env.cluster.ensemble.cache)
	}
	if (include.meme) {
		# repeat as above for meme data
		if (file.exists(env$file.env.meme.data.cache)) {
			cat("reading meme data cache...\n")
			env$meme.data <- readRDS(env$file.env.meme.data.cache)
		} else {
			env$meme.data <- env.meme.data.setup()
			saveRDS(env$meme.data, env$file.env.meme.data.cache)
		}
		# repeat as above for meme sites
		if (file.exists(env$file.env.meme.sites.cache)) {
		cat("reading meme sites cache...\n")
			env$meme.sites <- readRDS(env$file.env.meme.sites.cache)
		} else {
			env$meme.sites <- env.meme.sites.setup()
			saveRDS(env$meme.sites, env$file.env.meme.sites.cache)
		}
	}

	# sanity check the environment
	if (is.null(env$samples)) {
		warning()
		stop("env$samples is NULL! cannot continue")
	} else if (is.null(env$cluster.ensemble)) {
		warning()
		stop("env$cluster.ensemble is NULL! cannot continue")
	} else if (include.meme && is.null(env$meme.data)) {
		warning()
		stop("env$meme.data is NULL! cannot continue")
	} else if (include.meme && is.null(env$meme.sites)) {
		warning()
		stop("env$meme.sites is NULL! cannot continue")
	}
}
