# to be sourced

# data loading, computation for caching
source("env.samples.R")
source("env.genes.R")
source("env.cluster.ensemble.R")
source("env.meme.data.R")
source("env.meme.sites.R")

# utility functions
source("utilities.R")

# create a new empty environment to hold all of our data
env <- new.env(parent = emptyenv())

# load the environment configuration
source("env.config.R")

# include.meme: should the meme.data and meme.sites be loaded?
env.assemble <- function(include.meme = T) {
	# reconstruct the environment by loading each piece or executing code
	# if exists load
	if (file.exists(env$file.env.samples.cache)) {
		cat("reading samples cache...\n")
		env$samples <- readRDS(env$file.env.samples.cache)
	} else {
		# or compute
		env$samples <- env.samples.setup(env$file.sample.info, env$file.sample.ordering, env$file.sample.tracks, 
			env$file.rpkm, env$file.counts, env$reference.sample)
		saveRDS(env$samples, env$file.env.samples.cache);
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
		env$cluster.ensemble <- env.cluster.ensemble.setup(env$flexclust.seed, env$samples$log.ratio, env$flexclust.k.values, env$flexclust.nrep, env$flexclust.cores)
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
