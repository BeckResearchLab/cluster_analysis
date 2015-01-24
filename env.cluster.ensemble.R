library(flexclust)

# to be sourced
env.cluster.ensemble.setup <- function(flexclust.seed, input.data, flexclust.k.values, flexclust.nrep, flexclust.cores) {
	cat(paste("running flexclust on", flexclust.cores, "cores with", flexclust.nrep, "replicates\n"))
	options(mc.cores=flexclust.cores)
	cluster.ensemble <- stepFlexclust(input.data, flexclust.k.values, nrep=flexclust.nrep, save.data=TRUE, drop=FALSE, verbose=TRUE, seed=flexclust.seed, multicore=T)
	return(cluster.ensemble)
}
