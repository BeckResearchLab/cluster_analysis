library(flexclust)

# to be sourced
env.meme.sites.setup <- function() {
	cat("assembling meme.sites from meme.data\n");
	meme.sites <- list()
	for (i in 1:length(env$cluster.ensemble@k)) {
		# lookup the cluster data and meme data for this k run
		clusts <- clusters(env$cluster.ensemble[[i]])
		cmemes <- env$meme.data[[as.character(env$cluster.ensemble@k[i])]]
		# no meme data for this k?
		if (length(cmemes) <= 0) {
			stop(
				sprintf("no valid meme data for k = %d", env$cluster.ensemble@k[i])
			)
		}

		meme.sites.cluster <- list()
		for (j in 1:env$cluster.ensemble@k[i]) {
			# lookup the cluster specific ids and memes
			clust <- clusts[clusts==j]
	
			# get the list of memes for this cluster	
			memes <- cmemes[[j]]

			# empty list of motifs which can happen if the training set is too small
			if (is.na(memes) || length(memes) <= 0) {
				warning(
					sprintf("no valid motifs predicted for k = %d, cluster = %d", env$cluster.ensemble@k[i], j)
				)
				next
			}

			meme.sites.cluster[[as.character(j)]] <- meme.positions.to.sites(memes, env$genes$upstream.seqs, env$upstream.start)
		}
		meme.sites[[as.character(env$cluster.ensemble@k[i])]] <- meme.sites.cluster
	}
	return(meme.sites)
}
