library(flexclust)

# to be sourced
env.meme.sites.setup <- function() {
	cat("assembling meme.sites from meme.data\n");
	meme.sites <- list()
	for (i in 1:length(env$cluster.ensemble@k)) {
		# lookup the cluster data and meme data for this k run
		clusts <- clusters(env$cluster.ensemble[[i]])
		cmemes <- env$meme.data[[i]]
		meme.sites.cluster <- list()
		for (j in 1:env$cluster.ensemble@k[i]) {
			# lookup the cluster specific ids and memes
			clust <- clusts[clusts==j]
			memes <- cmemes[[j]]

			sites <- list()
			mi <- 1
			for (motif in memes) {
				# new data frame from position table, fill in motif and with with rep
				# compute xmin, xmax for easy ggplot viewing
				uplen <- env$genes$upstream.seqs[as.character(motif$posns$gene),"uplength"]
				df <- data.frame(
					cbind(motif$posns),
					motif=rep(mi, length(motif$posns$gene)),
					width=rep(motif$width, length(motif$posns$gene)),
					xmin=motif$posns$start-uplen,
					xmax=motif$posns$start+motif$width-uplen
				)
				sites[[as.character(mi)]] <- df
				mi <- mi + 1
			}
			msc <- do.call("rbind", sites)
			msc$motif <- as.factor(msc$motif)
			meme.sites.cluster[[as.character(j)]] <- msc
		}
		meme.sites[[as.character(env$cluster.ensemble@k[i])]] <- meme.sites.cluster
	}
	return(meme.sites)
}
