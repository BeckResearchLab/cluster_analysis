library(flexclust)
library(MASS)

source("memeIO.R")

load("cluster_analysis.3.RData")

# read meme data
meme.data <- list();
root.dir <- 'cluster_analysis.dir'
for (i in 1:length(env$cluster.ensemble@k)) {
	clusts <- clusters(env$cluster.ensemble[[i]])
	meme.data.cluster <- list()
	for (j in 1:env$cluster.ensemble@k[i]) {
		clust <- clusts[clusts==j]
		clust_seqs_upstream <- env$seqs.upstream[names(clust),]

		# setup pathes for input
		dir <- paste(root.dir, paste("k_", env$cluster.ensemble[[i]]@k, ".dir/cluster_", j, ".dir", sep=""), sep="/")
		meme_file <- paste(dir, "meme.txt", sep="/")
		
		# load the meme output file
		print(meme_file)
		if (file.exists(meme_file)) {
			meme_text <- readLines(meme_file)
			meme.data.cluster[[j]] <- memeParse(meme_text)
		} else {
			warning(paste("meme.txt missing for k =", env$cluster.ensemble[[i]]@k, "cluster =", j))
		}
	}
	meme.data[[i]] <- meme.data.cluster
}
warnings()

env$meme.data <- meme.data

save(env, file="cluster_analysis.4.RData")
