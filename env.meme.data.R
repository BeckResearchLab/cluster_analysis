library(flexclust)
library(MASS)

source("memeParse.R")

env.meme.data.setup <- function() {
	cat("scanning meme output files to construct meme.data\n")
	files.count <- 0
	# read meme data
	meme.data <- list()
	for (i in 1:length(env$cluster.ensemble@k)) {
		clusts <- clusters(env$cluster.ensemble[[i]])
		meme.data.cluster <- list()
		for (j in 1:env$cluster.ensemble@k[i]) {
			clust <- clusts[clusts==j]
	
			# setup pathes for input
			dir <- dir.k.cluster(env$dir.output, env$cluster.ensemble[[i]]@k, j)
			meme.file <- paste(dir, env$file.meme.txt, sep="/")

			# setup the training set data frame to be validated in memeParse
			training.set <- data.frame(length=env$genes$upstream.seqs[names(clust),"uplength"], row.names=names(clust))
			# remove any NA (i.e. the gene had no upstream sequence because of an overlap)
			training.set <- training.set[!is.na(training.set$length),"length", drop = F]
	
			if (length(training.set$length) > 1) {	
				# load the meme output file
				meme.data.cluster[[as.character(j)]] <- memeParse(meme.file, training.set)
				files.count <- files.count + 1
				if (files.count %% 100 == 0) {
					cat(paste("...processed", files.count, "files\n"))
				}
			} else {
				if (file.exists(meme.file)) {
					stop(
						sprintf("no valid training set for k = %d and cluster = %d, but meme output file was present: %s", env$cluster.ensemble[[i]]@k, j, meme.file)
					)
				}
				warning(
					sprintf("no valid training set for k = %d and cluster = %d", env$cluster.ensemble[[i]]@k, j)
				)
				meme.data.cluster[[as.character(j)]] <- NA
			}
		}
		meme.data[[as.character(env$cluster.ensemble@k[i])]] <- meme.data.cluster
	}
	warnings()
	cat(paste("...processed", files.count, "files total\n"))
	return(meme.data)
}
