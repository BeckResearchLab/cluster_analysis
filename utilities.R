# to be sourced

# return the directory path for a k / cluster pair
dir.k.cluster <- function(dir.output, k, cluster, make.dir = F) {
	dir <- sprintf("%s/k_%d.dir/cluster_%d.dir", dir.output, as.integer(k), as.integer(cluster))
	# equivalent to to mkdir -p
	if (identical(make.dir, T)) {
		dir.create(dir, recursive = T, showWarnings = F)
	}
	return(dir)
}

# the sum of member to centroid distances
get.distsum <- function() {
	x <- env$cluster.ensemble
	X <- x@k
	Y <- sapply(x@models, function(z) info(z, "distsum"))
	Z <- c(NA, diff(Y))
	df = data.frame(k=X, distsum=Y, distsum.delta=Z)
	return(df);
}

# take a set of motifs and process those lists into a single data frame
# with the sites that is easy to plot
meme.positions.to.sites <- function(motifs, upstream.seqs, upstream.start) {
	sites <- list()
	mi <- 1
	for (motif in motifs) {
		# new data frame from position table, fill in motif and with with rep
		# compute xmin, xmax for easy ggplot viewing
		uplen <- upstream.seqs[as.character(motif$positions$gene),"uplength"]
		df <- data.frame(
			cbind(motif$positions),
			motif = rep(mi, length(motif$positions$gene)),
			width = rep(motif$width, length(motif$positions$gene)),
			xmin = upstream.start - (uplen - motif$positions$start),
			xmax = upstream.start - (uplen - (motif$positions$start + motif$width)) - 1
		)
		sites[[as.character(mi)]] <- df
		mi <- mi + 1
	}
	msc <- do.call("rbind", sites)
	msc$motif <- as.factor(msc$motif)
	return(msc)
}

test.cluster.profile.plot <- function() {
	input <- list()
	input$k <- 10
	input$cluster <- 10
	clusts <- clusters(env$cluster.ensemble[[input$k]])
	clust <- clusts[clusts==input$cluster]
	rowFocus <- F
	cl <- clust
	profile.data <- env$samples$log.ratio[names(cl),]
	makeClusterProfilePlot(profile.data = profile.data,
		title = sprintf("K = %d : Cluster %d (%d genes)\nExpression profile",
			env$cluster.ensemble[[input$k]]@k, as.integer(input$cluster), length(names(cl))
		),
		focus = rowFocus,
		display.motif.gene.profile = c(1),
		motifs = env$meme.data[[input$k]][[input$cluster]],
		motif.colors = env$motif.colors,
		display.tracks = c("FAME"),
		tracks = env$samples$tracks
	)
}
