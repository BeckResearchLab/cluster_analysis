library(flexclust)
library(ggplot2)

load("cluster_analysis.4.RData")

motif.colors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")

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

		# setup output path
		dir <- paste(env$dir.output, paste("k_", env$cluster.ensemble[[i]]@k, ".dir/cluster_", j, ".dir/motif_plots.dir", sep=""), sep="/")
		dir.create(dir, recursive=T)

		sites <- list()
		mi <- 1
		for (motif in memes) {
			# new data frame from position table, fill in motif and with with rep
			# compute xmin, xmax for easy ggplot viewing
			uplen <- env$seqs.upstream[as.character(motif$posns$gene),"uplength"]
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
		names(motif.colors) <- levels(msc$motif)
		for (g in levels(msc$gene)) {
			mscg <- msc[msc$gene==g,]
			gene_file <- paste(dir, paste(g, "png", sep="."), sep="/")
			print(gene_file)
			png(filename=gene_file, width=180, height=18)
			ul <- env$seqs.upstream[as.character(g),"uplength"]
			sline <- data.frame(x=c(-ul+2, 2), y=c(.5,.5))
			#print(sline)
			print(
				ggplot(mscg) +
					scale_fill_manual(name="nmotif", values=motif.colors) +
					geom_rect(aes(xmin=xmin, xmax=xmax, ymin=0, ymax=1, fill=motif, group=1)) +
					geom_line(data=sline, aes(x=x, y=y, group=1), color="#000000") +
					theme(
						legend.position = "none",
						panel.background = element_blank(),
						panel.grid.major = element_blank(),
						panel.grid.minor = element_blank(),
						panel.margin = unit(0,"null"),
						plot.margin = rep(unit(0,"null"),4),
						axis.ticks = element_blank(),
						axis.text.x = element_blank(),
						axis.text.y = element_blank(),
						axis.title.x = element_blank(),
						axis.title.y = element_blank(),
						axis.ticks.length = unit(0,"null"),
						axis.ticks.margin = unit(0,"null")
					) +
					labs(x=NULL, y=NULL) + 
					xlim(-150,2) + ylim(0,1)
			)
			dev.off()
		}
	}
	meme.sites[[as.character(env$cluster.ensemble@k[i])]] <- meme.sites.cluster
}
env$meme.sites <- meme.sites

save(env, file="cluster_analysis.5.RData")
