# to be sourced

library(sROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# data and utilities for plotting cluster profiles
magenta <- "#FF00FF"
cyan <- "#00FFFF"
grey <- "#AAAAAA"
yellow <- "#FFD800"
# buffer to be added to the expression ratio plots
log.ratio.plot.adj <- 1.5
# buffer added to RPKM plots
rpkm.plot.adj <- 100
# get the 95% CI
cdf2cimin <- function(mycdf) {
        mycdf$x[which.min(abs(mycdf$Fhat - 0.05))]
}
cdf2cimax <- function(mycdf) {
        mycdf$x[which.min(abs(mycdf$Fhat - 0.95))]
}

distsum.plot <- function(kdsdf) {
	ggplot(kdsdf, aes(x=k, y=distSum)) +
		geom_bar(stat = "identity") +
		xlab("Input k") +
		ylab("Sum of within cluster distances") +
		theme_bw()
}

distsum.delta.plot <- function(kdsdf) {
	ggplot(kdsdf, aes(x=k, y=distSumDelta)) +
		geom_bar(stat = "identity") +
		xlab("Input k") +
		ylab("Difference of sum of within\ncluster distances to next smaller k") +
		theme_bw()
}

cluster.size.plot <- function(kclust) {
	csdf <- as.data.frame(table(kclust@cluster))
	
	ggplot(csdf, aes(x=Var1, y=Freq)) + 
		geom_bar(stat = "identity") + 
		ylab("Cluster size") +
		xlab("Cluster index") +
		theme_bw()
}

makeMyClusterProfilePlot <- function(log.ratio, genes, memes, displayMotifGeneProfile=F) {
	clustres <- log.ratio[genes,]
	genesdf <- data.frame(t(clustres), Sample=names(clustres))
	genesdf$Sample <- factor(genesdf$Sample, levels=names(clustres))
	mdf <- melt(genesdf, id.vars=c("Sample"))
	if (length(genes) > 4) {	
		cmin<-apply(clustres, 2, min)
		cmean<-apply(clustres, 2, mean)
		cmax<-apply(clustres, 2, max)
		# estimate CDF from data
		cis <- apply(clustres, 2, kCDF)
		# get min at 95% CI
		lwr <- sapply(cis, cdf2cimin)
		# get max at 95% CI
		upr <- sapply(cis, cdf2cimax)
		# cluster data frame
		clustdf <- data.frame(min=cmin, max=cmax, mean=cmean, lwr=lwr, upr=upr, Sample=names(cmean))
		clustdf$Sample <- factor(clustdf$Sample, levels=names(cmean))
		myplot <- ggplot(clustdf, aes(x=Sample)) + 
			geom_ribbon(aes(ymax=upr, ymin=lwr, group=1, alpha=0.7), colour=magenta)
	} else {
		myplot <- ggplot(mdf, aes(x=Sample, y=value, group=variable))
	}
	myplot <- myplot +
			geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow) +
			geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow) +
			ggtitle("Expression profile") +
			ylab("Log2(Sample / T0)\nNormalized counts") +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")
	if (!identical(displayMotifGeneProfile, F) && length(displayMotifGeneProfile) > 0) {
			for (m in 1:length(displayMotifGeneProfile)) {
				motif <- memes[[displayMotifGeneProfile[m]]]
				motifdf <- data.frame(t(clustres[as.character(motif$posns$gene),]), Sample=names(clustres))
				motifdf$Sample <- factor(genesdf$Sample, levels=names(clustres))
				mdf <- melt(motifdf, id.vars=c("Sample"))
				myplot <- myplot + 
					geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]]) +
					geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]])
			}
	}
	return(myplot)
}

makeClusterProfilePlot <- function(env, k, cluster, simple=F, focus=F, displayMotifGeneProfile=F) {
		mcpp_clusts <- clusters(env$cluster.ensemble[[k]])
		mcpp_clust <- mcpp_clusts[mcpp_clusts==cluster]
		clustres <- env$log.ratio[names(mcpp_clust),]
		cmin<-apply(clustres, 2, min)
		cmean<-apply(clustres, 2, mean)
		cmax<-apply(clustres, 2, max)
		if (!simple) {
			# estimate CDF from data
			cis <- apply(clustres, 2, kCDF)
			# get min at 95% CI
			lwr <- sapply(cis, cdf2cimin)
			# get max at 95% CI
			upr <- sapply(cis, cdf2cimax)
			# cluster data frame
			clustdf <- data.frame(min=cmin, max=cmax, mean=cmean, lwr=lwr, upr=upr, Sample=names(cmean))
			clustdf$Sample <- factor(clustdf$Sample, levels=names(cmean))
			myplot <- ggplot(clustdf, aes(x=Sample)) + 
				ggtitle(paste(paste("K =", k, ": Cluster", cluster, paste("(", length(mcpp_clust), " genes)", sep="")), "Expression profile", sep="\n")) +
				ylab("Log2(Sample / T0)\nNormalized counts") +
				geom_ribbon(aes(ymax=upr, ymin=lwr, group=1, alpha=0.7), colour=magenta) +
				theme_bw() +
				theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")
		} else {
			clustdf <- data.frame(min=cmin, max=cmax, mean=cmean, Sample=names(cmean))
			clustdf$Sample <- factor(clustdf$Sample, levels=names(cmean))

			# cluster label for center of plot
			clusterGrob <- grobTree(textGrob(cluster, x=.5, y=.5, gp=gpar(col="black", fontsize=30)))
			clusterGrobGenes <- grobTree(textGrob(paste(length(mcpp_clust), "genes"), x=.5, y=.1, gp=gpar(col="black", fontsize=15)))

			myplot <- ggplot(clustdf, aes(x=Sample)) +
				theme_bw() +
				theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank()) +
  				annotation_custom(clusterGrob) +
  				annotation_custom(clusterGrobGenes) +
				ylim(resmin - resPlotAdj, resmax + resPlotAdj)
		}
		if (!identical(focus, F)) {
			focusdf <- data.frame(t(clustres[focus,]), Sample=names(clustres))
			focusdf$Sample <- factor(focusdf$Sample, levels=names(clustres))
			mdf <- melt(focusdf, id.vars=c("Sample"))
			myplot <- myplot + 
				geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow) +
				geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow)
		}
		if (!identical(displayMotifGeneProfile, F) && length(displayMotifGeneProfile) > 0) {
				memes <- env$meme.data[[k]][[cluster]]
				for (m in 1:length(displayMotifGeneProfile)) {
					motif <- memes[[displayMotifGeneProfile[m]]]
					motifdf <- data.frame(t(clustres[as.character(motif$posns$gene),]), Sample=names(clustres))
					motifdf$Sample <- factor(names(clustres), levels=names(clustres))
					mdf <- melt(motifdf, id.vars=c("Sample"))
					myplot <- myplot + 
						geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]]) +
						geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]])
				}
		}
		myplot +
			geom_point(aes(y=mean), colour=cyan) + 
			geom_line(aes(x=Sample, y=mean, group=1), colour=cyan) + 
			geom_line(aes(x=Sample, y=min, group=1), colour=grey) + 
			geom_line(aes(x=Sample, y=max, group=1), colour=grey)
}

renderMotifPlots <- function(dir, genes, meme.data) {
	dir.create(dir, recursive=T)
	memes <- meme.data
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
	names(motif.colors) <- levels(msc$motif)
	for (g in levels(msc$gene)) {
		mscg <- msc[msc$gene==g,]
		gene_file <- paste(dir, paste(g, "png", sep="."), sep="/")
		png(filename=gene_file, width=180, height=18)
		ul <- env$seqs.upstream[as.character(g),"uplength"]
		sline <- data.frame(x=c(-ul+2, 2), y=c(.5,.5))
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
	return(msc)
}
