# to be sourced

library(sROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

source("utilities.R")

# data and utilities for plotting cluster profiles
magenta <- "#FF00FF"
cyan <- "#00FFFF"
grey <- "#AAAAAA"
yellow <- "#FFD800"
# get the 95% CI
cdf2cimin <- function(mycdf) {
        mycdf$x[which.min(abs(mycdf$Fhat - 0.05))]
}
cdf2cimax <- function(mycdf) {
        mycdf$x[which.min(abs(mycdf$Fhat - 0.95))]
}

distsum.plot <- function(kdsdf) {
	ggplot(kdsdf, aes(x = k, y = distsum)) +
		geom_bar(stat = "identity") +
		xlab("Input k") +
		ylab("Sum of within cluster distances") +
		theme_bw()
}

distsum.delta.plot <- function(kdsdf) {
	ggplot(kdsdf, aes(x = k, y = distsum.delta)) +
		geom_bar(stat = "identity") +
		xlab("Input k") +
		ylab("Difference of sum of within\ncluster distances to next smaller k") +
		theme_bw()
}

cluster.size.plot <- function(kclust) {
	csdf <- as.data.frame(table(kclust@cluster))
	
	ggplot(csdf, aes(x = Var1, y = Freq)) + 
		geom_bar(stat = "identity") + 
		ylab("Cluster size") +
		xlab("Cluster index") +
		theme_bw()
}

makeClusterProfilePlot <- function(profile.data, title, y.range.adj = 0, simple = F, focus = F, 
			display.motif.gene.profile = F, motifs = F, motif.colors = F, 
			display.tracks = F, tracks = F, alt.sample.names = F) {
		pd.min <- min(profile.data)
		pd.max <- max(profile.data)
		cmin<-apply(profile.data, 2, min)
		cmean<-apply(profile.data, 2, mean)
		cmax<-apply(profile.data, 2, max)
		if (identical(simple, F)) {
			if (dim(profile.data)[1] > 4) {
				# estimate CDF from data
				cis <- apply(profile.data, 2, kCDF)
				# get min at 95% CI
				lwr <- sapply(cis, cdf2cimin)
				# get max at 95% CI
				upr <- sapply(cis, cdf2cimax)
			} else {
				lwr <- NA
				upr <- NA
			}
			# cluster data frame
			clustdf <- data.frame(min = cmin, max = cmax, mean = cmean, lwr = lwr, upr = upr, Sample = names(cmean))
			clustdf$Sample <- factor(clustdf$Sample, levels = names(cmean))
			myplot <- ggplot(clustdf, aes(x = Sample))
			if (!is.null(title)) {
				myplot <- myplot + ggtitle(title)
			}
			myplot <- myplot +
				ylab("Log2(Sample / T0)\nNormalized counts")
			if (!is.na(lwr) && !is.na(upr)) {
				myplot <- myplot +
					geom_ribbon(aes(ymax = upr, ymin = lwr, group = 1, alpha=0.7), colour = magenta)
			}
			myplot <- myplot +
				theme_bw() +
				theme(axis.text.x = element_text(angle = 90, hjust = 1), 
					legend.position = "none",
					plot.margin = unit(c(0,0.5,0,0), "lines")
				)
			if (!identical(alt.sample.names, F)) {
				myplot <- myplot + scale_x_discrete(labels=alt.sample.names)
			}
		} else {
			clustdf <- data.frame(min = cmin, max = cmax, mean = cmean, Sample = names(cmean))
			clustdf$Sample <- factor(clustdf$Sample, levels = names(cmean))

			# cluster label for center of plot
			clusterGrob <- grobTree(textGrob(title, x = .5, y = .5, gp = gpar(col = "black", fontsize=30)))
			clusterGrobGenes <- grobTree(textGrob(paste(length(profile.data[,1]), "genes"), x = .5, y = .1, gp = gpar(col = "black", fontsize = 15)))

			myplot <- ggplot(clustdf, aes(x = Sample)) +
				theme_bw() +
				theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank()) +
  				annotation_custom(clusterGrob) +
  				annotation_custom(clusterGrobGenes) +
				ylim(pd.min - y.range.adj, pd.max + y.range.adj)
		}
		if (!identical(focus, F)) {
			focusdf <- data.frame(t(profile.data[focus,]), Sample = names(profile.data))
			focusdf$Sample <- factor(focusdf$Sample, levels = names(profile.data))
			mdf <- melt(focusdf, id.vars = c("Sample"))
			myplot <- myplot + 
				geom_point(data = mdf, aes(x = Sample, y = value, group = variable), color = yellow) +
				geom_line(data = mdf, aes(x = Sample, y = value, group = variable), color = yellow)
		}
	
		# for each motif in display.motif.gene.profile vector use the provided motif data to show profiles
		# using the colors in motif.colors
		# only attempt after some basic binary checks on variables
		if (!identical(display.motif.gene.profile, F) && length(display.motif.gene.profile) > 0 
					&& !identical(motifs, F) && length(motifs) > 0 
					&& !identical(motif.colors, F) && length(motif.colors) > 0) {
				for (m in 1:length(display.motif.gene.profile)) {
					motif <- motifs[[display.motif.gene.profile[m]]]
					motifdf <- data.frame(t(profile.data[as.character(motif$positions$gene),]), Sample = names(profile.data))
					motifdf$Sample <- factor(names(profile.data), levels = names(profile.data))
					mdf <- melt(motifdf, id.vars = c("Sample"))
					myplot <- myplot + 
						geom_point(data = mdf, 
							aes(x = Sample, y = value, group = variable), 
							color = motif.colors[display.motif.gene.profile[m]]
						) +
						geom_line(data = mdf, 
							aes(x = Sample, y = value, group = variable), 
							color = motif.colors[display.motif.gene.profile[m]]
						)
				}
		}

		# show min, mean and max
		myplot <- myplot +
			geom_point(aes(y = mean), colour = cyan) + 
			geom_line(aes(x = Sample, y = mean, group = 1), colour = cyan) + 
			geom_line(aes(x = Sample, y = min, group = 1), colour = grey) + 
			geom_line(aes(x = Sample, y = max, group = 1), colour = grey)

		# for each track in display.tracks plot the track above the profile
		if (identical(simple, F) && !identical(display.tracks, F) && length(display.tracks) > 0
				&& !identical(tracks, F) && length(tracks) > 0) {
			track.data <- data.frame(cbind(tracks[env$samples$ordering, display.tracks]), row.names=env$samples$ordering)
			names(track.data) <- display.tracks
			track.data$Sample <- factor(rownames(track.data), levels = rownames(track.data))
			#print(head(track.data))
			trackplot <- ggplot(track.data, aes(x = Sample))
			for (y in display.tracks) {
				trackplot <- trackplot + geom_bar(stat = "identity", aes_string(y = y))
			}
			trackplot <- trackplot +
				theme_bw() +
				theme(axis.title.x = element_blank(), 
					axis.ticks = element_blank(), 
					axis.text.x = element_blank(),
					plot.margin = unit(c(0,0.5,0,0.05), "lines")
				)
			ggplot_gtable(ggplot_build(trackplot))
			ggplot_gtable(ggplot_build(myplot))
			#max.width = unit.pmax(trackplot$widths[2:3], myplot$widths[2:3])
			#trackplot$widths[2:3] <- max.width
			#myplot$widths[2:3] <- max.width
			return(grid.arrange(trackplot, myplot, heights=c(1,5)))
		}
		return(myplot)
}

renderMotifPlots <- function(dir, genes, upstream.seqs, upstream.start, upstream.end, motifs, motif.colors, msc = F) {
	dir.create(dir, recursive = T, showWarnings = F)
	if (identical(msc, F)) {
		msc <- meme.positions.to.sites(motifs, upstream.seqs, upstream.start)
	}
	names(motif.colors) <- levels(msc$motif)
	for (g in levels(msc$gene)) {
		mscg <- msc[msc$gene == g,]
		gene_file <- paste(dir, paste(g, "png", sep = "."), sep = "/")
		png(filename = gene_file, width = 180, height = 18)
		ul <- upstream.seqs[as.character(g),"uplength"]
		sline <- data.frame(x = c(upstream.start - ul, upstream.start), y = c(.5,.5))
		print(
			ggplot(mscg) +
				scale_fill_manual(name = "nmotif", values = motif.colors) +
				geom_rect(aes(xmin = xmin, xmax = xmax, ymin=0, ymax = 1, fill = motif, group = 1)) +
				geom_line(data=sline, aes(x = x, y = y, group = 1), color = "#000000") +
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
				labs(x = NULL, y = NULL) + 
				xlim(upstream.end - 1, upstream.start + 1) + ylim(0,1)
		)
		dev.off()
	}
	return(msc)
}
