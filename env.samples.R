# to be sourced
env.samples.setup <- function(file.sample.info, file.sample.ordering, file.sample.tracks, file.rpkm, file.counts, reference.sample) {
	samples <- list()
	cat("reading sample info...\n")
	samples$info <- read.delim(file.sample.info, header=T, row.names=1, sep="\t")
	cat("reading sample ordering...\n")
	samples$ordering <- read.delim(file.sample.ordering, header=T, row.names=1, sep="\t")
	cat("reading sample tracks...\n")
	samples$tracks <- read.delim(file.sample.tracks, header=T, row.names=1, sep="\t")
	cat("reading sample rpkm...\n")
	samples$rpkm <- read.delim(file.rpkm, header=T, row.names=1, sep="\t")
	cat("reading sample counts...\n")
	samples$counts <- read.delim(file.counts, header=T, row.names=1, sep="\t")

	# integration test
	head(samples$rpkm[,names(samples$counts)])

	# preprocess counts
	library(DESeq2)
	cat("normalizing counts with DESeq2...\n")
	dexp <- data.frame(row.names=colnames(samples$counts), sample=c(names(samples$counts)), condition=samples$info[names(samples$counts),"shortd"])
	dds <- DESeqDataSetFromMatrix(countData = samples$counts, colData = dexp, design = ~ condition)
	# collapse techincal replicates??!?!
	# do log transform
	rld <- rlog(dds)
	# process into env
	samples$deseq2.dds <- dds
	samples$deseq2.rld <- rld
	#head(assay(rld))

	cat ("computing log ratios against refernce...\n")
	# normalize by reference sample log ratio
	res <- data.frame(assay(rld) - assay(rld)[,reference.sample])
	# remove reference sample from count matrix
	res <- res[, !names(res) %in% reference.sample]
	# remove corresponding QC for reference sample if it exists
	res <- res[, !names(res) %in% paste(reference.sample, "_QC", sep = "")]
	#head(res)

	# update our environment object
	samples$log.ratio <- res
	# create an ordering vector for just the samples in the log.ratio table
	samples$ordering <- rownames(samples$ordering)[rownames(samples$ordering) %in% names(samples$log.ratio)]
	# reorder table mostly for viewing pleasure
	samples$log.ratio <- samples$log.ratio[, samples$ordering]
	# derive fancy names for samples
	samples$info$fancy.names <- paste(samples$info$shortd, " (", rownames(samples$info), ")", sep="")

	cat("performing MDS...\n")
	library(MASS)
	# create MDS
	dlr <- dist(samples$log.ratio)
	# adjust identities to prevent failures
	dlr[dlr==0] <- 0.01
	mds.dim <- 2
	samples$cmds <- cmdscale(dlr, mds.dim)
	samples$mds <- isoMDS(dlr, samples$cmds, mds.dim)
	# create PCA
	cat("performing PCA...\n")
	samples$prcomp <- prcomp(samples$log.ratio)

	return(samples)
}
