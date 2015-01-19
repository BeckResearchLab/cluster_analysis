env <- list (
		study.title="5GB1",
		file.sample.info="sample_info.xls",
		file.sample.ordering="sample_ordering.xls",
		file.rpkm="5G_rpkm.xls",
		file.counts="5G_counts.xls",
		reference.sample="FM23_TR3"
		dir.output="cluster_analysis.dir"
	)

env$sample.info <- read.delim(env$file.sample.info, header=T, row.names=1, sep="\t")
head(env$sample.info)
names(env$sample.info)
env$sample.ordering <- read.delim(env$file.sample.ordering, header=T, row.names=1, sep="\t")

env$rpkm <- read.delim(env$file.rpkm, header=T, row.names=1, sep="\t")
head(env$rpkm)
names(env$rpkm)

env$cnts.raw <- read.delim(env$file.counts, header=T, row.names=1, sep="\t")
head(env$cnts.raw)
names(env$cnts.raw)

# integration test
head(env$rpkm[,names(env$cnts.raw)])

# preprocess counts
library(DESeq2)
dexp <- data.frame(row.names=colnames(env$cnts.raw), sample=c(names(env$cnts.raw)), condition=env$sample.info[names(env$cnts.raw),"shortd"])
dds <- DESeqDataSetFromMatrix(countData = env$cnts.raw, colData = dexp, design = ~ condition)
# collapse techincal replicates??!?!

rld <- rlog(dds)

env$deseq2.dds <- dds
env$deseq2.rld <- rld

head(assay(rld))

# normalize by reference sample log ratio
res <- data.frame(assay(rld) - assay(rld)[,env$reference.sample])
# remove reference sample from count matrix
res <- res[, !names(res) %in% env$reference.sample]

head(res)

# update our environment object
env$log.ratio <- res
# create an ordering vector for just the samples in the log.ratio table
env$sample.ordering <- rownames(env$sample.ordering)[rownames(env$sample.ordering) %in% names(env$log.ratio)]
# reorder table mostly for viewing pleasure
env$log.ratio <- env$log.ratio[, env$sample.ordering]

library(MASS)
# create MDS
dlr <- dist(env$log.ratio)
# adjust identities to prevent failures
dlr[dlr==0] <- 0.01
mds.dim <- 2
env$cmds <- cmdscale(dlr, mds.dim)
env$mds <- isoMDS(dlr, env$cmds, mds.dim)
# PCA
env$prcomp <- prcomp(env$log.ratio)

save(env, file = "cluster_analysis.1.RData", compress=T)
