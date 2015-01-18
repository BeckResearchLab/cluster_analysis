env <- list (
		study.title="5GB1",
		file.sample_info="sample_info.xls",
		file.rpkm="5G_rpkm.xls",
		file.counts="5G_counts.xls",
		reference.sample="FM23_TR3"
	)

sample_info <- read.delim(env$file.sample_info, header=T, row.names=1, sep="\t")
head(sample_info)
names(sample_info)

rpkm <- read.delim(env$file.rpkm, header=T, row.names=1, sep="\t")
head(rpkm)
names(rpkm)

d <- read.delim(env$file.counts, header=T, row.names=1, sep="\t")
head(d)
names(d)

# integration test
head(rpkm[,names(d)])

# preprocess counts
library(DESeq2)
dexp <- data.frame(row.names=colnames(d), sample=c(names(d)), condition=sample_info[names(d),"shortd"])
dds <- DESeqDataSetFromMatrix(countData = d, colData = dexp, design = ~ condition)
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
env$sample.info <- sample_info
env$rpkm <- rpkm
env$cnts.raw <- d
env$log.ratio <- res

save(env, file = "cluster_analysis.1.RData", compress=T)
