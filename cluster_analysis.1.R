dat <- list (
		file.sample_info="sample_info.xls",
		file.rpkm="5G_rpkm.xls",
		file.counts="5G_counts.xls",
		reference.sample="FM23_TR3"
	)

sample_info <- read.delim(dat[["file.sample_info"]], header=T, row.names=1, sep="\t")
head(sample_info)
names(sample_info)

rpkm <- read.delim(dat[["file.rpkm"]], header=T, row.names=1, sep="\t")
head(rpkm)
names(rpkm)

d <- read.delim(dat[["file.counts"]], header=T, row.names=1, sep="\t")
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

head(assay(rld))

# normalize by reference sample log ratio
res <- data.frame(assay(rld) - assay(rld)[,dat[["reference.sample"]]])
# remove reference sample from count matrix
res <- res[, !names(res) %in% dat[["reference.sample"]]]

head(res)

save.image(file = "cluster_analysis.1.RData", compress=T)
