
rpkm <- read.delim("5G_rpkm_Cu.xls", header=T, row.names=1, sep="\t")
head(rpkm)
names(rpkm)

d <- read.delim("5G_counts_Cu.xls", header=T, row.names=1, sep="\t")
head(d)
names(d)
d$product <- NULL
head(d)
split_names <- do.call(rbind, strsplit(names(d), '_'))[,2:3]
short_names <- paste(split_names[,1], split_names[,2], sep='_')
short_names

names(rpkm)
names(d)
head(rpkm[,names(d)])

library(DESeq2)
dexp <- data.frame(row.names=colnames(d), sample=c(names(d)), condition=c(short_names))
dds <- DESeqDataSetFromMatrix(countData = d, colData = dexp, design = ~ condition)
# collapse techincal replicates??!?!

rld <- rlog(dds)

head(assay(rld))

res <- data.frame(FM40_T10m_TR3 = assay(rld)[,2] - assay(rld)[,1],
	FM40_T40m_TR1 = assay(rld)[,3] - assay(rld)[,1],
	FM40_T60m_TR1 = assay(rld)[,4] - assay(rld)[,1],
	FM40_T90m_TR2 = assay(rld)[,5] - assay(rld)[,1],
	FM40_T150m_TR1 = assay(rld)[,6] - assay(rld)[,1],
	FM40_T180m_TR1 = assay(rld)[,7] - assay(rld)[,1],
	FM34_T3_TR3_QC = assay(rld)[,9] - assay(rld)[,8],
	FM34_T4_TR3_QC = assay(rld)[,10] - assay(rld)[,8],
	FM34_T5_TR2_QC = assay(rld)[,11] - assay(rld)[,8],
	FM34_T6_TR3_QC = assay(rld)[,12] - assay(rld)[,8],
	FM34_T7_TR3_QC = assay(rld)[,13] - assay(rld)[,8],
	FM34_T8_TR1_QC = assay(rld)[,14] - assay(rld)[,8])

head(res)

save.image(file = "cluster_analysis.1.RData", compress=T)
