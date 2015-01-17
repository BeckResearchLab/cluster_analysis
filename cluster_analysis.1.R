
rpkm <- read.delim("5G_rpkm.xls", header=T, row.names=1, sep="\t")
head(rpkm)
names(rpkm)

d <- read.delim("5G_counts.xls", header=T, row.names=1, sep="\t")
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

res <- data.frame(assay(rld) - assay(rld)[,"X5GB1_FM23_TR3"])

head(res)

save.image(file = "cluster_analysis.1.RData", compress=T)
