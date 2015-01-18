
load("cluster_analysis.2.RData")

# read in upstream sequences
seqs_upstream <- read.delim("5G.upstream.tab", skip = 3, header = T, row.names = 1, sep='\t');
head(seqs_upstream);

save.image("cluster_analysis.3.RData")

exit()
library(genomeIntervals)
library(IRanges)

# read in GFF
gff <- readGff3("5G.genbank.gff")
gff_genes <- RangedData( 
		IRanges( start=gff[,1], end=gff[,2]), 
		space=gff$seq_name,
		strand=gff$strand,
		type=as.factor(gff$type),
		parent=as.vector(
			getGffAttribute(gff,"Parent")),
		name=as.vector(
			getGffAttribute(gff,"Name")),
		ID=as.vector(
			getGffAttribute(gff,"ID"))
	)

