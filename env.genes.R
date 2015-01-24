# to be sourced

env.genes.setup <- function(file.genes.annotations, file.genes.upstream.seqs) {
	genes <- list()

	cat("reading annotations...\n")
	genes$annotations <- read.delim(file.genes.annotations, header = T, row.names = 1, sep='\t')
	cat("reading upstream sequences...\n")
	genes$upstream.seqs <- read.delim(file.genes.upstream.seqs, skip = 3, header = T, row.names = 1, sep='\t')

	return(genes)
}
