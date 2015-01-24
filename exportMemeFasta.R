source("env.R")

# populate env but skip the meme parsing (data don't exist yet!)
env.assemble(include.meme = F)

cat("making fasta files for meme for all clusters...\n")
# iterate over all k and clusters and produce a fasta for each
# then run meme on the fasta
if (file.exists(env$file.meme.jobs)) {
	file.remove(env$file.meme.jobs)
}

files.count <- 0
for (i in 1:length(env$cluster.ensemble@k)) {
	clusts <- clusters(env$cluster.ensemble[[i]])
	for (j in 1:env$cluster.ensemble@k[i]) {
		clust <- clusts[clusts==j]
		clust.upstream.seqs <- env$genes$upstream.seqs[names(clust),]
		# drop genes with NA for upstream seq
		clust.upstream.seqs <- clust.upstream.seqs[!is.na(clust.upstream.seqs$sequence),]

		# setup pathes for output
		dir <- dir.k.cluster(env$dir.output, env$cluster.ensemble[[i]]@k, j, make.dir = T)
		fasta.file <- paste(dir, env$file.upstream.fa, sep="/")
		if (file.exists(fasta.file)) {
			file.remove(fasta.file)
		}
		fasta.lines <- mapply(
			function (n, s) {
				sprintf(">%s\n%s\n", n, s)
			},
			n = rownames(clust.upstream.seqs),
			s = clust.upstream.seqs$sequence
		)
		if (length(fasta.lines) <= 1) {
			warning(
				sprintf("too few sequences for motif prediction, fasta file not written: %s", fasta.file)
			)
			next
		}
		lapply(fasta.lines, cat, "\n", file=fasta.file, append=TRUE)
		files.count <- files.count + 1
		if (files.count %% 100 == 0) {
			cat(sprintf("...processed %d files\n", files.count))
		}
		meme.cmd <- paste(env$path.to.meme, fasta.file, "-nmotifs", env$meme.nmotifs, env$meme.base.args, "-oc", dir, "-bfile", env$file.meme.bfile, ">&", paste(dir, env$file.meme.txt, sep="/"))
		cat(meme.cmd, "\n", file=env$file.meme.jobs, append=T)
	}
}
cat(sprintf("...processed %d files total\n", files.count))
cat(sprintf("meme commands saved to file: %s\n", env$file.meme.jobs))
