library(flexclust)

load("cluster_analysis.2.RData")

env$path.to.meme <- "meme"
env$meme.base.args <- "-dna -maxsize 1000000 -evt 1e10 -minw 6 -maxw 25 -mod zoops -nostatus -text"
env$meme.nmotifs <- 4
env$file.meme.bfile <- "5G.genome.bfile"
env$file.meme.jobs <- "meme.jobs"
env$file.upstream.seqs <- "5G.upstream.tab"

# read in upstream sequences
env$seqs.upstream <- read.delim(env$file.upstream.seqs, skip = 3, header = T, row.names = 1, sep='\t');
head(env$seqs.upstream);

# iterate over all k and clusters and produce a fasta for each
# then run meme on the fasta
if (file.exists(env$file.meme.jobs)) {
	file.remove(env$file.meme.jobs)
}
for (i in 1:length(env$cluster.ensemble@k)) {
	clusts <- clusters(env$cluster.ensemble[[i]])
	for (j in 1:env$cluster.ensemble@k[i]) {
		clust <- clusts[clusts==j]
		clust_seqs_upstream <- env$seqs.upstream[names(clust),]

		# setup pathes for output
		dir <- paste(env$dir.output, paste("k_", env$cluster.ensemble[[i]]@k, ".dir/cluster_", j, ".dir", sep=""), sep="/")
		dir.create(dir, recursive=T)
		fafile <- paste(dir, "upstream.fa", sep="/")
		if (file.exists(fafile)) {
			file.remove(fafile);
		}
		for (k in 1:length(rownames(clust_seqs_upstream))) {
			if (!is.na(clust_seqs_upstream$sequence[k])) {
				cat(paste(">", rownames(clust_seqs_upstream)[k], "\n", sep="") , file=fafile, append=T)
				cat(paste(clust_seqs_upstream$sequence[k], "\n", sep="") , file=fafile, append=T)
			}
		}
		meme.cmd <- paste(env$path.to.meme, fafile, "-nmotifs", env$meme.nmotifs, env$meme.base.args, "-oc", dir, "-bfile", env$file.meme.bfile, ">&", paste(dir, env$file.meme.txt, sep="/"))
		cat(meme.cmd, "\n", file=env$file.meme.jobs, append=T)
	}
}
warnings()

save(env, file="cluster_analysis.3.RData")
