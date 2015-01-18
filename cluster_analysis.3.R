library(flexclust)

load("cluster_analysis.2.RData")

env$path.to.meme <- "meme"
env$meme.base.args <- "-dna -maxsize 1000000 -evt 1e10 -minw 6 -maxw 25 -mod zoops -nostatus -text",
env$meme.nmotifs <- 4
env$file.meme.bfile <- "5G.genome.bfile"
env$file.meme.jobs <- "meme.jobs"
env$file.upstream.seqs <- "5G.upstream.tab"

# read in upstream sequences
seqs_upstream <- read.delim("5G.upstream.tab", skip = 3, header = T, row.names = 1, sep='\t');
head(seqs_upstream);

# iterate over all k and clusters and produce a fasta for each
# then run meme on the fasta
if (file.exists(env$meme.jobs)) {
	file.remove(env$meme.jobs)
}
root.dir <- 'cluster_analysis.dir'
for (i in 1:length(clustEnsemble@k)) {
	clusts <- clusters(clustEnsemble[[i]])
	for (j in 1:clustEnsemble@k[i]) {
		clust <- clusts[clusts==j]
		clust_seqs_upstream <- seqs_upstream[names(clust),]
		fafile <- 

		# setup pathes for output
		dir <- paste(root.dir, paste("k_", clustEnsemble[[i]]@k, ".dir/cluster_", j, ".dir", sep=""), sep="/")
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
		meme.cmd <- paste(env$path.to.meme, fafile, "-nmotifs", env$meme.nmotifs, env$meme.base.args, "-oc", dir, "-bfile", env$meme.bfile)
		cat(meme.cmd, "\n", file=env$meme.jobs, append=T)
	}
}
warnings()

save(env, file="cluster_analysis.3.RData")
