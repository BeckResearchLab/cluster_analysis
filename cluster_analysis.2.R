library(flexclust)

load("cluster_analysis.1.RData")

env$clustSeed <- 8576

#apply k-means with various cluster number
# values of k to try
cls<-c(seq(10,150,by=10))
#cls<-c(seq(4,40,by=2), seq(50,100,by=10))
#cls<-c(seq(5,100,by=1))
# cores to use
options(mc.cores=8)
env$cluster.ensemble <- stepFlexclust(env$log.ratio, cls, nrep=8, save.data=TRUE, drop=FALSE, verbose=TRUE, seed=env$clustSeed, multicore=T)

pdf("cluster_analysis.2.pdf")
plot(env$cluster.ensemble)
dev.off()

save(env, file="cluster_analysis.2.RData")
