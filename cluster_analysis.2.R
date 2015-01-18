library(flexclust)

load("cluster_analysis.1.RData")

env$clustSeed <- 8576

#apply K-Means with various cluster number
cls<-c(seq(10,40,by=5))
#cls<-c(seq(4,40,by=2), seq(50,100,by=10))
#cls<-c(seq(5,100,by=1))
env$cluster.ensemble <- stepFlexclust(env$cnts.norm, cls, nrep=3, save.data=TRUE, drop=FALSE, verbose=TRUE, seed=env$clustSeed)

pdf("cluster_analysis.2.pdf")
plot(env$cluster.ensemble)
dev.off()

save(env, file="cluster_analysis.2.RData")
