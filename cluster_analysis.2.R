library(flexclust)

load("cluster_analysis.1.RData")

#apply K-Means with various cluster number
#cls<-c(seq(4,40,by=2), seq(50,100,by=10))
cls<-c(seq(10,40,by=5))
clustEnsemble<-stepFlexclust(res,cls, nrep=3, save.data=TRUE, drop=FALSE, verbose=TRUE)

pdf("cluster_analysis.2.pdf")
plot(clustEnsemble)
dev.off()

save.image("cluster_analysis.2.RData")
