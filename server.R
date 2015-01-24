library(shiny)
library(flexclust)
library(sROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(seqLogo)
library(pdist)

source("env.R")

env.assemble()

getDistSum <- function() {
	x <- env$cluster.ensemble
	X <- x@k
	Y <- sapply(x@models, function(z) info(z, "distsum"))
	Z <- c(NA, diff(Y))
	df = data.frame(k=X, distSum=Y, distSumDelta=Z)
	return(df);
}

# utilities for plotting cluster profiles
magenta <- "#FF00FF"
cyan <- "#00FFFF"
grey <- "#AAAAAA"
yellow <- "#FFD800"
resmin <- min(env$log.ratio)
resmax <- max(env$log.ratio)
resPlotAdj <- 1.5
rpkmmin <- min(env$rpkm[,names(env$log.ratio)])
rpkmmax <- max(env$rpkm[,names(env$log.ratio)])
rpkmPlotAdj <- 100
cdf2cimin <- function(mycdf) {
        mycdf$x[which.min(abs(mycdf$Fhat - 0.05))]
}
cdf2cimax <- function(mycdf) {
        mycdf$x[which.min(abs(mycdf$Fhat - 0.95))]
}

addResourcePath("cluster_analysis.dir", "/Users/dacb/work/5G/cluster_analysis/cluster_analysis.dir")

# this should really be setup in env!
motif.colors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")

# demo genes, could also just fetch first 4 (4 gets the kCDF display)
# but these 3 make for a nice demo
default.my.cluster.genes <- c("MBURv2_160308", "MBURv2_160304", "MBURv2_160312")

shinyServer(
	function(input, output, session) {
		# all k tab
		kdsdf <- getDistSum()
		output$kDistSumPlot <- renderPlot({
			ggplot(kdsdf, aes(x=k, y=distSum)) +
				geom_bar(stat = "identity") +
				xlab("Input k") +
				ylab("Sum of within cluster distances") +
				theme_bw()
		})
		output$kDistSumDeltaPlot <- renderPlot({
			ggplot(kdsdf, aes(x=k, y=distSumDelta)) +
				geom_bar(stat = "identity") +
				xlab("Input k") +
				ylab("Difference of sum of within\ncluster distances to next smaller k") +
				theme_bw()
		})

		# choose k tab
		kclust <- reactive({
			env$cluster.ensemble[[input$k]]
		})
		output$k <- renderText({
			input$k
		})
		output$clusterSizePlot <- renderPlot({
			csdf <- as.data.frame(table(kclust()@cluster))
			ggplot(csdf, aes(x=Var1, y=Freq)) + 
				geom_bar(stat = "identity") + 
				ylab("Cluster size") +
				xlab("Cluster index") +
				theme_bw()
		})
		output$clusterOverviewPlot <- renderPlot({
			plot(kclust(), project=env$prcomp)
		})
		output$clusterProfileOverviewPlotArea <- renderPlot ({
			profilePlots <- lapply(1:input$k, k=input$k, simple=T, makeClusterProfilePlot)
			do.call(grid.arrange, c(profilePlots, list(ncol=3)))
		})

		# choose cluster tab
		clusts <- reactive({
			clusters(env$cluster.ensemble[[input$k]])
		})
		clust <- reactive({
			clusts()[clusts()==input$cluster]
		})
		output$cluster <- renderText({
			input$cluster
		})
		output$clusterSelection <- renderUI({
			# list of clusters
			clist <- 1:input$k
			# use a preselected cluster if available
			csr <- getClusterSearchResults(input$k, input$searchText)
			selectInput("cluster", "Choose cluster", clist, selected=csr[input$clusterSearchResultSelectedRow, "Cluster"])
		})
		output$clusterProfilePlot <- renderPlot({
			if (is.null(input$clusterSelectedRows)) {
				rowFocus <- F
			} else {
				rowFocus <- input$clusterSelectedRows
			}
			makeClusterProfilePlot(input$k, input$cluster, focus=rowFocus,
				displayMotifGeneProfile=c(1:env$meme.nmotifs)[
					c(
						input$displayMotif1GeneProfile,
						input$displayMotif2GeneProfile,
						input$displayMotif3GeneProfile,
						input$displayMotif4GeneProfile
					)
				]
			)
		})
		memes <- reactive({
			env$meme.data[[input$k]][[input$cluster]]
		})
		output$clusterMotif1Summary <- renderText({
			paste("E-value:", memes()[[1]]$e.value, "- genes: ", length(memes()[[1]]$posns$gene))
		})
		output$clusterMotif1Consensus <- renderText({
			paste("Consesus:", memes()[[1]]$consensus);
		})
		output$clusterMotif1Plot <- renderPlot({
			seqLogo(t(memes()[[1]]$pssm))
		})
		output$clusterMotif2Summary <- renderText({
			paste("E-value:", memes()[[2]]$e.value, "- genes: ", length(memes()[[2]]$posns$gene))
		})
		output$clusterMotif2Consensus <- renderText({
			paste("Consesus:", memes()[[2]]$consensus);
		})
		output$clusterMotif2Plot <- renderPlot({
			seqLogo(t(memes()[[2]]$pssm))
		})
		output$clusterMotif3Summary <- renderText({
			paste("E-value:", memes()[[3]]$e.value, "- genes: ", length(memes()[[3]]$posns$gene))
		})
		output$clusterMotif3Consensus <- renderText({
			paste("Consesus:", memes()[[3]]$consensus);
		})
		output$clusterMotif3Plot <- renderPlot({
			seqLogo(t(memes()[[3]]$pssm))
		})
		output$clusterMotif4Summary <- renderText({
			paste("E-value:", memes()[[4]]$e.value, "- genes: ", length(memes()[[4]]$posns$gene))
		})
		output$clusterMotif4Consensus <- renderText({
			paste("Consesus:", memes()[[4]]$consensus);
		})
		output$clusterMotif4Plot <- renderPlot({
			seqLogo(t(memes()[[4]]$pssm))
		})
		# old way, variable plots
		#	grid.newpage()
		#	for (i in 1:length(memes())) {
		#		pushViewport(viewport(x=(i*.25)-.125, y=0.5, width=0.25, height=1))
		#		print(memes()[[i]]$pssm)
		#		mySeqLogo(t(memes()[[i]]$pssm))
		#		popViewport()
		#	}
		output$clusterMembers <- renderDataTable({
			ns <- names(clust())

			# hierachical clustering of rows for row ordering
			# could this be precomputed?
			clustres <- env$log.ratio[ns,]
			hclustres <- hclust(dist(clustres), method="complete")
			ns <- ns[hclustres$order]

			dir <- paste(env$dir.output, paste("k_", env$cluster.ensemble[[input$k]]@k, ".dir/cluster_", input$cluster, ".dir/motif_plots.dir", sep=""), sep="/")
			motif_img <- paste("<img src='", paste("http://127.0.0.1:4202", dir, paste(ns, ".png", sep=""), sep="/"), "' alt=''></img>", sep="")
			ms <- env$meme.sites[[input$k]][[input$cluster]]
			for (n in 1:length(ns)) {
				if (dim(ms[ms$gene==ns[n],])[1] == 0) {
					motif_img[n] <- ""
				}
			}
			data.frame(locus_tag=ns, product=env$rpkm[ns,"product"], motifs=motif_img)
		}, options = list(paging = F),
			callback = "function(table) {
      				table.on('click.dt', 'tr', function() {
        				$(this).toggleClass('selected');
						var seldata = table.rows('.selected').indexes().toArray();
						var data = table.rows('.selected').data().data();
						var genes = [];
						for (sel in seldata) {
							genes.push(data[seldata[sel]][0])
						}
						console.log(genes);
        				Shiny.onInputChange('clusterSelectedRows', genes);
					});
				}"
		)
		output$downloadClusterData <- downloadHandler(
			filename = function() { paste("k", input$k, "_cluster", input$cluster, ".xls", sep='') },
			content = function(file) {
				ns <- names(clust())
				write.table(
					data.frame(locus_tag=ns, 
						product=env$rpkm[ns,"product"], 
						env$rpkm[ns,2:length(names(env$rpkm))],
						env$log.ratio[ns,]
					),
					file, quote=F, sep='\t', row.names=F)
			}
		)

		# search cluster tab
		output$clusterSearchResults <- renderDataTable({
			getClusterSearchResults(input$k, input$searchText)
		}, options = list(paging = F),
			callback = "function(table) {
      				table.on('click.dt', 'tr', function() {
						table.$('tr.selected').removeClass('selected');
        				$(this).toggleClass('selected');
						var seldata = table.rows('.selected').indexes().toArray();
						var data = table.rows('.selected').data().data();
						var genes = [];
						for (sel in seldata) {
							genes.push(data[seldata[sel]][0])
						}
						console.log(genes);
        				Shiny.onInputChange('clusterSearchResultSelectedRow', genes);
        				Shiny.onInputChange('clusterSelectedRows', genes);

						 tabs = $('.nav li')
					 	 tabs.each(function() {
							$(this).removeClass('active')
					 	 })
						 $(tabs[2]).addClass('active')
						
						 tabsContents = $('.tab-content .tab-pane')
					 	 tabsContents.each(function() {
							$(this).removeClass('active')
					 	 })
						 $(tabsContents[2]).addClass('active')

						 $('#cluster').trigger('change').trigger('shown');
						 
					});
				}"
		)
		output$clusterSearchResultSelectedRows <- renderText({
			csr <- getClusterSearchResults(input$k, input$searchText)
    		paste(c('Cluster:', csr[input$clusterSearchResultSelectedRow, "Cluster"]), collapse = ' ')
  		})
		# My cluster tab
		observe({
			if (input$myClusterRecruitButton != 0) {
				isolate({
					my.cluster.log.ratio <- env$log.ratio[myClusterGenes(),]
					other.log.ratio <- env$log.ratio[!rownames(env$log.ratio) %in% myClusterGenes(),]
					switch(input$myClusterRecruitBy,
						min2centroid = {
								cmean<-apply(my.cluster.log.ratio, 2, mean)
								other.log.ratio$dist <- sqrt(rowSums(t(t(other.log.ratio)-cmean)^2))
								new_genes <- rownames(other.log.ratio[order(other.log.ratio$dist),])[1:input$myClusterRecruitN]
							},
						min2member = {
								print("m2m")
								pdm <- as.matrix(pdist(other.log.ratio, my.cluster.log.ratio))
								# find the minimum for each row (gene to each member)
								rmin <- t(sapply(seq(nrow(pdm)), function(i) {
									j <- which.min(pdm[i,])
									pdm[i,j]
								}))
								other.log.ratio$dist <- t(rmin)
								new_genes <- rownames(other.log.ratio[order(other.log.ratio$dist),])[1:input$myClusterRecruitN]
							},
						random = {
								rrow <- sample(nrow(other.log.ratio), input$myClusterRecruitN)
								new_genes <- rownames(other.log.ratio[rrow,])
							},
						{	# default case, report a warning
							warning(paste("unhandled input$myClusterRecruitBy case:", input$myClusterRecruitBy))
						}
					)
					# send a client side message about the update to the textarea which will cascade the whole tab
					message <- list(
						value=paste(paste(myClusterGenes(), collapse="\n"), paste(new_genes, collapse="\n"), sep="\n")
					)
					session$sendInputMessage("myClusterGenes", message)
				})
			}
		})
		myClusterGenes <- reactive({
			if (!is.null(input$myClusterGenes)) {
				return(unlist(strsplit(input$myClusterGenes, "\n", fixed=T)))
			}
			return(default.my.cluster.genes)
		})
		myClusterMemes <- reactive({
			if (length(myClusterGenes()) > 1) {
				clust_seqs_upstream <- env$seqs.upstream[myClusterGenes(),]
				dir <- paste(env$dir.output, "my_cluster.dir", sep="/")
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
				meme_file <- paste(dir, env$file.meme.txt, sep="/")
                meme.cmd <- paste(env$path.to.meme, fafile, "-nmotifs", env$meme.nmotifs, env$meme.base.args, "-oc", dir, "-bfile", paste("..", env$file.meme.bfile, sep="/"), ">&", meme_file)
				print(meme.cmd)
				system(meme.cmd)
				meme_text <- readLines(meme_file)
				meme.data <- memeParse(meme_text)
				meme.sites <- renderMotifPlots(paste(dir, "motif_plots.dir", sep="/"), myClusterGenes, meme.data)
				return(list("meme.data"=meme.data, "meme.sites"=meme.sites))
			}
			return(NULL)
		})
		output$myClusterGenesUI <- renderUI({
			tags$textarea(id="myClusterGenes", rows=8, cols=32, paste(myClusterGenes(), collapse="\n"), style="display: block; margin-left: auto; margin-right: auto;")
		})
		output$myClusterProfilePlot <- renderPlot({
			makeMyClusterProfilePlot(myClusterGenes(), myClusterMemes()$meme.data,
				displayMotifGeneProfile=c(1:4)[
					c(
						input$displayMyMotif1GeneProfile,
						input$displayMyMotif2GeneProfile,
						input$displayMyMotif3GeneProfile,
						input$displayMyMotif4GeneProfile
					)
				]
			)
		})
		output$myClusterMembers <- renderDataTable({
			ns <- myClusterGenes()
			mcm <- myClusterMemes()
			# hierachical clustering of rows for row ordering
			# could this be precomputed?
			clustres <- env$log.ratio[ns,]
			if (length(ns) > 1) {
				hclustres <- hclust(dist(clustres), method="complete")
				ns <- ns[hclustres$order]
			}

			# put together path of the motif image for each gene
			dir <- paste(env$dir.output, "my_cluster.dir", "motif_plots.dir", sep="/")
			# use runif to append a random number to prevent all caching here
			motif_img <- paste("<img src='", paste("http://127.0.0.1:4202", dir, paste(ns, ".png?", runif(1, min=0, max=10), sep=""), sep="/"), "' alt=''></img>", sep="")
			# for genes with no sites, empty out the image url
			ms <- mcm$meme.sites
			for (n in 1:length(ns)) {
				if (dim(ms[ms$gene==ns[n],])[1] == 0) {
					motif_img[n] <- ""
				}
			}
			data.frame(locus_tag=ns, product=env$rpkm[ns,"product"], motifs=motif_img)
		}, options = list(paging = F))
		output$myClusterMotif1Summary <- renderText({
			paste("E-value:", myClusterMemes()$meme.data[[1]]$e.value, "- genes: ", length(myClusterMemes()$meme.data[[1]]$posns$gene))
		})
		output$myClusterMotif1Consensus <- renderText({
			paste("Consesus:", myClusterMemes()$meme.data[[1]]$consensus);
		})
		output$myClusterMotif1Plot <- renderPlot({
			seqLogo(t(myClusterMemes()$meme.data[[1]]$pssm))
		})
		output$myClusterMotif2Summary <- renderText({
			paste("E-value:", myClusterMemes()$meme.data[[2]]$e.value, "- genes: ", length(myClusterMemes()$meme.data[[2]]$posns$gene))
		})
		output$myClusterMotif2Consensus <- renderText({
			paste("Consesus:", myClusterMemes()$meme.data[[2]]$consensus);
		})
		output$myClusterMotif2Plot <- renderPlot({
			seqLogo(t(myClusterMemes()$meme.data[[2]]$pssm))
		})
		output$myClusterMotif3Summary <- renderText({
			paste("E-value:", myClusterMemes()$meme.data[[3]]$e.value, "- genes: ", length(myClusterMemes()$meme.data[[3]]$posns$gene))
		})
		output$myClusterMotif3Consensus <- renderText({
			paste("Consesus:", myClusterMemes()$meme.data[[3]]$consensus);
		})
		output$myClusterMotif3Plot <- renderPlot({
			seqLogo(t(myClusterMemes()$meme.data[[3]]$pssm))
		})
		output$myClusterMotif4Summary <- renderText({
			paste("E-value:", myClusterMemes()$meme.data[[4]]$e.value, "- genes: ", length(myClusterMemes()$meme.data[[4]]$posns$gene))
		})
		output$myClusterMotif4Consensus <- renderText({
			paste("Consesus:", myClusterMemes()$meme.data[[4]]$consensus);
		})
		output$myClusterMotif4Plot <- renderPlot({
			seqLogo(t(myClusterMemes()$meme.data[[4]]$pssm))
		})
		output$myClusterMemeLog <- renderText({
			# register reactivity with the gene list text area
			input$myClusterGenes
			dir <- paste(env$dir.output, "my_cluster.dir", sep="/")
			meme_file <- paste(dir, "meme.txt", sep="/")
			paste(readLines(meme_file), "\n")
		})

		# blastn
		output$blastnResults <- renderDataTable({
			data.frame(BLASTn=c("disabled"), reason=c("insuffecient resources"))
		}, options = list(paging=F))

		# blastp
		output$blastpResults <- renderDataTable({
			data.frame(BLASTp=c("disabled"), reason=c("insuffecient resources"))
		}, options = list(paging=F))

		# likes
		output$likesTable <- renderDataTable({
			data.frame(ID=c("1421879893"), description=c("demonstration workflow showing recruitment to My Cluster by BLASTn of motif consensus"))
		}, options = list(paging=F))
	}
)

getClusterSearchResults <- function(k, searchText) {
		rowSelect <- union(
				grep(searchText, rownames(env$rpkm), ignore.case=T),
				grep(searchText, env$rpkm$product, ignore.case=T)
			)
		gcsr_clusts <- clusters(env$cluster.ensemble[[k]])
		gcsr_clust <- gcsr_clusts[rowSelect]
		data.frame(locus_tag=names(gcsr_clust), product=env$rpkm[names(gcsr_clust),"product"], Cluster=gcsr_clust)
}

makeMyClusterProfilePlot <- function(genes, memes, displayMotifGeneProfile=F) {
	clustres <- env$log.ratio[genes,]
	genesdf <- data.frame(t(clustres), Sample=names(clustres))
	genesdf$Sample <- factor(genesdf$Sample, levels=names(clustres))
	mdf <- melt(genesdf, id.vars=c("Sample"))
	if (length(genes) > 4) {	
		cmin<-apply(clustres, 2, min)
		cmean<-apply(clustres, 2, mean)
		cmax<-apply(clustres, 2, max)
		# estimate CDF from data
		cis <- apply(clustres, 2, kCDF)
		# get min at 95% CI
		lwr <- sapply(cis, cdf2cimin)
		# get max at 95% CI
		upr <- sapply(cis, cdf2cimax)
		# cluster data frame
		clustdf <- data.frame(min=cmin, max=cmax, mean=cmean, lwr=lwr, upr=upr, Sample=names(cmean))
		clustdf$Sample <- factor(clustdf$Sample, levels=names(cmean))
		myplot <- ggplot(clustdf, aes(x=Sample)) + 
			geom_ribbon(aes(ymax=upr, ymin=lwr, group=1, alpha=0.7), colour=magenta)
	} else {
		myplot <- ggplot(mdf, aes(x=Sample, y=value, group=variable))
	}
	myplot <- myplot +
			geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow) +
			geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow) +
			ggtitle("Expression profile") +
			ylab("Log2(Sample / T0)\nNormalized counts") +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")
	if (!identical(displayMotifGeneProfile, F) && length(displayMotifGeneProfile) > 0) {
			for (m in 1:length(displayMotifGeneProfile)) {
				motif <- memes[[displayMotifGeneProfile[m]]]
				motifdf <- data.frame(t(clustres[as.character(motif$posns$gene),]), Sample=names(clustres))
				motifdf$Sample <- factor(genesdf$Sample, levels=names(clustres))
				mdf <- melt(motifdf, id.vars=c("Sample"))
				myplot <- myplot + 
					geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]]) +
					geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]])
			}
	}
	return(myplot)
}

makeClusterProfilePlot <- function(k, cluster, simple=F, focus=F, displayMotifGeneProfile=F) {
		mcpp_clusts <- clusters(env$cluster.ensemble[[k]])
		mcpp_clust <- mcpp_clusts[mcpp_clusts==cluster]
		clustres <- env$log.ratio[names(mcpp_clust),]
		cmin<-apply(clustres, 2, min)
		cmean<-apply(clustres, 2, mean)
		cmax<-apply(clustres, 2, max)
		if (!simple) {
			# estimate CDF from data
			cis <- apply(clustres, 2, kCDF)
			# get min at 95% CI
			lwr <- sapply(cis, cdf2cimin)
			# get max at 95% CI
			upr <- sapply(cis, cdf2cimax)
			# cluster data frame
			clustdf <- data.frame(min=cmin, max=cmax, mean=cmean, lwr=lwr, upr=upr, Sample=names(cmean))
			clustdf$Sample <- factor(clustdf$Sample, levels=names(cmean))
			myplot <- ggplot(clustdf, aes(x=Sample)) + 
				ggtitle(paste(paste("K =", k, ": Cluster", cluster, paste("(", length(mcpp_clust), " genes)", sep="")), "Expression profile", sep="\n")) +
				ylab("Log2(Sample / T0)\nNormalized counts") +
				geom_ribbon(aes(ymax=upr, ymin=lwr, group=1, alpha=0.7), colour=magenta) +
				theme_bw() +
				theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")
		} else {
			clustdf <- data.frame(min=cmin, max=cmax, mean=cmean, Sample=names(cmean))
			clustdf$Sample <- factor(clustdf$Sample, levels=names(cmean))

			# cluster label for center of plot
			clusterGrob <- grobTree(textGrob(cluster, x=.5, y=.5, gp=gpar(col="black", fontsize=30)))
			clusterGrobGenes <- grobTree(textGrob(paste(length(mcpp_clust), "genes"), x=.5, y=.1, gp=gpar(col="black", fontsize=15)))

			myplot <- ggplot(clustdf, aes(x=Sample)) +
				theme_bw() +
				theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank()) +
  				annotation_custom(clusterGrob) +
  				annotation_custom(clusterGrobGenes) +
				ylim(resmin - resPlotAdj, resmax + resPlotAdj)
		}
		if (!identical(focus, F)) {
			focusdf <- data.frame(t(clustres[focus,]), Sample=names(clustres))
			focusdf$Sample <- factor(focusdf$Sample, levels=names(clustres))
			mdf <- melt(focusdf, id.vars=c("Sample"))
			myplot <- myplot + 
				geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow) +
				geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=yellow)
		}
		if (!identical(displayMotifGeneProfile, F) && length(displayMotifGeneProfile) > 0) {
				memes <- env$meme.data[[k]][[cluster]]
				for (m in 1:length(displayMotifGeneProfile)) {
					motif <- memes[[displayMotifGeneProfile[m]]]
					motifdf <- data.frame(t(clustres[as.character(motif$posns$gene),]), Sample=names(clustres))
					motifdf$Sample <- factor(names(clustres), levels=names(clustres))
					mdf <- melt(motifdf, id.vars=c("Sample"))
					myplot <- myplot + 
						geom_point(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]]) +
						geom_line(data=mdf, aes(x=Sample, y=value, group=variable), color=motif.colors[displayMotifGeneProfile[m]])
				}
		}
		myplot +
			geom_point(aes(y=mean), colour=cyan) + 
			geom_line(aes(x=Sample, y=mean, group=1), colour=cyan) + 
			geom_line(aes(x=Sample, y=min, group=1), colour=grey) + 
			geom_line(aes(x=Sample, y=max, group=1), colour=grey)
}

renderMotifPlots <- function(dir, genes, meme.data) {
	dir.create(dir, recursive=T)
	memes <- meme.data
	sites <- list()
	mi <- 1
	for (motif in memes) {
		# new data frame from position table, fill in motif and with with rep
		# compute xmin, xmax for easy ggplot viewing
		uplen <- env$seqs.upstream[as.character(motif$posns$gene),"uplength"]
		df <- data.frame(
			cbind(motif$posns),
			motif=rep(mi, length(motif$posns$gene)),
			width=rep(motif$width, length(motif$posns$gene)),
			xmin=motif$posns$start-uplen,
			xmax=motif$posns$start+motif$width-uplen
		)
		sites[[as.character(mi)]] <- df
		mi <- mi + 1
	}
	msc <- do.call("rbind", sites)
	msc$motif <- as.factor(msc$motif)
	names(motif.colors) <- levels(msc$motif)
	for (g in levels(msc$gene)) {
		mscg <- msc[msc$gene==g,]
		gene_file <- paste(dir, paste(g, "png", sep="."), sep="/")
		png(filename=gene_file, width=180, height=18)
		ul <- env$seqs.upstream[as.character(g),"uplength"]
		sline <- data.frame(x=c(-ul+2, 2), y=c(.5,.5))
		print(
			ggplot(mscg) +
				scale_fill_manual(name="nmotif", values=motif.colors) +
				geom_rect(aes(xmin=xmin, xmax=xmax, ymin=0, ymax=1, fill=motif, group=1)) +
				geom_line(data=sline, aes(x=x, y=y, group=1), color="#000000") +
				theme(
					legend.position = "none",
					panel.background = element_blank(),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					panel.margin = unit(0,"null"),
					plot.margin = rep(unit(0,"null"),4),
					axis.ticks = element_blank(),
					axis.text.x = element_blank(),
					axis.text.y = element_blank(),
					axis.title.x = element_blank(),
					axis.title.y = element_blank(),
					axis.ticks.length = unit(0,"null"),
					axis.ticks.margin = unit(0,"null")
				) +
				labs(x=NULL, y=NULL) + 
				xlim(-150,2) + ylim(0,1)
		)
		dev.off()
	}
	return(msc)
}
