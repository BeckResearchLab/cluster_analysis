library(shiny)
library(flexclust)
library(sROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

load("../cluster_analysis.2.RData")

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

shinyServer(
	function(input, output) {
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
			selectInput("cluster", "Choose cluster", clist, selected=csr[input$clusterSearchResultSelectedRow + 1, "Cluster"])
		})
		output$clusterProfilePlot <- renderPlot({
			if (is.null(input$clusterSelectedRows)) {
				rowFocus <- F
			} else {
				rowFocus <- input$clusterSelectedRows + 1
			}
			makeClusterProfilePlot(input$k, input$cluster, focus=rowFocus)
		})
		output$clusterMembers <- renderDataTable({
			ns <- names(clust())
			data.frame(locus_tag=ns, product=env$rpkm[ns,"product"])
		}, options = list(paging = F),
			callback = "function(table) {
      				table.on('click.dt', 'tr', function() {
        				$(this).toggleClass('selected');
        				Shiny.onInputChange('clusterSelectedRows', table.rows('.selected').indexes().toArray());
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
        				Shiny.onInputChange('clusterSearchResultSelectedRow', table.rows('.selected').indexes().toArray());
        				Shiny.onInputChange('clusterSelectedRows', null);

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
    		paste(c('Cluster:', csr[input$clusterSearchResultSelectedRow + 1, "Cluster"]), collapse = ' ')
  		})
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

makeClusterProfilePlot <- function(k, cluster, simple=F, focus=F) {
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
		myplot +
			geom_point(aes(y=mean), colour=cyan) + 
			geom_line(aes(x=Sample, y=mean, group=1), colour=cyan) + 
			geom_line(aes(x=Sample, y=min, group=1), colour=grey) + 
			geom_line(aes(x=Sample, y=max, group=1), colour=grey)
}
