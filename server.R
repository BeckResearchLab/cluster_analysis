library(shiny)
library(flexclust)
library(sROC)
library(ggplot2)
library(grid)
library(gridExtra)

load("cluster_analysis.2.RData")

getDistSum <- function() {
	x <- clustEnsemble
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
resmin <- min(res)
resmax <- max(res)
resPlotAdj <- 1.5
rpkmmin <- min(rpkm[,names(d)])
rpkmmax <- max(rpkm[,names(d)])
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
		output$k <- renderText({
			input$k
		})
		kclust <- reactive({
			clustEnsemble[[input$k]]
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
			plot(kclust())
		})
		output$clusterProfileOverviewPlot <- renderPlot({
			profilePlots <- lapply(1:input$k, k=input$k, simple=T, makeClusterProfilePlot)
			do.call(grid.arrange, c(profilePlots, list(ncol=3)))
		})

		# choose cluster tab
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
			if (is.null(input$clusterSelectedRow)) {
				rowFocus <- F
			} else {
				rowFocus <- input$clusterSelectedRow + 1
			}
			makeClusterProfilePlot(input$k, input$cluster, focus=rowFocus)
		})
		output$clusterMembers <- renderDataTable({
			getClusterMembers(input$k, input$cluster)
		}, options = list(paging = F),
			callback = "function(table) {
      				table.on('click.dt', 'tr', function() {
						table.$('tr.selected').removeClass('selected');
        				$(this).toggleClass('selected');
        				Shiny.onInputChange('clusterSelectedRow', table.rows('.selected').indexes().toArray());
					});
				}"
		)
		output$downloadClusterData <- downloadHandler(
			filename = function() { paste("k", input$k, "_cluster", input$cluster, ".xls", sep='') },
			content = function(file) {
				write.table(
					getClusterMembersSpreadsheet(input$k, input$cluster),
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
				grep(searchText, rownames(rpkm), ignore.case=T),
				grep(searchText, rpkm$product, ignore.case=T)
			)
		clusts <- clusters(clustEnsemble[[k]])
		clust <- clusts[rowSelect]
		data.frame(locus_tag=names(clust), product=rpkm[names(clust),"product"], Cluster=clust)
}

getClusterMembers <- function(k, cluster) {
		clusts <- clusters(clustEnsemble[[k]])
		clust <- clusts[clusts==cluster]
		data.frame(locus_tag=names(clust), product=rpkm[names(clust),"product"])
}

getClusterMembersSpreadsheet <- function(k, cluster) {
		clusts <- clusters(clustEnsemble[[k]])
		clust <- clusts[clusts==cluster]
		data.frame(locus_tag=names(clust), 
			product=rpkm[names(clust),"product"], 
			rpkm[names(clust),2:length(names(rpkm))],
			res[names(clust),]
		)
}

makeClusterProfilePlot <- function(k, cluster, simple=F, focus=F) {
		clusts <- clusters(clustEnsemble[[k]])
		clust <- clusts[clusts==cluster]
		clustres <- res[names(clust),]
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
				ggtitle(paste(paste("K =", k, ": Cluster", cluster, paste("(", length(clust), " genes)", sep="")), "Expression profile", sep="\n")) +
				ylab("Log2(Sample / T0)\nNormalized counts") +
				geom_ribbon(aes(ymax=upr, ymin=lwr, group=1, alpha=0.7), colour=magenta) +
				theme_bw() +
				theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")
		} else {
			clustdf <- data.frame(min=cmin, max=cmax, mean=cmean, Sample=names(cmean))
			clustdf$Sample <- factor(clustdf$Sample, levels=names(cmean))

			# cluster label for center of plot
			clusterGrob <- grobTree(textGrob(cluster, x=.5, y=.5, gp=gpar(col="black", fontsize=30)))
			clusterGrobGenes <- grobTree(textGrob(paste(length(clust), "genes"), x=.5, y=.1, gp=gpar(col="black", fontsize=15)))

			myplot <- ggplot(clustdf, aes(x=Sample)) +
				theme_bw() +
				theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank()) +
  				annotation_custom(clusterGrob) +
  				annotation_custom(clusterGrobGenes) +
				ylim(resmin - resPlotAdj, resmax + resPlotAdj)
		}
		if (!identical(focus, F)) {
			print(clustres[focus,])
			print(names(clustres[focus,]))
			focusdf <- data.frame(t(clustres[focus,]), Sample=names(clustres))
			names(focusdf)[1] <- 'y'
			focusdf$Sample <- factor(focusdf$Sample, levels=names(clustres))
			print(focusdf)
			print(row.names(focusdf))
			print(names(focusdf))
			myplot <- myplot + 
				geom_point(data=focusdf, aes(x=Sample, y=y, group=1), color=yellow) +
				geom_line(data=focusdf, aes(x=Sample, y=y, group=1), color=yellow)
		}
		myplot +
			geom_point(aes(y=mean), colour=cyan) + 
			geom_line(aes(x=Sample, y=mean, group=1), colour=cyan) + 
			geom_line(aes(x=Sample, y=min, group=1), colour=grey) + 
			geom_line(aes(x=Sample, y=max, group=1), colour=grey)
}
