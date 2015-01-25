library(shiny)
library(flexclust)
library(sROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(seqLogo)
library(pdist)

options(error = recover)

source("env.R")
source("utilities.R")
source("plotting.R")

env.assemble()

addResourcePath("cluster_analysis.dir", sprintf("%s/%s", env$dir.root, env$dir.output))

# demo genes, could also just fetch first 4 (4 gets the kCDF display)
# but these 3 make for a nice demo
default.my.cluster.genes <- c("MBURv2_160308", "MBURv2_160304", "MBURv2_160312")

shinyServer(
	function(input, output, session) {
		# all k tab
		kdsdf <- get.distsum()
		output$kDistSumPlot <- renderPlot({
			distsum.plot(kdsdf)
		})
		output$kDistSumDeltaPlot <- renderPlot({
			distsum.delta.plot(kdsdf)
		})

		# choose k tab
		kclust <- reactive({
			env$cluster.ensemble[[input$k]]
		})
		output$k <- renderText({
			input$k
		})
		output$clusterSizePlot <- renderPlot({
			cluster.size.plot(kclust())
		})
		output$clusterOverviewPlot <- renderPlot({
			plot(kclust(), project=env$prcomp)
		})
		output$clusterProfileOverviewPlotArea <- renderPlot ({
			profilePlots <- lapply(1:input$k,
				function(cluster) {
					clust <- clusts()[clusts() == cluster]
					profile.data <- env$samples$log.ratio[names(clust),]
					makeClusterProfilePlot(profile.data, 
						title = cluster,
						y.range.adj = 1.5,
						simple = T
					)
				}
			)
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
			cl <- clust()
			if (!length(cl)) {
				return(NULL)
			}
			profile.data <- env$samples$log.ratio[names(cl),]
			makeClusterProfilePlot(profile.data,
				title = sprintf("K = %d : Cluster %d (%d genes)\nExpression profile",
					env$cluster.ensemble[[input$k]]@k, as.integer(input$cluster), length(names(cl))
				),
				focus = rowFocus,
				display.motif.gene.profile = c(1:env$meme.nmotifs)[
					c(
						input$displayMotif1GeneProfile,
						input$displayMotif2GeneProfile,
						input$displayMotif3GeneProfile,
						input$displayMotif4GeneProfile
					)
				],
				motifs <- env$meme.data[[input$k]][[input$cluster]]
			)
		})
		motifs <- reactive({
			if (is.null(input$k) || is.null(input$cluster)) {
				return(NULL)
			}
			return(env$meme.data[[input$k]][[input$cluster]])
		})
		output$clusterMotif1Summary <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 1) {
				return(NULL)
			}
			paste("E-value:", ms[[1]]$e.value, "- genes: ", length(ms[[1]]$positions$gene))
		})
		output$clusterMotif1Consensus <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 1) {
				return(NULL)
			}
			paste("Consesus:", ms[[1]]$consensus);
		})
		output$clusterMotif1Plot <- renderPlot({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 1) {
				return(NULL)
			}
			seqLogo(t(ms[[1]]$pssm))
		})
		output$clusterMotif2Summary <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 2) {
				return(NULL)
			}
			paste("E-value:", ms[[2]]$e.value, "- genes: ", length(ms[[2]]$positions$gene))
		})
		output$clusterMotif2Consensus <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 2) {
				return(NULL)
			}
			paste("Consesus:", ms[[2]]$consensus);
		})
		output$clusterMotif2Plot <- renderPlot({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 2) {
				return(NULL)
			}
			seqLogo(t(ms[[2]]$pssm))
		})
		output$clusterMotif3Summary <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 3) {
				return(NULL)
			}
			paste("E-value:", ms[[3]]$e.value, "- genes: ", length(ms[[3]]$positions$gene))
		})
		output$clusterMotif3Consensus <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 3) {
				return(NULL)
			}
			paste("Consesus:", ms[[3]]$consensus);
		})
		output$clusterMotif3Plot <- renderPlot({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 3) {
				return(NULL)
			}
			seqLogo(t(ms[[3]]$pssm))
		})
		output$clusterMotif4Summary <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 4) {
				return(NULL)
			}
			paste("E-value:", ms[[4]]$e.value, "- genes: ", length(ms[[4]]$positions$gene))
		})
		output$clusterMotif4Consensus <- renderText({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 4) {
				return(NULL)
			}
			paste("Consesus:", ms[[4]]$consensus);
		})
		output$clusterMotif4Plot <- renderPlot({
			ms <- motifs()
			if (is.null(ms) || length(ms) < 4) {
				return(NULL)
			}
			seqLogo(t(ms[[4]]$pssm))
		})
		output$clusterMembers <- renderDataTable({
			cl <- clust()
			if (!length(cl)) {
				return(NULL)
			}
			ns <- names(cl)

			# hierachical clustering of rows for row ordering
			# could this be precomputed?
			clustres <- env$samples$log.ratio[ns,]
			hclustres <- hclust(dist(clustres), method="complete")
			ns <- ns[hclustres$order]

			dir <- paste(
					dir.k.cluster(env$dir.output, env$cluster.ensemble[[input$k]]@k, input$cluster, make.dir = T),
					env$dir.motif.plots,
					sep = "/"
			)
			png.path = paste(dir, paste(ns, ".png", sep=""), sep="/")
			motif.img <- paste("<img src='", png.path, "' alt=''></img>", sep="")
			# get the list of sites
			msc <- env$meme.sites[[input$k]][[input$cluster]]
			# go through list and empty out image url for genes with no motif positions
			# if we need an image, check if it exists or set a flag to render all pngs
			render.pngs <- F
			for (n in 1:length(ns)) {
				if (dim(msc[msc$gene==ns[n], ])[1] == 0) {
					motif.img[n] <- ""
				} else if (identical(render.pngs, F) && !file.exists(png.path[n])) {
					render.pngs <- T
				}
			}
			# render the pngs if necessary
			if (identical(render.pngs, T)) {
				ms <- motifs()
				cat(sprintf("rendering %d pngs...", length(ns)))
				renderMotifPlots(dir,
					genes = ns,
					upstream.seqs = env$genes$upstream.seqs[ns,],
					upstream.start = env$upstream.start,
					upstream.end = env$upstream.end,
					motifs = ms,
					motif.colors = env$motif.colors,
					msc = msc
				)
				cat("done!\n")
			}
			data.frame("Locus tag" = ns, 
				"Product" = env$genes$annotations[ns, "product"],
				"Motif images" = motif.img,
				check.names = F
			)
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
					data.frame(locus.tag = ns, 
						product = env$genes$annotations[ns, "product"],
						env$samples$rpkm[ns,],
						env$samples$log.ratio[ns,]
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
					my.cluster.log.ratio <- env$samples$log.ratio[my.cluster.genes(),]
					other.log.ratio <- env$samples$log.ratio[!rownames(env$samples$log.ratio) %in% my.cluster.genes(),]
					switch(input$myClusterRecruitBy,
						min2centroid = {
								cmean<-apply(my.cluster.log.ratio, 2, mean)
								other.log.ratio$dist <- sqrt(rowSums(t(t(other.log.ratio)-cmean)^2))
								new.genes <- rownames(other.log.ratio[order(other.log.ratio$dist),])[1:input$myClusterRecruitN]
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
								new.genes <- rownames(other.log.ratio[order(other.log.ratio$dist),])[1:input$myClusterRecruitN]
							},
						random = {
								rrow <- sample(nrow(other.log.ratio), input$myClusterRecruitN)
								new.genes <- rownames(other.log.ratio[rrow,])
							},
						{	# default case, report a warning
							warning(paste("unhandled input$myClusterRecruitBy case:", input$myClusterRecruitBy))
						}
					)
					# send a client side message about the update to the textarea which will cascade the whole tab
					message <- list(
						value=paste(paste(my.cluster.genes(), collapse="\n"), paste(new.genes, collapse="\n"), sep="\n")
					)
					session$sendInputMessage("myClusterGenes", message)
				})
			}
		})
		my.cluster.genes <- reactive({
			input$myClusterGenesUpdateButton

			isolate({
			if (!is.null(input$myClusterGenes)) {
				genes <- unlist(strsplit(input$myClusterGenes, "\n", fixed=T))
				valid.genes <- genes %in% rownames(env$samples$log.ratio)
				return(genes[valid.genes])
			}
			return(default.my.cluster.genes)
			})
		})
		my.cluster.motifs <- reactive({
			mcg <- my.cluster.genes()

			# setup the training set data frame to be validated in memeParse
			training.set <- data.frame(length=env$genes$upstream.seqs[mcg, "uplength"], row.names = mcg)
			# remove any NA (i.e. the gene had no upstream sequence because of an overlap)
			training.set <- training.set[!is.na(training.set$length),"length", drop = F]

			if (length(training.set$length) > 1) {
				clust.seqs.upstream <- env$genes$upstream.seqs[mcg,]
				dir <- paste(env$dir.output, env$dir.my.cluster, sep="/")
                dir.create(dir, recursive = T, showWarnings = F)
                fasta.file <- paste(dir, env$file.upstream.fa, sep="/")
                if (file.exists(fasta.file)) {
                       file.remove(fasta.file);
                }
                for (k in 1:length(rownames(clust.seqs.upstream))) {
                        if (!is.na(clust.seqs.upstream$sequence[k])) {
                                cat(paste(">", rownames(clust.seqs.upstream)[k], "\n", sep="") , file=fasta.file, append=T)
                                cat(paste(clust.seqs.upstream$sequence[k], "\n", sep="") , file=fasta.file, append=T)
                        }
                }

				meme.file <- paste(dir, env$file.meme.txt, sep="/")
                meme.cmd <- paste(env$path.to.meme, fasta.file, "-nmotifs", env$meme.nmotifs, env$meme.base.args, "-oc", dir, "-bfile", 
					env$file.meme.bfile,
					">&", 
					meme.file
				)
				print(meme.cmd)
				system(meme.cmd)

				# load the meme output file
				motifs <- memeParse(meme.file, training.set)
				cat(sprintf("rendering %d pngs...", length(mcg)))
				meme.sites <- renderMotifPlots(paste(dir, env$dir.motif.plots, sep="/"),
					genes = mcg,
					upstream.seqs = env$genes$upstream.seqs[mcg,],
					upstream.start = env$upstream.start,
					upstream.end = env$upstream.end,
					motifs = motifs,
					motif.colors = env$motif.colors
				)
				cat("done!\n")
				return(list("meme.data" = motifs, "meme.sites" = meme.sites))
			}
			return(NULL)
		})
		output$myClusterGenesUI <- renderUI({
			tags$textarea(id="myClusterGenes", rows=6, cols=32, paste(my.cluster.genes(), collapse="\n"), style="display: block; margin-left: auto; margin-right: auto;")
		})
		output$myClusterProfilePlot <- renderPlot({
			if (is.null(input$myClusterSelectedRows)) {
				rowFocus <- F
			} else {
				rowFocus <- input$myClusterSelectedRows
			}
			profile.data <- env$samples$log.ratio[my.cluster.genes(),]
			makeClusterProfilePlot(profile.data = profile.data,
				title = "",
				y.range.adj = 1.5,
				simple = F,
				display.motif.gene.profile = c(1:4)[
					c(
					input$displayMyMotif1GeneProfile,
						input$displayMyMotif2GeneProfile,
						input$displayMyMotif3GeneProfile,
						input$displayMyMotif4GeneProfile
					)
				],
				motifs = my.cluster.motifs()$meme.data,
				motif.colors = env$motif.colors
			)
		})
		output$myClusterMembers <- renderDataTable({
			ns <- my.cluster.genes()
			mcm <- my.cluster.motifs()
			# could this be precomputed?
			if (length(ns) > 1) {
				clustres <- env$samples$log.ratio[ns,]
				hclustres <- hclust(dist(clustres), method="complete")
				ns <- ns[hclustres$order]
			}

			# put together path of the motif image for each gene
			dir <- paste(env$dir.output, env$dir.my.cluster, env$dir.motif.plots, sep="/")
			# use runif to append a random number to prevent all caching here
			motif.img <- paste("<img src='", 
					paste(env$url.prefix, dir, paste(ns, ".png?", runif(1, min=0, max=10), sep=""), sep="/"),
					"' alt=''></img>", sep="")
			# for genes with no sites, empty out the image url
			ms <- mcm$meme.sites
			if (length(ms)) {
				for (n in 1:length(ns)) {
					if (dim(ms[ms$gene==ns[n],])[1] == 0) {
						motif.img[n] <- ""
					}
				}
			} else {
				motif.img <- rep("", length(motif.img))
			}
			data.frame("Locus tag" = ns, 
				"Product" = env$genes$annotations[ns, "product"],
				"Motif images" = motif.img,
				check.names = F
			)
		}, options = list(paging = F))
		output$myClusterMotif1Summary <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 1) {
				return(NULL)
			}
			paste("E-value:", mcm$meme.data[[1]]$e.value, "- genes: ", length(mcm$meme.data[[1]]$positions$gene))
		})
		output$myClusterMotif1Consensus <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 1) {
				return(NULL)
			}
			paste("Consesus:", mcm$meme.data[[1]]$consensus);
		})
		output$myClusterMotif1Plot <- renderPlot({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 1) {
				return(NULL)
			}
			seqLogo(t(mcm$meme.data[[1]]$pssm))
		})
		output$myClusterMotif2Summary <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 2) {
				return(NULL)
			}
			paste("E-value:", mcm$meme.data[[2]]$e.value, "- genes: ", length(mcm$meme.data[[2]]$positions$gene))
		})
		output$myClusterMotif2Consensus <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 2) {
				return(NULL)
			}
			paste("Consesus:", mcm$meme.data[[2]]$consensus);
		})
		output$myClusterMotif2Plot <- renderPlot({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 2) {
				return(NULL)
			}
			seqLogo(t(mcm$meme.data[[2]]$pssm))
		})
		output$myClusterMotif3Summary <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 3) {
				return(NULL)
			}
			paste("E-value:", mcm$meme.data[[3]]$e.value, "- genes: ", length(mcm$meme.data[[3]]$positions$gene))
		})
		output$myClusterMotif3Consensus <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 3) {
				return(NULL)
			}
			paste("Consesus:", mcm$meme.data[[3]]$consensus);
		})
		output$myClusterMotif3Plot <- renderPlot({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 3) {
				return(NULL)
			}
			seqLogo(t(mcm$meme.data[[3]]$pssm))
		})
		output$myClusterMotif4Summary <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$mem.data) < 4) {
				return(NULL)
			}
			paste("E-value:", mcm$meme.data[[4]]$e.value, "- genes: ", length(mcm$meme.data[[4]]$positions$gene))
		})
		output$myClusterMotif4Consensus <- renderText({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 4) {
				return(NULL)
			}
			paste("Consesus:", mcm$meme.data[[4]]$consensus);
		})
		output$myClusterMotif4Plot <- renderPlot({
			mcm <- my.cluster.motifs()
			if (is.null(mcm) || length(mcm$meme.data) < 4) {
				return(NULL)
			}
			seqLogo(t(mcm$meme.data[[4]]$pssm))
		})
		output$myClusterMemeLog <- renderText({
			# register reactivity with the gene list text area and update button
			my.cluster.genes()
			dir <- paste(env$dir.output, env$dir.my.cluster, sep="/")
			meme.file <- paste(dir, env$file.meme.txt, sep="/")
			if (file.exists(meme.file)) {
				return(paste(readLines(meme.file), "\n"))
			}
			return(NULL)
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
		row.select <- union(
				grep(searchText, rownames(env$genes$annotations), ignore.case=T),
				grep(searchText, env$genes$annotations$product, ignore.case=T)
			)
		gcsr.clusts <- clusters(env$cluster.ensemble[[k]])
		gcsr.clust <- gcsr.clusts[row.select]
		data.frame("Locus tag" = names(gcsr.clust), 
			"Product" = env$genes$annotations[names(gcsr.clust), "product"], 
			"Cluster" = gcsr.clust,
			check.names = F
		)
}
