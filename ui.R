library(shiny)
library(shinyBS)

# demo genes, could also just fetch first 4 (4 gets the kCDF display)
# but these 3 make for a nice demo
default.my.cluster.genes <- c("MBURv2_160308", "MBURv2_160304", "MBURv2_160312")

myModal <- function (id, title, trigger, ..., href) {
	mo <- tags$div(class = "modal sbs-modal fade", id = id, 
		'data-trigger' = trigger, tabindex="-1", role="dialog", 
		'aria-labelledby'=sprintf("%sLabel", id), 'aria-hidden'="true",
		tags$div(class = "modal-dialog", 
			tags$div(class = "modal-content", 
				tags$div(class = "modal-header", 
					tags$button(type = "button", class = "close", 'data-dismiss' = "modal", HTML("&times;")), 
					tags$h3(class = "modal-title", id=sprintf("%sLabel", id), title)
				), 
				body <- tags$div(class = "modal-body"),
				tags$div(class = "modal-footer", 
					tags$button(type = "button", 
						class = "btn btn-warning btn-small",
						"data-dismiss" = "modal", "Cancel"),
					tags$button(type = "button", 
						id = paste(id, "Close", sep=""), 
						class = "btn btn-success action-button btn-small",
						"data-dismiss" = "modal", "Like")
				) 
			)
		)
	)

	if (!missing(href)) {
		mo <- addAttribs(mo, 'data-remote' = href)
	} else {
		mo$children[[1]]$children[[1]]$children[[2]] <- tagAppendChildren(mo$children[[1]]$children[[1]]$children[[2]],
			list = list(...)
		)
	}
	return(mo)
}

myLikeModal <- function(id_prefix, description) {
	myModal(
		id = sprintf("%sLikeReasonModal", id_prefix),
		title = sprintf("Reason for liking", id_prefix),
		trigger = "",
		tags$p(HTML(sprintf("Why did you like this %s?", description)),
			tags$div(class = "row-fluid",
				textInput(sprintf("%sLikeReason", id_prefix), "", ""),
				tags$head(tags$style(type="text/css", sprintf("#%sLikeReason {width: 510px}", id_prefix)))
			)
		),
		tags$script(paste("$('#", sprintf("%sLikeReasonModal", id_prefix), "').on('hidden.bs.modal', function (e) {
				$('#", sprintf("%sLikeReason", id_prefix), "').val('');
				Shiny.onInputChange('", sprintf("%sLikeReason", id_prefix), "', '');
			})",
			sep = "")
		)
	)
}

columnMotif <- function(width, offset = 0, n, id_prefix="cluster") {
	column(width = width, offset = offset,
		fluidRow(
			tags$b(sprintf("Motif %d", n)), style=paste("color:", env$motif.colors[n], sep=" "),
			bsButton(size = "sm", sprintf("%sMotif%dLikeButton", id_prefix, n), label = icon("thumbs-o-up")),
			bsTooltip(sprintf("%sMotif%dLikeButton", id_prefix, n), "Save the workflow that generated this motif", "left")
		),
		fluidRow(
			myLikeModal(sprintf("%sMotif%d", id_prefix, n), "motif"),
			align = "left"
		),
		fluidRow(textOutput(sprintf("%sMotif%dSummary", id_prefix, n), container = span), align = "center"),
		fluidRow(
			textOutput(sprintf("%sMotif%dConsensus", id_prefix, n), container = span), 
			bsTooltip(sprintf("%sMotif%dConsensus", id_prefix, n), "Top level of multilevel consensus", "left")
		),
		fluidRow(
			checkboxInput(sprintf("%sDisplayMotif%dGeneProfile", id_prefix, n), "Display profile", value=F)
		),
		fluidRow(
			plotOutput(sprintf("%sMotif%dPlot", id_prefix, n), height="180px")
		),
		align = "center"
	)
}

inputTextarea <- function(inputId, value="", nrows, ncols) {
	tagList(
		singleton(tags$head(tags$script(src = "textarea.js"))),
		tags$textarea(id = inputId,
					class = "inputtextarea",
					rows = nrows,
					cols = ncols,
					as.character(value),
					style="display: block; margin-left: auto; margin-right: auto;"
			)
	)
}

shinyUI(
	navbarPage(
		title = env$study.title,
		theme = "bootstrap.css",

		tabPanel("All k", fluidPage(
			fluidRow(
				column(12,
					plotOutput("kDistSumPlot"),
					bsTooltip("kDistSumPlot", "Sum across all clusters of the distance to cluster centroids", "top"),
					plotOutput("kDistSumDeltaPlot"),
					bsTooltip("kDistSumDeltaPlot", "Difference between sum for a k - the sum for next lower k", "top")
				)
			)
		)),

		tabPanel("Choose k", fluidPage(
			# K overview panes
			fluidRow(
				column(2, offset = 5,
					selectInput("k", "Choose k", env$cluster.ensemble@k),
					bsTooltip("k", "Select a value for the number of clusters", "right")
				)
			),
			hr(),
			fluidRow(
				column(12,
					h3("Results from clustering with k = ", textOutput("k", container = span), 
						bsButton(size = "sm", "kLikeButton", label = icon("thumbs-o-up")),
						bsTooltip("kLikeButton", "Save the workflow that got you to this k", "right"),
						align = "center"
					),
					myLikeModal("k", "value of k")
				)
			),
			fluidRow(
				h3("Cluster sizes", align = "center"),
				plotOutput("clusterSizePlot"),
				bsTooltip("clusterSizePlot", "Shows the number of genes in each cluster generated for this k", "top")
			),
			fluidRow(
				h3("Clustering overview", align = "center"),
				p(align = "center",
					plotOutput("clusterOverviewPlot", width="600px", height="600px"),
					bsTooltip("clusterOverviewPlot", "Identify clusters by color after projecting genes first two components of PCA", "right")
				)
			),
			fluidRow(
				h3("Cluster profile overview", align = "center"),
				uiOutput("clusterProfileOverviewPlotArea"),
				bsTooltip("clusterProfileOverviewPlotArea", "Expression profile for genes in each cluster", "top")
			)
		)),

		tabPanel("Choose cluster", fluidPage(
			# for K, specific cluster view panes
			fluidRow(
				column(3, offset = 3,
					uiOutput("clusterSelection"),
					bsTooltip("clusterSelection", "Select a cluster index to display", "left")
				),
				column(3, 
					p("Spreadsheet"),
					downloadButton('downloadClusterData', 'Download'),
					bsTooltip("downloadClusterData", "Download a spreadsheet for genes in cluster", "right"),
					align = "center"
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12, 
					h3("Cluster ", textOutput("cluster", container = span), 
						bsButton(size = "sm", "clusterLikeButton", label = icon("thumbs-o-up")),
						bsTooltip("clusterLikeButton", "Save the workflow that got you to this cluster", "top"),
						align = "center"
					),
					myLikeModal("cluster", "cluster")
				)
			),
			fluidRow(
				column(6,
					radioButtons("clusterProfilePlotTracks", "", names(env$samples$tracks), inline = T),
					bsTooltip("clusterProfilePlotTracks", "Which sample tracks should be displayed above the profile plot", "top"),
					align = "center"
				),
				column(6,
					radioButtons("clusterProfilePlotSampleNames", "",
							c("Short sample names" = "Short",
								"Full names" = "Full",
								"Sample ID" = "ID"), 
							inline = T
					),
					bsTooltip("clusterProfilePlotSampleNames", "Should samples be labeled by a short name, full name or ID", "top"),
					align = "center"
				)
			),
			fluidRow(
				column(12,
					plotOutput("clusterProfilePlot"),
					bsTooltip("clusterProfilePlot", "Depicts any tracks you have displayed and the expression profile for genes in cluster.  The band indicates a 95% CI", "top")
				)
			),
			fluidRow(
				column(12, 
					fluidRow(
						h3("Motifs", align = "center")
					),
					fluidRow(
						columnMotif(width = 3, n = 1),
						columnMotif(width = 3, n = 2),
						columnMotif(width = 3, n = 3),
						columnMotif(width = 3, n = 4)
					)
				)
			),
			fluidRow(
				column(12,
					dataTableOutput('clusterMembers'),
					bsTooltip("clusterMembers", "Searchable, orderable table of genes in this cluster", "top")
				)
			)
		)),

		tabPanel("Find cluster", fluidPage(
			# for K, search for clusters by gene or product
			fluidRow(
				column(6, offset = 3,
					textInput("searchText", "Search:", "methane monooxygenase"),
					bsTooltip("searchText", "Search for locus tags and annotations matching this text", "left")
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				dataTableOutput('clusterSearchResults'),
				bsTooltip("clusterSearchResults", "Genes matching your search, click to view cluster", "top")
			)
		)),

		tabPanel("My cluster", fluidPage(
			# allow user to create their own cluster and run meme
			fluidRow(
				column(4, offset = 1,
					fluidRow(
						column(12, 
							h4("Enter CDS locus tags:", 
								align = "center"
							)
						)
					),
					fluidRow(
						column(12, 
							inputTextarea("myClusterGenes", paste(default.my.cluster.genes, collapse="\n"), 6, 32),
							bsTooltip("myClusterGenes", "Enter a list of locus tags to include in the cluster, one per line", "right"),
							align = "center"
						)
					),
					fluidRow(
						column(12,
							bsButton(size = "sm", "myClusterGenesUpdateButton", "Update"),
							bsTooltip("myClusterGenesUpdateButton", "Submit the gene list for analyses"),
							align = "center"
						)
					)
				),
				column(6,
					fluidRow(
						column(12,
							h4("Recruit", 
								align = "center"
							)
						)
					),
					fluidRow(
						column(6,
							selectInput("myClusterRecruitN", "Choose number of genes to recruit", 1:10),
							bsTooltip("myClusterRecruitN", "Choose the number of new genes to recruit to this cluster", "left")
						),
						column(6,
							radioButtons("myClusterRecruitBy", "Recruiting metric:",
								c("Minimum distance to centroid" = "min2centroid",
									"Minimum distance to any member" = "min2member",
									"Random" = "random")),
							bsTooltip("myClusterRecruitBy", "Select the method that should be used to recruit new genes to this cluster", "left")
						)
					),
					fluidRow(
						column(12,
							bsButton(size = "sm", "myClusterRecruitButton", "Recruit"),
							bsTooltip("myClusterRecruitButton", "Search for additional transcripts"),
							align = "center"
						)
					)
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12, 
					h3("My cluster", 
						bsButton(size = "sm", "myClusterLikeButton", label = icon("thumbs-o-up")),
						bsTooltip("myClusterLikeButton", "Save the workflow that got you to this cluster", "top"),
						align = "center"
					),
					myLikeModal("myCluster", "recruited cluster")
				)
			),
			fluidRow(
				column(6,
					radioButtons("myClusterProfilePlotTracks", "", names(env$samples$tracks), inline = T),
					bsTooltip("myClusterProfilePlotTracks", "Which sample tracks should be displayed above the profile plot", "top"),
					align = "center"
				),
				column(6,
					radioButtons("myClusterProfilePlotSampleNames", "",
							c("Short sample names" = "Short",
								"Full names" = "Full",
								"Sample ID" = "ID"), 
							inline = T
					),
					bsTooltip("myClusterProfilePlotSampleNames", "Should samples be labeled by a short name, full name or ID", "top"),
					align = "center"
				)
			),
			fluidRow(
				column(12,
					plotOutput("myClusterProfilePlot"),
					bsTooltip("myClusterProfilePlot", "Expression profile for genes in this cluster", "top")
				)
			),
			fluidRow(
				column(12, h3("Motifs", align = "center"),
					fluidRow(
						columnMotif(width = 3, n = 1, id_prefix = "myCluster"),
						columnMotif(width = 3, n = 2, id_prefix = "myCluster"),
						columnMotif(width = 3, n = 3, id_prefix = "myCluster"),
						columnMotif(width = 3, n = 4, id_prefix = "myCluster")
					)
				)
			),
			fluidRow(
				column(12,
					dataTableOutput('myClusterMembers'),
					bsTooltip("myClusterMembers", "Searchable, orderable table of genes in this cluster", "top")
				)
			),
			fluidRow(
				column(12,
					h3("meme output", align = "center")
				),
				column(12,
					tags$code(
						textOutput("myClusterMemeLog"),
						style="white-space: pre;"
					)
				)
			)
		)),

		tabPanel("BLASTn", fluidPage(
			# allow user to run BLASTn
			fluidRow(
				column(12, h3("Enter nucleic acid sequences in FASTA format:", align = "center"))
			),
			fluidRow(
				column(12,
					tags$textarea(id="blastn.input", rows=8, cols=100, ">reverse transcriptase upstream motif\nACCGATGCGTGATACTGGGGCGGAG\n", style="width: 600px; height: 150px; display: block; margin-left: auto; margin-right: auto;")
				)
			),
			fluidRow(
				column(12,
					radioButtons("blastnDatabase", "BLASTn database:",
						c("Upstream sequences (-150:2)" = "upstream",
							"Features (CDS, tRNA, rRNA, etc.)" = "feature",
							"Genome" = "genome"),
						inline=T
					),
					align = "center"
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12,
					dataTableOutput("blastnResults")
				)
			)
		)),

		tabPanel("BLASTp", fluidPage(
			# allow user to run BLASTp
			fluidRow(
				column(12, h3("Enter amino acid sequences in FASTA format:", align = "center"))
			),
			fluidRow(
				column(12,
					tags$textarea(id="blastp.input", rows=8, cols=100, ">MBURv2_160304\nMARPLIQMALDSLDFEQTVALAEQVAPYVDIFEIGTPCIKYNGVGLVKELRQRFPDQLLL\nVDLKTMDAGEYEAAPFYAAGADICTVLGVSGLATIGGVIKAARAHNAEVQVDLINVPDKV\nECARESAKLGAQIVGVHTGLDAQAAGQTPFADLQAIADLGLNVRVSVAGGIKQATVQQVV\nASGASIIVVGAAIYGAPSPAEAAREIRQLVDAASA\n", style="width: 600px; height: 150px; display: block; margin-left: auto; margin-right: auto;")
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12,
					dataTableOutput("blastpResults")
				)
			)
		)),

		tabPanel("Likes", fluidPage(
			fluidRow(
				column(12, h4("Result sets that were liked", align = "center"))
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12,
					dataTableOutput("likesTable")
				),
				singleton(
					tags$head(tags$script(src = "message-handler.js"))
				)
			)
		))
	)
)
