library(shiny)
library(shinyBS)

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
						class = "btn btn-success sbs-action-button btn-small",
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

myLikeModal <- function(tall, short) {
	myModal(
		id = sprintf("like%sReasonModal", tall),
		title = sprintf("Like %s reason", short),
		trigger = "",
		tags$p(HTML(sprintf("Why did you like this %s?", short)),
			tags$div(class = "row-fluid",
				textInput(sprintf("like%sReason", tall), "", ""),
				tags$head(tags$style(type="text/css", sprintf("#like%sReason {width: 510px}", tall)))
			)
		)
	)
}

shinyUI(
	navbarPage(
		title = env$study.title,
		theme = "bootstrap.css",

		tabPanel("All k", fluidPage(fluidRow(column(12,
			plotOutput("kDistSumPlot"),
			bsTooltip("kDistSumPlot", "Sum across all clusters of the distance to cluster centroids", "top"),
			plotOutput("kDistSumDeltaPlot"),
			bsTooltip("kDistSumDeltaPlot", "Difference between sum for a k - the sum for next lower k", "top")
		)))),

		tabPanel("Choose k", fluidPage(
			# K overview panes
			fluidRow(
				column(2, offset=5,
					selectInput("k", "Choose k", env$cluster.ensemble@k),
					bsTooltip("k", "Select a value for the number of clusters", "right")
				)
			),
			hr(),
			fluidRow(
				column(12,
					h3("Results from clustering with k = ", textOutput("k", container = span), 
						bsButton("kLikeButton", label = icon("thumbs-o-up")),
						bsTooltip("kLikeButton", "Save the workflow that got you to this k", "right"),
						align="center"
					),
					myLikeModal("K", "k")
				)
			),
			fluidRow(
				h4("Cluster sizes", align="center"),
				plotOutput("clusterSizePlot"),
				bsTooltip("clusterSizePlot", "Shows the number of genes in each cluster generated for this k", "top")
			),
			fluidRow(
				h4("Clustering overview", align="center"),
				p(align="center",
					plotOutput("clusterOverviewPlot", width="600px", height="600px"),
					bsTooltip("clusterOverviewPlot", "Identify clusters by color after projecting genes first two components of PCA", "right")
				)
			),
			fluidRow(
				h4("Cluster profile overview", align="center"),
				uiOutput("clusterProfileOverviewPlotArea"),
				bsTooltip("clusterProfileOverviewPlotArea", "Expression profile for genes in each cluster", "top")
			)
		)),

		tabPanel("Choose cluster", fluidPage(
			# for K, specific cluster view panes
			fluidRow(
				column(3, offset=3,
					uiOutput("clusterSelection"),
					bsTooltip("clusterSelection", "Select a cluster index to display", "left")
				),
				column(3, 
					p("Spreadsheet"),
					downloadButton('downloadClusterData', 'Download'),
					bsTooltip("downloadClusterData", "Download a spreadsheet for genes in cluster", "right")
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12, 
					h4("Cluster ", textOutput("cluster", container = span), 
						bsButton("clusterLikeButton", label = icon("thumbs-o-up")),
						bsTooltip("clusterLikeButton", "Save the workflow that got you to this cluster", "top"),
						align = "center"
					),
					myLikeModal("Cluster", "cluster")
				)
			),
			fluidRow(
				column(6, offset = 3,
					radioButtons("clusterProfilePlotTracks", "", names(env$samples$tracks), inline = T),
					bsTooltip("clusterProfilePlotTracks", "Which sample tracks should be displayed above the profile plot", "top"),
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
				column(12, h4("Motifs", align="center"),
					fluidRow(
						column(width = 3,
							p(tags$b("Motif 1"), style=paste("color:", env$motif.colors[1], sep=" "), 
								bsButton("clusterMotif1LikeButton", label = icon("thumbs-o-up")),
								bsTooltip("clusterMotif1LikeButton", "Save the workflow that generated this motif", "right"),
								align="center"
							),
							p(textOutput("clusterMotif1Summary", container = span), 
								align="center"
							),
							p(textOutput("clusterMotif1Consensus", container = span), 
								bsTooltip("clusterMotif1Consensus", "Top level of multilevel consensus", "right"),
								align="center"
							),
							column(width=9, offset=3,
								checkboxInput("displayMotif1GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif1Plot", height="180px")
						),
						column(width = 3,
							p(tags$b("Motif 2"), style=paste("color:", env$motif.colors[2], sep=" "),
								bsButton("clusterMotif2LikeButton", label = icon("thumbs-o-up")),
								bsTooltip("clusterMotif2LikeButton", "Save the workflow that generated this motif", "right"),
								align="center"
							),
							p(textOutput("clusterMotif2Summary", container = span), align="center"),
							p(textOutput("clusterMotif2Consensus", container = span), 
								bsTooltip("clusterMotif2Consensus", "Top level of multilevel consensus", "right"),
								align="center"
							),
							column(width = 9, offset = 3,
								checkboxInput("displayMotif2GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif2Plot", height="180px")
						),
						column(width = 3,
							p(tags$b("Motif 3"), style=paste("color:", env$motif.colors[3], sep=" "),
								bsButton("clusterMotif3LikeButton", label = icon("thumbs-o-up")),
								bsTooltip("clusterMotif3LikeButton", "Save the workflow that generated this motif", "left"),
								align="center"
							),
							p(textOutput("clusterMotif3Summary", container = span), align="center"),
							p(textOutput("clusterMotif3Consensus", container = span), 
								bsTooltip("clusterMotif3Consensus", "Top level of multilevel consensus", "left"),
								align="center"
							),
							column(width = 9, offset = 3,
								checkboxInput("displayMotif3GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif3Plot", height="180px")
						),
						column(width = 3,
							p(tags$b("Motif 4"), style=paste("color:", env$motif.colors[4], sep=" "),
								bsButton("clusterMotif4LikeButton", label = icon("thumbs-o-up")),
								bsTooltip("clusterMotif4LikeButton", "Save the workflow that generated this motif", "left"),
								align="center"
							),
							p(textOutput("clusterMotif4Summary", container = span), align="center"),
							p(textOutput("clusterMotif4Consensus", container = span), 
								bsTooltip("clusterMotif4Consensus", "Top level of multilevel consensus", "left"),
								align="center"
							),
							column(width = 9, offset = 3,
								checkboxInput("displayMotif4GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif4Plot", height="180px")
						)
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
				column(2, offset=5,
					textInput("searchText", "Search:", "methane monooxygenase"),
					bsTooltip("searchText", "Search for locus tags and annotations matching this text", "left")
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				h4(textOutput("clusterSearchResultSelectedRows", container = span), align="center"),
				dataTableOutput('clusterSearchResults'),
				bsTooltip("clusterSearchResults", "Genes matching your search, click to view cluster", "top")
			)
		)),

		tabPanel("My cluster", fluidPage(
			# allow user to create their own cluster and run meme
			fluidRow(
				column(3, offset=1,
					column(12, 
						h4("Enter CDS locus tags:", 
							align="center"
						),
						uiOutput("myClusterGenesUI"),
						bsTooltip("myClusterGenesUI", "Enter a list of locus tags to include in the cluster, one per line", "right")
					),
					column(2, offset=3,
						actionButton("myClusterGenesUpdateButton", "Update..."),
						bsTooltip("myClusterGenesUpdateButton", "Submit the gene list for analyses")
					)
				),
				column(6, offset=1,
					column(12, 
						column(12,
							h4("Recruit", align="center"),
							column(12,
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
							column(2, offset=4,
								actionButton("myClusterRecruitButton", "Recruit..."),
								bsTooltip("myClusterRecruitButton", "Search for additional transcripts")
							)
						)
					)
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12, 
					h4("My cluster", 
						bsButton("myClusterLikeButton", label = icon("thumbs-o-up")),
						bsTooltip("myClusterLikeButton", "Save the workflow that got you to this cluster", "top"),
						align = "center"
					)
				)
			),
			myLikeModal("MyCluster", "cluster"),
			fluidRow(
				column(6, offset = 3,
					radioButtons("myClusterProfilePlotTracks", "", names(env$samples$tracks), inline = T),
					bsTooltip("myClusterProfilePlotTracks", "Which sample tracks should be displayed above the profile plot", "top"),
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
				column(12, h4("Motifs", align="center"),
					fluidRow(
						column(width = 3,
							p(tags$b("Motif 1"), style=paste("color:", env$motif.colors[1], sep=" "),
								bsButton("myClusterMotif1LikeButton", label = icon("thumbs-o-up")),
								align="center"
							),
							p(textOutput("myClusterMotif1Summary", container = span), align="center"),
							p(textOutput("myClusterMotif1Consensus", container = span), align="center"),
							column(width = 9, offset = 3,
								checkboxInput("displayMyMotif1GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif1Plot", height="180px")
						),
						column(width = 3,
							p(tags$b("Motif 2"), style=paste("color:", env$motif.colors[2], sep=" "),
								bsButton("myClusterMotif2LikeButton", label = icon("thumbs-o-up")),
								align="center"
							),
							p(textOutput("myClusterMotif2Summary", container = span), align="center"),
							p(textOutput("myClusterMotif2Consensus", container = span), align="center"),
							column(width = 9, offset = 3,
								checkboxInput("displayMyMotif2GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif2Plot", height="180px")
						),
						column(width = 3,
							p(tags$b("Motif 3"), style=paste("color:", env$motif.colors[3], sep=" "),
								bsButton("myClusterMotif3LikeButton", label = icon("thumbs-o-up")),
								align="center"
							),
							p(textOutput("myClusterMotif3Summary", container = span), align="center"),
							p(textOutput("myClusterMotif3Consensus", container = span), align="center"),
							column(width = 9, offset = 3,
								checkboxInput("displayMyMotif3GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif3Plot", height="180px")
						),
						column(width = 3,
							p(tags$b("Motif 4"), style=paste("color:", env$motif.colors[4], sep=" "),
								bsButton("myClusterMotif4LikeButton", label = icon("thumbs-o-up")),
								align="center"
							),
							p(textOutput("myClusterMotif4Summary", container = span), align="center"),
							p(textOutput("myClusterMotif4Consensus", container = span), align="center"),
							column(width = 9, offset = 3,
								checkboxInput("displayMyMotif4GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif4Plot", height="180px")
						)
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
					h4("meme output", align="center")
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
				column(12, h4("Enter nucleic acid sequences in FASTA format:", align="center"))
			),
			fluidRow(
				column(12,
					tags$textarea(id="blastn.input", rows=8, cols=100, ">reverse transcriptase upstream motif\nACCGATGCGTGATACTGGGGCGGAG\n", style="width: 600px; height: 150px; display: block; margin-left: auto; margin-right: auto;")
				)
			),
			fluidRow(
				column(8, offset=3,
					radioButtons("blastnDatabase", "BLASTn database:",
						c("Upstream sequences (-150:2)" = "upstream",
							"Features (CDS, tRNA, rRNA, etc.)" = "feature",
							"Genome" = "genome"),
						inline=T
					)
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
				column(12, h4("Enter amino acid sequences in FASTA format:", align="center"))
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
				column(12, h4("Collection of workflows for interesting results:", align="center"))
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12, h2("This is under development", align="center"))
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12,
					dataTableOutput("likesTable")
				)
			)
		))
	)
)
