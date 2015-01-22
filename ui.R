library(shiny)

load("../cluster_analysis.5.RData")

motif.colors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")

shinyUI(
	navbarPage(
		title = env$study.title,

		tabPanel("All k",
			plotOutput("kDistSumPlot"),
			plotOutput("kDistSumDeltaPlot")
		),

		tabPanel("Choose k",
			# K overview panes
			fluidRow(
				column(2, offset=5,
					selectInput("k", "Choose k", env$cluster.ensemble@k)
				)
			),
			hr(),
			fluidRow(
				h3("Results from clustering with k = ", textOutput("k", container = span), align="center")
			),
			fluidRow(
				h4("Cluster sizes", align="center"),
				plotOutput("clusterSizePlot")
			),
			fluidRow(
				h4("Clustering overview", align="center"),
				p(align="center",
					plotOutput("clusterOverviewPlot", width="600px", height="600px")
				)
			),
			fluidRow(
				h4("Cluster profile overview", align="center"),
				plotOutput("clusterProfileOverviewPlotArea", height="1600px")
			)
		),

		tabPanel("Choose cluster",
			# for K, specific cluster view panes
			fluidRow(
				column(2, offset=3,
					uiOutput("clusterSelection")
				),
				column(2, offset=1,
					p("Spreadsheet"),
					downloadButton('downloadClusterData', 'Download')
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				column(12, h4("Cluster ", textOutput("cluster", container = span), align="center"),
					fluidRow(
						column(12,
							plotOutput("clusterProfilePlot")
						)
					)
				)
			),
			fluidRow(
				column(12, h4("Motifs", align="center"),
					fluidRow(
						column(width=3,
							p(tags$b("Motif 1"), style=paste("color:", motif.colors[1], sep=" "), align="center"),
							p(textOutput("clusterMotif1Summary", container = span), align="center"),
							#p(textOutput("clusterMotif1Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMotif1GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif1Plot", height="180px")
						),
						column(width=3,
							p(tags$b("Motif 2"), style=paste("color:", motif.colors[2], sep=" "), align="center"),
							p(textOutput("clusterMotif2Summary", container = span), align="center"),
							#p(textOutput("clusterMotif2Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMotif2GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif2Plot", height="180px")
						),
						column(width=3,
							p(tags$b("Motif 3"), style=paste("color:", motif.colors[3], sep=" "), align="center"),
							p(textOutput("clusterMotif3Summary", container = span), align="center"),
							#p(textOutput("clusterMotif3Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMotif3GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif3Plot", height="180px")
						),
						column(width=3,
							p(tags$b("Motif 4"), style=paste("color:", motif.colors[4], sep=" "), align="center"),
							p(textOutput("clusterMotif4Summary", container = span), align="center"),
							#p(textOutput("clusterMotif4Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMotif4GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif4Plot", height="180px")
						)
					)
				)
			),
			fluidRow(
				column(12,
					dataTableOutput('clusterMembers')
				)
			)
		),

		tabPanel("Find cluster",
			# for K, search for clusters by gene or product
			fluidRow(
				column(2, offset=5,
					textInput("searchText", "Search:", "methane monooxygenase")
				)
			),
			fluidRow(
				column(12, hr())
			),
			fluidRow(
				h4(textOutput("clusterSearchResultSelectedRows", container = span), align="center"),
				dataTableOutput('clusterSearchResults')
			)
		),

		tabPanel("My cluster",
			# allow user to create their own cluster and run meme
			fluidRow(
				column(3, offset=1,
					column(12, 
						h4("Enter CDS locus tags:", align="center"),
						uiOutput("myClusterGenesUI")
						#tags$textarea(id="myClusterGenes", rows=8, cols=32, "MBURv2_160308\nMBURv2_160304\nMBURv2_160312", style="display: block; margin-left: auto; margin-right: auto;")
					)
				),
				column(6, offset=1,
					column(12, 
						column(12,
							h4("Recruit", align="center"),
							column(12,
								column(6,
									selectInput("myClusterRecruitN", "Choose number of genes to recruit", 1:10)
								),
								column(6,
									radioButtons("myClusterRecruitBy", "Recruiting metric:",
										c("Minimum distance to centroid" = "min2centroid",
											"Minimum distance to any member" = "min2member",
											"Random" = "random"))
								)
							),
							column(2, offset=5,
								actionButton("myClusterRecruitButton", "Recruit...")
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
					plotOutput("myClusterProfilePlot")
				)
			),
			fluidRow(
				column(12, h4("Motifs", align="center"),
					fluidRow(
						column(width=3,
							p(tags$b("Motif 1"), style=paste("color:", motif.colors[1], sep=" "), align="center"),
							p(textOutput("myClusterMotif1Summary", container = span), align="center"),
							#p(textOutput("myClusterMotif1Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMyMotif1GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif1Plot", height="180px")
						),
						column(width=3,
							p(tags$b("Motif 2"), style=paste("color:", motif.colors[2], sep=" "), align="center"),
							p(textOutput("myClusterMotif2Summary", container = span), align="center"),
							#p(textOutput("myClusterMotif2Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMyMotif2GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif2Plot", height="180px")
						),
						column(width=3,
							p(tags$b("Motif 3"), style=paste("color:", motif.colors[3], sep=" "), align="center"),
							p(textOutput("myClusterMotif3Summary", container = span), align="center"),
							#p(textOutput("myClusterMotif3Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMyMotif3GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif3Plot", height="180px")
						),
						column(width=3,
							p(tags$b("Motif 4"), style=paste("color:", motif.colors[4], sep=" "), align="center"),
							p(textOutput("myClusterMotif4Summary", container = span), align="center"),
							#p(textOutput("myClusterMotif4Consensus", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMyMotif4GeneProfile", "Display profile", value=F)
							),
							plotOutput("myClusterMotif4Plot", height="180px")
						)
					)
				)
			),
			fluidRow(
				column(12,
					p(tags$br(), tags$b("NOTE: "), "rows in this table are not selectable", align="center")
				),
				column(12,
					dataTableOutput('myClusterMembers')
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
		),

		tabPanel("BLASTn",
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
				column(12, hr())
			)
		),

		tabPanel("BLASTp",
			# allow user to run BLASTp
			fluidRow(
				column(12, h4("Enter protein sequences in FASTA format:", align="center"))
			),
			fluidRow(
				column(12,
					tags$textarea(id="blastp.input", rows=8, cols=100, ">MBURv2_160304\nMARPLIQMALDSLDFEQTVALAEQVAPYVDIFEIGTPCIKYNGVGLVKELRQRFPDQLLL\nVDLKTMDAGEYEAAPFYAAGADICTVLGVSGLATIGGVIKAARAHNAEVQVDLINVPDKV\nECARESAKLGAQIVGVHTGLDAQAAGQTPFADLQAIADLGLNVRVSVAGGIKQATVQQVV\nASGASIIVVGAAIYGAPSPAEAAREIRQLVDAASA\n", style="width: 600px; height: 150px; display: block; margin-left: auto; margin-right: auto;")
				)
			),
			fluidRow(
				column(12, hr())
			)
		)
	)
)
