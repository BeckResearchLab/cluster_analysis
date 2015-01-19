library(shiny)

load("../cluster_analysis.4.RData")

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
				plotOutput("clusterOverviewPlot")
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
							p("Motif 1", align="center"),
							p(textOutput("clusterMotif1Summary", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMotif1GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif1Plot", height="180px")
						),
						column(width=3,
							p("Motif 2", align="center"),
							p(textOutput("clusterMotif2Summary", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMotif2GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif2Plot", height="180px")
						),
						column(width=3,
							p("Motif 3", align="center"),
							p(textOutput("clusterMotif3Summary", container = span), align="center"),
							column(width=6, offset=3,
								checkboxInput("displayMotif3GeneProfile", "Display profile", value=F)
							),
							plotOutput("clusterMotif3Plot", height="180px")
						),
						column(width=3,
							p("Motif 4", align="center"),
							p(textOutput("clusterMotif4Summary", container = span), align="center"),
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
		)
	)
)
