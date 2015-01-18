library(shiny)

load("../cluster_analysis.2.RData")

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
				plotOutput("clusterProfileOverviewPlot", height="800px")
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
			hr(),
			fluidRow(
				h4("Cluster ", textOutput("cluster", container = span), align="center"),
				plotOutput("clusterProfilePlot")
			),
			fluidRow(
				dataTableOutput('clusterMembers')
			)
		),

		tabPanel("Find cluster",
			# for K, search for clusters by gene or product
			fluidRow(
				column(2, offset=5,
					textInput("searchText", "Search:", "methane monooxygenase")
				)
			),
			hr(),
			fluidRow(
				h4(textOutput("clusterSearchResultSelectedRows", container = span), align="center"),
				dataTableOutput('clusterSearchResults')
			)
		)
	)
)
