
	#removeResourcePath('Visualize')
	addResourcePath(prefix = 'Visualize', directoryPath = file.path(path.package("FuSViz"), "www"))
	#removeResourcePath('tmp')
	addResourcePath(prefix = 'tmp', directoryPath = tempdir())
	addResourcePath(prefix = 'Reference', directoryPath = file.path(path.package("FuSViz"), "tmp"))

	#---------------------------
	# Define dashboard manual 
	#---------------------------
	manual <- dashboardSidebar(
		width = 150,
		sidebarMenu(
			menuItem("Introduction", tabName = "main", icon = icon("mortar-board", verify_fa = FALSE)),
			menuItem("Table", tabName = "table", icon = icon("file-text-o", verify_fa = FALSE)),
			menuItem("Circle", tabName = "cir", icon = icon("venus-mars", verify_fa = FALSE)),
			menuItem("Linear", tabName = "igv", icon = icon("bar-chart", verify_fa = FALSE)),
			menuItem("Two-way", tabName = "fusion", icon = icon("line-chart", verify_fa = FALSE)),
			menuItem("Network", tabName = "net", icon = icon("th", verify_fa = FALSE)),
			menuItem("Document", tabName = "doc", icon = icon("map-o", verify_fa = FALSE))
		)
	)
	
	#----------------------
	# Introduction page 
	#----------------------
	bodyIntro <- tabItem(tabName = "main", value="main_panel",
					tags$head(tags$style(".box-header h3 {font-weight: bold; }")),
					fluidRow(
						box(
							tags$head(tags$link(rel="shortcut icon", href="Visualize/logo.ico")),
							title = "Welcome to FuSViz Shiny App!", width = 12, status = "primary",
							div(style="text-align: justify;", HTML("<b>FuSViz</b> is a web-based tool to visualize, interpret and prioritize genomic/transcriptomic structural variations (SVs) of multiple samples. It provides multiple solutions and an interactive interface to investigate the prevalence and recurrence of SVs and their relevant genes in a cohort of cases. This tool is used for combining SVs called from DNA-seq and RNA-seq data together to illustrate a biological consequence of SVs to cancer or genetic-disease related genes in context of various genomic annotations."))
						)
					),
					fluidRow(
						box(
							title = "Circular module", width = 3, status = "warning",
							p(style="text-align: justify;", "The circular module visualizes the distribution of SVs on a whole genome scale and is informative for viewing intra-chrom SVs with a large distance and inter-chrom SVs."),
							img(src='Visualize/Circle.png', align = "center", width="72%")
						),
						box(
							title = "Linear module", width = 3, status = "warning",
							p(style="text-align: justify;", "The linear module displays a reference genome and SVs horizontally with various types of customized annotation tracks added as parallel lines. It provides a high resolution for visualizing intra-chrom SVs."),
							img(src='Visualize/Linear.png', align = "center", width="102%")
						),
						box(
							title = "Two-way module", width = 3, status = "warning",
							p(style="text-align: justify;", "The two-way module is specific for displaying fusion genes or transcripts in a single panel, for which two distant genomic regions involved in a fusion event are shown."),
							img(src='Visualize/Two-way.png', align = "center", width="90%")
						),
						box(
							title = "Network module", width = 3, status = "warning",
							p(style="text-align: justify;", "The network module is used to identify hubs (i.e. a component with a high-degree connection) in SV interaction network and infer structural variation complexity of genes in a sub-network."),
							img(src='Visualize/Network.png', align = "center", width="89%")
						)
					),
					fluidRow(					 
						box(title = "Import genomic and transcriptomic annotations", width = 12, status = "info", solidHeader = TRUE,
							fluidRow(
								column(3, h5(icon("database", verify_fa = FALSE), "Annotation data", style = "font-weight: bold; font-size: 16px;")),
								column(4, selectizeInput("genome", label = "Genome version", choices = c("hg19", "hg38", "hg19_offline", "hg38_offline"), multiple = F, 
											options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
								column(3, style = "margin-top: 25px;", actionButton("Import_genome_data", "Import", icon = icon("bolt", verify_fa = FALSE), width = "100%"))
							)
						)
					),
					fluidRow(
						box(title = "Upload structural variation file", width = 12, status = "info", solidHeader = TRUE,
							 fluidRow(
								column(3, h5(icon("upload", verify_fa = FALSE), "SV callings from RNA-seq data", style = "font-weight: bold; font-size: 16px;")),
								column(4, fileInput('file_rna_data', label = "Upload file in txt format", accept = c('text/csv', 'text/comma-separated-values', 
										'text/tab-separated-values', '.csv', '.tsv', '.txt'), placeholder = "txt, csv or tsv format")),
								column(4, radioButtons('sep_rna_file', label = HTML('<p>Separator <a href="https://fusviz.s3.eu-north-1.amazonaws.com/RNA_SV_example.txt" target="_blank">(See and Download an exampale)</a></p>'), choices = c(Comma = ',', Tab = '\t', Semicolon = ';'), selected = '\t', inline = T))
							),
							hr(tags$style("border-top: 1px solid #000000; color: black; background-color: white")),
							fluidRow(
								column(3, h5(icon("upload", verify_fa = FALSE), "SV callings from DNA-seq data", style = "font-weight: bold; font-size: 16px;")),
								column(4, fileInput('file_dna_data', label = "Upload file in bedpe format", accept = c('text/csv', 'text/comma-separated-values', 
										'text/tab-separated-values', '.csv', '.tsv', '.txt'), placeholder = "bedpe format")),
								column(4, radioButtons('sep_dna_file', label = HTML('<p>Separator <a href="https://fusviz.s3.eu-north-1.amazonaws.com/DNA_SV_example.txt" target="_blank">(See and Download an exampale)</a></p>'), choices = c(Comma = ',', Tab = '\t', Semicolon = ';'), selected = '\t', inline = T))
							)
						)
					),
					fluidRow(
						box(title = "Upload mutation profile file (optional)", width = 12, status = "info", solidHeader = TRUE,
							fluidRow(
								column(3, h5(icon("upload", verify_fa = FALSE), "Mutation variants (SNVs and Indels)", style = "font-weight: bold; font-size: 16px;")),
								column(4, fileInput('file_maf_data', label = "Upload file in MAF format", accept = c('text/maf', '.maf'), placeholder = "maf format")),
								column(4, radioButtons('sep_maf_file', label = HTML('<p>Separator <a href="https://fusviz.s3.eu-north-1.amazonaws.com/TCGA.PRAD.mutect.somatic.maf" target="_blank">(See and Download an exampale)</a></p>'), choices = c(Tab = '\t'), selected = '\t', inline = T))
							)
						)
					)
				)
	
	#------------------
	# table vis
	#------------------
	bodyTable <- tabItem(tabName = "table", value="table_vis",
					h4("Table overview"),
					tags$head(tags$style("td[data-type='factor'] input { width: 100px !important; }")),
					tags$head(tags$style("table.dataTable .selectize-input { font-size: 11px; }")), # adjust font-size of select-input as 11px
					tags$head(tags$style("table.dataTable .selectize-dropdown { font-size: 11px; }")), # adjust font-size of select-dropdown as 11px
					tabsetPanel(id = "tabset-data",
						tabPanel("SV from RNA-seq",
							br(),
							fluidRow(
								# DT and dataTable library implementary here for RNA-seq
								column(12, shinycssloaders::withSpinner(DT::dataTableOutput("RNAcontents", height = 450), type = 5, size = 0.7)) 
							),
							br(),
							fluidRow(
								box(title = "Feature summary", status="warning", width=12, height = "100%", collapsible = TRUE, collapsed = TRUE,
									tabBox(id = "Feature_rna", height = "100%", width = "100%",
										tabPanel("SV relevant gene wordcloud",
											fluidRow(class = "rna_box",
												column(12, wordcloud2::wordcloud2Output("wordcloud_rna", width = "100%", height = 600)),
												tags$script(HTML(
												"$(document).on('mousemove', '#wordcloud_rna', function(){",
												"if ( document.getElementById('wordcloud_rna').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args') != null ) {",
												"var n = document.getElementById('wordcloud_rna').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args').lastIndexOf(':');",
												"var count = document.getElementById('wordcloud_rna').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args').substring(n + 1);",
												"var freq = count.slice(0, -1);",
												"var new_freq = Math.pow(freq, 2);",
												"var m = document.getElementById('wordcloud_rna').getElementsByClassName('wcSpan')[0].innerHTML.lastIndexOf(':');",
												"var job = document.getElementById('wordcloud_rna').getElementsByClassName('wcSpan')[0].innerHTML.substr(0, m);",
												"document.getElementById('wordcloud_rna').getElementsByClassName('wcSpan')[0].innerHTML = job + ':' + Math.round(new_freq);",
												"}});"
												))
											),
											fluidRow(class = "rna_control",
												column(4, sliderInput(inputId = "word_size_rna", label = "Word size", value=0.45, min=0,001, max=1, step=0.05)),
												column(4, sliderInput(inputId = "gene_freq_rna", label = "Gene freq", value=1, min=1, max=10, step=1)),
												column(4, selectInput(inputId="word_shape_rna", label="Word shape", choices=c("circle","square","cardioid"), selected="square", multiple=F))
											),
											tags$head(tags$style(HTML("
											.rna_box{height:620px;}
											.rna_control{height:-5px;}
											#word_size_rna{height:20px;}
											#gene_freq_rna{height:20px;}
											#word_shape_rna{height:20px;}"
											)))
										),
										tabPanel("SV distribution across samples",
											fluidRow(class = "rna_hist_sample",
												column(12, shinycssloaders::withSpinner(plotOutput("rna_hist_1", width = "100%", height = 150), type = 5, size = 0.7)),
												column(12, plotOutput("rna_hist_2", width = "100%", height = 150)),
												column(12, plotOutput("rna_hist_3", width = "100%", height = 150)),
												column(12, plotOutput("rna_hist_4", width = "100%", height = 150))
											),
											fluidRow(class = "rna_control_sample",
												column(3, sliderInput(inputId = "rna_hist_size", label = "Font size", value=0.6, min=0,001, max=1, step=0.05)),
												column(3, sliderInput(inputId = "rna_hist_pos", label = "X-axis pos", value=-20, min=-100, max=0, step=2)),
												column(3, sliderInput(inputId = "rna_hist_angle", label = "X-axis rotate", value=90, min=0, max=90, step=1)),
												column(3, sliderInput(inputId = "rna_hist_space", label = "Bar space", value=0.1, min=0, max=1, step=0.05))
											),
											tags$head(tags$style(HTML("
											.rna_hist_sample{height:620px;}
											.rna_control_sample{height:-5px;}
											#rna_hist_size{height:20px;}
											#rna_hist_pos{height:20px;}
											#rna_hist_angle{height:20px;}
											#rna_hist_space{height:20px;}"
											)))
										)
									)
								)
							)
						),
						tabPanel("SV from DNA-seq",
							br(),
							fluidRow(
								# DT and dataTable library implementary here for DNA-seq
								column(12, shinycssloaders::withSpinner(DT::dataTableOutput("DNAcontents", height = 450), type = 5, size = 0.7))
							),
							br(),
							fluidRow(
								box(title = "Feature summary", status="warning", width=12, height = "100%", collapsible = TRUE, collapsed = TRUE,
									tabBox(id = "Feature_dna", height = "100%", width = "100%",
										tabPanel("SV relevant gene wordcloud",
											fluidRow(class = "dna_box",
												column(12, wordcloud2::wordcloud2Output("wordcloud_dna", width = "100%", height = 800)),
												tags$script(HTML(
												"$(document).on('mousemove', '#wordcloud_dna', function(){",
												"if ( document.getElementById('wordcloud_dna').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args') != null ) {",
												"var n = document.getElementById('wordcloud_dna').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args').lastIndexOf(':');",
												"var count = document.getElementById('wordcloud_dna').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args').substring(n + 1);",
												"var freq = count.slice(0, -1);",
												"var new_freq = Math.pow(freq, 2);",
												"var m = document.getElementById('wordcloud_dna').getElementsByClassName('wcSpan')[0].innerHTML.lastIndexOf(':');",
												"var job = document.getElementById('wordcloud_dna').getElementsByClassName('wcSpan')[0].innerHTML.substr(0, m);",
												"document.getElementById('wordcloud_dna').getElementsByClassName('wcSpan')[0].innerHTML = job + ':' + Math.round(new_freq);",
												"}});"
												))
											),
											fluidRow(class = "dna_control",
												column(4, sliderInput(inputId = "word_size_dna", label = "Word size", value=0.2, min=0.001, max=1, step=0.05)),
												column(4, sliderInput(inputId = "gene_freq_dna", label = "Gene freq", value=1, min=1, max=10, step=1)),
												column(4, selectInput(inputId="word_shape_dna", label="Word shape", choices=c("circle","square","cardioid"), selected="square", multiple=F))
											),
											tags$head(tags$style(HTML("
											.dna_box{height:820px;}
											.dna_control{height:-5px;}
											#word_size_dna{height:20px;}
											#gene_freq_dna{height:20px;}
											#word_shape_dna{height:20px;}"
											)))
										),
										tabPanel("SV distribution across samples",
											fluidRow(class = "dna_hist_sample",
												column(12, plotOutput("category", width = "100%", height = 40)),
												column(12, shinycssloaders::withSpinner(plotOutput("dna_hist_1", width = "100%", height = 190), type = 5, size = 0.7)),
												column(12, plotOutput("dna_hist_2", width = "100%", height = 190)),
												column(12, plotOutput("dna_hist_3", width = "100%", height = 190)),
												column(12, plotOutput("dna_hist_4", width = "100%", height = 190))
											),
											fluidRow(class = "dna_control_sample",
												column(3, sliderInput(inputId = "dna_hist_size", label = "Font size", value=0.6, min=0,001, max=1, step=0.05)),
												column(3, sliderInput(inputId = "dna_hist_pos", label = "X-axis pos", value=-100, min=-500, max=0, step=5)),
												column(3, sliderInput(inputId = "dna_hist_angle", label = "X-axis rotate", value=90, min=0, max=90, step=1)),
												column(3, sliderInput(inputId = "dna_hist_space", label = "Bar space", value=0.1, min=0, max=1, step=0.05))
											),
											tags$head(tags$style(HTML("
											.dna_hist_sample{height:820px;}
											.dna_control_sample{height:-5px;}
											#dna_hist_size{height:20px;}
											#dna_hist_pos{height:20px;}
											#dna_hist_angle{height:20px;}
											#dna_hist_space{height:20px;}"
											)))
										)
									)
								)
							)
						),
						tabPanel("Mutation profile",
							br(),
							fluidRow(
								#DT and dataTable library implementary here for SNVs and indels
								column(12, shinycssloaders::withSpinner(DT::dataTableOutput("DNAMutations", height = 450), type = 5, size = 0.7)) 
							),
							br(),
							fluidRow(
								box(title = "Feature summary", status="warning", width=12, height = "100%", collapsible = TRUE, collapsed = TRUE,
									tabBox(id = "Feature_mut", height = "100%", width = "100%",
										tabPanel("Mutated gene wordcloud",
											fluidRow(class = "mut_box",
												column(12, wordcloud2::wordcloud2Output("wordcloud_mut", width = "100%", height = 800)),
												tags$script(HTML(
												"$(document).on('mousemove', '#wordcloud_mut', function(){",
												"if ( document.getElementById('wordcloud_mut').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args') != null ) {",
												"var n = document.getElementById('wordcloud_mut').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args').lastIndexOf(':');",
												"var count = document.getElementById('wordcloud_mut').getElementsByClassName('wcSpan')[0].getAttribute('data-l10n-args').substring(n + 1);",
												"var freq = count.slice(0, -1);",
												"var new_freq = Math.pow(freq, 2);",
												"var m = document.getElementById('wordcloud_mut').getElementsByClassName('wcSpan')[0].innerHTML.lastIndexOf(':');",
												"var job = document.getElementById('wordcloud_mut').getElementsByClassName('wcSpan')[0].innerHTML.substr(0, m);",
												"document.getElementById('wordcloud_mut').getElementsByClassName('wcSpan')[0].innerHTML = job + ':' + Math.round(new_freq);",
												"}});"
												))
											),
											fluidRow(class = "mut_control",
												column(4, sliderInput(inputId = "word_size_mut", label = "Word size", value=0.2, min=0.001, max=1, step=0.05)),
												column(4, sliderInput(inputId = "gene_freq_mut", label = "Gene freq", value=1, min=1, max=10, step=1)),
												column(4, selectInput(inputId="word_shape_mut", label="Word shape", choices=c("circle","square","cardioid"), selected="square", multiple=F))
											),
											tags$head(tags$style(HTML("
											.mut_box{height:820px;}
											.mut_control{height:-5px;}
											#word_size_mut{height:20px;}
											#gene_freq_mut{height:20px;}
											#word_shape_mut{height:20px;}"
											)))
										),
										tabPanel("Correlation of mutation and SV (RNA-seq) burden",
											column(12, shinycssloaders::withSpinner(plotOutput("rna_mut_cor", width="100%", height=600, brush=brushOpts(id="rna_brush")), type=5, size=0.7)),
											column(12, h4("Selected samples"), verbatimTextOutput("rna_mut_info"))
										),
										tabPanel("Correlation of mutation and SV (DNA-seq) burden",
											column(12, shinycssloaders::withSpinner(plotOutput("dna_mut_cor", width="100%", height=600, brush=brushOpts(id="dna_brush")), type=5, size=0.7)),
											column(12, h4("Selected samples"), verbatimTextOutput("dna_mut_info"))
										)
									)
								)
							)
						),
						tabPanel("Drug target info", 
							br(),
							fluidRow(
								#DT and dataTable library implementary here for Drug target info
								column(12, shinycssloaders::withSpinner(DT::dataTableOutput("drug", height = 600), type = 5, size = 0.7))
							)
						)
					)
				)
	#------------------
	# circular module
	#------------------
	bodyCircle <- tabItem(tabName = "cir", value="cir_plot",
					h4("Circular module"),
					fluidRow( 
						box(status="warning", width=3,
							tabBox(id = "circle_plot_tab", height = "100%", width = "100%",
								tabPanel("RNA_SV_panel",
									selectizeGroupUI(id = "circle_set", inline = FALSE,
										params = list(
											chr = list(inputId = "chr", title = "Chrom", placeholder = 'select'),
											gene = list(inputId = "gene", title = "Gene", placeholder = 'select'),
											name = list(inputId = "name", title = "Sample", placeholder = 'select')
									)),
									sliderInput(inputId = "circle_split", label = "Num of split reads", value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_span", label = "Num of span reads", value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_num", label = "Num of samples", value=0, min=0, max=50, step=1),
									checkboxInput("cirmutcheck", "Load mutation data", FALSE),
									conditionalPanel(condition = "input.cirmutcheck == 1", 
										selectInput("mut_type", label = "Mutation type", choices = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", 
											"Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", 
											"De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "3'UTR", "5'UTR", "3'Flank", "Targeted_Region", 
											"Silent", "Intron", "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA"), 
											selected = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", 
											"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", 
											"Start_Codon_SNP", "Stop_Codon_Del"), multiple = T)),
									conditionalPanel(condition = "input.cirmutcheck == 1",
										selectizeInput("mut_sample", label = "Mutation sample", choices = NULL, multiple = T, options = list(placeholder = 'select', 
											onInitialize = I('function() { this.setValue(""); }')))),
									br(),
									fluidRow(column(12, actionButton("circos_plot", "Plot / Refresh"))),
									br(),
									fluidRow(column(12, downloadButton(outputId = "CircleDownRNA", label = "Download circular plot")))
								),
								tabPanel("DNA_SV_panel",
									selectizeGroupUI(id = "circle_set_dna", inline = FALSE,
										params = list(
											chr = list(inputId = "chr", title = "Chrom", placeholder = 'select'),
											gene = list(inputId = "gene", title = "Gene", placeholder = 'select'),
											name = list(inputId = "name", title = "Sample", placeholder = 'select')
									)),
									sliderInput(inputId = "circle_split_dna", label = "Num of split reads", value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_span_dna", label = "Num of span reads", value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_num_dna", label = "Num of samples", value=0, min=0, max=50, step=1),
									sliderInput(inputId = "circle_dis_dna", label = "Dist intra-chrom", value=1000000, min=0, max=10000000, step=100000),
									checkboxInput("cirmutcheck_dna", "Load mutation data", FALSE),
									conditionalPanel(condition = "input.cirmutcheck_dna == 1", 
										selectInput("mut_type_dna", label = "Mutation type", choices = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", 
											"Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", 
											"De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "3'UTR", "5'UTR", "3'Flank", "Targeted_Region", 
											"Silent", "Intron", "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA"),
											selected = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", 
											"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", 
											"Start_Codon_SNP", "Stop_Codon_Del"), multiple = T)),
									conditionalPanel(condition = "input.cirmutcheck_dna == 1",
										selectizeInput("mut_sample_dna", label = "Mutation sample", choices = NULL, multiple = T, options = list(placeholder = 'select', 
											onInitialize = I('function() { this.setValue(""); }')))),
									br(),
									fluidRow(column(12, actionButton("circos_plot_dna", "Plot / Refresh"))),
									br(),
									fluidRow(column(12, downloadButton(outputId = "CircleDownDNA", label = "Download circular plot")))
								)
							)
						),
						box(status="warning", width=9,
							tabBox(id = "fusion_plot_page", height = "100%", width = "100%",
								tabPanel("RNA_SV_circular_plot",
									shinycssloaders::withSpinner(BioCircos::BioCircosOutput("circle_1", height=700), type = 5, size = 0.7)
								),
								tabPanel("DNA_SV_circular_plot",
									shinycssloaders::withSpinner(BioCircos::BioCircosOutput("circle_2", height=700), type = 5, size = 0.7)
								)
							)
						)
					)
				)
	
	#---------------------
	# Linear module 
	#---------------------
	bodyLinear <- tabItem(tabName = "igv", value="igv_plot",
					h4("Linear module"),
					sidebarLayout(
						sidebarPanel(width = 3,
							tabsetPanel(id = "igv-panel-control",
								tabPanel("DNA_SV",
									fluidRow(
										column(6, numericInput("dna_split_bedpe", label = "Num_split", value=0, min=0, max=50)),
										column(6, numericInput("dna_span_bedpe", label = "Num_span", value=0, min=0, max=50)),
										column(6, numericInput("dna_dismin_bedpe", label = "Min_dist", value=0, min=0, max=500000000)),
										column(6, numericInput("dna_dismax_bedpe", label = "Max_dist", value=500000000, min=0, max=500000000)),
										column(4, selectizeInput("dna_type_bedpe", label = "SV_type", choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(8, selectizeInput("dna_sample_bedpe", label = "Sample", choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(12, actionButton("addTrackButtonBedPe", "Load DNA SV track in bedpe"))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(6, numericInput("dna_split_seg", label = "Num_split", value=0, min=0, max=50)),
										column(6, numericInput("dna_span_seg", label = "Num_span", value=0, min=0, max=50)),
										column(3, div(class="no-spinners", textInput("dna_seg_chr", label = "Chrom", value = "", placeholder = "chr1..chrX,chrY,chrM"))),
										tags$style(HTML("#dna_seg_chr { width: 55px; font-size:12px; }")),
										column(4, div(class="no-spinners", numericInput("dna_seg_start", label = "Start", value=0, min=0, max=500000000))),
										tags$style(HTML("#dna_seg_start { width: 80px; font-size:12px; }")),
										column(5, div(class="no-spinners", numericInput("dna_seg_end", label = "End", value=0, min=0, max=500000000))),
										tags$style(HTML("#dna_seg_end { width: 80px; font-size:12px; }")),
										tags$style(HTML(".no-spinners { display: inline-block; vertical-align:top; } 
														.no-spinners input[type=number] { -moz -appearance:textfield; } 
														.no-spinners input[type=number]::{ -moz-appearance:textfield; }
														.no-spinners input[type=number]::-webkit-outer-spin-button, 
														.no-spinners input[type=number]::-webkit-inner-spin-button { -webkit-appearance: none; margin: 0; } ")),
										column(8, selectizeInput("dna_sample_seg", label = "Sample", choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(12, actionButton("addTrackButtonSeg", "Load DNA SV track in seg"))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(4, numericInput("dna_split_bed", label = "Num_split", value=0, min=0, max=50)),
										column(4, numericInput("dna_span_bed", label = "Num_span", value=0, min=0, max=50)),
										column(4, selectizeInput("dna_type_bed", label = "SV_type", choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										tags$style(HTML("#dna_type_bed+ div>.selectize-input { width: 70px; } 
														#dna_type_bed+ div>.selectize-input input { width: 70px; } 
														#dna_type_bed+ .selectize-dropdown-content { width: 70px; }")),
										column(6, numericInput("dna_dismin_bed", label = "Min_dist", value=0, min=0, max=500000000)),
										column(6, numericInput("dna_dismax_bed", label = "Max_dist", value=500000000, min=0, max=500000000)),
										column(12, actionButton("addTrackButtonBed", "Load DNA breakpoints in bed"))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(12, textInput(inputId = "genomicPos", label = "Genome coordinate", value = "")),
										column(6, actionButton("GetPos", label = "Show coordinate")),
										column(6, actionButton("MovePos", label = "Clear"))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
#										column(12, fileInput('load_file', label = "Upload compressed BED, VCF or GTF", accept = c('.gz', '.tbi'), 
#													multiple = T, placeholder = "*.gz with index *.tbi format")),
										column(12, tags$p(tags$b("Upload local alignments (indexed ", tags$u("BAM"), " and ", tags$u("CRAM"), ") / annotations (bgzipped and tabixed", 
															tags$u("VCF"), ", ", tags$u("BED"), " and ", tags$u("GTF"), ") / genome references (indexed ", tags$u("fasta"), " and cytoband in ", 
															tags$u("txt"), " format - only usage in offline mode):"), " choose file and index!")),
										br(),
										column(12, HTML("<input id='uploadfile' class='hidden' type='file' multiple='true' accept='.bam,.bai,.cram,.crai,.gz,.tbi,.fasta,.fai,.txt' onchange=\"load()\"/>
														<label for='uploadfile'><img src='Visualize/Button.png' alt='Upload img' height='100' width='200'></label>")),
										br(),
										column(12, HTML("<div id='filename'></div>")),
										tags$script(src = "Visualize/Upload_file.js")
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(12, tags$p(tags$b("* Only usage in offline mode - upload genome references before click"))),
										column(12, actionButton("offlineTrackButton", "Load gene track (offline)"))
									)
								),
								tabPanel("RNA_SV",
									fluidRow(
										column(6, numericInput("rna_split_bedpe", label = "Num_split", value=0, min=0, max=50)),
										column(6, numericInput("rna_span_bedpe", label = "Num_span", value=0, min=0, max=50)),
										column(6, numericInput("rna_dismin_bedpe", label = "Min_dist", value=0, min=0, max=500000000)),
										column(6, numericInput("rna_dismax_bedpe", label = "Max_dist", value=500000000, min=0, max=500000000)),
										column(8, selectizeInput("rna_sample_bedpe", label = "Sample", choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(12, actionButton("addTrackButtonRNABedpe", "Load RNA SV track in bedpe"))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(4, numericInput("rna_split_bed", label = "Num_split", value=0, min=0, max=50)),
										column(4, numericInput("rna_span_bed", label = "Num_span", value=0, min=0, max=50)),
										column(4, selectizeInput("rna_type_bed", label = "SV_type", choices = c("Inter", "Intra"), multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(6, numericInput("rna_dismin_bed", label = "Min_dist", value=0, min=0, max=500000000)),
										column(6, numericInput("rna_dismax_bed", label = "Max_dist", value=500000000, min=0, max=500000000)),
										column(12, actionButton("addTrackButtonRNABed", "Load RNA breakpoints in bed"))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(12, textInput(inputId = "BAM", label = "Alignment URL", value = NULL)),
										column(12, textInput(inputId = "BAMindex", label = "Alignment index URL", value = NULL)),
										column(12, actionButton("addTrackButtonBAM", label = NULL, style = "width: 115px; height: 71px;
											background: url('Visualize/cloud_upload.png');  background-size: cover; background-position: center;"))
									)
								),
								tabPanel("Mut",
									fluidRow(
										selectizeInput("mut_type_line", label = "Mut_type", choices = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
											"Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", 
											"De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "3'UTR", "5'UTR", 
											"3'Flank", "Targeted_Region", "Silent", "Intron", "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA"), multiple = T, 
											options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }'))),
										selectizeInput("mut_sample_line", label = "Mut_sample", choices = NULL, multiple = T, 
											options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }'))),
										actionButton("addTrackButtonMut", "Load mutation profile")
									)
								)
							)
						),
						mainPanel(width = 9,
							tags$style(HTML(".igv-center-line { display: none; pointer-events:none; position: absolute; top: 0; bottom: 0; transform: translateX(-50%); z-index: 8; user-select: none; -moz-user-select: none; -webkit-user-select: none; border-left-style: dashed; border-left-width: thin; border-right-style: dashed; border-right-width: thin; }")),
							tags$style(HTML(".igv-popover-name-value { cursor: default; text-wrap: none; white-space: nowrap; max-width: 384px; }")),
							tags$style(HTML(".igv-popover-name { font-weight: bold; padding-right: 4px; float: left; }")),
							tags$style(HTML(".igv-popover-value { padding-left: 4px; overflow: hidden; white-space: nowrap; text-overflow: ellipsis; 
											max-width: 256px; display: inline-block; }")),
							fluidRow(
								FuSViz::FuSVizOutput("FuSViz")
							),
							br()
						)
					)
	)

	#----------------------
	# Two-way module
	#----------------------
	bodyTwoway <- tabItem(tabName = "fusion", value="fusion_plot", 
					h4("Two-way module (RNA-seq)"),
					tags$head(tags$style(".two-way-select .selectize-input { font-size: 10px; width: 110px; }")), # adjust font-size of select-input as 9px
					tags$head(tags$style(".two-way-select .selectize-dropdown { font-size: 10px; }")), # adjust font-size of select-dropdown as 9px
					fluidRow( 
						box(status="warning", width=3,
							tabBox(id = "fusion_plot_tab", height = "100%", width = "100%",
								tabPanel("Overview",
									selectizeGroupUI(
										id = "overview1",
										params = list(
											gene1 = list(inputId = "gene1", title = div(style="font-size:12px", "GeneA (*)"), placeholder = 'select'),
											gene2 = list(inputId = "gene2", title = div(style="font-size:12px", "GeneB (*)"), placeholder = 'select')
										)
									),
									selectizeGroupUI(
										id = "overview2",
										params = list(
											pos1 = list(inputId = "pos1", title = div(style="font-size:12px", "Breakpoint A"), placeholder = 'select'),
											pos2 = list(inputId = "pos2", title = div(style="font-size:12px", "Breakpoint B"), placeholder = 'select')
										)
									),
									br(),
									fluidRow(
										column(6, div(class="two-way-select", selectizeInput(inputId = "transA_overview", label = div(style="font-size:12px", "GeneA transcript"), 
													choices = NULL, multiple = TRUE, options = list(placeholder = 'Ensembl transcript id')))),
										column(6, div(class="two-way-select", selectizeInput(inputId = "transB_overview", label = div(style="font-size:12px", "GeneB transcript"), 
													choices = NULL, multiple = TRUE, options = list(placeholder = 'Ensembl transcript id'))))
									),
									sliderInput(inputId = "overview_split", label = div(style="font-size:12px", "Num of split reads"), value=0, min=0, max=20, step=1),
									sliderInput(inputId = "overview_span", label = div(style="font-size:12px", "Num of span reads"), value=0, min=0, max=20, step=1),
									fluidRow(column(12, HTML("<label for='dashline'>Ruler line:&nbsp;&nbsp;&nbsp;&nbsp;</label><input id='dashline' type='checkbox'>"))),
									fluidRow(column(12, div(style="font-size:12px", HTML("<b>(*)</b> indicates required input")))),
									br(),
									actionButton("overview_on", "Plot / Refresh")
								),
								tabPanel("Sample",	
									selectizeGroupUI(
										id = "individual1", inline = FALSE,
										params = list(
											gene1 = list(inputId = "gene1", title = div(style="font-size:12px", "GeneA (*)"), placeholder = 'select'),
											gene2 = list(inputId = "gene2", title = div(style="font-size:12px", "GeneB (*)"), placeholder = 'select')
										)
									),
									selectizeGroupUI(
										id = "individual2", inline = FALSE,
										params = list(
											pos1 = list(inputId = "pos1", title = div(style="font-size:12px", "Breakpoint A(*)"), placeholder = 'select'),
											pos2 = list(inputId = "pos2", title = div(style="font-size:12px", "Breakpoint B(*)"), placeholder = 'select'),
											sample = list(inputId = "name", title = div(style="font-size:12px", "Sample (*)"), placeholder = 'select'),
											strand1 = list(inputId = "strand1", title = div(style="font-size:12px", "Strand(fusion) A"), placeholder = 'select'),
											strand2 = list(inputId = "strand2", title = div(style="font-size:12px", "Strand(fusion) B"), placeholder = 'select')
										)
									),
									selectizeInput(inputId = "transA_individual", label = div(style="font-size:12px", "GeneA transcript"), choices = NULL, multiple = TRUE, 
										options = list(placeholder = 'Ensembl transcript id')),
									selectizeInput(inputId = "transB_individual", label = div(style="font-size:12px", "GeneB transcript"), choices = NULL, multiple = TRUE, 
										options = list(placeholder = 'Ensembl transcript id')),
									fluidRow(column(12, div(style="font-size:12px", HTML("<b>(*)</b> indicates required input")))),
									br(),
									actionButton("individual_on", "Plot / Refresh")
								),
								tabPanel("Domain",
									selectizeGroupUI(
										id = "domain1",
										params = list(
											gene1 = list(inputId = "gene1", title = div(style="font-size:12px", "GeneA (*)"), placeholder = 'select'),
											gene2 = list(inputId = "gene2", title = div(style="font-size:12px", "GeneB (*)"), placeholder = 'select')
										)
									),
									hr(tags$style("border-top: 1px dotted; color: black; background-color: black;")),
									fluidRow(
										column(6, div(class="two-way-select", selectizeInput(inputId = "domainA", label = div(style="font-size:12px", "TranscriptA (*)"), choices = NULL, 
												multiple = FALSE, options = list(placeholder = "Ensembl transcript id")))),
										column(6, div(class="two-way-select", selectizeInput(inputId = "domainB", label = div(style="font-size:12px", "TranscriptB (*)"), choices = NULL, 
												multiple = FALSE, options = list(placeholder = "Ensembl transcript id"))))
									),
									selectizeGroupUI(
										id = "domain2",
										params = list(
											pos1 = list(inputId = "pos1", title = div(style="font-size:12px", "Breakpoint A"), placeholder = 'select'),
											pos2 = list(inputId = "pos2", title = div(style="font-size:12px", "Breakpoint B"), placeholder = 'select')
										)
									),
									fluidRow(column(12, div(style="font-size:12px", HTML("<b>(*)</b> indicates required input")))),
									br(),
									actionButton("domain_on", "Activate")
								)
							)
						),
						box(status="warning", width=9, height = "100%",
							tabBox(id = "fusion_plot_page", height = "100%", width = "100%",
								tabPanel("Overview_plot", 
									fluidRow(
										column(12, div(class="overviewlink", plotOutput("chimerics_up", width = "100%", height = "100%"))),
										column(12, div(class="overviewmain", id="overviewmain",
													shinycssloaders::withSpinner(plotOutput("chimerics_down", width = "100%", height = "100%"), type = 5, size = 0.7)
										)),
										tags$head(tags$style(".overviewmain {display: flex; align-items: flex-end; background-color: transparent; z-index: 0;}
											.overviewlink {top: 1px; position: relative; z-index: 1;}
											.straightLine, .hrLine {position: absolute; background-color: black;}
											html, body {height: 100%; width: 100%; margin: 0px;}"
										)),
										tags$script(src = "Visualize/coordinate.js")
									),
									box(width = NULL,  status = "warning",
										sliderInput("overview_size", label = "Zoom in/out", value = 1.5, min = 0.5, max = 4, ticks = FALSE, step = 0.05),
										fluidRow(
											column(6, radioButtons(inputId = "file_fusion1", label = "Choose file type to download:",
												inline = TRUE, choices = list("PDF", "PNG"))),
											column(6, downloadButton(outputId = "FusionDown1", label = "Download plot"))
										)
									)
								),
								tabPanel("Sample_plot", 
									fluidRow(
										column(12, shinycssloaders::withSpinner(plotOutput("chimerics2", width = "100%", height = "100%"), type = 5, size = 0.7))
									),
									box(width = NULL,  status = "warning",
										sliderInput("persample_size", label = "Zoom in/out", value = 1.5, min = 0.5, max = 4, ticks = FALSE, step = 0.05),
										fluidRow(
											column(6, radioButtons(inputId = "file_fusion2", label = "Choose file type to download:", 
												inline = TRUE, choices = list("PDF", "PNG"))),
											column(6, downloadButton(outputId = "FusionDown2", label = "Download plot"))
										)
									)
								),
								tabPanel("Domain_plot", 
									fluidRow(
										column(12, div(class="domainA", plotOutput("domain_up", height=220))),
										column(12, div(class="fuseline", shinycssloaders::withSpinner(plotOutput("linking", height=60), type = 5, size = 0.7))),
										column(12, div(class="domainB", plotOutput("domain_down", height=220))),
										tags$head(tags$style("
										.domainA{background-color: transparent; z-index: -1;}
										.fuseline{top: -10px; position: relative; z-index: 1;}
										.domainB{background-color: transparent; top: -20px; position: relative; z-index: 0;}"
										))
									),
									box(width = NULL,  status = "warning",
										fluidRow(
											column(6, radioButtons(inputId = "file_fusion3", label = "Choose file type to download:", inline = TRUE, choices = list("PDF", "PNG"))),
											column(6, downloadButton(outputId = "FusionDown3", label = "Download plot"))
										)
									)
								)
							)
						)
					)	
				)
	
	#-----------------------
	# Network module
	#-----------------------
	bodyNetwork <- tabItem(tabName = "net", value="net_plot",
						h4("Network module"),
						fluidRow(
							box(status="warning", width=3,
								tabBox(id = "network_setting_tab", height = "100%", width = "100%",
									tabPanel("RNA_SV_panel",
										fluidRow(column(12, actionButton("net_rna_plot", "Plot / Refresh"))),
										br(),
										fluidRow(column(12, tags$p(tags$h5("Parameter Setting:", 
														style="margin-left: 0em; text-align: left; font-weight: bold; text-decoration: underline")))),
										fluidRow(
											column(12, checkboxInput(inputId="navi_rna", label="Display Navigation", TRUE)),
											column(12, numericInput(inputId="grav_cons_rna", label="Gravitational Constant", value=-2000, min=-100000, max=100000)),
											column(12, numericInput(inputId="cent_grav_rna", label="Central Gravity", value=0.3, min=0, max=1)),
											column(12, numericInput(inputId="spri_cons_rna", label="Spring Constant", value=0.01, min=0, max=1)),
											column(12, numericInput(inputId="edge_leng_rna", label="Edge length", value=90, min=1, max=30)),
											column(12, selectizeInput(inputId="net_select_rna", label="Node search", choices=NULL, multiple=F, 
													options=list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
											column(12, numericInput(inputId="node_font_s_rna", label="Node font size", value=14, min=1, max=40)),
											column(12, numericInput(inputId="node_font_p_rna", label="Node font pos", value=-20, min=-100, max=100)),
											column(12, sliderInput(inputId="node_scale_rna", label="Node scaling", value=c(10, 20), min=1, max=50, step=1))
										),
										br()
									),
									tabPanel("DNA_SV_panel",
										fluidRow(column(12, actionButton("net_dna_plot", "Plot / Refresh"))),
										br(),
										fluidRow(column(12, tags$p(tags$h5("Parameter Setting:", 
														style="margin-left: 0em; text-align: left; font-weight: bold; text-decoration: underline")))),
										fluidRow(
											column(12, checkboxInput(inputId="navi_dna", label="Display Navigation", TRUE)),
											column(12, numericInput(inputId="grav_cons_dna", label="Gravitational Constant", value=-800, min=-100000, max=100000)),
											column(12, numericInput(inputId="cent_grav_dna", label="Central Gravity", value=0.1, min=0, max=1)),
											column(12, numericInput(inputId="spri_cons_dna", label="Spring Constant", value=0.02, min=0, max=1)),
											column(12, numericInput(inputId="edge_leng_dna", label="Edge length", value=10, min=1, max=30)),
											column(12, selectizeInput(inputId="net_select_dna", label="Node search", choices=NULL, multiple=F,
													options=list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
											column(12, numericInput(inputId="node_font_s_dna", label="Node font size", value=14, min=1, max=40)),
											column(12, numericInput(inputId="node_font_p_dna", label="Node font pos", value=-44, min=-100, max=100)),
											column(12, sliderInput(inputId="node_scale_dna", label="Node scaling", value=c(10, 20), min=1, max=50, step=1))
										)
									)
								)
							),
							box(status="warning", width=9,
								tabBox(id = "network_plot_page", height = "100%", width = "100%",
									tabPanel("RNA_SV_network_plot",
										shinycssloaders::withSpinner(visNetwork::visNetworkOutput("network_rna", width = "100%", height = "700px"), type = 5, size = 0.7)
									),
									tabPanel("RNA_SV_network_hub",
										fluidRow(
											# DT and dataTable library implementary
											column(12, shinycssloaders::withSpinner(DT::dataTableOutput("RNAhubs", height = 450), type = 5))
										)
									),
									tabPanel("DNA_SV_network_plot",
										shinycssloaders::withSpinner(visNetwork::visNetworkOutput("network_dna", width = "100%", height = "700px"), type = 5, size = 0.7)
									),
									tabPanel("DNA_SV_network_hub",
										fluidRow(
											# DT and dataTable library implementary
											column(12, shinycssloaders::withSpinner(DT::dataTableOutput("DNAhubs", height = 450), type = 5, size = 0.7))
										)
									)
								)
							)
						)
	)

	bodyDoc <- tabItem(tabName = "doc", value="doc_vis",	
				fluidPage(tags$iframe(src = 'Visualize/document.html', width = '100%', height = '1000px', frameborder = 0, scrolling = 'auto'))
	)

	#------------------------------------
	# Assemble all modules to main body 
	#-------------------------------------
	ui <- dashboardPage(
		title="FuSViz - an interactive Shiny app for visualization of structure variations and fusion genes/transcripts",
			skin = "black",
		dashboardHeader(title = tags$a(href='https://github.com/senzhaocode',
				tags$img(src='Visualize/logo.png', height = "65%")), dropdownMenuOutput("messageMenu"), titleWidth = 150),
		manual,
		dashboardBody( #// use_busy_spinner(spin = "fading-circle", position = "top-right", margins = c(200, 400), height = "100px", width = "100px"),
			tabItems(
				bodyIntro,
				#Table contents
				bodyTable,
				#Circle plots
				bodyCircle,
				#Linear plots
				bodyLinear,
				#Two-way plots
				bodyTwoway,
				#Network plots
				bodyNetwork,
				#Documentation
				bodyDoc
			)
		)
	)

