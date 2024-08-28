
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
			menuItem("Home", tabName = "main", icon = icon("mortar-board", verify_fa = FALSE)),
			menuItem("Table", tabName = "table", icon = icon("file-text", verify_fa = FALSE)),
			menuItem("Circle", tabName = "cir", icon = icon("venus-mars", verify_fa = FALSE)),
			menuItem("Linear", tabName = "igv", icon = icon("bar-chart", verify_fa = FALSE)),
			menuItem("Two-way", tabName = "fusion", icon = icon("line-chart", verify_fa = FALSE)),
			menuItem("Network", tabName = "net", icon = icon("th", verify_fa = FALSE)),
			menuItem("Document", tabName = "doc", icon = icon("map", verify_fa = FALSE))
		)
	)
	
	#----------------------
	# Introduction page 
	#----------------------
	bodyIntro <- tabItem(tabName = "main", value="main_panel",
					tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("disablebutton", function(message) {
          				console.log(message)
          				eval(message.code);
        			});'))),
					tags$head(tags$style(".box-header h3 {font-weight: bold; }")),
					fluidRow(
						box(
							tags$head(tags$link(rel="shortcut icon", href="Visualize/logo.ico")),
							title = "Welcome to FuSViz Shiny App!", width = 12, status = "primary",
							div(style="text-align: justify;", HTML("<b>FuSViz</b> is a web-based tool to visualize, interpret and prioritize genomic/transcriptomic structural variations (SVs) of multiple samples. It provides multiple data view approaches in an interactive interface to investigate the prevalence and recurrence of SVs and associated partner genes in a sample cohort. This tool combines SVs called from DNA-seq and RNA-seq data together to illustrate the biological impact of SVs on the implicated genes and their relevant genomic regions contextualized with various annotations."))
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
								column(4, selectizeInput("genome", label = span("Genome version", div(style = "display:inline-block;", title = "Choose a genome version for annotation data import", icon("info-circle", style="font-size: 12px"))), 
											choices = c("hg19", "hg38", "hg19_offline", "hg38_offline", "GRCm39", "GRCm39_offline"), multiple = F, 
											options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
								column(3, style = "margin-top: 25px;", actionButton("Import_genome_data", "Import", icon = icon("bolt", verify_fa = FALSE), width = "100%"))
							)
						)
					),
					fluidRow(
						box(title = "Upload structural variation file", width = 12, status = "info", solidHeader = TRUE,
							 fluidRow(
								column(3, h5(icon("upload", verify_fa = FALSE), "SV calls from RNA-seq data", style = "font-weight: bold; font-size: 16px;")),
								column(4, fileInput('file_rna_data', label = span("Upload file in txt format", div(style = "display:inline-block;", title = "File upload button will be enabled after importing annotation resource data", icon("info-circle", style="font-size: 12px"))), accept = c('text/csv', 'text/comma-separated-values', 
										'text/tab-separated-values', '.csv', '.tsv', '.txt'), placeholder = "txt, csv or tsv format")),
								column(4, radioButtons('sep_rna_file', label = HTML('<p>Separator <a href="https://fusviz.s3.eu-north-1.amazonaws.com/RNA_SV_example.txt" target="_blank">(See and download an example on basis of hg38 genome assembly)</a></p>'), choices = c(Comma = ',', Tab = '\t', Semicolon = ';'), selected = '\t', inline = T))
							),
							hr(tags$style("border-top: 1px solid #000000; color: black; background-color: white")),
							fluidRow(
								column(3, h5(icon("upload", verify_fa = FALSE), "SV calls from DNA-seq data", style = "font-weight: bold; font-size: 16px;")),
								column(4, fileInput('file_dna_data', label = span("Upload file in bedpe format", div(style = "display:inline-block;", title = "File upload button will be enabled after importing annotation resource data", icon("info-circle", style="font-size: 12px"))), accept = c('text/csv', 'text/comma-separated-values', 
										'text/tab-separated-values', '.csv', '.tsv', '.txt'), placeholder = "bedpe format")),
								column(4, radioButtons('sep_dna_file', label = HTML('<p>Separator <a href="https://fusviz.s3.eu-north-1.amazonaws.com/DNA_SV_example.txt" target="_blank">(See and download an example on basis of hg38 genome assembly)</a></p>'), choices = c(Comma = ',', Tab = '\t', Semicolon = ';'), selected = '\t', inline = T))
							)
						)
					),
					fluidRow(
						box(title = "Upload mutation profile file (optional)", width = 12, status = "info", solidHeader = TRUE,
							fluidRow(
								column(3, h5(icon("upload", verify_fa = FALSE), "Mutation variants (SNVs and Indels)", style = "font-weight: bold; font-size: 16px;")),
								column(4, fileInput('file_maf_data', label = span("Upload file in MAF format", div(style = "display:inline-block;", title = "File upload button will be enabled after importing annotation resource data", icon("info-circle", style="font-size: 12px"))), accept = c('text/maf', '.maf'), placeholder = "maf format")),
								column(4, radioButtons('sep_maf_file', label = HTML('<p>Separator <a href="https://fusviz.s3.eu-north-1.amazonaws.com/TCGA.PRAD.mutect.somatic.maf" target="_blank">(See and download an example on basis of hg38 genome assembly)</a></p>'), choices = c(Tab = '\t'), selected = '\t', inline = T))
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
					tags$script(HTML("function myFunction(element) {",
						"$('a[data-toggle=\"tab\"]').on('shown.bs.tab', function() {",
						"var target = this.href.split('#');",
						"$('.sidebar-menu a').filter('a[href=\"#' + target[1] + '\"]').tab('show');",
						"});",
						"if ( igv.browser != undefined ) {",
						"igv.browser.search(element)",
						"}",
						"};"
					)),
					tabsetPanel(id = "tabset-data",
						tabPanel("SV from RNA-seq",
							br(),
							fluidRow(
								# DT and dataTable library implementary here for RNA-seq
								column(12, shinycssloaders::withSpinner(DT::dataTableOutput("RNAcontents", height = 450), type = 5, size = 0.7)),
							),
							br(),
							fluidRow(
								box(title = "Feature summary", status="warning", width=12, height = "100%", collapsible = TRUE, collapsed = FALSE,
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
												column(4, sliderInput(inputId = "word_size_rna", label = span("Word size", div(style = "display:inline-block;", title = "Control the size of symbol and text relative to the default size", icon("info-circle", style="font-size: 12px"))), value=0.45, min=0,001, max=1, step=0.05)),
												column(4, sliderInput(inputId = "gene_freq_rna", label = span("Gene freq", div(style = "display:inline-block;", title = "Gene with recurrence above the cutoff value is visible in the plot", icon("info-circle", style="font-size: 12px"))), value=1, min=1, max=10, step=1)),
												column(4, selectInput(inputId = "word_shape_rna", label = span("Word shape", div(style = "display:inline-block;", title = "Control the shape of wordcloud", icon("info-circle", style="font-size: 12px"))), choices=c("circle","square","cardioid"), selected="square", multiple=F))
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
												column(3, sliderInput(inputId = "rna_hist_size", label = span("Font size", div(style = "display:inline-block;", title = "Control the font size of sample name and text relative to the default size", icon("info-circle", style="font-size: 12px"))), value=0.6, min=0,001, max=1, step=0.05)),
												column(3, sliderInput(inputId = "rna_hist_pos", label = span("X-axis pos", div(style = "display:inline-block;", title = "Set the position of sample name horizontally in the histogram plot", icon("info-circle", style="font-size: 12px"))), value=-20, min=-100, max=0, step=2)),
												column(3, sliderInput(inputId = "rna_hist_angle", label = span("X-axis rotate", div(style = "display:inline-block;", title = "Rotate the sample name in the histogram plot [0 - horizontal; 90 - vertical]", icon("info-circle", style="font-size: 12px"))), value=90, min=0, max=90, step=1)),
												column(3, sliderInput(inputId = "rna_hist_space", label = span("Bar space", div(style = "display:inline-block;", title = "Control the interval space between sample names [as a fraction of the average bar width] in the histogram plot", icon("info-circle", style="font-size: 12px"))), value=0.1, min=0, max=1, step=0.05))
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
								box(title = "Feature summary", status="warning", width=12, height = "100%", collapsible = TRUE, collapsed = FALSE,
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
												column(4, sliderInput(inputId = "word_size_dna", label = span("Word size", div(style = "display:inline-block;", title = "Control the size of symbol and text relative to the default size", icon("info-circle", style="font-size: 12px"))), value=0.2, min=0.001, max=1, step=0.05)),
												column(4, sliderInput(inputId = "gene_freq_dna", label = span("Gene freq", div(style = "display:inline-block;", title = "Gene with recurrence above the cutoff value is visible in the plot", icon("info-circle", style="font-size: 12px"))), value=1, min=1, max=10, step=1)),
												column(4, selectInput(inputId = "word_shape_dna", label = span("Word shape", div(style = "display:inline-block;", title = "Control the shape of wordcloud", icon("info-circle", style="font-size: 12px"))), choices=c("circle","square","cardioid"), selected="square", multiple=F))
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
												column(3, sliderInput(inputId = "dna_hist_size", label = span("Font size", div(style = "display:inline-block;", title = "Control the font size of sample name and text relative to the default size", icon("info-circle", style="font-size: 12px"))), value=0.6, min=0,001, max=1, step=0.05)),
												column(3, sliderInput(inputId = "dna_hist_pos", label = span("X-axis pos", div(style = "display:inline-block;", title = "Set the position of sample name horizontally in the histogram plot", icon("info-circle", style="font-size: 12px"))), value=-100, min=-500, max=0, step=5)),
												column(3, sliderInput(inputId = "dna_hist_angle", label = span("X-axis rotate", div(style = "display:inline-block;", title = "Rotate the sample name in the histogram plot [0 - horizontal; 90 - vertical]", icon("info-circle", style="font-size: 12px"))), value=90, min=0, max=90, step=1)),
												column(3, sliderInput(inputId = "dna_hist_space", label = span("Bar space", div(style = "display:inline-block;", title = "Control the interval space between sample names [as a fraction of the average bar width] in the histogram plot", icon("info-circle", style="font-size: 12px"))), value=0.1, min=0, max=1, step=0.05))
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
								box(title = "Feature summary", status="warning", width=12, height = "100%", collapsible = TRUE, collapsed = FALSE,
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
												column(4, sliderInput(inputId = "word_size_mut", label = span("Word size", div(style = "display:inline-block;", title = "Control the size of symbol and text relative to the default size", icon("info-circle",style="font-size: 12px"))), value=0.2, min=0.001, max=1, step=0.05)),
												column(4, sliderInput(inputId = "gene_freq_mut", label = span("Gene freq", div(style = "display:inline-block;", title = "Gene with recurrence above the cutoff value is visible in the plot", icon("info-circle", style="font-size: 12px"))), value=1, min=1, max=10, step=1)),
												column(4, selectInput(inputId = "word_shape_mut", label = span("Word shape", div(style = "display:inline-block;", title = "Control the shape of wordcloud", icon("info-circle", style="font-size: 12px"))), choices=c("circle","square","cardioid"), selected="square", multiple=F))
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
									select_group_ui(id = "circle_set", inline = FALSE,
										params = list(
											list(inputId = "chr", label = div("Chrom", style = "display:inline-block; font-size:13px", title = "Limit to plot SVs related to the selected chromosomes", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select'),
											list(inputId = "gene", label = div("Gene", style = "display:inline-block; font-size:13px", title = "Limit to plot SVs related to the selected genes", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select'),
											list(inputId = "name", label = div("Sample", style = "display:inline-block; font-size:13px", title = "Limit to plot SVs related to the selected samples", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select')
										), vs_args = list(search = TRUE)
									),
									sliderInput(inputId = "circle_split", label = div("Num of split reads", style = "display:inline-block; font-size:13px", title = "Filter out SVs with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_span", label = div("Num of span reads", style = "display:inline-block; font-size:13px", title = "Filter out SVs with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_num", label = div("Num of samples", style = "display:inline-block; font-size:13px", title = "Filter out recurrent SVs with the support number of samples no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=50, step=1),
									checkboxInput("cirmutcheck", "Load mutation data", FALSE),
									conditionalPanel(condition = "input.cirmutcheck == 1", 
										selectInput("mut_type", label = div("Mutation type", style = "display:inline-block; font-size:13px", title = "Limit to plot small variant profile of the choosen mutation types", icon("info-circle", style="font-size: 10px")), choices = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", 
											"Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", 
											"De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "3'UTR", "5'UTR", "3'Flank", "Targeted_Region", 
											"Silent", "Intron", "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA"), 
											selected = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", 
											"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", 
											"Start_Codon_SNP", "Stop_Codon_Del"), multiple = T)),
									conditionalPanel(condition = "input.cirmutcheck == 1",
										selectizeInput("mut_sample", label = div("Mutation sample", style = "display:inline-block; font-size:13px", title = "Limit to plot small variant profile related to the choosen samples", icon("info-circle", style="font-size: 10px")), choices = NULL, multiple = T, options = list(placeholder = 'select', 
											onInitialize = I('function() { this.setValue(""); }')))),
									br(),
									fluidRow(column(12, actionButton("circos_plot", label = div("Plot / Refresh", title = "Please note that changes made to the data will not be reflected in the plot until the ‘plot/refresh’ button is clicked on")))),
									br(),
									fluidRow(column(12, downloadButton(outputId = "CircleDownRNA", label = "Download circular plot")))
								),
								tabPanel("DNA_SV_panel",
									select_group_ui(id = "circle_set_dna", inline = FALSE,
										params = list(
											list(inputId = "chr", label = div("Chrom", style = "display:inline-block; font-size:13px", title = "Limit to plot SVs related to the selected chromosomes", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select'),
											list(inputId = "gene", label = div("Gene", style = "display:inline-block; font-size:13px", title = "Limit to plot SVs related to the selected genes", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select'),
											list(inputId = "name", label = div("Sample", style = "display:inline-block; font-size:13px", title = "Limit to plot SVs related to the selected samples", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select')
										), vs_args = list(search = TRUE)
									),
									sliderInput(inputId = "circle_split_dna", label = div("Num of split reads", style = "display:inline-block; font-size:13px", title = "Filter out SVs with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_span_dna", label = div("Num of span reads", style = "display:inline-block; font-size:13px", title = "Filter out SVs with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=20, step=1),
									sliderInput(inputId = "circle_num_dna", label = div("Num of samples", style = "display:inline-block; font-size:13px", title = "Filter out recurrent SVs with the support number of samples no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=50, step=1),
									sliderInput(inputId = "circle_dis_dna", label = div("Dist intra-chrom", style = "display:inline-block; font-size:13px", title = "Filter out intra-chromosome SVs with a distance between breakpoints below the cutoff value", icon("info-circle", style="font-size: 10px")), value=1000000, min=0, max=10000000, step=100000),
									checkboxInput("cirmutcheck_dna", "Load mutation data", FALSE),
									conditionalPanel(condition = "input.cirmutcheck_dna == 1", 
										selectInput("mut_type_dna", label = div("Mutation type", style = "display:inline-block; font-size:13px", title = "Limit to plot small variant profile of the choosen mutation types", icon("info-circle", style="font-size: 10px")), choices = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", 
											"Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", 
											"De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "3'UTR", "5'UTR", "3'Flank", "Targeted_Region", 
											"Silent", "Intron", "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA"),
											selected = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", 
											"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", 
											"Start_Codon_SNP", "Stop_Codon_Del"), multiple = T)),
									conditionalPanel(condition = "input.cirmutcheck_dna == 1",
										selectizeInput("mut_sample_dna", label = div("Mutation sample", style = "display:inline-block; font-size:13px", title = "Limit to plot small variant profile related to the choosen samples", icon("info-circle", style="font-size: 10px")), choices = NULL, multiple = T, options = list(placeholder = 'select', 
											onInitialize = I('function() { this.setValue(""); }')))),
									br(),
									fluidRow(column(12, actionButton("circos_plot_dna", label = div("Plot / Refresh", title = "Please note that changes made to the data will not be reflected in the plot until the ‘plot/refresh’ button is clicked on")))),
									br(),
									fluidRow(column(12, downloadButton(outputId = "CircleDownDNA", label = "Download circular plot")))
								)
							)
						),
						box(status="warning", width=9,
							tabBox(id = "circle_plot_page", height = "100%", width = "100%",
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
							fluidRow(
								column(12, tags$p(tags$b("Choose local files to upload: ", style = "font-size:13px"))),
								column(12, tags$p(tags$b("* Upload read alignment and its index files (", tags$u("BAM"), " or ", tags$u("CRAM"), " format)", style = "font-size:12px"))),
								column(12, tags$p(tags$b("* Upload annotation and its index files (bgzipped and tabixed", tags$u("VCF"), ", ", tags$u("BED"), " and ", tags$u("GTF"), " format)", style = "font-size:12px"))),
								column(12, tags$p(tags$b("* Upload genome reference and its index files (indexed ", tags$u("fasta"), ") and cytoband information in ", tags$u("txt"), " format - used to enable offline mode", style = "font-size:12px"))),
								br(),
								column(12, HTML("<input id='uploadfile' class='hidden' type='file' multiple='true' accept='.bam,.bai,.cram,.crai,.gz,.tbi,.fasta,.fai,.txt' onchange=\"load()\"/>
											<label for='uploadfile'><img src='Visualize/Button.png' alt='Upload img' height='80' width='160'></label>")),
								br(),
								column(12, HTML("<div id='filename'></div>")),
								tags$script(src = "Visualize/Upload_file.js")
							),
							hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
							fluidRow(
								column(12, tags$p(tags$b("Activate the gene annotation track for offline mode (click it after genome reference upload)", style = "font-size:13px"))),
								column(12, actionButton("offlineTrackButton", "Load gene track (offline)"))
							),
							hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
							fluidRow(
								column(12, textInput(inputId = "BAM", label = div("Alignment URL", style = "display:inline-block; font-size:12px", title = "Load a hosted BAM/CRAM file via URL", icon("info-circle", style="font-size: 9px")), value = NULL)),
								column(12, textInput(inputId = "BAMindex", label = div("Alignment index URL", style = "display:inline-block; font-size:12px", title = "Load a hosted BAM/CRAM index file via URL", icon("info-circle", style="font-size: 9px")), value = NULL)),
								column(12, actionButton("addTrackButtonBAM", label = NULL, style = "width: 89px; height: 55px;
											background: url('Visualize/cloud_upload.png');  background-size: cover; background-position: center;"))
							),
							hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
							tabsetPanel(id = "igv-panel-control",
								tabPanel("DNA_SV",
									fluidRow(
										column(6, numericInput("dna_split_bedpe", label = div("Num_split", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(6, numericInput("dna_span_bedpe", label = div("Num_span", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(6, numericInput("dna_dismin_bedpe", label = div("Min_dist", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with a distance between breakpoints below the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=500000000)),
										column(6, numericInput("dna_dismax_bedpe", label = div("Max_dist", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with a distance between breakpoints above the cutoff value", icon("info-circle", style="font-size: 9px")), value=500000000, min=0, max=500000000)),
										column(5, selectizeInput("dna_type_bedpe", label = div("SV_type", style = "display:inline-block; font-size:12px", title = "Limit to plot SVs of a specific category [e.g., INV-inversion, DEL-deletion, DUP-duplication, INS-insertion]", icon("info-circle", style="font-size: 9px")), choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(7, selectizeInput("dna_sample_bedpe", label = div("Sample", style = "display:inline-block; font-size:12px", title = "Limit to plot SVs of the choosen samples", icon("info-circle", style="font-size: 9px")), choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(12, actionButton("addTrackButtonBedPe", label = div("Load and refresh DNA SV track in bedpe", title = "Please note that changes made to the data will not be reflected in the plot until the ‘Load and refresh DNA SV track in bedpe’ button is clicked on")))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(6, numericInput("dna_split_seg", label = div("Num_split", style = "display:inline-block; font-size:12px", title = "Filter out copy number segments with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(6, numericInput("dna_span_seg", label = div("Num_span", style = "display:inline-block; font-size:12px", title = "Filter out copy number segments with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(4, div(class="no-spinners", textInput("dna_seg_chr", label = div("Chrom", style = "display:inline-block; font-size:12px", title = "Limit to plot copy number segments that overlap with a genomic region of the choosen chromosome [it is utilized together with textInput boxes Start and End]", icon("info-circle", style="font-size: 9px")), 
											value = "", placeholder = "chr1..chrX,chrY,chrM"))),
										tags$style(HTML("#dna_seg_chr { width: 55px; font-size:12px; }")),
										column(4, div(class="no-spinners", numericInput("dna_seg_start", label = div("Start", style = "display:inline-block; font-size:12px", title = "Limit to plot copy number segments that overlap with a genomic region [defined by textInput boxes start and end]", icon("info-circle", style="font-size: 9px")), 
											value=0, min=0, max=500000000))),
										tags$style(HTML("#dna_seg_start { width: 80px; font-size:12px; }")),
										column(4, div(class="no-spinners", numericInput("dna_seg_end", label = div("End", style = "display:inline-block; font-size:12px", title = "Limit to plot copy number segments that overlap with a genomic region [defined by textInput boxes start and end]", icon("info-circle", style="font-size: 9px")), 
											value=0, min=0, max=500000000))),
										tags$style(HTML("#dna_seg_end { width: 80px; font-size:12px; }")),
										tags$style(HTML(".no-spinners { display: inline-block; vertical-align:top; } 
														.no-spinners input[type=number] { -moz -appearance:textfield; } 
														.no-spinners input[type=number]::{ -moz-appearance:textfield; }
														.no-spinners input[type=number]::-webkit-outer-spin-button, 
														.no-spinners input[type=number]::-webkit-inner-spin-button { -webkit-appearance: none; margin: 0; } ")),
										column(8, selectizeInput("dna_sample_seg", label = div("Sample", style = "display:inline-block; font-size:12px", title = "Limit to plot SVs of the choosen samples", icon("info-circle", style="font-size: 9px")), choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(12, actionButton("addTrackButtonSeg", label = div("Load and refresh DNA SV track in seg", title = "Please note that changes made to the data will not be reflected in the plot until the ‘Load and refresh DNA SV track in seg’ button is clicked on")))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(4, numericInput("dna_split_bed", label = div("Num_split", style = "display:inline-block; font-size:12px", title = "Filter out SVs with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(4, numericInput("dna_span_bed", label = div("Num_span", style = "display:inline-block; font-size:12px", title = "Filter out SVs with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(4, selectizeInput("dna_type_bed", label = div("SV_type", style = "display:inline-block; font-size:12px", title = "Limit to plot SVs of a specific category [e.g., INV-inversion, DEL-deletion, DUP-duplication, INS-insertion]", icon("info-circle", style="font-size: 9px")), choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										tags$style(HTML("#dna_type_bed+ div>.selectize-input { width: 70px; } 
														#dna_type_bed+ div>.selectize-input input { width: 70px; } 
														#dna_type_bed+ .selectize-dropdown-content { width: 70px; }")),
										column(6, numericInput("dna_dismin_bed", label = div("Min_dist", style = "display:inline-block; font-size:12px", title = "Filter out SVs with a distance between breakpoints below the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=500000000)),
										column(6, numericInput("dna_dismax_bed", label = div("Max_dist", style = "display:inline-block; font-size:12px", title = "Filter out SVs with a distance between breakpoints above the cutoff value", icon("info-circle", style="font-size: 9px")), value=500000000, min=0, max=500000000)),
										column(12, actionButton("addTrackButtonBed", label = div("Load and refresh DNA breakpoints in bed", title = "Please note that changes made to the data will not be reflected in the plot until the ‘Load and refresh DNA breakpoints in bed’ button is clicked on")))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(12, textInput(inputId = "genomicPos", label = div("Genome coordinate", style = "display:inline-block; font-size:12px", title = "Show genomic start and end coordinates of current IGV session window", icon("info-circle", style="font-size: 9px")), value = "")),
										column(6, actionButton("GetPos", label = "Show coordinate")),
										column(6, actionButton("MovePos", label = "Clear"))
									)
								),
								tabPanel("RNA_SV",
									fluidRow(
										column(6, numericInput("rna_split_bedpe", label = div("Num_split", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(6, numericInput("rna_span_bedpe", label = div("Num_span", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(6, numericInput("rna_dismin_bedpe", label = div("Min_dist", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with a distance between breakpoints below the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=500000000)),
										column(6, numericInput("rna_dismax_bedpe", label = div("Max_dist", style = "display:inline-block; font-size:12px", title = "Filter out intra-chromosome SVs with a distance between breakpoints above the cutoff value", icon("info-circle", style="font-size: 9px")), value=500000000, min=0, max=500000000)),
										column(8, selectizeInput("rna_sample_bedpe", label = div("Sample", style = "display:inline-block; font-size:12px", title = "Limit to plot SVs of the choosen samples", icon("info-circle", style="font-size: 9px")), choices = NULL, multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(12, actionButton("addTrackButtonRNABedpe", label = div("Load and refresh RNA SV track in bedpe", title = "Please note that changes made to the data will not be reflected in the plot until the ‘Load and refresh RNA SV track in bedpe’ button is clicked on")))
									),
									hr(tags$style("border-top: 3px solid #000000; color: black; background-color: black")),
									fluidRow(
										column(4, numericInput("rna_split_bed", label = div("Num_split", style = "display:inline-block; font-size:12px", title = "Filter out SVs with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(4, numericInput("rna_span_bed", label = div("Num_span", style = "display:inline-block; font-size:12px", title = "Filter out SVs with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=50)),
										column(4, selectizeInput("rna_type_bed", label = div("SV_type", style = "display:inline-block; font-size:12px", title = "Limit to plot SVs of a specific category [e.g., inter - interchromosome event; intra - intrachromosome event]", icon("info-circle", style="font-size: 9px")), choices = c("Inter", "Intra"), multiple = T, 
													options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
										column(6, numericInput("rna_dismin_bed", label = div("Min_dist", style = "display:inline-block; font-size:12px", title = "Filter out SVs with a distance between breakpoints below the cutoff value", icon("info-circle", style="font-size: 9px")), value=0, min=0, max=500000000)),
										column(6, numericInput("rna_dismax_bed", label = div("Max_dist", style = "display:inline-block; font-size:12px", title = "Filter out SVs with a distance between breakpoints above the cutoff value", icon("info-circle", style="font-size: 9px")), value=500000000, min=0, max=500000000)),
										column(12, actionButton("addTrackButtonRNABed", label = div("Load and refresh RNA breakpoints in bed", title = "Please note that changes made to the data will not be reflected in the plot until the ‘Load and refresh RNA breakpoints in bed’ button is clicked on")))
									)
								),
								tabPanel("Mut",
									fluidRow(
										selectizeInput("mut_type_line", label = div("Mut_type", style = "display:inline-block; font-size:12px", title = "Limit to plot small variant profile of the choosen mutation types", icon("info-circle", style="font-size: 9px")), choices = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
											"Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", 
											"De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "3'UTR", "5'UTR", 
											"3'Flank", "Targeted_Region", "Silent", "Intron", "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA"), multiple = T, 
											options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }'))),
										selectizeInput("mut_sample_line", label = div("Mut_sample", style = "display:inline-block; font-size:12px", title = "Limit to plot small variant profile of the choosen samples", icon("info-circle", style="font-size: 9px")), choices = NULL, multiple = T, 
											options = list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }'))),
										actionButton("addTrackButtonMut", label = div("Load and refresh mutation profile", title = "Please note that changes made to the data will not be reflected in the plot until the ‘Load and refresh mutation profile’ button is clicked on"))
									)
								)
							)
						),
						mainPanel(width = 9,
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
					tags$head(tags$style(".two-way-select .selectize-input { font-size: 10px; width: 120px; }")), # adjust font-size of select-input as 9px
					tags$head(tags$style(".two-way-select .selectize-dropdown { font-size: 10px; }")), # adjust font-size of select-dropdown as 9px
					fluidRow( 
						box(status="warning", width=3,
							tabBox(id = "fusion_plot_tab", height = "100%", width = "100%",
								tabPanel("Overview",
									select_group_ui(
										id = "overview1",
										params = list(
											list(inputId = "gene1", label = div("GeneA (*)", style = "display:inline-block; font-size:12px", title = "Choose the gene name of upstream partner [gene symbol or Ensembl entry is acceptable]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "gene2", label = div("GeneB (*)", style = "display:inline-block; font-size:12px", title = "Choose the gene name of downstream partner [gene symbol or Ensembl entry is acceptable]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL)
										), vs_args = list(search = TRUE, maxValues=1, autoSelectFirstOption = F)
									),
									hr(tags$style("border-top: 1px dotted; color: black; background-color: black;")),
									select_group_ui(
										id = "overview2",
										params = list(
											list(inputId = "breakpoint_A", label = div("Breakpoint A", style = "display:inline-block; font-size:12px", title = "Choose breakpoint coordinates of upstream partner gene [multiple selection is allowed]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "breakpoint_B", label = div("Breakpoint B", style = "display:inline-block; font-size:12px", title = "Choose breakpoint coordinates of downstream partner gene [multiple selection is allowed]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL)
										), vs_args = list(search = TRUE, autoSelectFirstOption = F)
									),
									br(),
									fluidRow(
										column(6, div(class="two-way-select", selectizeInput(inputId = "transA_overview", label = div("GeneA transcript", style = "display:inline-block; font-size:12px", title = "Choose transcript isoforms of upstream partner gene and the canonical one is highlighted with underline [multiple selection is allowed]", icon("info-circle", style="font-size: 10px")),
											choices = NULL, multiple = TRUE, options = list(placeholder = 'Ensembl transcript id')))),
										column(6, div(class="two-way-select", selectizeInput(inputId = "transB_overview", label = div("GeneB transcript", style = "display:inline-block; font-size:12px", title = "Choose transcript isoforms of downstream partner gene and the canonical one is highlighted with underline [multiple selection is allowed]", icon("info-circle", style="font-size: 10px")),
											choices = NULL, multiple = TRUE, options = list(placeholder = 'Ensembl transcript id'))))
									),
									sliderInput(inputId = "overview_split", label = div("Num of split reads", style = "display:inline-block; font-size:12px", title = "Filter out fusion events with the support number of split reads no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=20, step=1),
									sliderInput(inputId = "overview_span", label = div("Num of span reads", style = "display:inline-block; font-size:12px", title = "Filter out fusion events with the support number of spanning reads no more than the cutoff value", icon("info-circle", style="font-size: 10px")), value=0, min=0, max=20, step=1),
									fluidRow(column(12, HTML("<label for='dashline' style='font-size:12px'>Ruler line&nbsp;&nbsp;<span title='Activate ruler line function after onclick'><i class='fas fa-info-circle' style='font-size:10px'></i></span>&nbsp;&nbsp;:&nbsp;&nbsp;&nbsp;&nbsp;</label><input id='dashline' type='checkbox'>"))),
									fluidRow(column(12, div(style="font-size:12px", HTML("<b>(*)</b> indicates required input")))),
									br(),
									actionButton("overview_on", label = div("Plot / Refresh", title = "Please note that changes made to the data will not be reflected in the plot until the ‘plot/refresh’ button is clicked on"))
								),
								tabPanel("Sample",	
									select_group_ui(
										id = "individual", inline = FALSE,
										params = list(
                            				list(inputId = "gene1", label = div("GeneA (*)", style = "display:inline-block; font-size:12px", title = "Choose the gene name of upstream partner [gene symbol or Ensembl entry is acceptable]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "gene2", label = div("GeneB (*)", style = "display:inline-block; font-size:12px", title = "Choose the gene name of downstream partner [gene symbol or Ensembl entry is acceptable]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "pos1", label = div("Breakpoint A(*)", style = "display:inline-block; font-size:12px", title = "Choose a breakpoint coordinate of upstream partner gene", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "pos2", label = div("Breakpoint B(*)", style = "display:inline-block; font-size:12px", title = "Choose a breakpoint coordinate of downstream partner gene", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "name", label = div("Sample (*)", style = "display:inline-block; font-size:12px", title = "Choose a sample name", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "strand1", label = div("Fusion A Strand (*)", style = "display:inline-block; font-size:12px", title = "Choose the transcription direction of upstream fusion sequence", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
											list(inputId = "strand2", label = div("Fusion B Strand (*)", style = "display:inline-block; font-size:12px", title = "Choose the transcription direction of downstream fusion sequence", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL)
										), vs_args = list(search = TRUE, maxValues=1, autoSelectFirstOption = F)
									),
									selectizeInput(inputId = "transA_individual", label = div("GeneA transcript", style = "display:inline-block; font-size:12px", title = "Choose transcript isoforms of upstream partner gene and the canonical one is highlighted with underline [multiple selection is allowed]", icon("info-circle", style="font-size: 10px")), choices = NULL, multiple = TRUE, 
						    			options = list(placeholder = 'Ensembl transcript id')),
									selectizeInput(inputId = "transB_individual", label = div("GeneB transcript", style = "display:inline-block; font-size:12px", title = "Choose transcript isoforms of downstream partner gene and the canonical one is highlighted with underline [multiple selection is allowed]", icon("info-circle", style="font-size: 10px")), choices = NULL, multiple = TRUE, 
										options = list(placeholder = 'Ensembl transcript id')),
									fluidRow(column(12, div(style="font-size:12px", HTML("<b>(*)</b> indicates required input")))),
									br(),
									actionButton("individual_on", label = div("Plot / Refresh", title = "Please note that changes made to the data will not be reflected in the plot until the ‘plot/refresh’ button is clicked on"))
								),
								tabPanel("Domain",
									select_group_ui(
			            				id = "domain1",
			            				params = list(
				            				list(inputId = "gene1", label = div("GeneA (*)", style = "display:inline-block; font-size:12px", title = "Choose the gene name of upstream partner [gene symbol or Ensembl entry is acceptable]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
				            				list(inputId = "gene2", label = div("GeneB (*)", style = "display:inline-block; font-size:12px", title = "Choose the gene name of downstream partner [gene symbol or Ensembl entry is acceptable]", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL)
			            				), vs_args = list(search = TRUE, maxValues=1, autoSelectFirstOption = F)
			        				),
			        				hr(tags$style("border-top: 1px dotted; color: black; background-color: black;")),
			        				fluidRow(
				        				column(6, div(class="two-way-select", selectizeInput(inputId = "domainA", label = div("TranscriptA (*)", style = "display:inline-block; font-size:12px", title = "Choose transcript isoforms of upstream partner gene [the canonical transcript and ones with domain/motif annotation are highlighted with underline and bold font, respectively]", icon("info-circle", style="font-size: 10px")), choices = NULL, 
						        			multiple = TRUE, options = list(placeholder = "Ensembl transcript id", maxItems = 1)))),
				        				column(6, div(class="two-way-select", selectizeInput(inputId = "domainB", label = div("TranscriptB (*)", style = "display:inline-block; font-size:12px", title = "Choose transcript isoforms of downstream partner gene [the canonical transcript and ones with domain/motif annotation are highlighted with underline and bold font, respectively]", icon("info-circle", style="font-size: 10px")), choices = NULL, 
						        			multiple = TRUE, options = list(placeholder = "Ensembl transcript id", maxItems = 1))))
			        				),
                    				select_group_ui(
				        				id = "domain2",
				        				params = list(
					        				list(inputId = "pos1", label = div("Breakpoint A", style = "display:inline-block; font-size:12px", title = "Choose a breakpoint coordinate of upstream partner gene", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL),
					        				list(inputId = "pos2", label = div("Breakpoint B", style = "display:inline-block; font-size:12px", title = "Choose a breakpoint coordinate of downstream partner gene", icon("info-circle", style="font-size: 10px")), multiple=TRUE, placeholder = 'select', choices = NULL)
				        				), vs_args = list(search = TRUE, autoSelectFirstOption = F)
			        				),
			        				br(),
			        				fluidRow(
				        				column(6, numericInput("offset_base_A", label = div("Offset base for breakpointA", style = "display:inline-block; font-size:12px", title = "Offset value is the adjusted number of bases for provisional breakpoint coordinate of upstream parnter gene to match the genomic coordinate of nearby exon boundary in the codon frame calculation [default value: 5 - do offset adjustment if a distance between the coordinates of provisional breakpoint and nearby exon boundary <= 5 bases]", icon("info-circle", style="font-size: 10px")), value=5, min=0, max=10)),
				        				column(6, numericInput("offset_base_B", label = div("Offset base for breakpointB", style = "display:inline-block; font-size:12px", title = "Offset value is the adjusted number of bases for provisional breakpoint coordinate of downstream parnter gene to match the genomic coordinate of nearby exon boundary in the codon frame calculation [default value: 5 - do offset adjustment if a distance between the coordinates of provisional breakpoint and nearby exon boundary <= 5 bases]", icon("info-circle", style="font-size: 10px")), value=5, min=0, max=10))
			        				),
			        				fluidRow(column(12, div(style="font-size:12px", HTML("<b>(*)</b> indicates required input for plotting")))),
			        				br()
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
											.overviewlink {top: 5px; position: relative; z-index: 1;}
											.straightLine, .hrLine {position: absolute; background-color: black;}
											html, body {height: 100%; width: 100%; margin: 0px;}"
										)),
										tags$script(src = "Visualize/coordinate.js")
									),
									box(width = NULL,  status = "warning",
										sliderInput("overview_size", label = "Zoom in/out", value = 1.5, min = 0.5, max = 4, ticks = FALSE, step = 0.05),
										fluidRow(
											column(3, radioButtons(inputId = "file_fusion1", label = "Choose file type to download:", inline = TRUE, choices = list("PDF", "PNG"))),
											column(2, numericInput("overview_width", label = div("Layout_width", style = "display:inline-block; font-size:12px", title = "Set the width of plot layout (unit: inches)", icon("info-circle", style="font-size: 9px")), value=10, min=1, max=50)),
											column(2, numericInput("overview_height", label = div("Layout_height", style = "display:inline-block; font-size:12px", title = "Set the height of plot layout (unit: inches)", icon("info-circle", style="font-size: 9px")), value=10, min=1, max=50)),
											column(3, downloadButton(outputId = "FusionDown1", label = "Download plot"))
										)
									)
               					),
								tabPanel("Sample_plot", 
									fluidRow(
										column(12, shinycssloaders::withSpinner(plotOutput("chimerics2", width = "100%", height = "100%"), type = 5, size = 0.7))
									),
									box(width = NULL, status = "warning",
										sliderInput("persample_size", label = "Zoom in/out", value = 1.5, min = 0.5, max = 4, ticks = FALSE, step = 0.05),
										fluidRow(
											column(3, radioButtons(inputId = "file_fusion2", label = "Choose file type to download:", inline = TRUE, choices = list("PDF", "PNG"))),
											column(2, numericInput("sample_width", label = div("Layout_width", style = "display:inline-block; font-size:12px", title = "Set the width of plot layout (unit: inches)", icon("info-circle", style="font-size: 9px")), value=10, min=1, max=50)),
											column(2, numericInput("sample_height", label = div("Layout_height", style = "display:inline-block; font-size:12px", title = "Set the height of plot layout (unit: inches)", icon("info-circle", style="font-size: 9px")), value=10, min=1, max=50)),
											column(3, downloadButton(outputId = "FusionDown2", label = "Download plot"))
										)
									)
								),
								tabPanel("Domain_plot", 
									fluidRow(
										column(12, div(class="domainA", shinycssloaders::withSpinner(plotOutput("domain_up", height=220), type = 5, size = 0.7))),
										column(12, div(class="fuseline", plotOutput("linking", height=30))),
										column(12, div(class="domainB", shinycssloaders::withSpinner(plotOutput("domain_down", height=220), type = 5, size = 0.7))),
										tags$head(tags$style("
											.domainA{background-color: transparent; z-index: -1;}
											.fuseline{top: 5px; position: relative; z-index: 1;}
											.domainB{background-color: transparent; z-index: 0;}"
										))
									),
                   					box(width = NULL,  status = "warning",
										fluidRow(
											column(3, radioButtons(inputId = "file_fusion3", label = "Choose file type to download:", inline = TRUE, choices = list("PDF", "PNG"))),
											column(2, numericInput("domain_width", label = div("Layout_width", style = "display:inline-block; font-size:12px", title = "Set the width of plot layout (unit: inches)", icon("info-circle", style="font-size: 9px")), value=20, min=1, max=50)),
											column(2, numericInput("domain_height", label = div("Layout_height", style = "display:inline-block; font-size:12px", title = "Set the height of plot layout (unit: inches)", icon("info-circle", style="font-size: 9px")), value=10, min=1, max=50)),
											column(3, downloadButton(outputId = "FusionDown3", label = "Download plot"))
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
										fluidRow(column(12, actionButton("net_rna_plot", label = div("Plot / Refresh", title = "Please note that changes made to the data will not be reflected in the plot until the ‘plot/refresh’ button is clicked on")))),
										br(),
										fluidRow(column(12, tags$p(tags$h5("Parameter Setting:", 
														style="margin-left: 0em; text-align: left; font-weight: bold; text-decoration: underline")))),
										fluidRow(
											column(12, checkboxInput(inputId="navi_rna", label = div("Display Navigation", style = "display:inline-block; font-size:13px", title = "Hide the navigation control panel if unchecked", icon("info-circle", style="font-size: 11px")), TRUE)),
											column(12, numericInput(inputId="grav_cons_rna", label = div("Gravitational Constant", style = "display:inline-block; font-size:13px", title = "Adjust the distance of a sub-graph where a decrease of this value will increase the repulsion of the sub-graph", icon("info-circle", style="font-size: 11px")), value=-2000, min=-100000, max=100000)),
											column(12, numericInput(inputId="cent_grav_rna", label = div("Central Gravity", style = "display:inline-block; font-size:13px", title = "Set a central gravity attraction to pull the entire network back to the center [a larger value shows a stronger central force]", icon("info-circle", style="font-size: 11px")), value=0.3, min=0, max=1)),
											column(12, numericInput(inputId="spri_cons_rna", label = div("Spring Constant", style = "display:inline-block; font-size:13px", title = "Control the size of the subgraph [a smaller value shows a shorter edge]", icon("info-circle", style="font-size: 11px")), value=0.01, min=0, max=1)),
											column(12, numericInput(inputId="edge_leng_rna", label = div("Edge length", style = "display:inline-block; font-size:13px", title = "Set the length of edge in the graph", icon("info-circle", style="font-size: 11px")), value=90, min=1, max=30)),
											column(12, selectizeInput(inputId="net_select_rna", label = div("Node search", style = "display:inline-block; font-size:13px", title = "Choose and centralize a node [i.e., gene name]", icon("info-circle", style="font-size: 11px")), choices=NULL, multiple=F, 
													options=list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
											column(12, numericInput(inputId="node_font_s_rna", label = div("Node font size", style = "display:inline-block; font-size:13px", title = "Set the font size of node name", icon("info-circle", style="font-size: 11px")), value=14, min=1, max=40)),
											column(12, numericInput(inputId="node_font_p_rna", label = div("Node font pos", style = "display:inline-block; font-size:13px", title = "Adjust the position of node name", icon("info-circle", style="font-size: 11px")), value=-20, min=-100, max=100)),
											column(12, sliderInput(inputId="node_scale_rna", label = div("Node scaling", style = "display:inline-block; font-size:13px", title = "Scale the size of node in proportion to the support number of samples [min: the minimum font-size used for labels when scaling; max: the maximum font-size used for labels when scaling]", icon("info-circle", style="font-size: 11px")), value=c(10, 20), min=1, max=50, step=1))
										),
										br()
									),
									tabPanel("DNA_SV_panel",
										fluidRow(column(12, actionButton("net_dna_plot", label = div("Plot / Refresh", title = "Please note that changes made to the data will not be reflected in the plot until the ‘plot/refresh’ button is clicked on")))),
										br(),
										fluidRow(column(12, tags$p(tags$h5("Parameter Setting:", 
														style="margin-left: 0em; text-align: left; font-weight: bold; text-decoration: underline")))),
										fluidRow(
											column(12, checkboxInput(inputId="navi_dna", label = div("Display Navigation", style = "display:inline-block; font-size:13px", title = "Hide the navigation control panel if unchecked", icon("info-circle", style="font-size: 11px")), TRUE)),
											column(12, numericInput(inputId="grav_cons_dna", label = div("Gravitational Constant", style = "display:inline-block; font-size:13px", title = "Adjust the distance of a sub-graph where a decrease of this value will increase the repulsion of the sub-graph", icon("info-circle", style="font-size: 11px")), value=-800, min=-100000, max=100000)),
											column(12, numericInput(inputId="cent_grav_dna", label = div("Central Gravity", style = "display:inline-block; font-size:13px", title = "Set a central gravity attraction to pull the entire network back to the center [a larger value shows a stronger central force]", icon("info-circle", style="font-size: 11px")), value=0.1, min=0, max=1)),
											column(12, numericInput(inputId="spri_cons_dna", label = div("Spring Constant", style = "display:inline-block; font-size:13px", title = "Control the size of the subgraph [a smaller value shows a shorter edge]", icon("info-circle", style="font-size: 11px")), value=0.02, min=0, max=1)),
											column(12, numericInput(inputId="edge_leng_dna", label = div("Edge length", style = "display:inline-block; font-size:13px", title = "Set the length of edge in the graph", icon("info-circle", style="font-size: 11px")), value=10, min=1, max=30)),
											column(12, selectizeInput(inputId="net_select_dna", label = div("Node search", style = "display:inline-block; font-size:13px", title = "Choose and centralize a node [i.e., gene name]", icon("info-circle", style="font-size: 11px")), choices=NULL, multiple=F,
													options=list(placeholder = 'select', onInitialize = I('function() { this.setValue(""); }')))),
											column(12, numericInput(inputId="node_font_s_dna", label = div("Node font size", style = "display:inline-block; font-size:13px", title = "Set the font size of node name", icon("info-circle", style="font-size: 11px")), value=14, min=1, max=40)),
											column(12, numericInput(inputId="node_font_p_dna", label = div("Node font pos", style = "display:inline-block; font-size:13px", title = "Adjust the position of node name", icon("info-circle", style="font-size: 11px")), value=-44, min=-100, max=100)),
											column(12, sliderInput(inputId="node_scale_dna", label = div("Node scaling", style = "display:inline-block; font-size:13px", title = "Scale the size of node in proportion to the support number of samples [min: the minimum font-size used for labels when scaling; max: the maximum font-size used for labels when scaling]", icon("info-circle", style="font-size: 11px")), value=c(10, 20), min=1, max=50, step=1))
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

