options(shiny.maxRequestSize=2*1024*1024^2) #// max upload file size 2Gb
options(ucscChromosomeNames=FALSE)

	server <- function(input, output, session) {
		#/////////////////////////////////
		# message menu report in header //
		#/////////////////////////////////
		icon_g <- icon("github", verify_fa = FALSE)
		icon_g[["attribs"]][["class"]] <- "fa fa-github"
		output$messageMenu <- renderMenu({
			dropdownMenu(type = "messages",
				#// 1st option
				messageItem(from = "FuSViz in Github", message = "Code, Source, Documentation", icon = icon_g, href = "https://github.com/senzhaocode/FuSViz"),
				#// 2nd option
				messageItem(from = "Issues", message = "Report issues and bugs", icon = icon("exclamation-circle", verify_fa = FALSE), href = "https://github.com/senzhaocode/FuSViz/issues"),
				badgeStatus = NULL,
				icon = icon("info-circle", verify_fa = FALSE),
				headerText = "App Information")
		})

		#/////////////////////////////////////////////////////////
		# Load genomics/transcriptomics database (hg19 or hg38) //
		#/////////////////////////////////////////////////////////
		database <- reactiveValues(txdb=NULL, whole_txdb=NULL, grTrack=NULL, chrTrack=NULL, domain=NULL, motif=NULL, genome_cir=NULL, gene_range=NULL, karyto=NULL, ensembl_id=NULL, symbol_ensem=NULL, canonical=NULL, cancergenes=NULL, chrom_cir=NULL, organism=NULL);
		observeEvent(input$Import_genome_data, {
			#// load annotation db
			if ( input$genome == "" ) {
				shiny::showModal(modalDialog(title = "Warning message", "Please choose genome version!"));	req(NULL);
			} else {
				version = gsub('_offline', '', input$genome);
				organism = NULL;
				chrom_cir = NULL;
				if ( version == 'hg19' || version == 'hg38' ) {
					organism = "Human";
					chrom_cir = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY");
				} else if ( version == 'GRCm39' ) {
					organism = "Mouse";
					chrom_cir = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"); 
				}
				shiny::withProgress(message='Load genomic/transcriptomic annotations', detail="Please wait for a while...", min=0, max=1, value=0.1, {
					txdb <- suppressWarnings(suppressPackageStartupMessages(AnnotationDbi::loadDb(file=file.path(extdata, paste(organism, ".gencode.annotation.", version, ".sqlite", sep="")))));
					shiny::incProgress(0.2);
					load(file=file.path(extdata, paste(organism, ".grTrack.", version, ".Rd", sep="")));
					shiny::incProgress(0.1);
					load(file=file.path(extdata, paste(organism, ".cytoband.", version, ".Rd", sep="")));
					shiny::incProgress(0.1);
					whole_txdb <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names=TRUE); # group exons by transcript_id
					shiny::incProgress(0.2);

					#// load domain annotation - one data.frame: domain
					load(file=file.path(extdata, paste(organism, ".Domain_interval.", version, ".Rd", sep="")));
					shiny::incProgress(0.1);
					# // load motif annotation - one data.frame: motif
					load(file=file.path(extdata, paste(organism, ".Motif_interval.", version, ".Rd", sep="")));
					shiny::incProgress(0.1);

					#// create genome length for circle plot
					genome_cir = NULL;
					if ( input$genome == "hg38" || input$genome == "hg38_offline" ) {
						genome_cir = list("chr1"=248956422, "chr2"=242193529, "chr3"=198295559, "chr4"=190214555, "chr5"=181538259, "chr6"=170805979, "chr7"=159345973, "chr8"=145138636,
							"chr9"=138394717, "chr10"=133797422, "chr11"=135086622, "chr12"=133275309, "chr13"=114364328, "chr14"=107043718, "chr15"=101991189, "chr16"=90338345,
							"chr17"=83257441, "chr18"=80373285, "chr19"=58617616, "chr20"=64444167, "chr21"=46709983, "chr22"=50818468, "chrX"=156040895, "chrY"=57227415);
					} else if (  input$genome == "hg19" || input$genome == "hg19_offline" ) {
						genome_cir = list("chr1"=249250621, "chr2"=243199373, "chr3"=198022430, "chr4"=191154276, "chr5"=180915260, "chr6"=171115067, "chr7"=159138663, "chr8"=146364022, 
							"chr9"=141213431, "chr10"=135534747, "chr11"=135006516, "chr12"=133851895, "chr13"=115169878, "chr14"=107349540, "chr15"=102531392, "chr16"=90354753,
							"chr17"=81195210, "chr18"=78077248, "chr19"=59128983, "chr20"=63025520, "chr21"=48129895, "chr22"=51304566, "chrX"=155270560, "chrY"=59373566);
					} else if (  input$genome == "GRCm39" || input$genome == "GRCm39_offline" ) {
						genome_cir = list("chr1"=195154279, "chr2"=181755017, "chr3"=159745316, "chr4"=156860686, "chr5"=151758149, "chr6"=149588044, "chr7"=144995196, "chr8"=130127694,
							"chr9"=124359700, "chr10"=130530862, "chr11"=121973369, "chr12"=120092757, "chr13"=120883175, "chr14"=125139656, "chr15"=104073951, "chr16"=98008968,
							"chr17"=95294699, "chr18"=90720763, "chr19"=61420004, "chrX"=169476592, "chrY"=91455967);
					}

					#// load canonical transcript id
					load(file=file.path(extdata, paste(organism, ".canonical.Rd", sep="")));
					#// Load gene name data: gene_id a data.frame class with two columns: 'ensembl_gene_id' and 'gene_symbol'
					load(file=file.path(extdata, paste(organism, ".ensembl_symbol.Rd", sep="")));
					#// assign the ensembl_id to gene symbol
					symbol_ensem = gene_id$Gene_name;
					names(symbol_ensem) = gene_id$Gene_ID;
					shiny::incProgress(0.1);

					#// Load predefined gene dataset for oncogenes, tumor-suppress genes and cancer-related genes
					load(file=file.path(extdata, paste(organism, ".cancergenes.Rd", sep="")));
					shiny::incProgress(0.1);

					#// create gene annotation for circle plot
					gene_range = as.data.frame(genes(txdb), stringsAsFactors=F);
					names(gene_range)[6] = "Gene_ID";
					ensembl_id = unique(gene_range$Gene_ID);
					#// union of x and y by "Gene_ID" - merge(gene_range, gene_id, by="Gene_ID", all.x=T, all.y=T)
					#// intersection of x and y by "Gene_ID" - merge(gene_range, gene_id, by="Gene_ID", all.x=F, all.y=F)
					#// only use "Gene_id" of y - merge(gene_range, gene_id, by="Gene_ID", all.x=F, all.y=T)
					gene_range=merge(gene_range, gene_id, by="Gene_ID", all.x=T, all.y=F); #// only use "Gene_id" of x
					gene_range$seqnames = as.character(gene_range$seqnames);
					gene_range = gene_range[gene_range$seqnames %in% chrom_cir, ];
					names(gene_range)[7] = "gene";	names(gene_range)[2] = "chr";
					gene_range$strand = as.character(gene_range$strand);
					#// if gene_symbol is NA, substitute with ensembl_id
					gene_update=apply(gene_range, 1, function(x){
        					if ( is.na(x[7]) ) { return(as.character(x[1])); } else { return(as.character(x[7])); }
					})
					gene_range$gene = gene_update
					shiny::incProgress(0.1);

					#// Create cytoband info for circule plot
					karyto = chrTrack@bandTable;
					karyto = karyto[karyto$chrom %in% chrom_cir, ];
					names(karyto)[1] = "chr"; names(karyto)[2] = "start"; names(karyto)[3] = "end";
					karyto$gieStain = as.character(karyto$gieStain);
					karyto$name = as.character(karyto$name);
					cytoband_color = data.frame(gieStain=c("acen", "gneg", "gpos100", "gpos25", "gpos50", "gpos75", "gvar", "stalk"),
								color=c("#8b2323", "#ffffff", "#000000", "#C1C1C1", "#808080", "#404040", "#000000", "#cd3333"), stringsAsFactors=F);
					karyto = merge(karyto, cytoband_color, by="gieStain", all.x=T, all.y=F); #// only use "gieStain" of x
					karyto = karyto[order(karyto$chr), ]; #// re-order by chromosome names
					shiny::incProgress(0.1);
				})
				database$txdb=txdb;	database$whole_txdb=whole_txdb;	database$grTrack=grTrack;	database$chrTrack=chrTrack;	database$domain=domain;	database$symbol_ensem=symbol_ensem;	
				database$cancergenes=cancergenes;	database$motif=motif;	database$genome_cir=genome_cir;	database$gene_range=gene_range;	database$karyto=karyto; 
				database$ensembl_id=ensembl_id;	database$canonical=canonical;	database$chrom_cir=chrom_cir;	database$organism=organism;
			}
		})
		observe({
        	if ( input$Import_genome_data == 0 || input$genome == "") {
            	session$sendCustomMessage(type="disablebutton", list(code=paste("$('.btn-file').attr('disabled', 'true').attr('onclick', 'return false;')",sep="")))
        	} else {
            	session$sendCustomMessage(type="disablebutton", list(code=paste("$('.btn-file').removeAttr('disabled').attr('onclick', 'return true;')",sep="")))
       		}
    	})
		observe({ req(database$txdb); })

		#// inputFile data structure stores upload rna SV, dna SV and mutation profile
		inputFile <- reactiveValues(rnadata=NULL, dnadata=NULL, mutationdata=NULL);
		#/////////////////////////////////////////////////////////////
		#// Load external RNA-seq sv files as the data.frame class  //
		#/////////////////////////////////////////////////////////////
		observeEvent(input$file_rna_data$datapath, {
			#// Avoid 'input$file_rna_data' infinite problem when no file loading; the req() basically aborts the rest of the block 
			req(input$file_rna_data);
			#// Upload and read inputfile of RNA SVs (13 columns)
			tumordata = read.csv(input$file_rna_data$datapath, header=TRUE, sep=input$sep_rna_file, quote="");
			col_num_rna = colnames(tumordata);
			if ( length(col_num_rna) < 13 ) { showModal(modalDialog(title = "Error message", "RNA SV input file does not meet requirement - column number MUST be 12!")); req(NULL); }
			#// NOTE: column1-11 (chrom1, pos1, gene1, chrom2, pos2, gene2, name, split, span, strand1, strand2, untemplated_insert, comment) - NA unaccepted.
			if ( col_num_rna[1] != "chrom1" || any(is.na(tumordata$chrom1)) == T || any(tumordata$chrom1 == "") == T ) { showModal(modalDialog(title = "Error message", "'chrom1' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[2] != "pos1" || any(is.na(tumordata$pos1)) == T || any(tumordata$pos1 <= 0) == T ) { showModal(modalDialog(title = "Error message", "'pos1' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[3] != "gene1" || any(is.na(tumordata$gene1)) == T || any(tumordata$gene1 == "") == T ) { showModal(modalDialog(title = "Error message", "'gene1' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[4] != "chrom2" || any(is.na(tumordata$chrom2)) == T || any(tumordata$chrom2 == "") == T ) { showModal(modalDialog(title = "Error message", "'chrom2' column has an incorrect header name or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[5] != "pos2" || any(is.na(tumordata$pos2)) == T || any(tumordata$pos2 <= 0) == T ) { showModal(modalDialog(title = "Error message", "'pos2' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[6] != "gene2" || any(is.na(tumordata$gene2)) == T || any(tumordata$gene2 == "") == T ) { showModal(modalDialog(title = "Error message", "'gene2' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[7] != "name" || any(is.na(tumordata$name)) == T || any(tumordata$name == "") == T ) { showModal(modalDialog(title = "Error message", "'name' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[8] != "split" || any(is.na(tumordata$split)) == T || any(tumordata$split < 0) == T ) { showModal(modalDialog(title = "Error message", "'split' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[9] != "span" || any(is.na(tumordata$span)) == T || any(tumordata$sapn < 0) == T ) { showModal(modalDialog(title = "Error message", "'span' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[10] != "strand1" || any(is.na(tumordata$strand1)) == T || any(tumordata$strand1 == "") == T ) { showModal(modalDialog(title = "Error message", "'strand1' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[11] != "strand2" || any(is.na(tumordata$strand2)) == T || any(tumordata$strand2 == "") == T ) { showModal(modalDialog(title = "Error message", "'strand2' column has an incorrect header or empty/NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[12] != "untemplated_insert" ) { showModal(modalDialog(title = "Error message", "'untemplated_insert' column has an incorrect header for RNA SVs!")); req(NULL); }
			tumordata$untemplated_insert[is.na(tumordata$untemplated_insert)] = ""; # if NA for untemplated_insert, fill in ""
			if ( any(is.na(tumordata$untemplated_insert)) == T ) { showModal(modalDialog(title = "Error message", "'untemplated_insert' column has NA value for RNA SVs!")); req(NULL); }
			if ( col_num_rna[13] != "comment" ) { showModal(modalDialog(title = "Error message", "'comment' column has an incorrect header for RNA SVs!")); req(NULL); }
			tumordata$comment[is.na(tumordata$comment)] = ""; # if NA for comment, fill in ""
			if ( any(is.na(tumordata$comment)) == T ) { showModal(modalDialog(title = "Error message", "'comment' column has NA value for RNA SVs!")); req(NULL); }
			if ( database$organism == 'Human' ) {
				if ( FALSE %in% grepl("^[A-Z][\\.A-Za-z0-9_-]+$", tumordata$gene1) ||  FALSE %in% grepl("^[A-Z][\\.A-Za-z0-9_-]+$", tumordata$gene2) ) { # lowercase present in symbol
					showModal(modalDialog(title = "Error message", "'gene1' and 'gene2' columns has invalid symbol names for human!")); req(NULL); 
				}
			} else if  ( database$organism == 'Mouse' ) {
				if ( FALSE %in% grepl("^[A-Z][\\.a-z0-9_-]+$", tumordata$gene1) ||  FALSE %in% grepl("^[A-Z][\\.a-z0-9_-]+$", tumordata$gene2) ) { # lowercase no present in symbol
					showModal(modalDialog(title = "Error message", "'gene1' and 'gene2' columns has invalid symbol names for mouse!")); req(NULL); 
				}
			}
			no_empty=tumordata$untemplated_insert[nzchar(tumordata$untemplated_insert)]
			if ( length(no_empty) > 0 ) {
				if ( TRUE %in% grepl("[^A|^T|^G|^C|^N|^a|^t|^g|^c|^n]", no_empty)  ) {
					showModal(modalDialog(title = "Error message", "Unexpected letter is present in the string in addition to A|a, C|c, G|g, T|t and N|n!")); req(NULL);
				}
			}
			inputFile$rnadata = tumordata;
		})
		#// edit and update cell value of rna SV table
		observeEvent(input[["RNAcontents_cell_edit"]], {
			rnarow = input$RNAcontents_cell_edit$row;
			rnaclmn = input$RNAcontents_cell_edit$col;
			info = input[["RNAcontents_cell_edit"]];
			rnavalue = inputFile$rnadata[rnarow, rnaclmn];
			if ( rnaclmn == 1 || rnaclmn == 4 ) {
				if ( info[["value"]] %in% database$chrom_cir ) {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
				} else {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					info[["value"]] <- rnavalue;
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$RNAcontents_cell_edit$value, " is not a valid chromosome name! NOTE: chromosome name should be a 'chr1-22-X-Y-M'!", sep="")));
				}
			} else if ( rnaclmn == 2 || rnaclmn == 5 || rnaclmn == 8 || rnaclmn == 9) {
				if ( is(info[["value"]], "numeric") ) {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
				} else {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					info[["value"]] <- rnavalue;
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$RNAcontents_cell_edit$value, " is not a numeric feature! Please input a number!", sep="")));
				}
			} else if ( rnaclmn == 3 || rnaclmn == 6 ) {
				genename = NULL;
				if ( grepl("ENSG00", info[["value"]]) || grepl("ENSMUSG0", info[["value"]]) ) {
					genename  = database$ensembl_id[database$ensembl_id == input$RNAcontents_cell_edit$value];
				} else {
					genename = names(database$symbol_ensem[database$symbol_ensem == input$RNAcontents_cell_edit$value]);
				}
				if ( length(genename) > 0 ) {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
				} else {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					info[["value"]] <- rnavalue;
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$RNAcontents_cell_edit$value, " is not a valid gene name! Please input valid gene name.", sep="")));
				}
			} else if ( rnaclmn == 7 || rnaclmn == 12 || rnaclmn == 13 ) {
				if ( is(info[["value"]], "character") ) {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
				} else {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					info[["value"]] <- rnavalue;
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$RNAcontents_cell_edit$value, " is not a valid sample name! NOTE: sample name should be a character!", sep="")));
				}
			} else if ( rnaclmn == 10 ||  rnaclmn == 11 ) {
				if ( info[["value"]] == '+' || info[["value"]] == '-' ) {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
				} else {
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					info[["value"]] <- rnavalue;
					inputFile$rnadata <<- editData(inputFile$rnadata, info, 'RNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$RNAcontents_cell_edit$value, " is not a valid strand direction! NOTE: strand direction should be +/-!", sep="")));
				}
			}
    	})

		#/////////////////////////////////////////////////////////////
		#// Load external DNA-seq sv files as the data.frame class  //
		#/////////////////////////////////////////////////////////////
		observeEvent(input$file_dna_data$datapath, {
			#// Avoid 'input$file_dna_data' infinite problem when no file loading; the req() basically aborts the rest of the block
			req(input$file_dna_data);
			#// Upload and read inputfile of DNA SVs (13 columns)
			dnadata = data.table::fread(input$file_dna_data$datapath, sep=input$sep_dna_file, stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE,
						header = TRUE, fill = TRUE, quote = "");
			col_num_dna = colnames(dnadata);
			if ( length(col_num_dna) < 13 ) { showModal(modalDialog(title = "Error message", "DNA SV input file does not meet requirement - column number MUST be 13!")); req(NULL); }
			#// NOTE: column1-13 (chrom1, start1, end1, chrom2, start2, end2, name, type, split, span, gene1, gene2, comment) - no NA accepted.
			if ( col_num_dna[1] != "chrom1" || any(is.na(dnadata$chrom1)) == T || any(dnadata$chrom1 == "") == T ) { showModal(modalDialog(title = "Error message", "'chrom1' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[2] != "start1" || any(is.na(dnadata$start1)) == T || any(dnadata$start1 <= 0) == T ) { showModal(modalDialog(title = "Error message", "'start1' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[3] != "end1" || any(is.na(dnadata$end1)) == T || any(dnadata$end1 <= 0) == T ) { showModal(modalDialog(title = "Error message", "'end1' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[4] != "chrom2" || any(is.na(dnadata$chrom2)) == T || any(dnadata$chrom2 == "") == T ) { showModal(modalDialog(title = "Error message", "'chrom2' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[5] != "start2" || any(is.na(dnadata$start2)) == T || any(dnadata$start2 <= 0) == T ) { showModal(modalDialog(title = "Error message", "'start2' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[6] != "end2" || any(is.na(dnadata$end2)) == T || any(dnadata$end1 <= 0) == T ) { showModal(modalDialog(title = "Error message", "'end2' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[7] != "name" || any(is.na(dnadata$name)) == T || any(dnadata$name == "") == T ) { showModal(modalDialog(title = "Error message", "'name' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[8] != "type" || any(is.na(dnadata$type)) == T || any(dnadata$type == "") == T ) { showModal(modalDialog(title = "Error message", "'type' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[9] != "split" || any(is.na(dnadata$split)) == T || any(dnadata$split < 0) == T ) { showModal(modalDialog(title = "Error message", "'split' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[10] != "span" || any(is.na(dnadata$span)) == T || any(dnadata$sapn < 0) == T ) { showModal(modalDialog(title = "Error message", "'span' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[11] != "gene1" || any(is.na(dnadata$gene1)) == T || any(dnadata$gene1 == "") == T ) { showModal(modalDialog(title = "Error message", "'gene1' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[12] != "gene2" || any(is.na(dnadata$gene2)) == T || any(dnadata$gene2 == "") == T ) { showModal(modalDialog(title = "Error message", "'gene2' column has an incorrect header or empty/NA value for DNA SVs!")); req(NULL); }
			if ( col_num_dna[13] != "comment" ) { showModal(modalDialog(title = "Error message", "'comment' column has an incorrect header for DNA SVs!")); req(NULL); }			
			dnadata$comment[is.na(dnadata$comment)] = ""; # if NA for comment, fill in ""
			if ( any(is.na(dnadata$comment)) == T ) { showModal(modalDialog(title = "Error message", "'comment' column has NA value for DNA SVs!")); req(NULL); }			
			if ( database$organism == 'Human' ) {
				if ( FALSE %in% grepl("^[A-Z][\\.A-Za-z0-9_-]+$|^\\*$", dnadata$gene1) ||  FALSE %in% grepl("^[A-Z][\\.A-Za-z0-9_-]+$|^\\*$", dnadata$gene2) ) { # lowercase present in symbol
					showModal(modalDialog(title = "Error message", "'gene1' and 'gene2' columns has invalid symbol names for human!")); req(NULL); 
				}
			} else if  ( database$organism == 'Mouse' ) {
				if ( FALSE %in% grepl("^[A-Z][\\.a-z0-9_-]+$|^\\*$", dnadata$gene1) ||  FALSE %in% grepl("^[A-Z][\\.a-z0-9_-]+$|^\\*$", dnadata$gene2) ) { # lowercase no present in symbol
					showModal(modalDialog(title = "Error message", "'gene1' and 'gene2' columns has invalid symbol names for mouse!")); req(NULL); 
				}
			}
			inputFile$dnadata = dnadata;
		})
		#// edit and update cell value of dna SV table
		observeEvent(input$DNAcontents_cell_edit, {
			dnarow = input$DNAcontents_cell_edit$row;
			dnaclmn = input$DNAcontents_cell_edit$col;
			info = input[["DNAcontents_cell_edit"]];
			dnavalue = inputFile$dnadata[dnarow, dnaclmn, with=FALSE];
			inputFile$dnadata <<- editData(inputFile$dnadata, input$DNAcontents_cell_edit, 'DNAcontents')
			if ( dnaclmn == 1 || dnaclmn == 4 ) {
				if ( info[["value"]] %in% database$chrom_cir ) {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
				} else {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					info[["value"]] <- dnavalue;
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$DNAcontents_cell_edit$value, " is not a valid chromosome name! NOTE: chromosome name should be a 'chr1-22-X-Y-M'!", sep="")));
				}
			} else if ( dnaclmn == 2 || dnaclmn == 3 || dnaclmn == 5 || dnaclmn == 6 || dnaclmn == 9 || dnaclmn == 10) {
				if ( is(info[["value"]], "numeric") ) {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
				} else {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					info[["value"]] <- dnavalue;
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$DNAcontents_cell_edit$value, " is not a numeric feature! NOTE: input a number!", sep="")));
				}
			} else if ( dnaclmn == 11 || dnaclmn == 12 ) {
				genename = NULL;
				if ( info[["value"]] == '*' ) {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
				} else {
					if ( grepl("ENSG00", info[["value"]]) || grepl("ENSMUSG0", info[["value"]]) ) {
						genename = database$ensembl_id[database$ensembl_id == input$DNAcontents_cell_edit$value];
					} else {
						genename = names(database$symbol_ensem[database$symbol_ensem == input$DNAcontents_cell_edit$value]);
					}
					if ( length(genename) > 0 ) {
						inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					} else {
						inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
						info[["value"]] <- dnavalue;
						inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
						showModal(modalDialog(title = "Warning message", paste(input$DNAcontents_cell_edit$value, " is not a valid gene name! Please input valid gene name!", sep="")));
					}
				}
			} else if ( dnaclmn == 7 || dnaclmn == 13 ) {
				if ( is(info[["value"]], "character") ) {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
				} else {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					info[["value"]] <- dnavalue;
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$DNAcontents_cell_edit$value, " is not a valid sample name! NOTE: sample name should be a character!", sep="")));
				}
			} else if ( dnaclmn == 8 ) {
				if ( info[["value"]] %in% c("BND", "DEL", "DUP", "INS", "INV") ) {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
				} else {
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					info[["value"]] <- dnavalue;
					inputFile$dnadata <<- editData(inputFile$dnadata, info, 'DNAcontents')
					showModal(modalDialog(title = "Warning message", paste(input$DNAcontents_cell_edit$value, " is not a valid SV type! NOTE: SV type should be 'BND, DEL, DUP, INS or INV'!", sep="")));
				}
			}
    	})

		#//////////////////////////////////
		#// load drug target information //
		#//////////////////////////////////
		drug_target_match <- reactive({
			tmp = NULL;
			if ( input$genome == "GRCm39" ) { return(tmp); }
			if ( is.null(input$file_rna_data) ) { # dna_sv file not load
				if ( is.null(input$file_dna_data) ) { # rna_sv file not load
					return(NULL);
				} else {
					tmp = drug_target[drug_target$symbol %in% unique(c(as.character(inputFile$dnadata$gene1), as.character(inputFile$dnadata$gene2))), ]
				}
			} else {
				if ( is.null(input$file_dna_data) ) { # rna_sv file load
					tmp = drug_target[drug_target$symbol %in% unique(c(as.character(inputFile$rnadata$gene1), as.character(inputFile$rnadata$gene2))), ];
				} else {
					tmp = drug_target[drug_target$symbol %in% unique(c(as.character(inputFile$dnadata$gene1), as.character(inputFile$dnadata$gene2), as.character(inputFile$rnadata$gene1), as.character(inputFile$rnadata$gene2))), ];
				}
			}
			if ( nrow(tmp) > 0 ) {
				tmp$molecule_chembl_id = paste0('<a href="https://www.ebi.ac.uk/chembl/g/#search_results/all/query=', tmp$molecule_chembl_id, '" target="_blank">', tmp$molecule_chembl_id, '</a>')
				tmp$target_chembl_id = paste0('<a href="https://www.ebi.ac.uk/chembl/g/#search_results/all/query=', tmp$target_chembl_id, '" target="_blank">', tmp$target_chembl_id, '</a>')
			}
			return(tmp)
		})

		#/////////////////////////////////////////////////////////////////////////////
		#// Load external Mutation (SNV+indel) input files as the data.frame class  //
		#/////////////////////////////////////////////////////////////////////////////
		observeEvent(input$file_maf_data$datapath, { #// inputFile$mutationdata is a data.table class
			#// Avoid 'input$file_maf_data' infinite problem when no file loading; the req() basically aborts the rest of the block
			req(input$file_maf_data); 
			#// Upload and read inputfile of DNA mutations (maf format), subset samples matched to SVs from RNA-seq and DNA-seq
			mutdata = data.table::fread(input$file_maf_data$datapath, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, 
					showProgress = TRUE, header = TRUE, fill = TRUE, skip = "Hugo_Symbol", quote = "", 
					select = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Tumor_Sample_Barcode'));
			col_num_mut = colnames(mutdata);
			if ( length(col_num_mut) != 8 ) { showModal(modalDialog(title = "Error message", "Eight required columns (i.e. 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification' and 'Tumor_Sample_Barcode') of Mutation MAF file fail to be met!")); req(NULL); }
			#// NOTE: column1-8 (Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Tumor_Sample_Barcode) - no NA accepted.
			if ( any(is.na(mutdata$Chromosome)) == T || any(mutdata$Chromosome == "") == T ) { showModal(modalDialog(title = "Error message", "'Chromosome' column has an incorrect header or empty/NA value for mutation profile!")); req(NULL); }
			if ( any(is.na(mutdata$Start_Position)) == T || any(mutdata$Start_Position <= 0) == T ) { showModal(modalDialog(title = "Error message", "'Start_Position' column has an incorrect header or empty/NA value for mutation profile!")); req(NULL); }
			if ( any(is.na(mutdata$End_Position)) == T || any(mutdata$End_Position <= 0) == T ) { showModal(modalDialog(title = "Error message", "'End_Position' column has an incorrect header or empty/NA value for mutation profile!")); req(NULL); }
			if ( any(is.na(mutdata$Reference_Allele)) == T || any(mutdata$Reference_Allele == "") == T ) { showModal(modalDialog(title = "Error message", "'Reference_Allele' column has an incorrect header or empty/NA value for mutation profile!")); req(NULL); }
			if ( any(is.na(mutdata$Tumor_Seq_Allele2)) == T || any(mutdata$Tumor_Seq_Allele2 == "") == T ) { showModal(modalDialog(title = "Error message", "'Tumor_Seq_Allele2' column has an incorrect header or empty/NA value for mutation profile!")); req(NULL); }
			if ( any(is.na(mutdata$Variant_Classification)) == T || any(mutdata$Variant_Classification == "") == T ) { showModal(modalDialog(title = "Error message", "'Variant_Classification' column has an incorrect header or empty/NA value for mutation profile!")); req(NULL); }
			if ( any(is.na(mutdata$Tumor_Sample_Barcode)) == T || any(mutdata$Tumor_Sample_Barcode == "") == T ) { showModal(modalDialog(title = "Error message", "'Tumor_Sample_Barcode' column has an incorrect header or empty/NA value for mutation profile!")); req(NULL); }

			if ( is.null(input$file_rna_data) ) { # dna_sv file not load
				if ( is.null(input$file_dna_data) ) { # rna_sv file not load
					shiny::showModal(modalDialog(title = "Warning message", "Before load mutation data, please load SV data firstly!!!"));	req(NULL);
				} else { # dna_sv file load 
					mutdata = mutdata[Tumor_Sample_Barcode %in% unique(as.character(inputFile$dnadata$name))];
				}
			} else {
				if ( is.null(input$file_dna_data) ) { # rna_sv file load
					mutdata = mutdata[Tumor_Sample_Barcode %in% unique(as.character(inputFile$rnadata$name))];
				} else { # both dna and rna files load
					mutdata = mutdata[Tumor_Sample_Barcode %in% unique(c(as.character(inputFile$dnadata$name), as.character(inputFile$rnadata$name)))];
				}
			}
			#// remove duplication in mutation data
			mutdata = mutdata[, variantId := paste(Chromosome, Start_Position, Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2, sep = ':')];
			if(nrow(mutdata[duplicated(variantId)]) > 0) {
				cat("--Removed",  nrow(mutdata[duplicated(variantId)]) ,"duplicated variants\n");
				mutdata = mutdata[!duplicated(variantId)];
			}
			mutdata[,variantId := NULL];
			#// check empty value present in Hugo_Symbol (if yes assign 'UnknownGene')
			if(nrow(mutdata[Hugo_Symbol %in% ""]) > 0) {
				cat('--Found ', nrow(mutdata[Hugo_Symbol %in% ""]), ' variants with no Gene Symbols\n');
				mutdata$Hugo_Symbol = ifelse(test = mutdata$Hugo_Symbol == "", yes = 'UnknownGene', no = mutdata$Hugo_Symbol);
			}
			#// check NA present in Hugo_Symbol (if yes assign 'UnknownGene')
			if(nrow(mutdata[is.na(Hugo_Symbol)]) > 0) {
				cat('--Found ', nrow(mutdata[is.na(Hugo_Symbol) > 0]), ' variants with no Gene Symbols\n');
				mutdata$Hugo_Symbol = ifelse(test = is.na(mutdata$Hugo_Symbol), yes = 'UnknownGene', no = mutdata$Hugo_Symbol);
			}
			#// combine allele to Variant_Class
			mutdata = mutdata[, anno := as.character(paste(Reference_Allele, '>', Tumor_Seq_Allele2, '(', Variant_Classification, ')', sep = ''))];
			#// convert data.table to data.frame, and count freq of 'Hugo_Symbol-anno-Turmo_Sample_Barcode' per 'Chrom-Position'
			mutdata = mutdata[, c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "anno")];
			#// chromosome control not including 'chrM'
			mutdata = mutdata[mutdata$Chromosome %in% database$chrom_cir, ]; # chromosome control, not including 'chrM'
			inputFile$mutationdata = mutdata
		})
		mutation_bedgraph <- reactiveValues(value=NULL); #// summarize mutation freq 
		observe({ #// data.table class
			if ( is.null(inputFile$mutationdata) ) { req(NULL); }
			tmp = as.data.table(inputFile$mutationdata);
			tmp = tmp[, keyby = .(Hugo_Symbol, Chromosome, Start_Position, anno), label := paste(Tumor_Sample_Barcode, collapse= ', ')];
			tmp = tmp[, keyby = .(Hugo_Symbol, Chromosome, Start_Position, anno, label), .(value = uniqueN(Tumor_Sample_Barcode))];
			tmp$anno_label = paste(tmp$anno, '[', tmp$label, ']', sep = '');
			tmp = tmp[, keyby = .(Hugo_Symbol, Chromosome, Start_Position), anno_new := paste(anno_label, collapse= ' | ')];
			tmp = tmp[, keyby = .(Hugo_Symbol, Chromosome, Start_Position, anno_new), .(freq = sum(value))];
			tmp$anno_new = paste(tmp$Hugo_Symbol, ": ", tmp$anno_new, sep="");
			mutation_bedgraph$value = tmp;
		})
		#///////////////////////////////////////////////////////////
		# create 'wordcloud_data' vector class for wordcloud plot //
		#///////////////////////////////////////////////////////////
		#// wordcloud_svrna => for partner genes from RNA-seq
		wordcloud_svrna <- reactive({ # list with one data.frame class and one vector class
			if ( is.null(inputFile$rnadata) ) { req(NULL); }
			word = NULL;	word = inputFile$rnadata[,c("name","gene1","gene2")]; #// select 'name', 'gene1' and 'gene2'
			word = FuSViz::wordcloud_processs(word, input$gene_freq_rna, "RNA", database$cancergenes, onco_color, supp_color, rela_color);
			return(word);
		})
		#// wordcloud_svdna => for partner genes from DNA-seq
		wordcloud_svdna <- reactive({ # list with one data.frame class and one vector class
			if ( is.null(inputFile$dnadata) ) { req(NULL); }
			word = NULL;    word = inputFile$dnadata[,c("name","gene1","gene2")]; #// select 'name', 'gene1' and 'gene2'
			word = word[! (word$gene1 == '*' & word$gene2 == '*'), ]; #// delete row with both gene1 and gene1 == '*'
			word = FuSViz::wordcloud_processs(word, input$gene_freq_dna, "DNA", database$cancergenes, onco_color, supp_color, rela_color);
			return(word);
		})
		#// wordcloud_mutat => for partner genes from mutation profile
		wordcloud_mutat <- reactive({ # list with one data.frame class and one vector class
			if ( is.null(inputFile$mutationdata) ) { req(NULL); }
			word = NULL; #// define a data.frame
			#// only keep the gene with significant change mutations
			change = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", 
			  "Missense_Mutation", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del");
			reg_exp = paste(change, collapse="|");
			word = inputFile$mutationdata[anno %like% reg_exp];
			word = word[,c("Tumor_Sample_Barcode","Hugo_Symbol")]; #// select 'name', 'gene'
			names(word) = c("name", "gene");
			word = FuSViz::wordcloud_processs(word, input$gene_freq_mut, "Mut", database$cancergenes, onco_color, supp_color, rela_color);
			return(word);
		})

		#/////////////////////////////////////////////////////////
		# create 'histgram' vector for histgram plot per sample //
		#/////////////////////////////////////////////////////////
		#// => for plotting rna SVs distribution per sample
		hist_svrna <- reactive({
			if ( is.null(inputFile$rnadata) ) { req(NULL); }
			hist = table(inputFile$rnadata$name);	
			hist = FuSViz::histgram_process(hist);
			return(hist);
		})
		#// => for plotting dna SVs distribution per sample
		hist_svdna <- reactive({
			if ( is.null(inputFile$dnadata) ) { req(NULL); }
			hist = inputFile$dnadata[,c("name", "type")];
			hist = as.data.frame(hist[, .N, by = c("name", "type")]);
			sample = unique(hist$name);	
			sample = FuSViz::histgram_process(sample);
			section1 = hist[hist$name %in% sample$section1, ];	section2 = hist[hist$name %in% sample$section2, ];
			section3 = hist[hist$name %in% sample$section3, ];	section4 = hist[hist$name %in% sample$section4, ];
			return(list(section1=section1, section2=section2, section3=section3, section4=section4));
		})
		#// => for ploting mutation load
		cor_mut <- reactive({
			if ( is.null(inputFile$mutationdata) ) { req(NULL); }
			hist = NULL;	hist = table(inputFile$mutationdata$Tumor_Sample_Barcode);
			rna_data_frame = NULL;	dna_data_frame = NULL;
			if (! is.null(inputFile$rnadata) && length(hist) > 0 ) {
				tmp_sv_rna = c(hist_svrna()$section1, hist_svrna()$section2, hist_svrna()$section3, hist_svrna()$section4);
				tmp_sv_rna = tmp_sv_rna[names(tmp_sv_rna) %in% names(hist)];
				rna_data_frame = data.frame(log2_sv_num=as.vector(tmp_sv_rna), log2_mut_num=as.vector(hist), row.names=names(hist), stringsAsFactors=F);
				rna_data_frame = log2(rna_data_frame);
			}
			if (! is.null(inputFile$dnadata) && length(hist) > 0 ) {
				tmp_sv_dna = rbind(hist_svdna()$section1, hist_svdna()$section2, hist_svdna()$section3, hist_svdna()$section4); #// data.frame
				tmp_sv_dna = aggregate(tmp_sv_dna$N, by=list(name=tmp_sv_dna$name), FUN=sum);	
				names(tmp_sv_dna) = c("name", "sum_type");
				tmp_sv_dna = tmp_sv_dna[tmp_sv_dna$name %in% names(hist), ];
				tmp_sv_dna = tmp_sv_dna[order(tmp_sv_dna$name), ];
				dna_data_frame = data.frame(log2_sv_num=as.numeric(tmp_sv_dna$sum_type), log2_mut_num=as.vector(hist), row.names=names(hist), stringsAsFactors=F);
				dna_data_frame = log2(dna_data_frame);
			}
			return(list(rna_data_frame=rna_data_frame, dna_data_frame=dna_data_frame));
		})

		#///////////////////////////////////////////////////////////////////////////
		#// Create 'circle_twocol' data.frame class for Circle RNA and DNA module //
		#///////////////////////////////////////////////////////////////////////////
		circle_twocol <- reactive({ # a data.frame class
			if ( is.null(inputFile$rnadata) ) { req(NULL); }
			twocol = inputFile$rnadata[inputFile$rnadata$chrom1 %in% database$chrom_cir & inputFile$rnadata$chrom2 %in% database$chrom_cir, ]; # chromosome control for circular plot
			#// 'twocol' has three columns "chrom", "gene" and "name"
			twocol = data.frame(chr=c(as.character(twocol$chrom1), as.character(twocol$chrom2)), 
					gene=c(as.character(twocol$gene1), as.character(twocol$gene2)),
					name=c(as.character(twocol$name), as.character(twocol$name)), stringsAsFactors=FALSE);
			twocol = unique(twocol);
			return(twocol);
		})
		circle_twocol_dna <- reactive({ # a data.frame class
			if ( is.null(inputFile$dnadata) ) { req(NULL); }
			twocol_dna = inputFile$dnadata[inputFile$dnadata$chrom1 %in% database$chrom_cir & inputFile$dnadata$chrom2 %in% database$chrom_cir, ]; # chromosome control for circular plot
			#// 'twocol' has two columns "chrom", "gene" and "name"
			twocol_dna = data.frame(chr=c(as.character(twocol_dna$chrom1), as.character(twocol_dna$chrom2)), 
					gene=c(as.character(twocol_dna$gene1), as.character(twocol_dna$gene2)),
					name=c(as.character(twocol_dna$name), as.character(twocol_dna$name)), stringsAsFactors=FALSE);
			## twocol_dna = twocol_dna[twocol_dna$gene != '*', ]; # move row without gene name (like *)
			twocol_dna = unique(twocol_dna);
			return(twocol_dna);
		})

		#////////////////////////////////////////////////////////
		#// Create a class for Linear module for Mutation data //
		#////////////////////////////////////////////////////////
		input_twoway_mut_bed <- reactive({ # a data.frame class
			if ( is.null(inputFile$mutationdata) ) { return(NULL); }
			#// Initialize data.frame class 'mutbed' and 'mutbedgraph' for RNA SVs
			mutbed = data.frame(chr=as.character(inputFile$mutationdata$Chromosome), start=as.numeric(inputFile$mutationdata$Start_Position)-1, 
								end=as.numeric(inputFile$mutationdata$End_Position), gene_id=as.character(inputFile$mutationdata$Hugo_Symbol), 
								name=as.character(inputFile$mutationdata$Tumor_Sample_Barcode), anno=as.character(inputFile$mutationdata$anno), stringsAsFactors=FALSE);
			#// filter by mutation type
			if (! is.null(input$mut_type_line) ) { 
				if ( input$mut_type_line[1] != "" ) {
					reg_exp = paste(input$mut_type_line, collapse="|");
					mutbed = mutbed[mutbed$anno %like% reg_exp, ];
				}
			}
			#// filter by mutation sample
			if (! is.null(input$mut_sample_line) ) { 
				if ( input$mut_sample_line[1] != "" ) {
					mutbed = mutbed[which(mutbed$name %in% input$mut_sample_line), ];
				}
			}
			return(mutbed);
		})

		#////////////////////////////////////////////////////////////////////////////////
		#// Create 'bed' and 'bedgraph' list class for Linear module using RNA-seq SVs //
		#////////////////////////////////////////////////////////////////////////////////
		#// 'input_twoway' has two elements: 'rnabed' - a data.frame class for bed format and 'rnabedgraph' - a data.frame class for tumorbedgraph format
		input_twoway_rna_bed <- reactive({
			#// Defined data.frame class 'rnabed' and 'rnabedgraph' for RNA SVs
			rnabed = NULL;	rnabedgraph = NULL;
			if ( is.null(inputFile$rnadata) ) { return(list(rnabed=NULL, rnabedgraph=NULL)) }
			#// Set num of split and span reads for filtering control
			rnatmp = inputFile$rnadata[inputFile$rnadata$split >= input$rna_split_bed & inputFile$rnadata$span >= input$rna_span_bed, ];
			#// Assign SV type either 'Intra' or 'Inter''
			type = apply(rnatmp, 1, function(x){
				if ( as.character(x[1]) == as.character(x[4]) ) { return("Intra"); } else { return("Inter"); }
			})
			rnatmp = cbind(rnatmp, type); #// add new column named type
			if (! is.null(input$rna_type_bed) ) {
				if ( input$rna_type_bed[1] != "" ) {
					rnatmp = rnatmp[which(rnatmp$type %in% input$rna_type_bed), ];
				}
			}
			dis_bed = apply(rnatmp, 1, function(x){
				y=as.numeric(x[5]) - as.numeric(x[2]);
				if ( as.character(x[12]) == "Inter" ) { return(TRUE); } else { abs(y) >= input$rna_dismin_bed && abs(y) <= input$rna_dismax_bed; }
			})
			rnatmp = cbind(rnatmp, dis_bed); rnatmp = rnatmp[rnatmp$dis_bed == T,]; rnatmp = rnatmp[, c(1,2,3,4,5,6,7,8,9,10,11)]; # only keep 1-9 column

			#// Set 'rnabed' as a bed format for IGV input - note bed format start coordinate [start - 1]
			rnabed = data.frame(chr=c(as.character(rnatmp$chrom1), as.character(rnatmp$chrom2)), start=c(as.numeric(rnatmp$pos1) - 1, as.numeric(rnatmp$pos2) - 1),
					end=c(as.numeric(rnatmp$pos1), as.numeric(rnatmp$pos2)), gene_id=c(as.character(rnatmp$gene1), as.character(rnatmp$gene2)),
					name=c(as.character(rnatmp$name), as.character(rnatmp$name)), split=c(as.numeric(rnatmp$split), as.numeric(rnatmp$split)),
					span=c(as.numeric(rnatmp$span), as.numeric(rnatmp$span)), tag=c(as.character(rnatmp$strand1), as.character(rnatmp$strand2)),
					partner_chr=c(as.character(rnatmp$chrom2), as.character(rnatmp$chrom1)), partner_start=c(as.numeric(rnatmp$pos2), as.numeric(rnatmp$pos1)),
					partner_end=c(as.numeric(rnatmp$pos2), as.numeric(rnatmp$pos1)), partner_gene_id=c(as.character(rnatmp$gene2), as.character(rnatmp$gene1)), stringsAsFactors=FALSE)
			#// Set 'rnabedgraph' as a bedgraph format (only four columns accepted) for IGV input
			rnabedgraph = data.table(chr=c(as.character(rnatmp$chrom1), as.character(rnatmp$chrom2)), start=c(as.numeric(rnatmp$pos1) - 1, as.numeric(rnatmp$pos2) - 1),
					end=c(as.numeric(rnatmp$pos1), as.numeric(rnatmp$pos2)), count=c(rep(1, length(rnatmp[,4])*2)), stringsAsFactors=FALSE)
			#// Count number of unique breakpoints for 'tumorbedgraph'
			rnabedgraph = rnabedgraph[, by=.(chr, start, end), .(value = sum(count))];
			return(list(rnabed=rnabed, rnabedgraph=rnabedgraph))
		})
		#///////////////////////////////////////////////////////////////////
		#// Create 'bedpe' list class for Linear module using RNA-seq SVs //
		#///////////////////////////////////////////////////////////////////
		#// 'input_twoway_rna_bedpe' has one element: 'rnabedpe' - a data.frame class for bedpe format 
		rnabedpe_tmp <- reactive({
			if ( is.null(inputFile$rnadata) ) { return(NULL); }
			tmp = inputFile$rnadata[inputFile$rnadata$split >= input$rna_split_bedpe & inputFile$rnadata$span >= input$rna_span_bedpe, ];
			#// remove translocs events / intra-events > a given distance
			dis = apply(tmp, 1, function(x){
				if ( as.character(x[1]) == as.character(x[4]) ) { #// intra-chromosome
					y=as.numeric(x[5]) - as.numeric(x[2]); abs(y) >= input$rna_dismin_bedpe && abs(y) <= input$rna_dismax_bedpe;
				} else {
					return(FALSE);
				}
			})
			tmp = cbind(tmp, dis);		tmp = tmp[tmp$dis == T,]; #// rnabedpe = rnabedpe[, c(1,2,3,4,5,6,7,8,9,10)];
			return(tmp);
		})
		observe({
			name_bedpe_rna = as.character(rnabedpe_tmp()$name);
			updateSelectizeInput(session = session, inputId = "rna_sample_bedpe", choices = unique(name_bedpe_rna), selected = "", options = list(placeholder = 'select'));
		})

		input_twoway_rna_bedpe <- reactive({
			#// Define data.frame class 'rnabedpe' for RNA SVs
			if ( is.null(rnabedpe_tmp()) ) { return(NULL); }
			rnabedpe = rnabedpe_tmp();
			#// choose selected 'name'(samples) for control
			if (! is.null(input$rna_sample_bedpe) ) {
				if ( input$rna_sample_bedpe[1] != "" ) {
					rnabedpe = rnabedpe[which(rnabedpe$name %in% input$rna_sample_bedpe), ];
				}
			}
			#// Set 'rnabedpe' as a bedpe format for IGV input
			rnabedpe = data.table(chrom1=as.character(rnabedpe$chrom1), start1=as.numeric(rnabedpe$pos1) - 1, end1=as.numeric(rnabedpe$pos1), 
					chrom2=as.character(rnabedpe$chrom2), start2=as.numeric(rnabedpe$pos2) - 1, end2=as.numeric(rnabedpe$pos2), name=as.character(rnabedpe$name), 
					type=rep("Intra", length(rnabedpe[,1])), split=as.numeric(rnabedpe$split), span=as.numeric(rnabedpe$span), stringsAsFactors=FALSE);
			#// Adding a new column 'label' by concatenating 'names' grouped by 'chrom1'-'start1'-'end1'-'chrom2'-'start2'-'end2'
			rnabedpe = rnabedpe[, by=.(chrom1, start1, end1, chrom2, start2, end2), label := paste(name, collapse= ',')];
			#// Create a new data.frame class after removing duplicate record 'chrom1'-'start1'-'end1'-'chrom2'-'start2'-'end2'-'label'-'type'
			rnabedpe = rnabedpe[, by=.(chrom1, start1, end1, chrom2, start2, end2, label, type), .(score = uniqueN(name))];
			names(rnabedpe)[7] = "name";	names(rnabedpe)[1] = "chr1";	names(rnabedpe)[4] = "chr2"; #// change 'label' to 'name', 'chrom1' to 'chr1' and change 'chrom1' to 'chr2'
			rnabedpe = rnabedpe[, c(1,2,3,4,5,6,7,9,8)]; #// chr1, start1, end1, chr2, start2, end2, name(mutiple samples), score (num of samples) and type
			return(rnabedpe);
		})

		#/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		#// Create 'bed' and 'bedgraph' list class for Linear module using DNA-seq SVs (ignore 'gene1' and 'gene2') //
		#/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		#// 'input_twoway_dna_bed' has two element: 'dnabed' - a data.frame class for bed format and 'dnabedgraph' - a data.frame class for bedgraph format
		dnabed_tmp <- reactive({
			if ( is.null(inputFile$dnadata) ) { return(NULL); }
			tmp = inputFile$dnadata[inputFile$dnadata$split >= input$dna_split_bed & inputFile$dnadata$span >= input$dna_span_bed, ];
			#// use distance of two SVs as control for filtering (translocs PASS)
			dis_bed = apply(tmp, 1, function(x){
				y=as.numeric(x[6]) - as.numeric(x[2]); 
				if ( as.character(x[8]) == "BND") { return(TRUE); } else { abs(y) >= input$dna_dismin_bed && abs(y) <= input$dna_dismax_bed; }
			})
			tmp = cbind(tmp, dis_bed);	tmp = tmp[tmp$dis_bed == T,]; #// tmp = tmp[, c(1,2,3,4,5,6,7,8,9,10)];
			return(tmp);
		})
		observe({
			type_bed_dna = as.character(dnabed_tmp()$type);
			updateSelectizeInput(session = session, inputId = "dna_type_bed", choices = unique(type_bed_dna), selected = "", options = list(placeholder = 'select'));
		})

		input_twoway_dna_bed <- reactive({
			#// Defined data.frame class 'dnabed' and 'dnabedgraph' for DNA SVs
			dnabed = NULL;        dnabedgraph = NULL;
			if ( is.null(dnabed_tmp()) ) { return(list(dnabed=dnabed, dnabedgraph=dnabedgraph)); }
			tmpbed = dnabed_tmp();
			if (! is.null(input$dna_type_bed) ) {
				if ( input$dna_type_bed[1] != "" ) {
					tmpbed = tmpbed[which(tmpbed$type %in% input$dna_type_bed), ];
				}
			}
			#// Set 'dnabed' as a bed format for IGV input
			dnabed = data.frame(chr=c(as.character(tmpbed$chrom1), as.character(tmpbed$chrom2)), start=c(as.numeric(tmpbed$start1) - 1, as.numeric(tmpbed$start2) - 1),
					end=c(as.numeric(tmpbed$end1), as.numeric(tmpbed$end2)), name=c(as.character(tmpbed$name), as.character(tmpbed$name)),
					split=c(as.numeric(tmpbed$split), as.numeric(tmpbed$split)), span=c(as.numeric(tmpbed$span), as.numeric(tmpbed$span)),
					type=c(as.character(tmpbed$type), as.character(tmpbed$type)), partner_chr=c(as.character(tmpbed$chrom2), as.character(tmpbed$chrom1)),
					partner_start=c(as.numeric(tmpbed$start2), as.numeric(tmpbed$start1)), partner_end=c(as.numeric(tmpbed$end2), as.numeric(tmpbed$end1)), stringsAsFactors=FALSE);
			#// Set 'dnabedgraph' as a bedgraph format (only four columns accepted) for IGV input
			dnabedgraph = data.table(chr=c(as.character(tmpbed$chrom1), as.character(tmpbed$chrom2)), start=c(as.numeric(tmpbed$start1) - 1, as.numeric(tmpbed$start2) - 1),
					end=c(as.numeric(tmpbed$end1), as.numeric(tmpbed$end2)), count=c(rep(1, length(tmpbed[,4])*2)), stringsAsFactors=FALSE);
			dnabedgraph = dnabedgraph[, by=.(chr, start, end), .(value = sum(count))];
			return(list(dnabed=dnabed, dnabedgraph=dnabedgraph));
		})

		#/////////////////////////////////////////////////////////////////////////////////////////////////////
		#// Create 'bedpe' list class for Linear module in bedpe using DNA-seq SVs (ignore gene1 and gene2) //
		#/////////////////////////////////////////////////////////////////////////////////////////////////////
		#// 'input_twoway_dna_bedpe' has one element: 'dnabedpe' - a data.frame class for bedpe format
		dnabedpe_tmp <- reactive({
			if ( is.null(inputFile$dnadata) ) { return(NULL); }
			tmp = inputFile$dnadata[inputFile$dnadata$split >= input$dna_split_bedpe & inputFile$dnadata$span >= input$dna_span_bedpe, ];
			tmp = tmp[tmp$type != 'BND', ]; # remove translocs due to no support for IGV
			#// use distance of two SVs as control for filtering (no translocs events included)
			dis = apply(tmp, 1, function(x){
				y=as.numeric(x[6]) - as.numeric(x[2]); 
				if ( as.character(x[8]) == "BND") { return(TRUE); } else { abs(y) >= input$dna_dismin_bedpe && abs(y) <= input$dna_dismax_bedpe; }
			})
			tmp = cbind(tmp, dis);	tmp = tmp[tmp$dis == T,];
			return(tmp);
		})
		observe({
			type_bedpe_dna = as.character(dnabedpe_tmp()$type);	name_bedpe_dna = as.character(dnabedpe_tmp()$name);
			updateSelectizeInput(session = session, inputId = "dna_type_bedpe", choices = unique(type_bedpe_dna), selected = "", options = list(placeholder = 'select'));
			updateSelectizeInput(session = session, inputId = "dna_sample_bedpe", choices = unique(name_bedpe_dna), selected = "", options = list(placeholder = 'select'));
		})

		input_twoway_dna_bedpe <- reactive({
			#// defined data.frame class 'dnabedpe' for DNA SVs
			if ( is.null(dnabedpe_tmp()) ) { return(NULL); }
			dnabedpe = dnabedpe_tmp();
			#// choose selected 'name'(samples) for control
			if (! is.null(input$dna_sample_bedpe) ) {
				if ( input$dna_sample_bedpe[1] != "" ) {
					dnabedpe = dnabedpe[which(dnabedpe$name %in% input$dna_sample_bedpe),];
				}
			}
			#// adjust genomic coordinate of start1 and start2
			dnabedpe$start1 = dnabedpe$start1 - 1;	dnabedpe$start2 = dnabedpe$start2 - 1;
			#// Adding a new column 'label' by concatenating 'names' grouped by 'chrom1'-'start1'-'end1'-'chrom2'-'start2'-'end2'
			dnabedpe = dnabedpe[, by=.(chrom1, start1, end1, chrom2, start2, end2), label := paste(name, collapse= ',')];
			#// Create a new data.frame class after removing duplicate record 'chrom1'-'start1'-'end1'-'chrom2'-'start2'-'end2'-'label'-'type'
			dnabedpe = dnabedpe[, by=.(chrom1, start1, end1, chrom2, start2, end2, label, type), .(score = uniqueN(name))];
			names(dnabedpe)[7] = "name";	names(dnabedpe)[1] = "chr1";	names(dnabedpe)[4] = "chr2"; #// change 'label' to 'name', 'chrom1' to 'chr1' and change 'chrom1' to 'chr2'
			dnabedpe = dnabedpe[, c(1,2,3,4,5,6,7,9,8)]; #// chr1, start1, end1, chr2, start2, end2, name(mutiple samples), score (num of samples) and type

			if (! is.null(input$dna_type_bedpe) ) {
				if ( input$dna_type_bedpe[1] != "" ) {
					dnabedpe = dnabedpe[which(dnabedpe$type %in% input$dna_type_bedpe), ];
				}
			}
			return(dnabedpe);
		})

		#//////////////////////////////////////////////////////////////////////////////////////////////
		#// Create 'segment' list class for Linear module using DNA-seq SVs (ignore gene1 and gene2) //
		#//////////////////////////////////////////////////////////////////////////////////////////////
		#// 'input_twoway_dna_seg' has one element: 'dnaseg' - a data.frame class for segment format 
		dnaseg_tmp <- reactive({
			if ( is.null(inputFile$dnadata) ) { return(NULL); }
			tmp = inputFile$dnadata[inputFile$dnadata$split >= input$dna_split_seg & inputFile$dnadata$span >= input$dna_span_seg, ];
			tmp = tmp[which(tmp$type %in% c('DUP','DEL')), ];
			return(tmp);
		})
		observe({
			name_seg_dna = as.character(dnaseg_tmp()$name);
			updateSelectizeInput(session = session, inputId = "dna_sample_seg", choices = unique(name_seg_dna), selected = "", options = list(placeholder = 'select'));
		})

		input_twoway_dna_seg <- reactive({
			#// Defined data.frame class 'seg' for DNA SVs
			if ( is.null(dnaseg_tmp()) ) { return(list(dup=NULL, del=NULL)); }
			dnaseg_del = NULL;		dnaseg_dup = dnaseg_tmp();
			dnaseg_del = dnaseg_dup[dnaseg_dup$type == "DEL", ];
			dnaseg_dup = dnaseg_dup[dnaseg_dup$type == "DUP", ];
			#// choose selected 'name'(samples) for control
			if (! is.null(input$dna_sample_seg) ) {
				if ( input$dna_sample_seg[1] != "" ) {
					dnaseg_dup = dnaseg_dup[which(dnaseg_dup$name %in% input$dna_sample_seg), ];
					dnaseg_del = dnaseg_del[which(dnaseg_del$name %in% input$dna_sample_seg), ];
				}
			}
			#// Set ampification in a bed format, adjust genomic coordinate of start1
			dnaseg_del = data.frame(chr=dnaseg_del$chrom1, start=dnaseg_del$start1-1, end=dnaseg_del$end2, value=rep(-1,length(dnaseg_del$chrom1)), 
									sample=dnaseg_del$name, stringsAsFactors=FALSE);
			#// Set deletion in a bed format, adjust genomic coordinate of start1
			dnaseg_dup = data.frame(chr=dnaseg_dup$chrom1, start=dnaseg_dup$start1-1, end=dnaseg_dup$end2, value=rep(1,length(dnaseg_dup$chrom1)), 
									sample=dnaseg_dup$name, stringsAsFactors=FALSE);
			if (! is.null(input$dna_seg_chr) ) {
				if ( grepl("chr[1-9XYMT]", input$dna_seg_chr) && input$dna_seg_start > 0 && input$dna_seg_end > 0) {
					if ( input$dna_seg_end > input$dna_seg_start ) {
						dnaseg_del = dnaseg_del[dnaseg_del$chr == input$dna_seg_chr,];
						dnaseg_dup = dnaseg_dup[dnaseg_dup$chr == input$dna_seg_chr,];
						#// judge whether seg interval overlaps with given region
						tag_del = apply(dnaseg_del, 1, function(x){ y<-c(sort(as.numeric(x[2:3]))); max(y[1], input$dna_seg_start) <= min(y[2], input$dna_seg_end)});
						tag_dup = apply(dnaseg_dup, 1, function(x){ y<-c(sort(as.numeric(x[2:3]))); max(y[1], input$dna_seg_start) <= min(y[2], input$dna_seg_end)});
						dnaseg_del = cbind(dnaseg_del, tag_del);	dnaseg_del = dnaseg_del[dnaseg_del$tag == T,];
						dnaseg_dup = cbind(dnaseg_dup, tag_dup);	dnaseg_dup = dnaseg_dup[dnaseg_dup$tag == T,];
						dnaseg_del = dnaseg_del[,c(1,2,3,4,5)];
						dnaseg_dup = dnaseg_dup[,c(1,2,3,4,5)];
					} else {
						shiny::showModal(modalDialog(title = "Error message", "Start coordinate > End coordinate!"));	req(NULL);
					}
				}
			}
			return(list(dup=dnaseg_dup, del=dnaseg_del));
		})

#-------------------------------------------------
# Circle module plot (for both RNA and DNA SVs)
#-------------------------------------------------
		#//----- For DNA SVs circular plot (only chr1..chrY)-----//#
		observe({
			if ( is.null(input$file_dna_data) ) { output$circle_2 <- NULL; } #// if no file loads, stop busy indicator bar
			if ( is.null(input$file_maf_data) ) { req(NULL); }
			sample_name = as.character(inputFile$mutationdata$Tumor_Sample_Barcode);
			updateSelectizeInput(session = session, inputId = "mut_sample_dna", choices = unique(sample_name), selected = "", options = list(placeholder = 'select'));
		})
		genome_select_dna = reactiveValues(value=NULL); #// use selected chromosomes for visualization
		circle_two_dna <- select_group_server(id = "circle_set_dna", data = circle_twocol_dna, vars = reactive(c("chr", "gene", "name")));
		circle_data_dna <- reactive({
			if ( is.null(inputFile$dnadata) ) { return(NULL); }
			#// select 'chr' (multiple choices accepted)
			tmp = inputFile$dnadata[inputFile$dnadata$chrom1 %in% unique(circle_two_dna()$chr) | inputFile$dnadata$chrom2 %in% unique(circle_two_dna()$chr), ];
			#// filter using selected 'gene' (multiple choices accepted)
			tmp = tmp[tmp$gene1 %in% unique(as.character(circle_two_dna()$gene)) | tmp$gene2 %in% unique(as.character(circle_two_dna()$gene)), ];
			tmp = tmp[tmp$name %in% unique(as.character(circle_two_dna()$name)), ]; #// select 'name' (multiple accepted)
			return(tmp);
		})
		circle_data_final_dna <- reactive({
			if ( is.null(circle_data_dna()) ) { return(NULL); }
			tmp = circle_data_dna()[circle_data_dna()$split >= input$circle_split_dna & circle_data_dna()$span >= input$circle_span_dna, ];
			#// Adding a new column 'label' by concatenating 'name' names grouped by 'chrom1'-'start1'-'end1'-'chrom2'-'start2'-'end2'-'gene1'-'gene2'
			tmp = tmp[, by=.(chrom1, start1, end1, chrom2, start2, end2, gene1, gene2), label := paste(name, collapse= ', ')];
			#// Create a new data.frame class after removing duplicate record 'chrom1'-'start1'-'end1'-'chrom2'-'start2'-'end2'-'gene1'-'gene2'-'label'
			tmp = tmp[, by=.(chrom1, start1, end1, chrom2, start2, end2, gene1, gene2, label), .(value = uniqueN(name))];
			#// if intra-chrom SV, only keep when distance > input$circle_dis_dna
			dis_tag = apply(tmp, 1, function(x){
				if ( as.character(x[1]) == as.character(x[4]) ) { #// intra-chromosome
					y=as.numeric(x[6]) - as.numeric(x[2]);
					abs(y) >= input$circle_dis_dna; #// more -> TRUE; less -> FALSE
				} else {
					return(TRUE);
				}
			})
			tmp = cbind(tmp, dis_tag);	tmp = tmp[tmp$dis_tag == T,];	tmp$dis_tag = NULL;
			tmp = tmp[tmp$value >= input$circle_num_dna, ]; #// 'input$circle_num_dna' (sliding bar)
			return(tmp);
		})
		#// create gene tracks for plotting (only genes related to SV events were selected)
		circle_gene_dna <- reactive({
			tmp = database$gene_range;
			tmp = tmp[tmp$gene %in% unique(c(circle_data_dna()$gene1, circle_data_dna()$gene2)), ];
			return(tmp);
		})
		#// create cytoband tracks for plotting (including all 24 chromosomes)
		circle_cytoband_dna <- reactive({
			tmp = database$karyto;
			tmp = tmp[tmp$chr %in% unique(c(circle_data_dna()$chrom1, circle_data_dna()$chrom2)), ];
			return(tmp);
		})
		#// create mutation tracks for plotting (SNVs+indels)
		circle_mutation_dna <- reactive({
			if ( is.null(mutation_bedgraph$value) ) { return(NULL); }
			#// if any of chr, sample were selected
			tmp = mutation_bedgraph$value;
			if ( is.null(input[['circle_set_dna-chr']]) && is.null(input[['circle_set_dna-gene']]) && is.null(input[['circle_set_dna-sample']]) ) {
			} else { #// if any of 'circle_set_dna-chr', 'circle_set_dna-gene' and 'circle_set_dna-sample' are selected
				tmp = tmp[tmp$Chromosome %in% unique(c(circle_data_dna()$chrom1, circle_data_dna()$chrom2)), ]; #// filter using chromosome
			}

			if (! is.null(input$mut_type_dna) ) {
				if ( input$mut_type_dna[1] != "" ) {
					reg_exp = paste(input$mut_type_dna, collapse="|");
					tmp = tmp[anno_new %like% reg_exp];
				}
			}
			if (! is.null(input$mut_sample_dna) ) {
				if ( input$mut_sample_dna[1] != "" ) {
					reg_exp = paste(input$mut_sample_dna, collapse="|");
					tmp = tmp[anno_new %like% reg_exp];
				}
			}
			return(tmp)
		})
		
		observeEvent(input$circos_plot_dna, {
			shinycssloaders::showSpinner("circle_2")
			circle_plot_object=NULL;
			if ( is.null(circle_data_final_dna()) || nrow(circle_data_final_dna()) == 0 || is.null(circle_gene_dna()) || is.null(circle_cytoband_dna()) ) {
				shiny::showModal(modalDialog(title = "Warning message", "No DNA SVs are available!"));
			} else {
				#// set SVs tracks
				chrom_A = circle_data_final_dna()$chrom1;	chrom_B = circle_data_final_dna()$chrom2;
				pos_A_s = circle_data_final_dna()$start1;	pos_A_e = circle_data_final_dna()$end1;
				pos_B_s = circle_data_final_dna()$start2;   pos_B_e = circle_data_final_dna()$end2;
				gene1 = circle_data_final_dna()$gene1;		gene2 = circle_data_final_dna()$gene2;
				labels = paste(chrom_A, ":", pos_A_s, "-", pos_A_e, ":", gene1, " | ", chrom_B, ":", pos_B_s, "-", pos_B_e, ":", gene2, 
						   " (", circle_data_final_dna()$value, ")", " [",circle_data_final_dna()$label, "]", sep="");
				#// assemble tracklist for plotting
				tracklist = BioCircos::BioCircosBackgroundTrack("BackgroundTrack", minRadius = 0, maxRadius = 0.73, borderSize = 0, fillColors = "#EEFFEE");
				tracklist = tracklist + BioCircos::BioCircosBackgroundTrack("BackgroundTrack", minRadius = 0.75, maxRadius = 0.95, borderSize = 0, fillColors = "#FFFFFF");
				tracklist = tracklist + BioCircos::BioCircosArcTrack("GeneTrack", minRadius = 1.3, maxRadius = 1.4, chromosomes = circle_gene_dna()$chr, starts = circle_gene_dna()$start,
									ends = circle_gene_dna()$end, labels = circle_gene_dna()$gene, colors = "red");
				tracklist = tracklist + BioCircos::BioCircosLinkTrack('LinkTrack', chrom_A, pos_A_s, pos_A_e, chrom_B, pos_B_s, pos_B_e, maxRadius=0.73, width="0.15em", 
									labels=labels, displayLabel=F, color="red");
				#// checkbox is true
				if ( input$cirmutcheck_dna ) {
					if (! is.null(circle_mutation_dna()) ) {
						tracklist = tracklist + BioCircos::BioCircosSNPTrack('MutationTrack', circle_mutation_dna()$Chromosome, circle_mutation_dna()$Start_Position, circle_mutation_dna()$freq, 
							labels = circle_mutation_dna()$anno_new, size = 1.4, maxRadius = 0.95, minRadius = 0.75, colors = "orange");
					} else {
						shiny::showModal(modalDialog(title = "Warning message", "No mutation data in MAF format is loaded!"));
					}
				}
				#// if no selection, all chroms for visual
				if ( is.null(input[['circle_set_dna-chr']]) && is.null(input[['circle_set_dna-gene']]) && is.null(input[['circle_set_dna-sample']]) ) { 
					genome_select_dna$value = database$genome_cir;
					tracklist = tracklist + BioCircos::BioCircosArcTrack("CytobanTrack", minRadius=1, maxRadius=1.15, chromosomes=database$karyto$chr, starts=database$karyto$start, 
												ends=database$karyto$end, labels=database$karyto$name, colors=database$karyto$color);
				} else {
					genome_select_dna$value = database$genome_cir[unique(circle_cytoband_dna()$chr)];
					tracklist = tracklist + BioCircos::BioCircosArcTrack("CytobanTrack", minRadius=1, maxRadius=1.15, chromosomes=circle_cytoband_dna()$chr, starts=circle_cytoband_dna()$start,
						ends=circle_cytoband_dna()$end, labels=circle_cytoband_dna()$name, colors=circle_cytoband_dna()$color);
				}
				circle_plot_object <- BioCircos::BioCircos(tracklist, genome=genome_select_dna$value, genomeFillColor=rep("white", length(genome_select_dna$value)), displayGenomeBorder=T,
					SNPMouseOverTooltipsBorderWidth = "0.5px", SNPMouseOverCircleSize = 1, SNPMouseOverCircleOpacity = 0.5, 
					ARCMouseOverTooltipsBorderWidth = "0.5px", ARCMouseOverArcOpacity = 0.5, 
					LINKMouseOverTooltipsBorderWidth = "0.5px", LINKMouseOverOpacity = 0.5, LINKMouseOverStrokeWidth = 0.05,
					genomeBorderSize=0.4, genomeBorderColor="black", genomeLabelTextSize="7pt", chrPad=0.02, genomeTicksLen=2, 
					genomeTicksTextSize="0.5em", SNPMouseOverTooltipsHtml04="<br/>Anno: ", LINKMouseOverTooltipsHtml01="Link: ", 
					ARCMouseOverTooltipsHtml04="<br/>Name: ", TEXTModuleDragEvent=T);
			}
			shinycssloaders::hideSpinner("circle_2");
			output$circle_2 <- BioCircos::renderBioCircos({ circle_plot_object });
			output$CircleDownDNA <- downloadHandler(filename=function(){ paste("circle_plot_DNA", "html", sep =".") }, content=function(file){ saveWidget(circle_plot_object, file) });
		})

		#//----- For RNA (only for chr1..Y)-----//#
		observe({
			if ( is.null(input$file_rna_data) ) { output$circle_1 <- NULL; } #// if no file loads, stop busy indicator bar
			if ( is.null(input$file_maf_data) ) { req(NULL); }
			sample_name = as.character(inputFile$mutationdata$Tumor_Sample_Barcode);
			updateSelectizeInput(session = session, inputId = "mut_sample", choices = unique(sample_name), selected = "", options = list(placeholder = 'select'));
		})
		genome_select = reactiveValues(value=NULL); #// use selected chromosomes for visualization
		circle_two <- select_group_server(id = "circle_set", data = circle_twocol, vars = reactive(c("chr", "gene", "name")));
		circle_data <- reactive({
			if ( is.null(inputFile$rnadata) ) { return(NULL); }
			#// filter using selected 'chr' (multiple accepted)
			tmp = as.data.table(inputFile$rnadata[inputFile$rnadata$chrom1 %in% unique(circle_two()$chr) | inputFile$rnadata$chrom2 %in% unique(circle_two()$chr), ]); 
			#// filter using selected 'gene' (multiple accepted)
			tmp = tmp[tmp$gene1 %in% unique(as.character(circle_two()$gene)) | tmp$gene2 %in% unique(as.character(circle_two()$gene)), ]; 
			tmp = tmp[tmp$name %in% unique(as.character(circle_two()$name)), ]; #// filter using selected 'name' (multiple accepted)
			return(tmp);
		})
		circle_data_final <- reactive({
			if ( is.null(circle_data()) ) { return(NULL); }
			#// filter using 'input$circle_split' and 'input$circle_span' (sliding bar)
			tmp = circle_data()[circle_data()$split >= input$circle_split & circle_data()$span >= input$circle_span, ];
			#// Adding a new column 'label' by concatenating 'name' names grouped by 'chrom1'-'pos1'-'gene1'-'chrom2'-'pos2'-'gene2'
			tmp = tmp[, by=.(chrom1, pos1, gene1, chrom2, pos2, gene2), label := paste(name, collapse= ', ')];
			#// Create a new data.frame class after removing duplicate record 'chrom1'-'pos1'-'gene1'-'chrom2'-'pos2'-'gene2'-'label'
			tmp = tmp[, by=.(chrom1, pos1, gene1, chrom2, pos2, gene2, label), .(value = uniqueN(name))];
			tmp = tmp[tmp$value >= input$circle_num, ]; #// filter using 'input$circle_num' (sliding bar)
			return(tmp);
		})
		#// create gene tracks for plotting (only genes related to SV events were selected)
		circle_gene <- reactive({
			tmp = database$gene_range;
			tmp = tmp[tmp$gene %in% unique(c(circle_data()$gene1, circle_data()$gene2)), ];
			return(tmp);
		})
		#// create cytoband tracks for plotting (including all 24 chromosomes)
		circle_cytoband <- reactive({
			tmp = database$karyto;
			tmp = tmp[tmp$chr %in% unique(c(circle_data()$chrom1, circle_data()$chrom2)), ];
			return(tmp);
		})
		#// create mutation tracks for plotting (SNVs+indels)
		circle_mutation <- reactive({
			#// if any of chr, sample were selected
			tmp = mutation_bedgraph$value;
			if ( is.null(tmp) ) { return(NULL); }
			if ( is.null(input[['circle_set-chr']]) && is.null(input[['circle_set-gene']]) && is.null(input[['circle_set-sample']]) ) { 
			} else { #// if any of 'circle_set-chr', 'circle_set-gene' and 'circle_set-sample' are selected
				tmp = tmp[tmp$Chromosome %in% unique(c(circle_data()$chrom1, circle_data()$chrom2)), ]; #// filter using chromosome
			}
			if (! is.null(input$mut_type) ) {
				if ( input$mut_type[1] != "" ) {
					reg_exp = paste(input$mut_type, collapse="|");
					tmp = tmp[anno_new %like% reg_exp];
				}
			}
			if (! is.null(input$mut_sample) ) {
				if ( input$mut_sample[1] != "" ) {
					reg_exp = paste(input$mut_sample, collapse="|");
					tmp = tmp[anno_new %like% reg_exp];
				}
			}
			return(tmp);
		})

		observeEvent(input$circos_plot, {
			shinycssloaders::showSpinner("circle_1")
			circle_plot_object=NULL;
			if ( is.null(circle_data_final()) || nrow(circle_data_final()) == 0 || is.null(circle_gene()) || is.null(circle_cytoband()) ) { 
				shiny::showModal(modalDialog(title = "Warning message", "No RNA SVs are available!"));
			} else {
			#// set SVs tracks
				chrom_A = circle_data_final()$chrom1;	chrom_B = circle_data_final()$chrom2;
				pos_A_s = circle_data_final()$pos1;		pos_A_e = circle_data_final()$pos1 + 1;
				pos_B_s = circle_data_final()$pos2;		pos_B_e = circle_data_final()$pos2 + 1;
				gene1 = circle_data_final()$gene1;		gene2 = circle_data_final()$gene2;
				labels = paste(chrom_A, ":", pos_A_s, ":", gene1, " | ", chrom_B, ":", pos_B_s, ":", gene2, 
					" (", circle_data_final()$value, ")", " [",circle_data_final()$label, "]", sep="");
				#// assemble tracklist for plotting
				tracklist = BioCircos::BioCircosBackgroundTrack("BackgroundTrack", minRadius = 0, maxRadius = 0.73, borderSize = 0, fillColors = "#EEFFEE");
				tracklist = tracklist + BioCircos::BioCircosBackgroundTrack("BackgroundTrack", minRadius = 0.75, maxRadius = 0.95, borderSize = 0, fillColors = "#FFFFFF");
				tracklist = tracklist + BioCircos::BioCircosArcTrack("GeneTrack", minRadius = 1.3, maxRadius = 1.4, chromosomes = circle_gene()$chr, starts = circle_gene()$start,
									ends = circle_gene()$end, labels = circle_gene()$gene, colors = "red");
				tracklist = tracklist + BioCircos::BioCircosLinkTrack('LinkTrack', chrom_A, pos_A_s, pos_A_e, chrom_B, pos_B_s, pos_B_e, maxRadius=0.73, width="0.15em",
									labels=labels, displayLabel=F, color="red");
				#// checkbox is true
				if ( input$cirmutcheck ) {
					if (! is.null(circle_mutation()) ) {
						tracklist = tracklist + BioCircos::BioCircosSNPTrack('MutationTrack', circle_mutation()$Chromosome, circle_mutation()$Start_Position, circle_mutation()$freq, 
							labels = circle_mutation()$anno_new, size = 1.4, maxRadius = 0.95, minRadius = 0.75, colors = "orange");
					} else {
						shiny::showModal(modalDialog(title = "Warning message", "No mutation data in MAF format is loaded!"));
					}
				}
				#// if no selection, use all chroms for visualization
				if ( is.null(input[['circle_set-chr']]) && is.null(input[['circle_set-gene']]) && is.null(input[['circle_set-sample']]) ) { 
					genome_select$value = database$genome_cir;
					tracklist = tracklist + BioCircos::BioCircosArcTrack("CytobanTrack", minRadius=1, maxRadius=1.15, chromosomes=database$karyto$chr, starts=database$karyto$start, 
							ends=database$karyto$end, labels=database$karyto$name, colors=database$karyto$color);
				} else { #// use selected chromosomes for visualization
					genome_select$value = database$genome_cir[unique(circle_cytoband()$chr)];
					tracklist = tracklist + BioCircos::BioCircosArcTrack("CytobanTrack", minRadius=1, maxRadius=1.15, chromosomes=circle_cytoband()$chr, starts=circle_cytoband()$start,
							ends=circle_cytoband()$end, labels=circle_cytoband()$name, colors=circle_cytoband()$color);
				}
				circle_plot_object <- BioCircos::BioCircos(tracklist, genome=genome_select$value, genomeFillColor=rep("white", length(genome_select$value)), displayGenomeBorder=T,
					SNPMouseOverTooltipsBorderWidth = "0.5px", SNPMouseOverCircleSize = 1, SNPMouseOverCircleOpacity = 0.5,
					ARCMouseOverTooltipsBorderWidth = "0.5px", ARCMouseOverArcOpacity = 0.5,
					LINKMouseOverTooltipsBorderWidth = "0.5px", LINKMouseOverOpacity = 0.5, LINKMouseOverStrokeWidth = 0.05,
					genomeBorderSize=0.4, genomeBorderColor="black", genomeLabelTextSize="7pt", chrPad=0.02, genomeTicksLen=2, 
					genomeTicksTextSize="0.5em", SNPMouseOverTooltipsHtml04="<br/>Anno: ", LINKMouseOverTooltipsHtml01="Link: ", 
					ARCMouseOverTooltipsHtml04="<br/>Name: ", TEXTModuleDragEvent=T);
			}
			shinycssloaders::hideSpinner("circle_1");
			output$circle_1 <- BioCircos::renderBioCircos({ circle_plot_object });
			output$CircleDownRNA <- downloadHandler(filename=function(){ paste("circle_plot_RNA", "html", sep =".") }, content=function(file){ saveWidget(circle_plot_object, file) });
		})

#-------------------------------------------
# Two-way module plot (only for RNA SVs)
#-------------------------------------------
	#------ For fusion overview -------#
	overview_size_up <- reactive({ return(20*input$overview_size); })
	overview_size_down <- reactive({ return(480*input$overview_size); })
	#// 'object_over_A$value' and 'object_over_B$value' are list classes with five elements: 'select_region_f' is a data.frame class, 'dataset' is a data.frame class
	object_over_A <- reactiveValues(value=NULL)
	object_over_B <- reactiveValues(value=NULL) 
	#// 'plotbody' is list class with two elements: A1_xy - exon coordinates of geneA, B1_xy - exon coordinates of geneB)
	plotbody <- reactiveValues(collect=NULL)
	overviewAB <- reactiveValues(geneA=NULL, symbol_A=NULL, geneB=NULL, symbol_B=NULL)
    status_check <- reactiveValues(value=0) #// used for plotting
    filter_overview <- reactiveValues(value=NULL)
    overview_data_static <- reactiveValues(value=NULL)

	overview_tmp <- select_group_server(id = "overview1", data = reactive(inputFile$rnadata), vars = reactive(c("gene1", "gene2")))
    observe({
        if ( is.null(inputFile$rnadata) ) { req(NULL); }
        filter_origin = overview_tmp()[overview_tmp()$split >= input$overview_split & overview_tmp()$span >= input$overview_span, ]
		filter_origin = data.table(breakpoint_A=as.character(filter_origin$pos1), breakpoint_B=as.character(filter_origin$pos2), 
					count=rep(1, length(filter_origin$pos1)), stringsAsFactors=FALSE)
		filter_overview$value = as.data.frame(filter_origin[, by=.(breakpoint_A, breakpoint_B), .(value = sum(count))]);
    })
    overview_data <- select_group_server(id = "overview2", data = reactive(filter_overview$value), vars = reactive(c("breakpoint_A", "breakpoint_B")))

    #// when input[['overview1-gene1']] NOT null, assign 'object_over_A$value' using FUNCTION "get_annotation_db_extend"
	observe({
        if ( is.null(input[['overview1-gene1']]) || input[['overview1-gene1']] == "" ) {
            object_over_A$value = NULL
            updateSelectizeInput(session = session, inputId = "transA_overview", choices = "")
        } else {
			if ( is.null(database$txdb) || is.null(database$grTrack) || is.null(database$whole_txdb) ) { req(NULL); }
            status_check$value = 0
			if ( length(unique(as.character(overview_tmp()$gene1))) == 1 ) {
				choice_geneA = unique(as.character(overview_tmp()$gene1))
				ens_A = NULL;
				if ( grepl("ENSG00", choice_geneA) || grepl("ENSMUSG0", choice_geneA) ) {
					ens_A = database$ensembl_id[database$ensembl_id == choice_geneA];
				} else {
					ens_A = names(database$symbol_ensem[database$symbol_ensem == choice_geneA])
				}
				if ( length(ens_A) == 1 ) { #// if gene symbol matches to ensembl_id
					object_over_A$value <- FuSViz::get_annotation_db_extend(ens_A, database$txdb, database$grTrack, database$whole_txdb)
					if ( is.null(object_over_A$value) ) {
						updateSelectizeInput(session = session, inputId = "transA_overview", choices = "") #// NOTE: set choices = "" not (NULL)
						showModal(modalDialog(title = "Warning message", paste("GeneA(", choice_geneA, ") is not available in annotation database, and process stops here!", sep="")));	req(NULL)
					} else {
						#// update the choice values in selectInput
                        canonical_transA = database$canonical[database$canonical$ensembl_gene_id == ens_A,]$ensembl_transcript_id
                        selectinput_A = data.frame(TXNAME=object_over_A$value$dataset$TXNAME, TAG=object_over_A$value$dataset$TXNAME, stringsAsFactors=FALSE)
                        if (! is.na(canonical_transA[1]) ) { #// transcripts with canonical are tagged
                            if ( length(selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG) > 0 ) {
                                selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG = paste("@", selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG, sep="");
                            }
                        }
                        selectinput_A = as.list(setNames(selectinput_A$TAG, selectinput_A$TXNAME)); #// key is TAG; value is TXNAME
                        updateSelectizeInput(session = session, inputId = "transA_overview", choices = selectinput_A, selected = NULL,
						options = list(placeholder = 'Ensembl transcript id', render = I('{ option: function(item, escape) {
								const regExpStr = "@";
								return RegExp(regExpStr, "g").test(item.value) ? "<div><u>" + escape(item.label) + "</u></div>" : "<div>" + escape(item.label) + "</div>"
						}}')))
					}
				} else {
					object_over_A$value = NULL
					updateSelectizeInput(session = session, inputId = "transA_overview", choices = "") #// NOTE: set choices = "" not (NULL)
					showModal(modalDialog(title = "Warning message", paste("No or multiple geneA(", choice_geneA, ") is present (check gene symbol!)", sep="")));	req(NULL);
				}
			}
        }
	})
    #// when input[['overview1-gene2']] NOT null, assign 'object_over_B$value' using FUNCTION "get_annotation_db_extend"
	observe({
        if ( is.null(input[['overview1-gene2']]) || input[['overview1-gene2']] == "" ) {
            object_over_B$value = NULL
            updateSelectizeInput(session = session, inputId = "transB_overview", choices = "")
        } else {
			if ( is.null(database$txdb) || is.null(database$grTrack) || is.null(database$whole_txdb) ) { req(NULL); }
            status_check$value = 0
			if ( length(unique(as.character(overview_tmp()$gene2))) == 1 ) {
				choice_geneB = unique(as.character(overview_tmp()$gene2))
				ens_B = NULL;
				if ( grepl("ENSG00", choice_geneB) || grepl("ENSMUSG0", choice_geneB) ) {
					ens_B = database$ensembl_id[database$ensembl_id == choice_geneB];
				} else {
					ens_B = names(database$symbol_ensem[database$symbol_ensem == choice_geneB])
				}
				if ( length(ens_B) == 1 ) {
					object_over_B$value <- FuSViz::get_annotation_db_extend(ens_B, database$txdb, database$grTrack, database$whole_txdb)
					if ( is.null(object_over_B$value) ) {
						updateSelectizeInput(session = session, inputId = "transB_overview", choices = "") #// NOTE: set choices = "" not (NULL)
						showModal(modalDialog(title = "Warning message", paste("GeneB(", choice_geneB, ") is not available in annotation database, and process stops here!", sep=""))); req(NULL);
					} else {
						#// update the choice values in selectInput
                        canonical_transB = database$canonical[database$canonical$ensembl_gene_id == ens_B,]$ensembl_transcript_id
						selectinput_B = data.frame(TXNAME=object_over_B$value$dataset$TXNAME, TAG=object_over_B$value$dataset$TXNAME, stringsAsFactors=FALSE)
                        if (! is.na(canonical_transB[1]) ) { #// transcripts with canonical are tagged
                            if ( length(selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG) > 0 ) {
                                selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG = paste("@", selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG, sep="");
                            }
                        }
                        selectinput_B = as.list(setNames(selectinput_B$TAG, selectinput_B$TXNAME)); #// key is TAG; value is TXNAME
                        updateSelectizeInput(session = session, inputId = "transB_overview", choices = selectinput_B, selected = NULL,
						options = list(placeholder = 'Ensembl transcript id', render = I('{ option: function(item, escape) {
								const regExpStr = "@";
								return RegExp(regExpStr, "g").test(item.value) ? "<div><u>" + escape(item.label) + "</u></div>" : "<div>" + escape(item.label) + "</div>"
                        }}')))
					}
				} else {
					object_over_B$value = NULL
					updateSelectizeInput(session = session, inputId = "transB_overview", choices = "") #// NOTE: set choices = "" not (NULL)
					showModal(modalDialog(title = "Warning message", paste("No or multiple geneB(", choice_geneB, ") is present (check gene symbol!)", sep="")));	req(NULL);
				}
			}
        }
	})
    #// click button 'input$overview_on' for plotting
	observeEvent(input$overview_on, {
		shinycssloaders::showSpinner("chimerics_down");
		overviewAB$geneA=NULL;	overviewAB$symbol_A=NULL;	overviewAB$geneB=NULL;	overviewAB$symbol_B=NULL;
		if ( is.null(database$whole_txdb) || is.null(database$chrTrack) ) { showModal(modalDialog(title = "Error message", "Genome annotation database is not loaded!")); }
		if ( !is.null(input[['overview1-gene1']]) && !is.null(input[['overview1-gene2']]) ) {
			if ( !is.null(object_over_A$value) && !is.null(object_over_B$value) ) { #// make sure gene symbols match to valid ensembl_id for geneA and geneB
				symbol_A = overview_tmp()[1,]$gene1
				symbol_B = overview_tmp()[1,]$gene2
				breakpoint_set = unique(data.frame(breakpoint_A=overview_tmp()$pos1, breakpoint_B=overview_tmp()$pos2, stringsAsFactors=FALSE))
				#// a subset of 'object_over_A' as input for FUNCTION gene_trans_ex_overview
				tmp_control_transA = list();
				tmp_control_transA$value$chr = object_over_A$value$chr
				tmp_control_transA$value$start = object_over_A$value$start
				tmp_control_transA$value$end = object_over_A$value$end
				tmp_control_transA$value$strand = object_over_A$value$strand
				if (! is.null(input$transA_overview) ) { #// if selected transcript of geneA present
                    name_transA = input$transA_overview;	name_transA = gsub('@', '', name_transA);
					tmp_control_transA$value$dataset = object_over_A$value$dataset[object_over_A$value$dataset$TXNAME %in% name_transA, ][1,] #// select the first row represent geneA
					tmp_control_transA$value$txTr_f = object_over_A$value$txTr_f[object_over_A$value$txTr_f$transcript %in% name_transA, ]
				} else {
					tmp_control_transA$value$dataset = object_over_A$value$dataset[1,] #// select the first row represent geneA
					tmp_control_transA$value$txTr_f = object_over_A$value$txTr_f
				}
				#// a subset of 'object_over_B' as input for FUNCTION gene_trans_ex_overview
				tmp_control_transB = list(); 
				tmp_control_transB$value$chr = object_over_B$value$chr
				tmp_control_transB$value$start = object_over_B$value$start
				tmp_control_transB$value$end = object_over_B$value$end
				tmp_control_transB$value$strand = object_over_B$value$strand
				if (! is.null(input$transB_overview) ) { #// if selected transcript of geneB present
                    name_transB = input$transB_overview;	name_transB = gsub('@', '', name_transB);
					tmp_control_transB$value$dataset = object_over_B$value$dataset[object_over_B$value$dataset$TXNAME %in% name_transB, ][1,] #// select the first row represent geneB
					tmp_control_transB$value$txTr_f = object_over_B$value$txTr_f[object_over_B$value$txTr_f$transcript %in% name_transB, ]
				} else {
					tmp_control_transB$value$dataset = object_over_B$value$dataset[1,] #// select the first row represent geneB
					tmp_control_transB$value$txTr_f = object_over_B$value$txTr_f
				}
				
				#// annotate and evaluate breakpoints of geneA and geneB
				geneA = list(); geneB = list(); #// geneA and geneB are list class
				for (i in 1:length(breakpoint_set[,1])) {
					print(paste("start overview loop: ", i))
					if (! exists(as.character(breakpoint_set[i,1]), where=geneA) ) {
						geneA[[as.character(breakpoint_set[i,1])]] <- FuSViz::gene_trans_ex_overview(breakpoint_set[i,1], tmp_control_transA$value, database$whole_txdb) #// For geneA
					}
					if (! exists(as.character(breakpoint_set[i,2]), where=geneB) ) {
						geneB[[as.character(breakpoint_set[i,2])]] <- FuSViz::gene_trans_ex_overview(breakpoint_set[i,2], tmp_control_transB$value, database$whole_txdb) #// For geneB
					}
				}
				if ( length(geneA) > 0 ) { 
					geneA[sapply(geneA, is.null)] <- NULL; #// remove element with NULL in the list geneA (if all elements are NULL, geneA = NULL)
				} else {
					showModal(modalDialog(title = "Warning message", "No breakpoints within geneA or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!"));
				}
				if ( length(geneB) > 0 ) { 
					geneB[sapply(geneB, is.null)] <- NULL; #// remove element with NULL in the list geneB (if all elements are NULL, geneB = NULL)
				} else {
					showModal(modalDialog(title = "Warning message", "No breakpoints within geneB or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!)"));
				}
				overviewAB$geneA=geneA; overviewAB$symbol_A=symbol_A;
				overviewAB$geneB=geneB; overviewAB$symbol_B=symbol_B;
                overview_data_static$value = overview_data();	# set fusion breakpoint in static status
                status_check$value = 1
			} else {
				showModal(modalDialog(title = "Warning message", paste("No or multiple geneA(", choice_geneA, ") / geneB(", choice_geneB, ") is present (check gene symbol!)", sep="")));
				shinycssloaders::hideSpinner("chimerics_down")
			}
		} else {
			showModal(modalDialog(title = "Warning message", "No gene symbols are chosen for fusion partners (choose gene symbol first!)"));
			shinycssloaders::hideSpinner("chimerics_down")
		}
	})
    
    #// start overview plotting
	output$chimerics_down <- renderPlot({
		if ( length(overviewAB$geneA) > 0 && length(overviewAB$geneB) > 0 ) {
			print("start to plotting down overview fusion")
			reactive_control = overview_size_down();
			plotbody$collect = tryCatch(FuSViz::plot_separate_overview(overviewAB$geneA, overviewAB$symbol_A, overviewAB$geneB, overviewAB$symbol_B, database$chrTrack), error = function(e) {});
			if ( is.null(plotbody$collect) ) {
				showModal(modalDialog(title = "Warning message", paste("No enough space to draw all ", overviewAB$symbol_A, " / ", overviewAB$symbol_B, " transcript isoforms in the current layout device.", " Please increase the layout device size by pulling 'Zoom in/out' slide bar!", sep="")));
			}
		} else {
			return(NULL);
		}
		shinycssloaders::hideSpinner("chimerics_down")
	}, height = overview_size_down, width = overview_size_down)

	observe({
		if ( !is.null(object_over_A$value) && !is.null(object_over_B$value) ) {
			print("start to plotting up overview fusion")
			output$chimerics_up <- renderPlot({
				if ( length(overviewAB$geneA) > 0 && length(overviewAB$geneB) > 0 && length(plotbody$collect) > 0 ) {
					if ( status_check$value == 0 ) {
						if ( length(rownames(overview_data_static$value)) > 0 ) {
							reactive_control = overview_size_up();
							FuSViz::plot_arrow_overview(plotbody$collect, overviewAB$geneA, overviewAB$geneB, overview_data_static$value, reactive_control)
						}
					} else {
						if ( length(rownames(overview_data())) > 0 ) {
							reactive_control = overview_size_up();
							FuSViz::plot_arrow_overview(plotbody$collect, overviewAB$geneA, overviewAB$geneB, overview_data(), reactive_control)
						}
					}
				}
			}, height = overview_size_up, width = overview_size_down)
		} else {
			output$chimerics_up <- renderPlot({
				if ( length(overviewAB$geneA) > 0 && length(overviewAB$geneB) > 0 && length(plotbody$collect) > 0 ) {
					print("start to plotting up overview fusion static")
					if ( length(rownames(overview_data_static$value)) > 0 ) {
						reactive_control = overview_size_up();
						FuSViz::plot_arrow_overview(plotbody$collect, overviewAB$geneA, overviewAB$geneB, overview_data_static$value, reactive_control)
					}
				}
			}, height = overview_size_up, width = overview_size_down)
		}
	})

	#// download overview plot
	output$FusionDown1 <- downloadHandler(
	filename = function(){ paste("Two-way_overview_plot", tolower(input$file_fusion1), sep =".") },
		content = function(file) {
			pixelratio = 2
			if( input$file_fusion1 == "PNG" ) {
				png(file, width = input$overview_width, height = input$overview_height, res=72*pixelratio, units = "in")
			} else {
				pdf(file, width = input$overview_width, height = input$overview_height)
			}
			if ( length(overviewAB$geneA) > 0 && length(overviewAB$geneB) > 0 && length(rownames(overview_data_static$value)) > 0 ) {
				overview_log = tryCatch(FuSViz::plot_separate_overview_download(overviewAB$geneA, overviewAB$symbol_A, overviewAB$geneB, overviewAB$symbol_B, database$chrTrack, overview_data()), error = function(e) {});
				if ( is.null(overview_log) ) {
					showModal(modalDialog(title = "Warning message", paste("Download an incomplete plot due to a small layout size. ", "Please increase height and width values in 'Layout_width' and 'Layout_height' boxes to get a full plot!", sep="")));
				}
			} else {
				plot.new()
			}
			dev.off()
		}
	)

    #------ For individual visualisation -----#
	persample_size_full <- reactive({ return(500*input$persample_size); })
	object_individual_A <- reactiveValues(value=NULL) #// 'object_individual_A$value' is list class with five elements:
	object_individual_B <- reactiveValues(value=NULL) #// 'object_individual_B$value' is list class with five elements:
	individualAB <- reactiveValues(geneA=NULL, symbol_A=NULL, geneB=NULL, symbol_B=NULL, breakpoint_set=NULL)
	individual_data <- select_group_server(id = "individual", data = reactive(inputFile$rnadata), vars = reactive(c("gene1", "gene2", "pos1", "pos2", "name", "strand1", "strand2")))
    #// when input[['individual-gene1']] NOT null, assign 'object_individual_A$value' using FUNCTION "get_annotation_db"
    observeEvent(input[['individual-gene1']], {
        if ( is.null(input[['individual-gene1']]) || input[['individual-gene1']] == "" ) {
            object_individual_A = NULL
            updateSelectizeInput(session = session, inputId = "transA_individual", choices = "")
        } else {
			if ( is.null(database$txdb) || is.null(database$grTrack) ) { req(NULL);}
			if ( length(unique(as.character(individual_data()$gene1))) == 1 ) {
				choice_geneA = unique(as.character(individual_data()$gene1))
				ens_A = NULL;
				if ( grepl("ENSG00", choice_geneA) ||  grepl("ENSMUSG0", choice_geneA) ) {
					ens_A = database$ensembl_id[database$ensembl_id == choice_geneA];
				} else {
					ens_A = names(database$symbol_ensem[database$symbol_ensem == choice_geneA])
				}
				if ( length(ens_A) == 1 ) { #// if gene symbol matches to ensembl_id
					object_individual_A$value <- FuSViz::get_annotation_db(ens_A, database$txdb, database$grTrack)
					if ( is.null(object_individual_A$value) ) {
						updateSelectizeInput(session = session, inputId = "transA_individual", choices = "") #// NOTE: set choices = "" not (NULL)
						howModal(modalDialog(title = "Warning message", paste("GeneA(", choice_geneA, ") is not available in annotation database, and process stops here!", sep=""))); req(NULL)
					} else {
						#// update the choice values in selectInput
                        canonical_transA = database$canonical[database$canonical$ensembl_gene_id == ens_A,]$ensembl_transcript_id
                        selectinput_A = data.frame(TXNAME=object_individual_A$value$dataset$TXNAME, TAG=object_individual_A$value$dataset$TXNAME, stringsAsFactors=FALSE)
                        if (! is.na(canonical_transA[1]) ) { #// transcripts with canonical are tagged
                            if ( length(selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG) > 0 ) {
                                selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG = paste("@", selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG, sep="");
                            }
                        }
                        selectinput_A = as.list(setNames(selectinput_A$TAG, selectinput_A$TXNAME)); #// key is TAG; value is TXNAME
                        updateSelectizeInput(session = session, inputId = "transA_individual", choices = selectinput_A, selected = NULL,
						options = list(placeholder = 'Ensembl transcript id', render = I('{ option: function(item, escape) {
								const regExpStr = "@";
								return RegExp(regExpStr, "g").test(item.value) ? "<div><u>" + escape(item.label) + "</u></div>" : "<div>" + escape(item.label) + "</div>"
						}}')))
                    }
				} else {
					object_individual_A$value = NULL
					updateSelectizeInput(session = session, inputId = "transA_individual", choices = "") #// NOTE: set choices = "" not (NULL)
					showModal(modalDialog(title = "Warning message", paste("No or multiple geneA(", choice_geneA, ") is present (check gene symbol!)", sep="")));	req(NULL)
				}
			}
        }
	})
    #// when input[['individual-gene2']] NOT null, assign 'object_individual_B$value' using FUNCTION "get_annotation_db"
	observeEvent(input[['individual-gene2']], {
        if ( is.null(input[['individual-gene2']]) || input[['individual-gene2']] == "" ) {
            object_individual_B = NULL
            updateSelectizeInput(session = session, inputId = "transB_individual", choices = "")
        } else {
			if ( is.null(database$txdb) || is.null(database$grTrack) ) { req(NULL); }
			if ( length(unique(as.character(individual_data()$gene2))) == 1 ) {
				choice_geneB = unique(as.character(individual_data()$gene2))
				ens_B = NULL;
				if ( grepl("ENSG00", choice_geneB) || grepl("ENSMUSG0", choice_geneB) ) {
					ens_B = database$ensembl_id[database$ensembl_id == choice_geneB];
				} else {
					ens_B = names(database$symbol_ensem[database$symbol_ensem == choice_geneB])
				}
				if ( length(ens_B) == 1 ) { #// if gene symbol matches to ensembl_id
					object_individual_B$value <- FuSViz::get_annotation_db(ens_B, database$txdb, database$grTrack)
					if ( is.null(object_individual_B$value) ) {
						updateSelectizeInput(session = session, inputId = "transB_individual", choices = "") #// NOTE: set choices = "" not (NULL)
						howModal(modalDialog(title = "Warning message", paste("GeneB(", choice_geneB, ") is not available in annotation database, and process stops here!", sep="")));	req(NULL)
					} else {
						#// update the choice values in selectInput
                        canonical_transB = database$canonical[database$canonical$ensembl_gene_id == ens_B,]$ensembl_transcript_id
						selectinput_B = data.frame(TXNAME=object_individual_B$value$dataset$TXNAME, TAG=object_individual_B$value$dataset$TXNAME, stringsAsFactors=FALSE)
                        if (! is.na(canonical_transB[1]) ) { #// transcripts with canonical are tagged
                            if ( length(selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG) > 0 ) {
                                selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG = paste("@", selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG, sep="");
                            }
                        }
                        selectinput_B = as.list(setNames(selectinput_B$TAG, selectinput_B$TXNAME)); #// key is TAG; value is TXNAME
                        updateSelectizeInput(session = session, inputId = "transB_individual", choices = selectinput_B, selected = NULL,
						options = list(placeholder = 'Ensembl transcript id', render = I('{ option: function(item, escape) {
								const regExpStr = "@";
								return RegExp(regExpStr, "g").test(item.value) ? "<div><u>" + escape(item.label) + "</u></div>" : "<div>" + escape(item.label) + "</div>"
                        }}')))

					}
				} else {
					object_individual_B$value = NULL
					updateSelectizeInput(session = session, inputId = "transB_individual", choices = "") #// NOTE: set choices = "" not (NULL)
					showModal(modalDialog(title = "Warning message", paste("No or multiple geneB(", choice_geneB, ") is present (check gene symbol!)", sep="")));	req(NULL)
				}
			}
        }
	})
    	#// click button 'input$individual_on' for plotting
		observeEvent(input$individual_on, {
			shinycssloaders::showSpinner("chimerics2");
			individualAB$geneA=NULL; individualAB$symbol_A=NULL; individualAB$geneB=NULL; individualAB$symbol_B=NULL; individualAB$breakpoint_set=NULL;
			if ( !is.null(object_individual_A$value) && !is.null(object_individual_B$value) ) { #// make sure gene symbols match to valid ensembl_id for geneA and geneB
				if ( length(individual_data()[,1]) == 1 ) {
					geneA = list(); geneB = list(); # // geneA and geneB are list structure
					symbol_A = individual_data()[1,]$gene1
					symbol_B = individual_data()[1,]$gene2
					breakpoint_set = data.frame(breakpoint_A=individual_data()[1,]$pos1, breakpoint_B=individual_data()[1,]$pos2, 
							split=individual_data()[1,]$split, span=individual_data()[1,]$span, strand_A=individual_data()[1,]$strand1, strand_B=individual_data()[1,]$strand2, stringsAsFactors=FALSE)
					
					mytmp_A = list(); #// a subset of 'object_individual_A' as input for FUNCTION gene_trans_ex
					mytmp_A$value$txTr_f = object_individual_A$value$txTr_f; 
					mytmp_A$value$chr = object_individual_A$value$chr; 
					mytmp_A$value$start = object_individual_A$value$start; 
					mytmp_A$value$end = object_individual_A$value$end; 
					mytmp_A$value$strand = object_individual_A$value$strand; 
					if (! is.null(input$transA_individual) ) { #// if selected transcript of geneA present
                        name_transA = input$transA_individual;	name_transA = gsub('@', '', name_transA);
						mytmp_A$value$dataset = object_individual_A$value$dataset[object_individual_A$value$dataset$TXNAME %in% name_transA, ]
					} else {
						mytmp_A$value$dataset = object_individual_A$value$dataset
					}				
					geneA[[as.character(breakpoint_set[1,1])]] <- FuSViz::gene_trans_ex(breakpoint_set[1,1], mytmp_A$value, database$whole_txdb) #// For geneA
					mytmp_B = list(); #// a subset of 'object_individual_B' as input for FUNCTION gene_trans_ex
					mytmp_B$value$txTr_f = object_individual_B$value$txTr_f; 
					mytmp_B$value$chr = object_individual_B$value$chr;
					mytmp_B$value$start = object_individual_B$value$start; 
					mytmp_B$value$end = object_individual_B$value$end; 
					mytmp_B$value$strand = object_individual_B$value$strand; 
					if (! is.null(input$transB_individual) ) { #// if selected transcript of geneB present
                        name_transB = input$transB_individual;	name_transB = gsub('@', '', name_transB);
						mytmp_B$value$dataset = object_individual_B$value$dataset[object_individual_B$value$dataset$TXNAME %in% name_transB, ]
					} else {
						mytmp_B$value$dataset = object_individual_B$value$dataset
					}
					geneB[[as.character(breakpoint_set[1,2])]] <- FuSViz::gene_trans_ex(breakpoint_set[1,2], mytmp_B$value, database$whole_txdb) #// For geneB

					if ( length(geneA) > 0 ) {
						geneA[sapply(geneA, is.null)] <- NULL; #// remove element with NULL
					} else {
						showModal(modalDialog(title = "Warning message", "No breakpoints within geneA or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!)"));
					}
					if ( length(geneB) > 0 ) {
						geneB[sapply(geneB, is.null)] <- NULL; #// remove element with NULL
					} else {
						showModal(modalDialog(title = "Warning message", "No breakpoints within geneB or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!)"));
					}
					individualAB$geneA=geneA;	individualAB$symbol_A=symbol_A;
					individualAB$geneB=geneB;	individualAB$symbol_B=symbol_B;
					individualAB$breakpoint_set=breakpoint_set;
				} else {
					showModal(modalDialog(title = "Warning message", "Multiple entries are not allowed!"));
					shinycssloaders::hideSpinner("chimerics2")
				}
			} else {
				showModal(modalDialog(title = "Warning message", "No gene symbols are chosen for fusion partners (choose gene symbol first!)"));
				shinycssloaders::hideSpinner("chimerics2");
			}
		})
    #// start persample_size
	output$chimerics2 <- renderPlot({
		if ( length(individualAB$geneA) > 0 && length(individualAB$geneB) > 0 ) {
			print("start to plotting individual fusion");
			individual_control = persample_size_full();
			individual_layout = tryCatch(FuSViz::plot_separate_individual(individualAB$geneA, individualAB$symbol_A, individualAB$geneB, individualAB$symbol_B, individualAB$breakpoint_set, database$chrTrack), error = function(e) {});
			if ( is.null(individual_layout) ) {
				showModal(modalDialog(title = "Warning message", paste("No enough space to draw all ", individualAB$symbol_A, " / ", individualAB$symbol_B, " transcript isoforms in the current layout device.", " Please increase the layout device size by pulling 'Zoom in/out' slide bar!", sep="")));
			}
		} else {
			return(NULL);
		}
		shinycssloaders::hideSpinner("chimerics2");
	}, height = persample_size_full, width = persample_size_full)
	#// download persample plot
	output$FusionDown2 <- downloadHandler(
		filename = function(){ paste("Two-way_persample_plot", tolower(input$file_fusion2), sep =".") },
		content = function(file){
			pixelratio = 2
			if( input$file_fusion2 == "PNG" ) {
				png(file, width = input$sample_width, height = input$sample_height, res=72*pixelratio, units = "in")
			} else {
				pdf(file, width = input$sample_width, height = input$sample_height)
			}
			if ( length(individualAB$geneA) > 0 && length(individualAB$geneB) > 0 ) {
				individual_output = tryCatch(FuSViz::plot_separate_individual(individualAB$geneA, individualAB$symbol_A, individualAB$geneB, individualAB$symbol_B, individualAB$breakpoint_set, database$chrTrack), error = function(e) {});
				if ( is.null(individual_output) ) {
					showModal(modalDialog(title = "Warning message", paste("Download an incomplete plot due to a small layout size. ", "Please increase height and width values in 'Layout_width' and 'Layout_height' boxes to get a full plot!", sep="")));
				}
            } else {
				plot.new()
			}
			dev.off()
		}
	)

    #------ For domain visualization -----#
	observe({
		 #// if no file loads, stop busy indicator bar
		if ( is.null(inputFile$rnadata) ) { output$linking <- NULL; }
    })
    object_domain_A <- reactiveValues(value=NULL) #// 'object_domain_A$value' is list class with 6 elements:
	object_domain_B <- reactiveValues(value=NULL) #// 'object_domain_B$value' is list class with 6 elements:
	control_break_AB <- reactiveValues(A=NULL, B=NULL) #// 'control_break_AB$A' and 'control_break_AB$B' are vector of breakpoints for geneA and geneB
	domain_plotA <- reactiveValues(geneA=NULL, symbol_A=NULL, domainA=NULL, motifA=NULL)
	domain_plotB <- reactiveValues(geneB=NULL, symbol_B=NULL, domainB=NULL, motifB=NULL)
	domain_plot_link <- reactiveValues(A1_xy=NULL, B1_xy=NULL)
    domain_filter <- reactiveValues(value=NULL)

    domain_1 <- select_group_server(id = "domain1", data = reactive(inputFile$rnadata), vars = reactive(c("gene1", "gene2")))
    #// filter breakpoints (i.e. keep breakpoints inside of the selected transcript_id of geneA and geneB)
	observe({
        if ( is.null(inputFile$rnadata) ) { req(NULL); }
		domain_break = unique(data.frame(pos1=domain_1()$pos1, pos2=domain_1()$pos2, strand1=domain_1()$strand1, strand2=domain_1()$strand2, comment=domain_1()$untemplated_insert, stringsAsFactors=FALSE))
		if ( length(control_break_AB$A) > 0 ) { #// control_break_AB$A not empty
			if ( length(control_break_AB$B) > 0 ) { #// control_break_AB$B not empty
				domain_filter$value = domain_break[domain_break$pos1 %in% control_break_AB$A & domain_break$pos2 %in% control_break_AB$B, ]
			} else {
				domain_filter$value = domain_break[domain_break$pos1 %in% control_break_AB$A, ]
			}
		} else {
			if ( length(control_break_AB$B) > 0 ) {
				domain_filter$value = domain_break[domain_break$pos2 %in% control_break_AB$B, ]
			} else {
				domain_filter$value = domain_break
			}
		}
	})
    domain_data <- select_group_server(id = "domain2", data = reactive(domain_filter$value), vars = reactive(c("pos1", "pos2")))

    #// if gene symbol is selected in input[['domain1-gene1']], set'control_break_AB$A', 'control_break_AB$B' and 'domain_plotA$geneA' as NULL when only one symbol selected
    observe({
        if ( is.null(input[['domain1-gene1']]) || input[['domain1-gene1']] == "" ) {
            object_domain_A$value = NULL
            domain_plot_link$A1_xy = NULL
			updateSelectizeInput(session = session, inputId = "domainA", choices = "")
            updateNumericInput(session = session, inputId = "offset_base_A", value = 5)
            output$linking <- renderPlot({ return(NULL); })
            output$domain_up <- renderPlot({ return(NULL); })
        } else {
		    if ( is.null(database$txdb) || is.null(database$grTrack) || is.null(database$domain) || is.null(database$motif) ) { req(NULL); }
            if ( length(unique(as.character(domain_1()$gene1))) == 1 ) { #// make sure only one symbol of geneA selected
            	choice_geneA = unique(as.character(domain_1()$gene1))
				ens_A = NULL;
				if ( grepl("ENSG00", choice_geneA) || grepl("ENSMUSG0", choice_geneA) ) {
					ens_A = database$ensembl_id[database$ensembl_id == choice_geneA];
				} else {
					ens_A = names(database$symbol_ensem[database$symbol_ensem == choice_geneA])
				}
				if ( length(ens_A) == 1 ) { #// if gene symbol matches to ensembl_id
					object_domain_A$value <- FuSViz::get_annotation_db(ens_A, database$txdb, database$grTrack)
					if ( is.null(object_domain_A$value) ) {
						updateSelectizeInput(session = session, inputId = "domainA", choices = "") 
						showModal(modalDialog(title = "Warning message", paste("GeneA(", choice_geneA, ") is not available in annotation database, and process stops here!", sep=""))); req(NULL)
					} else {
						domain_geneA = database$domain[database$domain$Gene_id == ens_A, ]
						motif_geneA = database$motif[database$motif$Gene_id == ens_A, ]
						#// get transcript_id with domain or motif annotation
						DM_transA = unique(c(as.character(domain_geneA$Transcript_id), as.character(motif_geneA$Transcript_id)))
                   		canonical_transA = database$canonical[database$canonical$ensembl_gene_id == ens_A,]$ensembl_transcript_id
						#// only keep transcripts of geneA that harbor breakpoints
						break_max_A = max(domain_1()$pos1)
						break_min_A = min(domain_1()$pos1)
						tmp_AA = object_domain_A$value$dataset
						tmp_AA = object_domain_A$value$dataset[break_max_A >= object_domain_A$value$dataset$TXSTART & break_min_A <= object_domain_A$value$dataset$TXEND, ]
						if ( length(dim(tmp_AA)[1]) == 0 ) { #// judge whether 0-row in tmp_AA
							showModal(modalDialog(title = "Warning message", "No breakpoints within geneA or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!)")); req(NULL);
						}
						selectinput_A = data.frame(TXNAME=tmp_AA$TXNAME, TAG=tmp_AA$TXNAME, stringsAsFactors=FALSE)
						if (! is.na(DM_transA[1]) ) { #// transcripts harboring domain/motif are present
							if ( length(selectinput_A[selectinput_A$TXNAME %in% DM_transA, ]$TAG) > 0 ) {
								selectinput_A[selectinput_A$TXNAME %in% DM_transA, ]$TAG = paste("-", selectinput_A[selectinput_A$TXNAME %in% DM_transA, ]$TAG, sep="");
							}
						}
                    	if (! is.na(canonical_transA[1]) ) { #// transcripts with canonical are tagged
                        	if ( length(selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG) > 0 ) {
                            	selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG = paste("@", selectinput_A[selectinput_A$TXNAME %in% canonical_transA, ]$TAG, sep="");
                        	}
                    	}
						selectinput_A = as.list(setNames(selectinput_A$TAG, selectinput_A$TXNAME)); #// key is TAG; value is TXNAME
						#// font of transcript_id harboring domain/motif is bold
						updateSelectizeInput(session = session, inputId = "domainA", choices = selectinput_A, selected = NULL,
							options = list(placeholder = 'Ensembl transcript id', maxItems = 1, render = I('{ option: function(item, escape) {
                                const regExpStr_dm = "-";
                                const regExpStr_can = "@";
								if ( RegExp(regExpStr_dm, "g").test(item.value) === true ) {
                                    if ( RegExp(regExpStr_can, "g").test(item.value) === true ) {
                                        return("<div><b><u>" + escape(item.label) + "</u></b></div>")
                                    } else {
                                        return("<div><b>" + escape(item.label) + "</b></div>")
                                    }
                                } else {
                                    if ( RegExp(regExpStr_can, "g").test(item.value) ) {
                                        return("<div><u>" + escape(item.label) + "</u></div>")
                                    } else {
                                        return("<div>" + escape(item.label) + "</div>")
                                    }
                            	}
						}}')))
					}
				} else {
					object_domain_A$value = NULL
					updateSelectizeInput(session = session, inputId = "domainA", choices = "")
					showModal(modalDialog(title = "Warning message", paste("No or multiple geneA(", choice_geneA, ") is present (check gene symbol!)", sep="")));	req(NULL);
				}
			}
        }
	})
    #// if gene symbol is selected in input[['domain1-gene2']], set 'control_break_AB$A', 'control_break_AB$B' and 'domain_plotB$geneB' as NULL when only one symbol selected
	observe({
        if ( is.null(input[['domain1-gene2']]) || input[['domain1-gene2']] == "" ) {
            object_domain_B$value = NULL
            domain_plot_link$B1_xy = NULL
			updateSelectizeInput(session = session, inputId = "domainB", choices = "")
            updateNumericInput(session = session, inputId = "offset_base_B", value = 5)
            output$linking <- renderPlot({ return(NULL); })
            output$domain_down <- renderPlot({ return(NULL); })
        } else {
		    if ( is.null(database$txdb) || is.null(database$grTrack) || is.null(database$domain) || is.null(database$motif) ) { req(NULL); }
            if ( length(unique(as.character(domain_1()$gene2))) == 1 ) { #// make sure only one symbol of geneB selected
				choice_geneB = unique(as.character(domain_1()$gene2))
				ens_B = NULL;
				if ( grepl("ENSG00", choice_geneB) || grepl("ENSMUSG0", choice_geneB) ) {
					ens_B = database$ensembl_id[database$ensembl_id == choice_geneB];
				} else {
					ens_B = names(database$symbol_ensem[database$symbol_ensem == choice_geneB])
				}
				if ( length(ens_B) == 1 ) { #// if gene symbol matches to ensembl_id
					object_domain_B$value <- FuSViz::get_annotation_db(ens_B, database$txdb, database$grTrack)
					if ( is.null(object_domain_B$value) ) {
						updateSelectizeInput(session = session, inputId = "domainB", choices = "")
						showModal(modalDialog(title = "Warning message", paste("GeneB ", choice_geneB, " is not available in annotation database, and process stops here!", sep=""))); req(NULL)
					} else {
						domain_geneB = database$domain[database$domain$Gene_id == ens_B, ]
						motif_geneB = database$motif[database$motif$Gene_id == ens_B, ]
						DM_transB = unique(c(as.character(domain_geneB$Transcript_id), as.character(motif_geneB$Transcript_id))) #// get transcript_id with domain or motif annotation
                    	canonical_transB = database$canonical[database$canonical$ensembl_gene_id == ens_B,]$ensembl_transcript_id
						#// only keep transcripts of geneB that harbor breakpoints
						break_max_B = max(domain_1()$pos2)
						break_min_B = min(domain_1()$pos2)
						tmp_BB = object_domain_B$value$dataset
						tmp_BB = object_domain_B$value$dataset[break_max_B >= object_domain_B$value$dataset$TXSTART & break_min_B <= object_domain_B$value$dataset$TXEND, ]
						if ( length(dim(tmp_BB)[1]) == 0 ) { #// judge whether 0-row in tmp_BB 
							showModal(modalDialog(title = "Warning message", "No breakpoints within geneB or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!)")); req(NULL);
						}
						selectinput_B = data.frame(TXNAME=tmp_BB$TXNAME, TAG=tmp_BB$TXNAME, stringsAsFactors=FALSE)
						if (! is.na(DM_transB[1]) ) { #// transcript harboring domain/motif are present
							if ( length(selectinput_B[selectinput_B$TXNAME %in% DM_transB, ]$TAG) > 0 ) {
								selectinput_B[selectinput_B$TXNAME %in% DM_transB, ]$TAG = paste("-", selectinput_B[selectinput_B$TXNAME %in% DM_transB, ]$TAG, sep="");
							}
						}
                    	if (! is.na(canonical_transB[1]) ) { #// transcripts with canonical are tagged
                        	if ( length(selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG) > 0 ) {
                        		selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG = paste("@", selectinput_B[selectinput_B$TXNAME %in% canonical_transB, ]$TAG, sep="");
                        	}
                    	}
						selectinput_B = as.list(setNames(selectinput_B$TAG, selectinput_B$TXNAME))
						#// font of transcript_id harboring domain/motif is bold
						updateSelectizeInput(session = session, inputId = "domainB", choices = selectinput_B, selected = NULL,
							options = list(placeholder = 'Ensembl transcript id', maxItems = 1, render = I('{ option: function(item, escape) {
								const regExpStr_dm = "-";
                                const regExpStr_can = "@";
								if ( RegExp(regExpStr_dm, "g").test(item.value) === true ) {
                                    if ( RegExp(regExpStr_can, "g").test(item.value) === true ) {
                                        return("<div><b><u>" + escape(item.label) + "</u></b></div>")
                                    } else {
                                        return("<div><b>" + escape(item.label) + "</b></div>")
                                    }
                                } else {
                                    if ( RegExp(regExpStr_can, "g").test(item.value) ) {
                                        return("<div><u>" + escape(item.label) + "</u></div>")
                                    } else {
                                        return("<div>" + escape(item.label) + "</div>")
                                    }
                                }
							}}')))
					}
				} else {
					object_domain_B$value = NULL
					updateSelectizeInput(session = session, inputId = "domainB", choices = "")
					showModal(modalDialog(title = "Warning message", paste("No or multiple geneB(", choice_geneB, ") is present (check gene symbol!)", sep="")));	req(NULL);
				}
			}
        }
	})
    #// select transcript_id of geneA
	observeEvent(input$domainA, {
		if ( is.null(object_domain_A$value) || is.null(input$domainA) || input$domainA == "" ) { #// if 'object_domain_A$value' is NULL, set 'domain_plotA$geneA, $symbol_A, $domainA, $motifA' as NULL
			domain_plotA$geneA = NULL; domain_plotA$symbol_A = NULL; domain_plotA$domainA = NULL; domain_plotA$motifA = NULL;
            domain_plot_link$A1_xy = NULL;
            output$domain_up <- renderPlot({ return(); })
            output$linking <- renderPlot({ return(); })
		} else {
			name_domainA = input$domainA;	name_domainA = gsub('-', '', name_domainA); name_domainA = gsub('@', '', name_domainA);
			#// a subset transcript of geneA after update the choice values in selectInput
			mytmp_A = list(); #// a subset of 'object_domain_A' (for selected transcripts of geneA)
			mytmp_A$value$dataset = object_domain_A$value$dataset[object_domain_A$value$dataset$TXNAME %in% name_domainA, ]; 
			mytmp_A$value$txTr_f = object_domain_A$value$txTr_f; 
			mytmp_A$value$chr = object_domain_A$value$chr; 
			mytmp_A$value$start = object_domain_A$value$start;
			mytmp_A$value$end = object_domain_A$value$end; 
			mytmp_A$value$strand = object_domain_A$value$strand;

			geneA = list(); #// geneA is a list class
			symbol_A = domain_1()[1,]$gene1
			breakpoint_set = unique(domain_1()$pos1)
			for (i in 1:length(breakpoint_set)) { # print(paste("start domain A loop: ", i))
				if (! exists(as.character(breakpoint_set[i]), where=geneA) ) {
					geneA[[as.character(breakpoint_set[i])]] <- FuSViz::gene_trans_ex_reduce(breakpoint_set[i], mytmp_A$value, database$whole_txdb, "upstream", input$offset_base_A) #// For geneA; 
				}
			}
			#// remove the elements in the list where breakpoint outside transcript region of geneA
			if ( length(geneA) == 0 ) { showModal(modalDialog(title = "Warning message", "No breakpoints within geneA or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!")); req(NULL); }
			geneA = Filter(Negate(function(x) is.null(unlist(x))), geneA)
			control_break_AB$A = as.numeric(names(geneA)); # print(paste("selected breakpoint for geneA: ", as.numeric(names(geneA))));
			domain_plotA$geneA=geneA
			domain_plotA$symbol_A=symbol_A
			domain_plotA$domainA=database$domain[database$domain$Transcript_id %in% name_domainA, ]
			domain_plotA$motifA=database$motif[database$motif$Transcript_id %in% name_domainA, ]
        	if ( length(domain_plotA$geneA) > 0 ) {
				output$domain_up <- renderPlot({
					print("start to plotting domain of upstream gene")
					domain_plot_link$A1_xy = FuSViz::plot_separate_domain_geneA(domain_plotA$geneA, domain_plotA$symbol_A, domain_plotA$domainA, domain_plotA$motifA)
				})
        	}
        }
	})
    observeEvent(input$offset_base_A, {
		if ( is.null(object_domain_A$value) || is.null(input$domainA) || input$domainA == "" ) { #// if 'object_domain_A$value' is NULL, set 'domain_plotA$geneA, $symbol_A, $domainA, $motifA' as NULL
			domain_plotA$geneA = NULL; domain_plotA$symbol_A = NULL; domain_plotA$domainA = NULL; domain_plotA$motifA = NULL;
            domain_plot_link$A1_xy = NULL;
            output$domain_up <- renderPlot({ return(); })
            output$linking <- renderPlot({ return(); })
		} else {
			name_domainA = input$domainA;	name_domainA = gsub('-', '', name_domainA); name_domainA = gsub('@', '', name_domainA);
			#// a subset transcript of geneA after update the choice values in selectInput
			mytmp_A = list(); #// a subset of 'object_domain_A' (for selected transcripts of geneA)
			mytmp_A$value$dataset = object_domain_A$value$dataset[object_domain_A$value$dataset$TXNAME %in% name_domainA, ]; 
			mytmp_A$value$txTr_f = object_domain_A$value$txTr_f; 
			mytmp_A$value$chr = object_domain_A$value$chr; 
			mytmp_A$value$start = object_domain_A$value$start;
			mytmp_A$value$end = object_domain_A$value$end; 
			mytmp_A$value$strand = object_domain_A$value$strand;

			geneA = list(); #// geneA is a list class
			symbol_A = domain_1()[1,]$gene1
			breakpoint_set = unique(domain_1()$pos1)
			for (i in 1:length(breakpoint_set)) { # print(paste("start domain A loop: ", i))
				if (! exists(as.character(breakpoint_set[i]), where=geneA) ) {
					geneA[[as.character(breakpoint_set[i])]] <- FuSViz::gene_trans_ex_reduce(breakpoint_set[i], mytmp_A$value, database$whole_txdb, "upstream", input$offset_base_A) #// For geneA; 
				}
			}
			#// remove the elements in the list where breakpoint outside transcript region of geneA
			if ( length(geneA) == 0 ) { showModal(modalDialog(title = "Warning message", "No breakpoints within geneA or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!")); req(NULL); }
			geneA = Filter(Negate(function(x) is.null(unlist(x))), geneA)
			control_break_AB$A = as.numeric(names(geneA)); # print(paste("selected breakpoint for geneA: ", as.numeric(names(geneA))));
			domain_plotA$geneA=geneA
			domain_plotA$symbol_A=symbol_A
			domain_plotA$domainA=database$domain[database$domain$Transcript_id %in% name_domainA, ]
			domain_plotA$motifA=database$motif[database$motif$Transcript_id %in% name_domainA, ]
       		if ( length(domain_plotA$geneA) > 0 ) {
				output$domain_up <- renderPlot({
					print("start to plotting domain of upstream gene")
					domain_plot_link$A1_xy = FuSViz::plot_separate_domain_geneA(domain_plotA$geneA, domain_plotA$symbol_A, domain_plotA$domainA, domain_plotA$motifA)
				})
        	}
        }
	})
    #// select transcript_id of geneB
    observeEvent(input$domainB, {
		if ( is.null(object_domain_B$value) || is.null(input$domainB) || input$domainB == "" ) { #// if 'object_domain_B$value' is NULL, set 'domain_plotB$geneB, $symbol_B, $domainB, $motifB' as NULL
			domain_plotB$geneB = NULL; domain_plotB$symbol_B = NULL; domain_plotB$domainB = NULL; domain_plotB$motifB = NULL;
            domain_plot_link$B1_xy = NULL;
            output$domain_down <- renderPlot({ return(); })
            output$linking <- renderPlot({ return(); })
		} else {
			name_domainB = input$domainB;	name_domainB = gsub('-', '', name_domainB); name_domainB = gsub('@', '', name_domainB);
			#// a subset transcript of geneB after update the choice values in selectInput
			mytmp_B = list(); #// a subset of 'object_domain_B' (for selected transcripts of geneB)
			mytmp_B$value$dataset = object_domain_B$value$dataset[object_domain_B$value$dataset$TXNAME %in% name_domainB, ];
			mytmp_B$value$txTr_f = object_domain_B$value$txTr_f; 
			mytmp_B$value$chr = object_domain_B$value$chr; 
			mytmp_B$value$start = object_domain_B$value$start;
			mytmp_B$value$end = object_domain_B$value$end; 
			mytmp_B$value$strand = object_domain_B$value$strand;

			geneB = list(); #// geneB is a list class
			symbol_B = domain_1()[1,]$gene2
			breakpoint_set = unique(domain_1()$pos2)
			for (i in 1:length(breakpoint_set)) { # print(paste("start domain B loop: ", i))
				if (! exists(as.character(breakpoint_set[i]), where=geneB) ) {
					geneB[[as.character(breakpoint_set[i])]] <- FuSViz::gene_trans_ex_reduce(breakpoint_set[i], mytmp_B$value, database$whole_txdb, "downstream", input$offset_base_B) #// For geneB
				}
			}
			#// remove the elements in the list where breakpoint outside transcript region of geneB
			if ( length(geneB) == 0 ) { showModal(modalDialog(title = "Warning message", "No breakpoints within geneB or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!")); req(NULL); }
			geneB = Filter(Negate(function(x) is.null(unlist(x))), geneB)
			control_break_AB$B = as.numeric(names(geneB)); # print(paste("selected breakpoint for gene2: ", as.numeric(names(geneB)))); 
			domain_plotB$geneB=geneB
			domain_plotB$symbol_B=symbol_B
			domain_plotB$domainB=database$domain[database$domain$Transcript_id %in% name_domainB, ]
			domain_plotB$motifB=database$motif[database$motif$Transcript_id %in% name_domainB, ]
        	if ( length(domain_plotB$geneB) > 0 ) {
				output$domain_down <- renderPlot({
					print("start to plotting domain of downstream gene")
					domain_plot_link$B1_xy = FuSViz::plot_separate_domain_geneB(domain_plotB$geneB, domain_plotB$symbol_B, domain_plotB$domainB, domain_plotB$motifB)
				})
        	}
        }
	})
    observeEvent(input$offset_base_B, {
		if ( is.null(object_domain_B$value) || is.null(input$domainB) || input$domainB == "" ) { #// if 'object_domain_B$value' is NULL, set 'domain_plotB$geneB, $symbol_B, $domainB, $motifB' as NULL
			domain_plotB$geneB = NULL; domain_plotB$symbol_B = NULL; domain_plotB$domainB = NULL; domain_plotB$motifB = NULL;
            domain_plot_link$B1_xy = NULL;
            output$domain_down <- renderPlot({ return(); })
            output$linking <- renderPlot({ return(); })
		} else {
			name_domainB = input$domainB;	name_domainB = gsub('-', '', name_domainB); name_domainB = gsub('@', '', name_domainB);
			#// a subset transcript of geneB after update the choice values in selectInput
			mytmp_B = list(); #// a subset of 'object_domain_B' (for selected transcripts of geneB)
			mytmp_B$value$dataset = object_domain_B$value$dataset[object_domain_B$value$dataset$TXNAME %in% name_domainB, ];
			mytmp_B$value$txTr_f = object_domain_B$value$txTr_f; 
			mytmp_B$value$chr = object_domain_B$value$chr; 
			mytmp_B$value$start = object_domain_B$value$start;
			mytmp_B$value$end = object_domain_B$value$end; 
			mytmp_B$value$strand = object_domain_B$value$strand;

			geneB = list(); #// geneB is a list class
			symbol_B = domain_1()[1,]$gene2
			breakpoint_set = unique(domain_1()$pos2)
			for (i in 1:length(breakpoint_set)) { # print(paste("start domain B loop: ", i))
				if (! exists(as.character(breakpoint_set[i]), where=geneB) ) {
					geneB[[as.character(breakpoint_set[i])]] <- FuSViz::gene_trans_ex_reduce(breakpoint_set[i], mytmp_B$value, database$whole_txdb, "downstream", input$offset_base_B) #// For geneB
				}
			}
			#// remove the elements in the list where breakpoint outside transcript region of geneB
			if ( length(geneB) == 0 ) { showModal(modalDialog(title = "Warning message", "No breakpoints within geneB or breakpoint coordinates do not match the genome version of imported annotation resource. Plot stops and please check coordinates!")); req(NULL); }
			geneB = Filter(Negate(function(x) is.null(unlist(x))), geneB)
			control_break_AB$B = as.numeric(names(geneB)); # print(paste("selected breakpoint for gene2: ", as.numeric(names(geneB)))); 
			domain_plotB$geneB=geneB
			domain_plotB$symbol_B=symbol_B
			domain_plotB$domainB=database$domain[database$domain$Transcript_id %in% name_domainB, ]
			domain_plotB$motifB=database$motif[database$motif$Transcript_id %in% name_domainB, ]
        	if ( length(domain_plotB$geneB) > 0 ) {
				output$domain_down <- renderPlot({
					print("start to plotting domain of downstream gene")
					domain_plot_link$B1_xy = FuSViz::plot_separate_domain_geneB(domain_plotB$geneB, domain_plotB$symbol_B, domain_plotB$domainB, domain_plotB$motifB)
				})
        	}
        }
	})

    observe({
        if ( is.null(inputFile$rnadata) ) { req(NULL) }
		if ( length(domain_plotA$geneA) > 0 && length(domain_plotB$geneB) > 0 && !is.null(domain_plot_link$A1_xy) && !is.null(domain_plot_link$B1_xy) ) {
            if ( length(rownames(domain_data())) > 0 ) {
                output$linking <- renderPlot({
					print("start to plotting rows to connect partner genes")
					FuSViz::plot_separate_domain_arrow(domain_plot_link$A1_xy, domain_plot_link$B1_xy, domain_plotA$geneA, domain_plotB$geneB, domain_data())
                })
            }
		}
	})
    #// download domain plot
	output$FusionDown3 <- downloadHandler(
		filename = function(){ paste("Two-way_domain_plot", tolower(input$file_fusion3), sep =".") },
		content = function(file){
			pixelratio = 2
			if( input$file_fusion3 == "PNG" ) {
				png(file, width=input$domain_width, height=input$domain_height, res=72*pixelratio, units = "in")
			} else {
				pdf(file, width=input$domain_width, height=input$domain_height)
			}
            if ( !is.null(object_domain_A$value) && !is.null(object_domain_B$value) ) {
			    FuSViz::plot_separate_domain_download(domain_plotA$geneA, domain_plotA$symbol_A, domain_plotA$domainA, domain_plotA$motifA,
							domain_plotB$geneB, domain_plotB$symbol_B, domain_plotB$domainB, domain_plotB$motifB, domain_data())
            } else {
                plot.new()
            }
			dev.off()
		}
	)

#--------------------
# Vis in network
#--------------------
		observe({ #// if no files load, stop busy indicator
			if ( is.null(input$file_dna_data) ) {
				output$DNAhubs <- NULL;
				output$network_dna <- NULL;
			}
			if ( is.null(input$file_rna_data) ) {
				output$RNAhubs <- NULL;
				output$network_rna <- NULL;
			}
		})

		network_rna <- reactive({
			if ( is.null(inputFile$rnadata) ) { return(NULL); }
			tmp = inputFile$rnadata[, c("gene1", "gene2", "name")];	tmp = as.data.frame(tmp);
			assembly = FuSViz::network_process(tmp, "RNA", database$cancergenes, onco_color, supp_color, rela_color, intergenic_color, other_color);
			if ( is.null(assembly$degree_score$nodes) ) {
				showModal(modalDialog(title = "Warning message", "No oncogenes, tumor suppressors and cancer-related genes are available in RNA SV data OR gene symbols between RNA SV data and database resource do not match!")); req(NULL); 
				return(NULL);
			} else {
				#// assemble final data.frame
				nodes <- data.frame(id=as.character(assembly$nodes$nodes), value=as.numeric(assembly$nodes$Freq), label=rep('', length(assembly$nodes$nodes)),
					color=list(background=as.character(assembly$nodes$tag_color), border="#0d0d0d"),
					title=paste0("<p><b>", as.character(assembly$nodes$nodes),"</b></p>"), stringsAsFactors = F);
				edges <- data.frame(from=assembly$edges$gene1, to=assembly$edges$gene2, value=as.numeric(assembly$edges$value),
					title=paste0("<p><br>Num_sample: <b>", as.numeric(assembly$edges$value),"</b></p>"), stringsAsFactors = F);
				return(list(nodes=nodes, edges=edges, degree_score=assembly$degree_score));
			}
		})
		network_dna <- reactive({
			if ( is.null(inputFile$dnadata) ) { return(NULL); }
			tmp = inputFile$dnadata[, c("gene1", "gene2", "name")];	tmp = as.data.frame(tmp);
			assembly = FuSViz::network_process(tmp, "DNA", database$cancergenes, onco_color, supp_color, rela_color, intergenic_color, other_color);
			if ( is.null(assembly$degree_score$nodes) ) {
				showModal(modalDialog(title = "Warning message", "No oncogenes, tumor suppressors and cancer-related genes are available in DNA SV data OR gene symbols between DNA SV data and database resource do not match!")); req(NULL); 
				return(NULL)
			} else {
				#// assemble final data.frame
				nodes <- data.frame(id=as.character(assembly$nodes$nodes), value=as.numeric(assembly$nodes$Freq), label=rep('', length(assembly$nodes$nodes)),
					color=list(background=as.character(assembly$nodes$tag_color), border="#0d0d0d"),
					title=paste0("<p><b>", as.character(assembly$nodes$nodes),"</b></p>"), stringsAsFactors = F);
				edges <- data.frame(from=assembly$edges$gene1, to=assembly$edges$gene2, value=as.numeric(assembly$edges$value),
					title=paste0("<p><br>Num_sample: <b>", as.numeric(assembly$edges$value),"</b></p>"), stringsAsFactors = F);
				return(list(nodes=nodes, edges=edges, degree_score=assembly$degree_score));
			}
		})

		#// Plot network for RNA SVs
		observe({
			name_network_rna = as.character(network_rna()$nodes$id);
			updateSelectizeInput(session = session, inputId="net_select_rna", choices = unique(name_network_rna), selected = "", options = list(placeholder = 'select'), server = TRUE);
		})
		observeEvent(input$net_rna_plot, { 
			#// create table for node hubs
			output$RNAhubs <- DT::renderDataTable({
				if ( is.null(network_rna()) ) { showModal(modalDialog(title = "Warning message", "No RNA SV data is available!"));	return(NULL); }
				DT::datatable(network_rna()$degree_score, options = list(autoWidth = TRUE, initComplete = JS("function(settings, json) {",
					"$(this.api().table().header()).css({'background-color': '#34495E', 'color': '#AEB6BF'});}"))) %>%
                    DT::formatStyle(c('nodes'), backgroundColor = DT::styleEqual(names(database$cancergenes$oncogene), rep(onco_color, length(database$cancergenes$oncogene)))) %>%
					DT::formatStyle(c('nodes'), backgroundColor = DT::styleEqual(names(database$cancergenes$tumorsuppress), rep(supp_color, length(database$cancergenes$tumorsuppress)))) %>%
					DT::formatStyle(c('nodes'), backgroundColor = DT::styleEqual(names(database$cancergenes$related), rep(rela_color, length(database$cancergenes$related))))
			})

			#// create the network for RNA SVs
			output$network_rna <- visNetwork::renderVisNetwork({
				if ( is.null(network_rna()) ) { showModal(modalDialog(title = "Warning message", "No RNA SV data is available!"));	return(NULL); }
				vis_rna_object <- visNetwork::visNetwork(network_rna()$nodes, network_rna()$edges) %>%
					visNetwork::visPhysics(solver="barnesHut", barnesHut=list(gravitationalConstant=-2000, centralGravity=0.3, springConstant=0.01), stabilization=FALSE) %>%
					visNetwork::visIgraphLayout(layout="layout_nicely") %>% #// using igraph layout for visualization to reduce complexity
					visNetwork::visEdges(shadow=FALSE, smooth=FALSE, length=90, arrows=list(to=list(enabled=TRUE, scaleFactor=1)), physics=TRUE) %>%
					visNetwork::visNodes(shadow=FALSE, scaling=list(min=10, max=20), font=list(size=14, vadjust=-20), physics=TRUE) %>%
					visNetwork::visInteraction(navigationButtons=TRUE, multiselect=TRUE) %>%
					visNetwork::visEvents(selectNode = "function(properties) { // console.log(properties); console.log(this);
						var node_name = this.body.data.nodes.get(properties.nodes[0]).label; 
						// console.log(node_name);  console.log(this.body.data.nodes);
						if ( node_name === '' ) {
							if ( properties.nodes[0].match(/\\*/) ) {
								this.body.data.nodes.update({id: properties.nodes[0], label: 'intergenic'});
							} else {
								this.body.data.nodes.update({id: properties.nodes[0], label: properties.nodes[0]});
							}
						} else {
							this.body.data.nodes.update({id: properties.nodes[0], label: ''});
						}}"
					) %>% visNetwork::visLayout(randomSeed = 123)
				visNetwork::visExport(vis_rna_object, type = "png");
			})
		})
		#// modify the display of network for RNA SVs
		observe({
			if (! is.null(network_rna()) ) {
				visNetwork::visNetworkProxy("network_rna") %>% visNetwork::visInteraction(navigationButtons=input$navi_rna) %>%
				visNetwork::visPhysics(solver="barnesHut", barnesHut=list(gravitationalConstant=input$grav_cons_rna, centralGravity=input$cent_grav_rna,
						springConstant=input$spri_cons_rna)) %>% visNetwork::visEdges(length=input$edge_leng_rna) %>%
				visNetwork::visNodes(scaling=list(min=input$node_scale_rna[1], max=input$node_scale_rna[2]), font=list(size=input$node_font_s_rna, vadjust=input$node_font_p_rna))
			}
		})
		observe({
			if (! is.null(network_rna()) ) {
				if (! is.null(input$net_select_rna) ) {
					if ( input$net_select_rna != "" ) {
						visNetwork::visNetworkProxy("network_rna") %>% visNetwork::visFocus(id=input$net_select_rna, scale=2)
					}
				}
			}                                                                                                                                                                                
		})

		#// Plot network for DNA SVs
		observe({
			name_network_dna = as.character(network_dna()$nodes$id);
			updateSelectizeInput(session = session, inputId="net_select_dna", choices = unique(name_network_dna), selected = "", options = list(placeholder = 'select'), server = TRUE);
		})
		observeEvent(input$net_dna_plot, { 
			#// create table for node hubs
			output$DNAhubs <- DT::renderDataTable({
				if ( is.null(network_dna()) ) { showModal(modalDialog(title = "Warning message", "No DNA SV data is available!"));	return(NULL); }
				DT::datatable(network_dna()$degree_score, options = list(autoWidth = TRUE, initComplete = JS("function(settings, json) {",
					"$(this.api().table().header()).css({'background-color': '#34495E', 'color': '#AEB6BF'});}"))) %>%
					DT::formatStyle(c('nodes'), backgroundColor = DT::styleEqual(names(database$cancergenes$oncogene), rep(onco_color, length(database$cancergenes$oncogene)))) %>%
					DT::formatStyle(c('nodes'), backgroundColor = DT::styleEqual(names(database$cancergenes$tumorsuppress), rep(supp_color, length(database$cancergenes$tumorsuppress)))) %>%
					DT::formatStyle(c('nodes'), backgroundColor = DT::styleEqual(names(database$cancergenes$related), rep(rela_color, length(database$cancergenes$related))))			})

			#// create the network for DNA SVs
			output$network_dna <- visNetwork::renderVisNetwork({
				if ( is.null(network_dna()) ) { showModal(modalDialog(title = "Warning message", "No DNA SV data is available!"));	return(NULL); }
				vis_dna_object <- visNetwork::visNetwork(network_dna()$nodes, network_dna()$edges) %>%
					visNetwork::visPhysics(solver="barnesHut", barnesHut=list(gravitationalConstant=-800, centralGravity=0.1, springConstant=0.02), stabilization = FALSE) %>%
					visNetwork::visIgraphLayout(layout="layout_nicely") %>% #// using igraph layout for visualization to reduce complexity
					visNetwork::visEdges(shadow=FALSE, smooth=FALSE, length=10) %>% #// set a fixed length for edges to avoid too crowd
					visNetwork::visNodes(shadow=FALSE, scaling=list(min=10, max=20), font=list(size=14, vadjust=-44), physics=TRUE) %>%
					visNetwork::visInteraction(navigationButtons=TRUE, multiselect=TRUE) %>%
					visNetwork::visEvents(selectNode = "function(properties) { // console.log(properties); console.log(this);
						var node_name = this.body.data.nodes.get(properties.nodes[0]).label; 
						// console.log(node_name);  console.log(this.body.data.nodes);
						if ( node_name === '' ) {
							if ( properties.nodes[0].match(/\\*/) ) {
								this.body.data.nodes.update({id: properties.nodes[0], label: 'intergenic'});
							} else {
								this.body.data.nodes.update({id: properties.nodes[0], label: properties.nodes[0]});
							}
						} else {
							this.body.data.nodes.update({id: properties.nodes[0], label: ''});
						}}"
					) %>% visNetwork::visLayout(randomSeed = 123)
				visNetwork::visExport(vis_dna_object, type = "png");
			})
		})
		#// modify the display of network for DNA SVs
		observe({
			if (! is.null(network_dna()) ) {
				visNetwork::visNetworkProxy("network_dna") %>% visNetwork::visInteraction(navigationButtons=input$navi_dna) %>% 
				visNetwork::visPhysics(solver="barnesHut", barnesHut=list(gravitationalConstant=input$grav_cons_dna, centralGravity=input$cent_grav_dna, 
						springConstant=input$spri_cons_dna)) %>% visNetwork::visEdges(length=input$edge_leng_dna) %>%
				visNetwork::visNodes(scaling=list(min=input$node_scale_dna[1], max=input$node_scale_dna[2]), font=list(size=input$node_font_s_dna, vadjust=input$node_font_p_dna))
			}
		})
		observe({
			if (! is.null(network_dna()) ) {
				if (! is.null(input$net_select_dna) ) {
					if ( input$net_select_dna != "" ) {
						visNetwork::visNetworkProxy("network_dna") %>% visNetwork::visFocus(id=input$net_select_dna, scale=2)
					}
				}
			}
		})

#--------------------
# Vis in bodyTable
#--------------------
		#// summarize of RNA SVs in upload data
		output$RNAcontents <- DT::renderDataTable({
			if ( is.null(inputFile$rnadata) ) { req(NULL);  }
			tmp = inputFile$rnadata;	
			tag1 = apply(tmp, 1, function(x){
				if ( as.character(x[3]) %in% names(database$cancergenes$oncogene) ) { 
					return(1); 
				} else if ( as.character(x[3]) %in% names(database$cancergenes$tumorsuppress) ) { 
					return(3); 
				} else if ( as.character(x[3]) %in% names(database$cancergenes$related) ) { 
					return(5);
				} else {
					return(7);
				} 
			})
			tag2 = apply(tmp, 1, function(x){
				if ( as.character(x[6]) %in% names(database$cancergenes$oncogene) ) { 
					return(1); 
				} else if ( as.character(x[6]) %in% names(database$cancergenes$tumorsuppress) ) { 
					return(3); 
				} else if ( as.character(x[6]) %in% names(database$cancergenes$related) ) { 
					return(5);
				} else {
					return(7);
				} 
			})
			tmp$tag1=tag1;	tmp$tag2=tag2;
			tmp$name = as.factor(tmp$name);	tmp$gene1 = as.factor(tmp$gene1);	tmp$gene2 = as.factor(tmp$gene2);
			if ( nrow(tmp) > 0 ) {
				attr(tmp, "row.names") = paste0('<a href=\'#shiny-tab-igv\' data-toggle=\'tab\' class=\'tablesenz\' onclick=\'myFunction(\"', tmp$chrom1, ':', tmp$pos1-50, '-', tmp$pos1+50, ' ', tmp$chrom2, ':', tmp$pos2-50, '-', tmp$pos2+50, '\");\'>', rownames(tmp), '</a>');
				DT::datatable(tmp, escape = FALSE, editable = 'cell', extensions = 'Buttons', options = list(dom = 'Blfrtip', buttons = list(
					list(extend = "collection", text = 'Download Data', action = DT::JS("function(e, dt, node, config) { Shiny.setInputValue('rnadownload', true, {priority: 'event'}); }")),
					list(extend = "collection", text = 'Delete Row', action = DT::JS("function(e, dt, node, config) { Shiny.setInputValue('rnaDeleteRow', true, {priority: 'event'}); }"))),
					autoWidth = TRUE, columnDefs = list(list(targets = c(14,15), visible = FALSE)), initComplete = JS("function(settings, json) { $(this.api().table().header()).css({'background-color': '#34495E', 'color': '#AEB6BF'});}")),
					filter = list(position = 'top', clear = FALSE, plain = TRUE)) %>% DT::formatStyle(names(tmp), fontSize = '11px') %>%
					DT::formatStyle('gene1','tag1', backgroundColor = DT::styleInterval(c(2, 4, 6), c("#ff8566", "#00ccff", "#ffcc33", "white"))) %>%
					DT::formatStyle('gene2','tag2', backgroundColor = DT::styleInterval(c(2, 4, 6), c("#ff8566", "#00ccff", "#ffcc33", "white")))
			}
		})
		observeEvent(input$rnadownload, {
			showModal(div(id = "rnadownload", modalDialog(downloadButton("rnacsvdownload","Download as csv format"), br(), br(),
				downloadButton("rnatsvdownload","Download as tsv format"), easyClose = TRUE, title = "Download Panel for RNA SVs")))
		})
		output$rnacsvdownload <- downloadHandler(
			filename = function() { paste("rna-sv-data-", Sys.Date(), ".csv", sep="") },
			content = function(file) { write.table(inputFile$rnadata, file, sep=",", col.names=TRUE, row.names=FALSE, quote=F) }
		)
		output$rnatsvdownload <- downloadHandler(
			filename = function() { paste("rna-sv-data-", Sys.Date(), ".tsv", sep="") },
			content = function(file) { write.table(inputFile$rnadata, file, sep="\t", col.names=TRUE, row.names=FALSE, quote=F) }
		)
		observeEvent(input$rnaDeleteRow, {
			if ( is.null(inputFile$rnadata) ) { req(NULL); }
			rrr = inputFile$rnadata
			if (!is.null(input$RNAcontents_rows_selected)) {
				rrr <- rrr[-as.numeric(input$RNAcontents_rows_selected),]
			}
			inputFile$rnadata = rrr
		})

		#// partner gene wordcloud of RNA SVs
		output$wordcloud_rna <- wordcloud2::renderWordcloud2({
			if ( is.null(wordcloud_svrna()$freq)	) { return(NULL); }
			wordcloud2::wordcloud2(wordcloud_svrna()$freq, color=wordcloud_svrna()$colorlist, size=input$word_size_rna, shuffle=T, shape=input$word_shape_rna,
				fontFamily = 'Calibri', fontWeight='normal', minRotation= -pi/2, maxRotation= -pi/2, ellipticity=1, rotateRatio=0)
		})
		#// RNA SVs distribution per sample
		output$rna_hist_1 <- renderPlot({ #// plot 1st section of canvas for RNA SVs distribution
			if ( is.null(hist_svrna()$section1)	) { return(NULL); }
			ymax = max(hist_svrna()$section1);
			xxx = barplot(hist_svrna()$section1, xaxt="n", space=input$rna_hist_space, ylim=c(0, ymax), cex.axis=0.5, tcl=-0.2, mgp=c(3,0.5,0));
			labs = names(hist_svrna()$section1);
			text(cex=input$rna_hist_size, x=xxx, y=input$rna_hist_pos, labs, xpd=TRUE, srt=input$rna_hist_angle);
		})
		output$rna_hist_2 <- renderPlot({ #// plot 2nd section of canvas for RNA SVs distribution
			if ( is.null(hist_svrna()$section2)	) { return(NULL); }
			ymax = max(hist_svrna()$section2);
			xxx = barplot(hist_svrna()$section2, xaxt="n", space=input$rna_hist_space, ylim=c(0, ymax), cex.axis=0.5, tcl=-0.2, mgp=c(3,0.5,0));
			labs = names(hist_svrna()$section2);
			text(cex=input$rna_hist_size, x=xxx, y=input$rna_hist_pos, labs, xpd=TRUE, srt=input$rna_hist_angle);
		})
		output$rna_hist_3 <- renderPlot({ #// plot 3rd section of canvas for RNA SVs distribution
			if ( is.null(hist_svrna()$section3)	) { return(NULL); }
			ymax = max(hist_svrna()$section3);
			xxx = barplot(hist_svrna()$section3, xaxt="n", space=input$rna_hist_space, ylim=c(0, ymax), cex.axis=0.5, tcl=-0.2, mgp=c(3,0.5,0));
			labs = names(hist_svrna()$section3);
			text(cex=input$rna_hist_size, x=xxx, y=input$rna_hist_pos, labs, xpd=TRUE, srt=input$rna_hist_angle);
		})
		output$rna_hist_4 <- renderPlot({ #// plot 4th section of canvas for RNA SVs distribution
			if ( is.null(hist_svrna()$section4)	) { return(NULL); }
			ymax = max(hist_svrna()$section4);
			xxx = barplot(hist_svrna()$section4, xaxt="n", space=input$rna_hist_space, ylim=c(0, ymax), cex.axis=0.5, tcl=-0.2, mgp=c(3,0.5,0));
			labs = names(hist_svrna()$section4);
			text(cex=input$rna_hist_size, x=xxx, y=input$rna_hist_pos, labs, xpd=TRUE, srt=input$rna_hist_angle);
		})

		#// summarize DNA SVs in upload data
		output$DNAcontents <- DT::renderDataTable({ 
			if ( is.null(inputFile$dnadata) ) { req(NULL);  }
			tmp = inputFile$dnadata;	
			tag1 = apply(tmp, 1, function(x){
				if ( as.character(x[11]) %in% names(database$cancergenes$oncogene) ) { 
					return(1); 
				} else if ( as.character(x[11]) %in% names(database$cancergenes$tumorsuppress) ) { 
					return(3); 
				} else if ( as.character(x[11]) %in% names(database$cancergenes$related) ) { 
					return(5);
				} else {
					return(7);
				} 
			})
			tag2 = apply(tmp, 1, function(x){
				if ( as.character(x[12]) %in% names(database$cancergenes$oncogene) ) { 
					return(1); 
				} else if ( as.character(x[12]) %in% names(database$cancergenes$tumorsuppress) ) { 
					return(3); 
				} else if ( as.character(x[12]) %in% names(database$cancergenes$related) ) { 
					return(5);
				} else {
					return(7);
				} 
			})
			tmp$tag1=tag1;	tmp$tag2=tag2;
			tmp$name = as.factor(tmp$name);	tmp$gene1 = as.factor(tmp$gene1);	tmp$gene2 = as.factor(tmp$gene2);	tmp$type = as.factor(tmp$type);
			if ( nrow(tmp) > 0 ) {
				attr(tmp, "row.names") = paste0('<a href=\'#shiny-tab-igv\' data-toggle=\'tab\' onclick=\'myFunction(\"', tmp$chrom1, ':', tmp$start1-50, '-', tmp$end1+50, ' ', tmp$chrom2, ':', tmp$start2-50, '-', tmp$end2+50, '\");\'>', rownames(tmp), '</a>');
				DT::datatable(tmp, escape = FALSE, editable = 'cell', extensions = 'Buttons', options = list(dom = 'Blfrtip', buttons = list(
					list(extend = "collection", text = 'Download Data', action = DT::JS("function(e, dt, node, config) { Shiny.setInputValue('dnadownload', true, {priority: 'event'}); }")),
					list(extend = "collection", text = 'Delete Row', action = DT::JS("function(e, dt, node, config) { Shiny.setInputValue('dnaDeleteRow', true, {priority: 'event'}); }"))),
					autoWidth = TRUE, columnDefs = list(list(targets = c(14,15), visible = FALSE)), initComplete = JS("function(settings, json) { $(this.api().table().header()).css({'background-color': '#34495E', 'color': '#AEB6BF'});}")),
					filter = list(position = 'top', clear = FALSE, plain = TRUE)) %>% DT::formatStyle(names(tmp), fontSize = '11px') %>%
					DT::formatStyle('gene1','tag1', backgroundColor = DT::styleInterval(c(2, 4, 6), c("#ff8566", "#00ccff", "#ffcc33", "white"))) %>%
					DT::formatStyle('gene2','tag2', backgroundColor = DT::styleInterval(c(2, 4, 6), c("#ff8566", "#00ccff", "#ffcc33", "white")))
			}
		})
		observeEvent(input$dnadownload, {
			showModal(div(id = "dnadownload", modalDialog(downloadButton("dnacsvdownload","Download as csv format"), br(), br(),
				downloadButton("dnatsvdownload","Download as tsv format"), easyClose = TRUE, title = "Download Panel for DNA SVs")))
		})
		output$dnacsvdownload <- downloadHandler(
			filename = function() { paste("dna-sv-data-", Sys.Date(), ".csv", sep="") },
			content = function(file) { write.table(inputFile$dnadata, file, sep=",", col.names=TRUE, row.names=FALSE, quote=F) }
		)
		output$dnatsvdownload <- downloadHandler(
			filename = function() { paste("dna-sv-data-", Sys.Date(), ".tsv", sep="") },
			content = function(file) { write.table(inputFile$dnadata, file, sep="\t", col.names=TRUE, row.names=FALSE, quote=F) }
		)
		observeEvent(input$dnaDeleteRow, {
			if ( is.null(inputFile$dnadata) ) { req(NULL); }
			ddd = inputFile$dnadata
			if (!is.null(input$DNAcontents_rows_selected)) {
				ddd <- ddd[-as.numeric(input$DNAcontents_rows_selected),]
			}
			inputFile$dnadata = ddd
		})

		#// partner gene wordcloud of DNA SVs
		output$wordcloud_dna <- wordcloud2::renderWordcloud2({
			if ( is.null(wordcloud_svdna()$freq)	) { return(NULL); }
			wordcloud2::wordcloud2(wordcloud_svdna()$freq, color=wordcloud_svdna()$colorlist, size=input$word_size_dna, shuffle=T, shape=input$word_shape_dna,
				fontFamily = 'Calibri', fontWeight='normal', minRotation= -pi/2, maxRotation= -pi/2, ellipticity=1, rotateRatio=0)
		})
		#// DNA SVs distribution per sample
		output$category <- renderPlot({
			if ( is.null(inputFile$dnadata) ) { req(NULL); }
			legend("center", fill=c("red", "blue", "yellow", "green", "pink"), legend=c("BND", "DEL", "DUP", "INS", "INV"), ncol = 5, cex=0.8, adj=c(0, 0.2))
		})
		output$dna_hist_1 <- renderPlot({ #// plot 1st section of canvas for DNA SVs distribution
			if ( is.null(hist_svdna()$section1) || nrow(hist_svdna()$section1) == 0 ) { return(NULL); }
			xxx = barplot(N ~ type + name, data = hist_svdna()$section1, legend = F, col=c("red", "blue", "yellow", "green", "pink"), xlab = "", ylab = "", xaxt="n", 
						  cex.axis=0.5, space=input$dna_hist_space, tcl=-0.2, mgp=c(3,0.5,0));
			labs = unique(hist_svdna()$section1$name);
			text(cex=input$dna_hist_size, x=xxx, y=input$dna_hist_pos, labs, xpd=TRUE, srt=input$dna_hist_angle);
		})
		output$dna_hist_2 <- renderPlot({ #// plot 2nd section of canvas for DNA SVs distribution
			if ( is.null(hist_svdna()$section2) || nrow(hist_svdna()$section2) == 0 ) { return(NULL); }
			xxx = barplot(N ~ type + name, data = hist_svdna()$section2, legend = F, col=c("red", "blue", "yellow", "green", "pink"), xlab = "", ylab = "", xaxt="n", 
						  cex.axis=0.5, space=input$dna_hist_space, tcl=-0.2, mgp=c(3,0.5,0));
			labs = unique(hist_svdna()$section2$name);
			text(cex=input$dna_hist_size, x=xxx, y=input$dna_hist_pos, labs, xpd=TRUE, srt=input$dna_hist_angle);
		})
		output$dna_hist_3 <- renderPlot({ #// plot 3rd section of canvas for DNA SVs distribution
			if ( is.null(hist_svdna()$section3) || nrow(hist_svdna()$section3) == 0 ) { return(NULL); }
			xxx = barplot(N ~ type + name, data = hist_svdna()$section3, legend = F, col=c("red", "blue", "yellow", "green", "pink"), xlab = "", ylab = "", xaxt="n", 
						  cex.axis=0.5, space=input$dna_hist_space, tcl=-0.2, mgp=c(3,0.5,0));
			labs = unique(hist_svdna()$section3$name);
			text(cex=input$dna_hist_size, x=xxx, y=input$dna_hist_pos, labs, xpd=TRUE, srt=input$dna_hist_angle);
		})
		output$dna_hist_4 <- renderPlot({ #// plot 4th section of canvas for DNA SVs distribution
			if ( is.null(hist_svdna()$section4) || nrow(hist_svdna()$section4) == 0 ) { return(NULL); }
			xxx = barplot(N ~ type + name, data = hist_svdna()$section4, legend = F, col=c("red", "blue", "yellow", "green", "pink"), xlab = "", ylab = "", xaxt="n", 
						  cex.axis=0.5, space=input$dna_hist_space, tcl=-0.2, mgp=c(3,0.5,0));
			labs = unique(hist_svdna()$section4$name);
			text(cex=input$dna_hist_size, x=xxx, y=input$dna_hist_pos, labs, xpd=TRUE, srt=input$dna_hist_angle);
		})

		#// summarize DNA mutations in upload data
		output$DNAMutations <- DT::renderDataTable({
			if ( is.null(inputFile$mutationdata) ) { req(NULL); }
			tmp = inputFile$mutationdata;	
			tag = apply(tmp, 1, function(x){
				if ( as.character(x[1]) %in% names(database$cancergenes$oncogene) ) { 
					return(1); 
				} else if ( as.character(x[1]) %in% names(database$cancergenes$tumorsuppress) ) { 
					return(3); 
				} else if ( as.character(x[1]) %in% names(database$cancergenes$related) ) { 
					return(5);
				} else {
					return(7);
				} 
			})
			tmp$tag = tag;
			tmp$Hugo_Symbol = as.factor(tmp$Hugo_Symbol);	tmp$Tumor_Sample_Barcode = as.factor(tmp$Tumor_Sample_Barcode);
			if ( nrow(tmp) > 0 ) {
				attr(tmp, "row.names") = paste0('<a href=\'#shiny-tab-igv\' data-toggle=\'tab\' onclick=\'myFunction(\"', tmp$Chromosome, ':', tmp$Start_Position-50, '-', tmp$Start_Position+50, '\");\'>', rownames(tmp), '</a>');
				DT::datatable(tmp, escape = FALSE, options = list(autoWidth = TRUE, columnDefs = list(list(targets = c(7), visible = FALSE)), initComplete = JS("function(settings, json) {",
					"$(this.api().table().header()).css({'background-color': '#34495E', 'color': '#AEB6BF'});}")),
					filter = list(position = 'top', clear = FALSE, plain = TRUE)) %>% DT::formatStyle(names(tmp), fontSize = '11px') %>%
					DT::formatStyle('Hugo_Symbol', 'tag', backgroundColor = DT::styleInterval(c(2, 4, 6), c("#ff8566", "#00ccff", "#ffcc33", "white")))
			}
		})
		#// wordcloud of mutated genes
		output$wordcloud_mut <- wordcloud2::renderWordcloud2({
			wordcloud2::wordcloud2(wordcloud_mutat()$freq, color=wordcloud_mutat()$colorlist, size=input$word_size_mut, shuffle=T, shape=input$word_shape_mut,
				fontFamily = 'Calibri', fontWeight='normal', minRotation= -pi/2, maxRotation= -pi/2, ellipticity=1, rotateRatio=0)
		})
		#// scatter plot of mutation load and SV load
		output$rna_mut_cor <- renderPlot({
			if ( is.null(cor_mut()$rna_data_frame) ) { return(NULL); }
			plot(cor_mut()$rna_data_frame, xlab = "log2 transformation of RNA SVs num", ylab = "log2 transformation of mutation num");
		})
		output$rna_mut_info <- renderPrint({
			if ( is.null(cor_mut()$rna_data_frame) ) { return(NULL); }
			brushedPoints(cor_mut()$rna_data_frame, input$rna_brush, xvar="log2_sv_num", yvar="log2_mut_num");
		})
		output$dna_mut_cor <- renderPlot({
			if ( is.null(cor_mut()$dna_data_frame) ) { return(NULL); }
			plot(cor_mut()$dna_data_frame, xlab = "log2 transformation of DNA SVs num", ylab = "log2 transformation of mutation num");
		})
		output$dna_mut_info <- renderPrint({
			if ( is.null(cor_mut()$dna_data_frame) ) { return(NULL); }
			brushedPoints(cor_mut()$dna_data_frame, input$dna_brush, xvar="log2_sv_num", yvar="log2_mut_num");
		})

		#// drug information
		output$drug = DT::renderDataTable({
			if ( is.null(drug_target_match()) ) { return(NULL); }
			if ( nrow(drug_target_match()) > 0 ) {
				escape_vector = -(which(colnames(drug_target_match()) %in% list("molecule_chembl_id","target_chembl_id"))  + 1)
				DT::datatable(drug_target_match(), extensions = c("Buttons","Responsive"), escape = escape_vector,
					options = list(buttons = c('csv','excel'), pageLength = 20, autoWidth = T, width = "100%", fixedColumns = T, fixedHeader = T,
					scrollCollapse = T, columnDefs = list(list(width = '10%', targets = c(1,2,3,4,5))), dom = 'Bfrtip',
					initComplete = JS("function(settings, json) {",
						"$(this.api().table().header()).css({'background-color': '#34495E', 'color': '#AEB6BF'});}"))) %>% DT::formatStyle(names(drug_target_match()), fontSize = '11px')
			}
		})

#---------------------
# Vis in Linear plot
#---------------------
		observe({
			if ( is.null(input$genome) || input$genome == "" ) {
			} else {
				if ( input$genome == "hg19" || input$genome == "hg38" || input$genome == "GRCm39" ) {
					output$FuSViz <- FuSViz::renderFuSViz(
						FuSViz::FuSViz(genomeName=input$genome, displayMode="EXPANDED", trackHeight=200)
					)
				}
			}
		})
		#// load gene annotation file in offline mode
		observeEvent(input$offlineTrackButton, {
			if ( input$genome == "hg19_offline" ) {
				FuSViz::Trackoffline(session, "hg19", name="RefSeq [hg19 offline]")
			} else if ( input$genome == "hg38_offline" ) {
				FuSViz::Trackoffline(session, "hg38", name="RefSeq [hg38 offline]")
			} else if ( input$genome == "GRCm39_offline" ) {
				FuSViz::Trackoffline(session, "GRCm39", name="RefSeq [GRCm39 offline]")
			} else {
				showModal(modalDialog(title = "Warning message", "Load gene track in offline mode not available!")); req(NULL);
			}
		})

		#// Get genomic coordinate
		observeEvent(input$GetPos, {
			FuSViz::Coordinate(session, "GetPosCoord");
		})
		observeEvent(input[["GetPosCoord"]], {
			genomic_region = input[["GetPosCoord"]];
			updateTextInput(session, "genomicPos", value = genomic_region);
		})
		observeEvent(input$MovePos, {
			updateTextInput(session, "genomicPos", value = "");
		})

		#// Load DNA SV events
		observeEvent(input$addTrackButtonBed, {
			print("@ Add Bed-Track and Bedgraph-Track for DNA-seq data@")
			if ( is.null(input_twoway_dna_bed()[["dnabed"]]) ) { showModal(modalDialog(title = "Warning message", "No DNA-SV data is loaded!")); req(NULL); } 
			if ( nrow(input_twoway_dna_bed()[["dnabed"]]) > 0 ) {
				FuSViz::TrackinBed(session, input_twoway_dna_bed()[["dnabed"]], name="DNA SV breakpoint", color="green")
				FuSViz::TrackinBedGraph(session, input_twoway_dna_bed()[["dnabedgraph"]], name="DNA SV breakpoint freq", color="blue", autoscale=TRUE)
			} else {
				showModal(modalDialog(title = "Warning message", "No DNA-SV data available for bed and bedgraph track after filtering control!")); req(NULL);
			}
		})
		observeEvent(input$addTrackButtonBedPe, {
			print("@ Add Bedpe-Track for DNA-seq data@")
			if ( is.null(input_twoway_dna_bedpe()) ) { showModal(modalDialog(title = "Warning message", "No DNA-SV data is loaded!")); req(NULL); }
			if ( nrow(input_twoway_dna_bedpe()) > 0 ) {
				FuSViz::TrackinBedPe(session, input_twoway_dna_bedpe(), name="DNA SV distribution", color="blue", logScale=TRUE)
			} else {
				showModal(modalDialog(title = "Warning message", "No DNA-SV data available for bedpe track after filtering control!")); req(NULL);
			}
		})
		observeEvent(input$addTrackButtonSeg, {
			print("@ Add Seg-Track for DNA-seq data@")
			if ( is.null(input_twoway_dna_seg()[['dup']]) && is.null(input_twoway_dna_seg()[['del']]) ) { showModal(modalDialog(title = "Warning message", "No DNA-SV data is loaded!")); req(NULL); }
			if ( nrow(input_twoway_dna_seg()[['dup']]) > 0 ) {
				if ( nrow(input_twoway_dna_seg()[['del']]) > 0) {
					FuSViz::TrackinSeg(session, input_twoway_dna_seg()[['dup']], name="DNA small CNV - DUP", trackHeight=100)
					FuSViz::TrackinSeg(session, input_twoway_dna_seg()[['del']], name="DNA small CNV - DEL", trackHeight=100)
				} else {
					FuSViz::TrackinSeg(session, input_twoway_dna_seg()[['dup']], name="DNA small CNV - DUP", trackHeight=100)
					showModal(modalDialog(title = "Warning message", "No DEL SVs available for seg track after filtering control!")); req(NULL);
				}
			} else {
				if ( nrow(input_twoway_dna_seg()[['del']]) > 0 ) {
					FuSViz::TrackinSeg(session, input_twoway_dna_seg()[['del']], name="DNA small CNV - DEL", trackHeight=100)
					showModal(modalDialog(title = "Warning message", "No DUP SVs available for seg track after filtering control!")); req(NULL);
				} else {
					showModal(modalDialog(title = "Warning message", "No DUP and DEL SVs available for seg track after filtering control!")); req(NULL);
				}
			}
		})
		#// Load annotation files
		observeEvent(input$load_file, {
			print("@ Load userdefined and indexed *.gz file with bed, vcf and gtf format@")
			if (! is.null(input$load_file) ) {
				FuSViz::TrackinFile(session, input$load_file)
			}
		})
		#// Load alignment url
		observeEvent(input$addTrackButtonBAM, {
			print("@ Add alignment tracks via URL@");
			if (! is.null(input$BAM) && ! is.null(input$BAMindex) ) {
				FuSViz::TrackinBAM(session, input$BAM, input$BAMindex)
			} else {
				showModal(modalDialog(title = "Warning message", "Either hosted BAM/CRAM alignment or its index file is not available!"));
			}
		})
		#// Load RNA SV events
		observeEvent(input$addTrackButtonRNABed, {
			print("@ Add Bed-Track and Bedgraph-Track for RNA-seq data@")
			if ( is.null(input_twoway_rna_bed()[["rnabed"]]) ) { showModal(modalDialog(title = "Warning message", "No RNA-SV data is loaded!")); req(NULL); }
			if ( nrow(input_twoway_rna_bed()[["rnabed"]]) > 0 ) {
				FuSViz::TrackinBed(session, input_twoway_rna_bed()[["rnabed"]], name="RNA SV breakpoint", color="green")
				FuSViz::TrackinBedGraph(session, input_twoway_rna_bed()[["rnabedgraph"]], name="RNA SV breakpoint freq", color="blue", autoscale=TRUE)
			} else {
				showModal(modalDialog(title = "Warning message", "No RNA-SV data available for bed and bedgraph track after filtering control!")); req(NULL);
			}
		})
		observeEvent(input$addTrackButtonRNABedpe, {
			print("@ Add Bedpe-Track for RNA-seq data@")
			if ( is.null(input_twoway_rna_bedpe()) ) { showModal(modalDialog(title = "Warning message", "No RNA-SV data is loaded!")); req(NULL); }
			if ( nrow(input_twoway_rna_bedpe()) > 0 ) {
				FuSViz::TrackinBedPe(session, input_twoway_rna_bedpe(), name="RNA SV distribution", color="blue", logScale=TRUE)
			} else {
				showModal(modalDialog(title = "Warning message", "No RNA-SV data available for bedpe track after filtering control!")); req(NULL);
			}
		})
		#// Load mutation data
		observeEvent(input$addTrackButtonMut, {
			print("@ Add Bed-Track mutation data@")
			if ( is.null(input_twoway_mut_bed()) ) { showModal(modalDialog(title = "Warning message", "No small mutation data is loaded!")); req(NULL); }
			if ( nrow(input_twoway_mut_bed()) > 0 ) {
				FuSViz::TrackinBed(session, input_twoway_mut_bed(), name="Mutation profile", color="orange")
			} else {
				showModal(modalDialog(title = "Warning message", "No mutation data available for bed after filtering control!")); req(NULL);
			}
		})

#----------------------
# End FuSViz app
#----------------------
		shiny::onStop(function() {
		#	unlink(file.path('tmp', '*'))
		#	removeResourcePath('tmp')
			cat("FuSViz application quit and cleanup\n")
		}, NULL)
	}

