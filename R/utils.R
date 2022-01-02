#' Start FuSViz web app
#'
#' @description The app is started with system default browser
#'
#' @param data Data to be loaded with app; the default value as NULL
#'
#' @export
FuSViz_app <- function(data = NULL) {
	DIR = system.file("app", package = "FuSViz")
	if (DIR == "") { stop("Could not find app directory. Try re-installing `FuSViz`.", call. = FALSE) }
	source(file.path(DIR, "ui.R"), local = TRUE, chdir = TRUE)
	source(file.path(DIR, "server.R"), local = TRUE, chdir = TRUE)
	shiny_app = shiny::shinyApp(ui = ui, server = server)
	shiny::runApp(shiny_app, launch.browser = TRUE, display.mode = "normal")
}

#' Input format control on uploading files
#'
#' @description Check the input format of upload files.
#'
#' @param inputdf A data.table object for RNA-seq SVs, DNA-seq SVs or mutation profile.
#' @param type Uploading file source (e.g. RNA SVs, DNA SVs, mutations).
#'
#' @export
check_input_format <- function(inputdf, type) {
	stopifnot(is.character(type));

	#// check whether there are NA or "" for specific column for RNA SV input
	if ( type == "RNA" ) {
		#// NOTE: column1-11 (chrom1, pos1, gene1, chrom2, pos2, gene2, name, split, span, strand1, strand2) - no NA accepted.
		if ( any(is.na(inputdf$chrom1)) == T || any(inputdf$chrom1 == "") == T ) { stop("chrom1 has empty or NA value for RNA SVs!"); }
		if ( any(is.na(inputdf$pos1)) == T || any(inputdf$pos1 <= 0) == T ) { stop("pos1 has empty, NA value or not numeric type for RNA SVs!"); }
		if ( any(is.na(inputdf$gene1)) == T || any(inputdf$gene1 == "") == T ) { stop("gene1 has empty or NA value for RNA SVs!"); }
		if ( any(is.na(inputdf$chrom2)) == T || any(inputdf$chrom2 == "") == T ) { stop("chrom2 has empty or NA value for RNA SVs!"); }
		if ( any(is.na(inputdf$pos2)) == T || any(inputdf$pos2 <= 0) == T ) { stop("pos2 has empty, NA value or not numeric type for RNA SVs!"); }
		if ( any(is.na(inputdf$gene2)) == T || any(inputdf$gene2 == "") == T ) { stop("gene2 has empty or NA value for RNA SVs!"); }
		if ( any(is.na(inputdf$name)) == T || any(inputdf$name == "") == T ) { stop("name has empty or NA value for RNA SVs!"); }
		if ( any(is.na(inputdf$split)) == T || any(inputdf$split < 0) == T ) { stop("split has empty, NA value or not numeric type for RNA SVs!"); }
		if ( any(is.na(inputdf$span)) == T || any(inputdf$sapn < 0) == T ) { stop("span has empty, NA value or not numeric type for RNA SVs!"); }
		if ( any(is.na(inputdf$strand1)) == T || any(inputdf$strand1 == "") == T ) { stop("strand1 has empty or NA value for RNA SVs!"); }
		if ( any(is.na(inputdf$strand2)) == T || any(inputdf$strand2 == "") == T ) { stop("strand2 has empty or NA value for RNA SVs!"); }
	} else if ( type == "DNA" ) {
		#// NOTE: column1-12 (chrom1, start1, end1, chrom2, start2, end2, name, type, split, span, gene1 and gene2) - no NA accepted.
		if ( any(is.na(inputdf$chrom1)) == T || any(inputdf$chrom1 == "") == T ) { stop("chrom1 has empty or NA value for DNA SVs!"); }
		if ( any(is.na(inputdf$start1)) == T || any(inputdf$start1 <= 0) == T ) { stop("start1 has empty, NA value or not numeric type for DNA SVs!"); }
		if ( any(is.na(inputdf$end1)) == T || any(inputdf$end1 <= 0) == T ) { stop("end1 has empty, NA value or numeric type for DNA SVs!"); }
		if ( any(is.na(inputdf$chrom2)) == T || any(inputdf$chrom2 == "") == T ) { stop("chrom2 has empty or NA value for DNA SVs!"); }
		if ( any(is.na(inputdf$start2)) == T || any(inputdf$start2 <= 0) == T ) { stop("start2 has empty, NA value or not numeric type for DNA SVs!"); }
		if ( any(is.na(inputdf$end2)) == T || any(inputdf$end1 <= 0) == T ) { stop("end2 has empty, NA value or not numeric type for DNA SVs!"); }
		if ( any(is.na(inputdf$name)) == T || any(inputdf$name == "") == T ) { stop("name has empty or NA value for DNA SVs!"); }
		if ( any(is.na(inputdf$type)) == T || any(inputdf$type == "") == T ) { stop("type has empty or NA value for DNA SVs!"); }
		if ( any(is.na(inputdf$split)) == T || any(inputdf$split < 0) == T ) { stop("split has empty, NA value or not numeric type for DNA SVs!"); }
		if ( any(is.na(inputdf$span)) == T || any(inputdf$sapn < 0) == T ) { stop("span has empty, NA value or not numeric type for DNA SVs!"); }
		if ( any(is.na(inputdf$gene1)) == T || any(inputdf$gene1 == "") == T ) { stop("gene1 has empty or NA value for DNA SVs!"); }
		if ( any(is.na(inputdf$gene2)) == T || any(inputdf$gene2 == "") == T ) { stop("gene2 has empty or NA value for DNA SVs!"); }
	} else if ( type == "Mut" ) {
		#// remove duplication in mutation data
		inputdf = inputdf[, variantId := paste(Chromosome, Start_Position, Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2, sep = ':')];
		if(nrow(inputdf[duplicated(variantId)]) > 0) {
			cat("--Removed",  nrow(inputdf[duplicated(variantId)]) ,"duplicated variants\n");
			inputdf = inputdf[!duplicated(variantId)];
		}
		inputdf[,variantId := NULL];
		#// check empty value present in Hugo_Symbol (if yes assign 'UnknownGene')
		if(nrow(inputdf[Hugo_Symbol %in% ""]) > 0) {
			cat('--Found ', nrow(inputdf[Hugo_Symbol %in% ""]), ' variants with no Gene Symbols\n');
			inputdf$Hugo_Symbol = ifelse(test = inputdf$Hugo_Symbol == "", yes = 'UnknownGene', no = inputdf$Hugo_Symbol);
		}
		#// check NA present in Hugo_Symbol (if yes assign 'UnknownGene')
		if(nrow(inputdf[is.na(Hugo_Symbol)]) > 0) {
			cat('--Found ', nrow(inputdf[is.na(Hugo_Symbol) > 0]), ' variants with no Gene Symbols\n');
			inputdf$Hugo_Symbol = ifelse(test = is.na(inputdf$Hugo_Symbol), yes = 'UnknownGene', no = inputdf$Hugo_Symbol);
		}
		#// combine allele to Variant_Class
		inputdf = inputdf[, anno := as.character(paste(Reference_Allele, '>', Tumor_Seq_Allele2, '(', Variant_Classification, ')', sep = ''))];
		#// convert data.table to data.frame, and count freq of 'Hugo_Symbol-anno-Turmo_Sample_Barcode' per 'Chrom-Position'
		inputdf = inputdf[, c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "anno")];
	}
	return(inputdf);
}

