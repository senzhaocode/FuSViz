#' Convert chromome notation (e.g. '1' to 'chr1') 
#'
#' This function it to convert chromosome notation in import bam file, which meets requirement of Gviz::AlignmentsTrack function
#'
#' @param file path of read alignment (BAM format).
#' @param selection A GRanges object defines a region for coverage plotting using the bamfile.
#'
#' @return A GRanges object with read alignment in the selection region.
chrom_name_function <- function(file, selection) {
	index=paste(file, ".bai", sep="")
	flag = Rsamtools::scanBamFlag(isUnmappedQuery=F)	

	# Rename chromosome names from ensembl to UCSC notation
	chrname = sub("chr", "", as.character(selection@seqnames@values));	chrname = sub("M", "MT", chrname)
	selection@seqnames@values = factor(chrname)
	chrname = sub("chr", "", as.character(selection@seqinfo@seqnames));	chrname = sub("M", "MT", chrname)
	selection@seqinfo@seqnames = chrname

	#// extract read alignment -- please use 'GenomicAlignments::readGAlignmentsList'
	param = Rsamtools::ScanBamParam(flag = flag, tag = "MD", what = Rsamtools::scanBamWhat(), which = selection)
	read_align = unlist(GenomicAlignments::readGAlignmentsList(file=file, index=index, use.names=FALSE, param=param))

	# change @elementMetadata$rname
	tmp = as.character(read_align@elementMetadata$rname)
	tmp = sapply(tmp, simplify=TRUE, USE.NAMES=FALSE, function(x){
		if (! is.na(as.character(x)) ) { 
			if ( as.character(x) == "MT") { return("chrM") } else { return(paste("chr", x, sep="")) } 
		} else {
			return(NA)
		}})
	# shift mapped sequence width to avoid out of [start, end]
	if ( length(read_align@elementMetadata$seq) ) {
		seq_shift = GenomicAlignments::sequenceLayer(read_align@elementMetadata$seq, read_align@elementMetadata$cigar)
		seq_new = Biostrings::stackStrings(seq_shift, start(selection@ranges), end(selection@ranges),
			shift = read_align@elementMetadata$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
		names(seq_new) = seq_along(read_align@elementMetadata$qname)
	} else {
		seq_new = Biostrings::DNAStringSet()
	}

	return(
		GenomicRanges::GRanges(seqnames = tmp, strand = read_align@elementMetadata$strand, ranges = IRanges::IRanges(start = read_align@elementMetadata$pos, 
			width = read_align@elementMetadata$qwidth), id = read_align@elementMetadata$qname, cigar = read_align@elementMetadata$cigar, 
			mapq = read_align@elementMetadata$mapq, flag = read_align@elementMetadata$flag, md = read_align@elementMetadata$MD, seq = seq_new, 
			isize = read_align@elementMetadata$isize, groupid = read_align@elementMetadata$groupid, status = read_align@elementMetadata$mate_status))
}

#' Plot fusion transcript in context of chromosome ideogram, gene annotation and read coverage
#' 
#' @description Plot curved lines to link breakpoints of fusion transcripts in context of gene annotation and read coverage
#'
#' @param first_name A string - symbol name of geneA (e.g. 'TMPRSS2').
#' @param second_name A string - symbol name of geneB (e.g. 'ERG').
#' @param breakpoint_A A number - genomic coordinate of breakpoint in geneA.
#' @param breakpoint_B A number - genomic coordinate of breakpoint in geneB.
#' @param rna_bam_path A string - path of alignment file from RNA-seq (BAM format).
#' @param dna_bam_path A string - path of alignment file from DNA-seq (BAM format), default value: NULL.
#' @param chrom_notation_rna A boolean type - if TRUE (default value) chromosome name in RNA-seq alignment denotes like 'chrX'; if FALSE chromosome denotes like 'X'.
#' @param chrom_notation_dna A boolean type - if TRUE (default value) chromosome name in DNA-seq alignment denotes like 'chrX'; if FALSE chromosome denotes like 'X'.
#' @param coverage_plot_trans A boolean type - if TRUE, RNA read coverage calculated by reads mapped to selected transcripts of geneA and geneB; otherwise (default: FALSE) RNA read coverage calculated by reads mapped to genomic regions of geneA and geneB. coverage_plot_trans is only valid for RNA-seq alignment.
#' @param offset A number - the value (default: 5) indicates how many bases are extended from 5'-end and 3'-end of geneA and geneB in plotting.
#' @param split A number - split reads support (default: NULL).
#' @param span A number - discordant spanning read pair support (default: NULL).
#' @param transcriptA A vector of string - selected transcripts (ensembl id) of geneA for plotting (default: NULL).
#' @param transcriptB A vector of string - selected transcripts (ensembl id) of geneB for plotting (default: NULL).
#' @param fusion_strandA A string - the strand direction of geneA in transcribed fusion (default: NULL; options: '+' or '-').
#' @param fusion_strandB A string - the strand direction of geneB in transcribed fusion (default: NULL; options: '+' or '-').
#' @param coverage_max_A A number - the range of y-axis (RNA-seq read coverage) scale for geneA (default: NULL).
#' @param coverage_max_B A number - the range of y-axis (RNA-seq read coverage) scale for geneB (default: NULL).
#' 
#' @return NULL
#' 
#' @export
plot_separate_individual_bam <- function(first_name, second_name, breakpoint_A, breakpoint_B, rna_bam_path, dna_bam_path=NULL, coverage_plot_trans=F, offset=5, chrom_notation_rna=T,
	chrom_notation_dna=T, split=NULL, span=NULL, transcriptA=NULL, transcriptB=NULL, fusion_strandA=NULL, fusion_strandB=NULL, coverage_max_A=NULL, coverage_max_B=NULL) {

	stopifnot(is(first_name, "character")); 
	stopifnot(is(second_name, "character"));
	stopifnot(is(breakpoint_A, "numeric"));
	stopifnot(is(breakpoint_B, "numeric"));
	stopifnot(is(rna_bam_path, "character"));
	stopifnot(is(offset, "numeric"));

	if (! is.null(dna_bam_path) ) { stopifnot(is(dna_bam_path, "character")); }
	if (! is.null(split) ) { stopifnot(is(split, "numeric")); }
	if (! is.null(span) ) { stopifnot(is(span, "numeric")); }
	if (! is.null(transcriptA) ) { stopifnot(is(transcriptA, "character")); }
	if (! is.null(transcriptB) ) { stopifnot(is(transcriptB, "character")); }
	if (! is.null(fusion_strandA) ) { stopifnot( fusion_strandA == '+' || fusion_strandA == '-' ); }
	if (! is.null(fusion_strandB) ) { stopifnot( fusion_strandB == '+' || fusion_strandB == '-' ); }
	if (! is.null(coverage_max_A) ) { stopifnot(is(coverage_max_A, "numeric")); }
	if (! is.null(coverage_max_B) ) { stopifnot(is(coverage_max_B, "numeric")); }

	ens_A = names(symbol_ensem[symbol_ensem==first_name])
	ens_B = names(symbol_ensem[symbol_ensem==second_name])
	if ( length(ens_A) > 0 ) { 
		object_individual_A <- get_annotation_db(ens_A, txdb, grTrack)
	} else { #// if symbol_A invalid
		stop(paste(ens_A, " is invalid or matched ensembl_id in database!"));
	}
	if ( length(ens_B) > 0 ) { 
		object_individual_B <- get_annotation_db(ens_B, txdb, grTrack)
	} else { #// if symbol_B invalid
		stop(paste(ens_B, " is invalid or matched ensembl_id in database!"));
	}

	first = list(); second = list(); # // geneA and geneB are list structure
	# select transcript for geneA
	if (! is.null(transcriptA) ) {
		object_individual_A$dataset = object_individual_A$dataset[object_individual_A$dataset$TXNAME %in% transcriptA, ]
		if ( nrow(object_individual_A$dataset) == 0 ) {
			stop(paste(transcriptA, " not present in annotation database!"));
		}
	}
	first[[as.character(breakpoint_A)]] <- gene_trans_ex(breakpoint_A, object_individual_A, whole_txdb)
	# select transcript for geneB
	if (! is.null(transcriptB) ) {
		object_individual_B$dataset = object_individual_B$dataset[object_individual_B$dataset$TXNAME %in% transcriptB, ]
		if ( nrow(object_individual_B$dataset) == 0 ) {
			stop(paste(transcriptB, " not present in annotation database!"));
		}
	}
	second[[as.character(breakpoint_B)]] <- gene_trans_ex(breakpoint_B, object_individual_B, whole_txdb)

	#// get chromosome name for 'first' and 'second'
	chrom_f = first[[1]]$transcript$Chrom[1];
	chrom_s = second[[1]]$transcript$Chrom[1];

	#// Set plotting tracks for geneA
	grTrack_f = Gviz::GeneRegionTrack(first[[1]]$select_region, chromosome=chrom_f, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", 
		just.group="below", fontsize=5, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="green", col.symbol="black", lwd=0.1, 
		lwd.border=1, lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);	
	first_vis_s = min(first[[1]]$transcript$GStart) #// start position of geneA for visualization
	first_vis_e = max(first[[1]]$transcript$GEnd) #// end position of geneA for visualization

	#// Set plotting tracks for geneB
	grTrack_s = Gviz::GeneRegionTrack(second[[1]]$select_region, chromosome=chrom_s, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", 
		just.group="below", fontsize=5, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="orange", col.symbol="black", lwd=0.1, 
		lwd.border=1, lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);
	second_vis_s = min(second[[1]]$transcript$GStart)	#// visualization start position of geneB
	second_vis_e = max(second[[1]]$transcript$GEnd) #// visualization start position of geneB

	#// load bam file from RNA-seq: the AlignmentsTrack can only read the bam file with chrom name 'chr' (Ensembl reference not recognised correctly)
	if ( coverage_plot_trans == F ) { #// coverage calculated by all reads mapped to regions of upstream and downstream genes
		if ( chrom_notation_rna == T ) {
			alTrack_f = Gviz::AlignmentsTrack(rna_bam_path, isPaired=TRUE, name="RNA read coverage")
		} else {  # ensembl chromosome name like '1'
			alTrack_f = Gviz::AlignmentsTrack(rna_bam_path, isPaired=TRUE, name="RNA read coverage", importFunction = chrom_name_function)
		}
		Gviz::displayPars(alTrack_f) = list(col.coverage="black", fill.coverage="green", col.axis="black", type="coverage", strand = '*', coverageHeight=0.08, 
			minCoverageHeight=5, col.title="black", background.title="transparent", background.panel="transparent", lty=1, lwd=0.2, cex=0.6, cex.axis=0.6, 
			cex.title=0.6, lwd.coverage=0.2, lty.coverage=1, fontsize=9, stackHeight=0.7, stacking="squish", min.width=0.1, min.height=3, min.distance=0)	

		if ( chrom_notation_rna == T ) {
			alTrack_s = Gviz::AlignmentsTrack(rna_bam_path, isPaired=TRUE, name="RNA read coverage")
		} else { # ensembl chromosome name like '1'
			alTrack_s = Gviz::AlignmentsTrack(rna_bam_path, isPaired=TRUE, name="RNA read coverage", importFunction = chrom_name_function)
		}
		Gviz::displayPars(alTrack_s) = list(col.coverage="black", fill.coverage="orange", col.axis="black", type="coverage", strand = '*', coverageHeight=0.08, 
			minCoverageHeight=5, col.title="black", background.title="transparent", background.panel="transparent", lty=1, lwd=0.2, cex=0.6, cex.axis=0.6, 
			cex.title=0.6, lwd.coverage=0.2, lty.coverage=1, fontsize=9, stackHeight=0.7, stacking="squish", min.width=0.1, min.height=3, min.distance=0)
	} else { #// coverage calculated by reads mapped to selected transcripts of upstream and downstrean genes
		flag = Rsamtools::scanBamFlag(isPaired=T, isProperPair=T, isUnmappedQuery=F, isDuplicate=F, isSecondaryAlignment=F)
		# Get first coverage
		df_A = GenomicRanges::makeGRangesFromDataFrame(first[[1]]$select_region, keep.extra.columns=TRUE)
		first_trans = unique(df_A@elementMetadata$transcript); # get all transcript id of geneA
		collapse_trans_A = df_A[df_A@elementMetadata$transcript=="ZZZZZZZZ",]; # a GRange object - collapse exon intervals of all transcripts
		#// collapse exon interval - using union mode
		for ( id in first_trans ) { collapse_trans_A = GenomicRanges::union(collapse_trans_A, df_A[df_A@elementMetadata$transcript==id, ]); }
		if ( chrom_notation_rna != T ) { # ensembl chromosome name like '1'
			chrname = as.character(collapse_trans_A@seqnames@values)
			chrname = sub("chr", "", chrname)
			collapse_trans_A@seqnames@values = factor(chrname)
		}
		#// calculate the coverage
		cov_A = GenomicAlignments::coverage(rna_bam_path, param = Rsamtools::ScanBamParam(flag = flag, tag = "MD", which = collapse_trans_A))
		cov_A = cov_A[collapse_trans_A]
		#// convert object collapse_trans_A from a GRange to a data.frame
		collapse_trans_A = as.data.frame(collapse_trans_A)
		bedgraph_A = NULL;
		if ( nrow(collapse_trans_A) == length(cov_A) ) {
			for ( i in 1:nrow(collapse_trans_A) ) {
				value = as.numeric(cov_A[[i]])
				end = seq(collapse_trans_A[i,]$start, collapse_trans_A[i,]$end, by = 1)
				start = end - 1;
				tmp = data.frame(start, end, value, stringsAsFactors = F)
				bedgraph_A = rbind(bedgraph_A, tmp) 
			}
		}
		bedgraph_A=unique(bedgraph_A)
		alTrack_f <- Gviz::DataTrack(range = bedgraph_A, chromosome = as.character(first[[1]]$select_region[1,]$seqnames), 
			strand = as.character(first[[1]]$select_region[1,]$strand), genome="NA", name="RNA read coverage", type="h", col="green", fill="green", 
			col.axis="black", col.title="black", background.title="transparent", background.panel="transparent", 
			lty=1, lwd=0.2, cex=0.6, cex.axis=0.6, cex.title=0.6, fontsize=9)

		# Get second coverage
		df_B = GenomicRanges::makeGRangesFromDataFrame(second[[1]]$select_region, keep.extra.columns=TRUE)
		second_trans = unique(df_B@elementMetadata$transcript); # get all transcript id of geneA
		collapse_trans_B = df_B[df_B@elementMetadata$transcript=="ZZZZZZZZ",]; # a GRange object - collapse exon intervals of all transcripts
		#// collapse exon interval - using union mode
		for ( id in second_trans ) { collapse_trans_B = GenomicRanges::union(collapse_trans_B, df_B[df_B@elementMetadata$transcript==id, ]); }
		if ( chrom_notation_rna != T ) { # ensembl chromosome name like '1'
			chrname = as.character(collapse_trans_B@seqnames@values)
			chrname = sub("chr", "", chrname)
			collapse_trans_B@seqnames@values = factor(chrname)
		}
		#// calculate the coverage
		cov_B = GenomicAlignments::coverage(rna_bam_path, param = Rsamtools::ScanBamParam(flag = flag, tag = "MD", which = collapse_trans_B))
		cov_B = cov_B[collapse_trans_B]
		#// convert object collapse_trans_A from a GRange to a data.frame
		collapse_trans_B = as.data.frame(collapse_trans_B)
		bedgraph_B = NULL;
		if ( nrow(collapse_trans_B) == length(cov_B) ) {
			for ( i in 1:nrow(collapse_trans_B) ) {
				value = as.numeric(cov_B[[i]])
				end = seq(collapse_trans_B[i,]$start, collapse_trans_B[i,]$end, by = 1)
				start = end - 1;
				tmp = data.frame(start, end, value, stringsAsFactors = F)
				bedgraph_B = rbind(bedgraph_B, tmp) 
			}
		}
		bedgraph_B=unique(bedgraph_B)
		alTrack_s <- Gviz::DataTrack(range = bedgraph_B, chromosome = as.character(second[[1]]$select_region[1,]$seqnames), 
			strand = as.character(second[[1]]$select_region[1,]$strand), genome="NA", name="RNA read coverage", type="h", col="orange", fill="orange",
			col.axis="black", col.title="black", background.title="transparent", background.panel="transparent", 
			lty=1, lwd=0.2, cex=0.6, cex.axis=0.6, cex.title=0.6, fontsize=9)
	}

	#// load bam file from DNA-seq: the AlignmentsTrack can only read the bam file with chrom name 'chr' (Ensembl reference not recognised correctly)
	if (! is.null(dna_bam_path) ) {
		if ( chrom_notation_dna == T ) {
			alTrack_f_dna = Gviz::AlignmentsTrack(dna_bam_path, isPaired=TRUE, name="DNA read coverage")
		} else {
			alTrack_f_dna = Gviz::AlignmentsTrack(dna_bam_path, isPaired=TRUE, name="DNA read coverage", importFunction = chrom_name_function)
		}
		Gviz::displayPars(alTrack_f_dna) = list(col.coverage="black", fill.coverage="green", col.axis="black", type="coverage", strand = '*', coverageHeight=0.08, 
				minCoverageHeight=5, col.title="black", background.title="transparent", background.panel="transparent", lty=1, lwd=0.2, cex=0.6, cex.axis=0.6, 
				cex.title=0.6, lwd.coverage=0.2, lty.coverage=1, fontsize=9, stackHeight=0.7, stacking="squish", min.width=0.1, min.height=3, min.distance=0)

		if ( chrom_notation_dna == T ) {
			alTrack_s_dna = Gviz::AlignmentsTrack(dna_bam_path, isPaired=TRUE, name="DNA read coverage")
		} else {
			alTrack_s_dna = Gviz::AlignmentsTrack(dna_bam_path, isPaired=TRUE, name="DNA read coverage", importFunction = chrom_name_function)
		}
		Gviz::displayPars(alTrack_s_dna) = list(col.coverage="black", fill.coverage="orange", col.axis="black", type="coverage", strand = '*', coverageHeight=0.08, 
				minCoverageHeight=5, col.title="black", background.title="transparent", background.panel="transparent", lty=1, lwd=0.2, cex=0.6, cex.axis=0.6, 
				cex.title=0.6, lwd.coverage=0.2, lty.coverage=1, fontsize=9, stackHeight=0.7, stacking="squish", min.width=0.1, min.height=3, min.distance=0)
	}

	#// set read coverage ylim scale
	if (! is.null(coverage_max_A) ) { Gviz::displayPars(alTrack_f) = list(ylim = c(0, coverage_max_A[1])) }
	if (! is.null(coverage_max_B) ) { Gviz::displayPars(alTrack_s) = list(ylim = c(0, coverage_max_B[1])) }

	#// adjust font size of transcript_id (It can be adjusted by customers in UI in further version)
	if ( length(first[[1]]$transcript$TXNAME) > 0 &&  length(first[[1]]$transcript$TXNAME) <= 5 ) { 
		grTrack_f@dp@pars$fontsize.group = 11; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 5 &&	 length(first[[1]]$transcript$TXNAME) <= 10 ) { 
		grTrack_f@dp@pars$fontsize.group = 9; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 10 &&  length(first[[1]]$transcript$TXNAME) <= 15 ) {
		grTrack_f@dp@pars$fontsize.group = 8; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 15 &&  length(first[[1]]$transcript$TXNAME) <= 20 ) {
		grTrack_f@dp@pars$fontsize.group = 7;
	} else if ( length(first[[1]]$transcript$TXNAME) > 20 &&  length(first[[1]]$transcript$TXNAME) <= 25 ) {
		grTrack_f@dp@pars$fontsize.group = 6;
	} else {
		grTrack_f@dp@pars$fontsize.group = 4;
	}
	if ( length(second[[1]]$transcript$TXNAME) > 0 &&  length(second[[1]]$transcript$TXNAME) <= 5 ) { 
		grTrack_s@dp@pars$fontsize.group = 11; 
	} else if ( length(second[[1]]$transcript$TXNAME) > 5 &&  length(second[[1]]$transcript$TXNAME) <= 10 ) { 
		grTrack_s@dp@pars$fontsize.group = 9; 
	} else if ( length(second[[1]]$transcript$TXNAME) > 10 &&  length(second[[1]]$transcript$TXNAME) <= 15 ) {
		grTrack_s@dp@pars$fontsize.group = 8; 
	} else if ( length(second[[1]]$transcript$TXNAME) > 15 &&  length(second[[1]]$transcript$TXNAME) <= 20 ) {
		grTrack_s@dp@pars$fontsize.group = 7;
	} else if ( length(second[[1]]$transcript$TXNAME) > 20 &&  length(second[[1]]$transcript$TXNAME) <= 25 ) {
		grTrack_s@dp@pars$fontsize.group = 6;
	} else {
		grTrack_s@dp@pars$fontsize.group = 4;
	}

	#// set the chromosome axis coordinate
	axis = Gviz::GenomeAxisTrack(name="Axis", labelPos="alternating", fontcolor="black", littleTicks=F, cex=0.5, cex.id=0.5, cex.axis=0.5, 
		fontsize=12, col="black", distFromAxis=0.4, lwd=1.4, lwd.border=1, min.width=0.1, min.height=3, min.distance=0); 

	#//////////////////////////////////////////////////////////////////////
	#// perform the plot ('chromosome', 'title' and 'gene+axis' tracks)	 //
	#//////////////////////////////////////////////////////////////////////
	# offset = 10; #// define plot region [start-offset, end+offset]
	#// plot 'chromosome', 'title' and 'gene+axis' tracks in proportion as 1 : 0.5 : 8 (can be adjusted)
	nf = grid::grid.layout(nrow=4, ncol=2, widths=grid::unit(c(1, 1), "null"), heights = grid::unit(c(0.8, 0.2, 0.4, 10), "null"))
	grid::grid.newpage(); #// start a new page
	grid::pushViewport(grid::viewport(layout = nf)); #// divide the page to different section (nf specifies the section composite)

	if ( chrom_f == chrom_s ) { #// if both geneA and geneB at the same chromosome
		#// plot cytoband of the chromosome
		names(chrTrack) = chrom_f;
		grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = NULL)); #// load on the section (1, null)
		if ( first[[1]]$transcript$GStart[1] < second[[1]]$transcript$GStart[1] ) {
			Gviz::plotTracks(chrTrack, chromosome=chrom_f, from=first_vis_s, to=second_vis_e, add=T, main=paste("Fusion: ", first_name, "-", second_name, sep=""),
				cex.main=0.6, fontsize=13, fontface.main=2, margin=10, innerMargin=6)
		} else {
			Gviz::plotTracks(chrTrack, chromosome=chrom_s, from=second_vis_s, to=first_vis_e, add=T, main=paste("Fusion: ", first_name, "-", second_name, sep=""),
				cex.main=0.6, fontsize=13, fontface.main=2, margin=10, innerMargin=6)
		}
		grid::popViewport(1); #// close the section (1, null)
	} else { #// if geneA and geneB at the different chromosomes
		chrTrack_A = chrTrack;	names(chrTrack_A) = chrom_f;
		chrTrack_B = chrTrack;	names(chrTrack_B) = chrom_s;
		grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = NULL)); #// load on the section (1, null)
		grid::grid.text(label=paste("Fusion: ", first_name, "-", second_name, sep=""), x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), 
				gp=grid::gpar(fontsize=10, fontface="bold", col="black"));
		grid::popViewport(1); #// close the section (1, null)
		#// plot cytoband of geneA chromosome	
		grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1)); #// load on the section (1, 1)
		Gviz::plotTracks(chrTrack_A, chromosome=chrom_f, from=first_vis_s, to=first_vis_e, add=T, cex.main=0.6, fontsize=13, fontface.main=2, margin=10, innerMargin=6);
		grid::popViewport(1); #// close the section (1, 1)
		#// plot cytoband of geneB chromosome 
		grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2)); #// load on the section (1, 2)
		Gviz::plotTracks(chrTrack_B, chromosome=chrom_s, from=second_vis_s, to=second_vis_e, add=T, cex.main=0.6, fontsize=13, fontface.main=2, margin=10, innerMargin=6);
		grid::popViewport(1); #// close the section (1, 2)
	}
	#// add sub-title for geneA
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1)); #// load on the section (2, 1)
	grid::grid.text(label=paste("Upstream gene: ", first_name, " (", first[[1]]$transcript$Strand[1], ")", sep=""), 
			x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), gp=grid::gpar(fontsize=9, fontface="bold", col="black"));
	grid::popViewport(1); #// close the section (2, 1)
	#// add sub-title for geneB
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2)); #// load on the section (2, 2)
	grid::grid.text(label=paste("Downstream gene: ", second_name, " (", second[[1]]$transcript$Strand[1], ")", sep=""), 
			x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), gp=grid::gpar(fontsize=9, fontface="bold", col="black"));
	grid::popViewport(1); #// close the section (2, 2)

	#// plot 'geneA' and 'axis' tracks in proportion as 1 : 0.2 (can be adjusted)
	grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = 1)) #// load on the section (3, 1)
	if ( first[[1]]$transcript$Strand[1] == "+" ) {
		if ( nrow(first[[1]]$transcript) < 3 && nrow(second[[1]]$transcript) < 3 ) {
			if ( is.null(dna_bam_path) ) {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 3, 0.6), add=T);
			} else {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, alTrack_f_dna, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 1.5, 1.5, 0.6), add=T);
			}
		} else {
			if ( is.null(dna_bam_path) ) {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 0.6), add=T);
			} else {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, alTrack_f_dna, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 1, 0.6), add=T);
			}
		}
	} else {
		if ( nrow(first[[1]]$transcript) < 3 && nrow(second[[1]]$transcript) < 3 ) {
			if ( is.null(dna_bam_path) ) {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 3, 0.6), reverseStrand=TRUE, add=T);
			} else {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, alTrack_f_dna, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 1.5, 1.5, 0.6), reverseStrand=TRUE, add=T);
			}
		} else {
			if ( is.null(dna_bam_path) ) {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 0.6), reverseStrand=TRUE, add=T);
			} else {
				A1_xy = Gviz::plotTracks(c(grTrack_f, alTrack_f, alTrack_f_dna, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 1, 0.6), reverseStrand=TRUE, add=T);
			}
		}
	}
	grid::popViewport(1); #// close the section (3, 1)
	#// plot 'geneB' and 'axis' tracks in proportion as 1 : 0.2 (can be adjusted)
	grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = 2)) #// load on the section (3, 2)
	if ( second[[1]]$transcript$Strand[1] == "+" ) {
		if ( nrow(first[[1]]$transcript) < 3 && nrow(second[[1]]$transcript) < 3 ) {
			if ( is.null(dna_bam_path) ) {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 3, 0.6), add=T);
			} else {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, alTrack_s_dna, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 1.5, 1.5, 0.6), add=T);
			}
		} else {
			if ( is.null(dna_bam_path) ) {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 0.6), add=T);
			} else {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, alTrack_s_dna, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 1, 0.6), add=T);
			}
		}
	} else {
		if ( nrow(first[[1]]$transcript) < 3 && nrow(second[[1]]$transcript) < 3 ) {
			if ( is.null(dna_bam_path) ) {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 3, 0.6), reverseStrand=TRUE, add=T);
			} else {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, alTrack_s_dna, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(1, 1.5, 1.5, 0.6), reverseStrand=TRUE, add=T);
			}
		} else {
			if ( is.null(dna_bam_path) ) {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 0.6), reverseStrand=TRUE, add=T);
			} else {
				B1_xy = Gviz::plotTracks(c(grTrack_s, alTrack_s, alTrack_s_dna, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
					margin=1, innerMargin=1, sizes=c(3, 1, 1, 0.6), reverseStrand=TRUE, add=T);
			}
		}
	}
	grid::popViewport(1); #// close the section (3, 2)

	#////////////////////////////////////////
	#// perform the plot of one curve line //
	#////////////////////////////////////////
	#// Control the degree of bezier curve (can be adjusted)
	bezier_control_point_offset = 32
	#// 'A1_xy$GeneRegionTrack' and 'B1_xy$GeneRegionTrack': coordinates(x1, y1, x2, y2) of all exons of transcripts for geneA and geneB
	breakpoint = data.frame(breakpoint_A=breakpoint_A, breakpoint_B=breakpoint_B, stringsAsFactors=FALSE)
	curve_coordinate = .plot_curve(A1_xy, first, B1_xy, second, breakpoint);
	#// 'curve_coordinate' is a vector class with three elements (i.e. x_pos_gene_upstream: x-coordinate of geneA breakpoint,
	#//  x_pos_gene_downstream: x-coordinate of geneB breakpoint, y_pos: y-coordinate of geneA and geneB breakpoint)

	#// Draw the bezier curve between transcripts
	grid::pushViewport(grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)));
	grid::grid.bezier(rep(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), each = 2), #// x coordinates
			c(curve_coordinate["y_pos"], curve_coordinate["y_pos"] - bezier_control_point_offset, 
			curve_coordinate["y_pos"] - bezier_control_point_offset, curve_coordinate["y_pos"]),
			default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), gp = grid::gpar(col="red", fill="red", lwd=1));

	number_position_x = curve_coordinate["x_pos_gene_upstream"] + (curve_coordinate["x_pos_gene_downstream"] - curve_coordinate["x_pos_gene_upstream"]) / 2
	number_position_y = curve_coordinate["y_pos"] - bezier_control_point_offset + 4
	#// 'supporting_reads_text' is (num of split reads, num of span reads)
	supporting_reads_text = paste("(", split, ", ", span, ")", sep="")
	grid::grid.text(supporting_reads_text, x = number_position_x / grDevices::dev.size(units = "px")[1], y = 1 - number_position_y / grDevices::dev.size(units = "px")[2],
			vp = grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)),
			gp = grid::gpar(fontsize = 7));
	
	#// Add arrow line showing transcription direction of upstream part
	if (! is.null(fusion_strandA) ) { #// breakpoint strand of geneA available
		if ( first[[1]]$transcript$Strand[1] == fusion_strandA ) {
			grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"] - 12, curve_coordinate["x_pos_gene_upstream"] - 2), #// x positions - left margin 15
    			c(curve_coordinate["y_pos"] - 2, curve_coordinate["y_pos"] - 2), #// y positions
        		default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), gp = grid::gpar(col="red", fill="red", lwd=1));
		} else {
			grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"] + 12, curve_coordinate["x_pos_gene_upstream"] + 2), #// x positions - right margin 255
    			c(curve_coordinate["y_pos"] - 2, curve_coordinate["y_pos"] - 2), #// y positions
        		default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), gp = grid::gpar(col="red", fill="red", lwd=1));
		}
	}
	#// Add arrow line showing transcription direction of downstream part
	if (! is.null(fusion_strandB) ) { #// breakpoint strand of geneA available
		if ( second[[1]]$transcript$Strand[1] == fusion_strandB ) {
			grid::grid.lines(c(curve_coordinate["x_pos_gene_downstream"] + 3, curve_coordinate["x_pos_gene_downstream"] + 13), #// x positions -right margin 501
    			c(curve_coordinate["y_pos"] - 2, curve_coordinate["y_pos"] - 2), #// y positions
        		default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), gp = grid::gpar(col="red", fill="red", lwd=1));
		} else {
			grid::grid.lines(c(curve_coordinate["x_pos_gene_downstream"] - 3, curve_coordinate["x_pos_gene_downstream"] - 13), #// x positions - legt margin 285
    			c(curve_coordinate["y_pos"] - 2, curve_coordinate["y_pos"] - 2), #// y positions
        		default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), gp = grid::gpar(col="red", fill="red", lwd=1));
		}
	}
	#// Add dash lines showing the breakpoint of upstream partner
	grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_upstream"]), #// x positions
    	c(curve_coordinate["y_pos"] - 5, 450), #// y positions
        default.units = "native", arrow=NULL, gp = grid::gpar(col="black", lty=3, lwd=0.3));
	#// Add dash lines showing the breakpoint of downstream partner
	grid::grid.lines(c(curve_coordinate["x_pos_gene_downstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
    	c(curve_coordinate["y_pos"] - 5, 450), #// y positions
        default.units = "native", arrow=NULL, gp = grid::gpar(col="black", lty=3, lwd=0.3));
}

